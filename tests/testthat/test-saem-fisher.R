test_that("Fisher jointly updates the mean and residual variance on fixed draws", {
  mod <- function() {
    ini({
      tbase <- 1.5
      tslope <- 0.6
      eta.base ~ 0.2
      add.sd <- 0.4
    })
    model({
      muBase <- tbase + eta.base
      pred <- muBase + tslope * TREND
      pred ~ add(add.sd)
    })
  }

  data <- expand.grid(TIME = 1:5, ID = 1:12)
  data$EVID <- 0L
  data$TREND <- data$TIME
  data$DV <- with(data, 1.5 + 0.75 * TIME + 0.2 * sin(ID) + 0.3 * cos(TIME * ID))
  fit <- function(opt, observations = data) {
    suppressWarnings(nlmixr2(
      mod, observations, "saem",
      saemControl(
        nBurn = 1L, nEm = 1L, nmc = 2L, seed = 1042L, print = 0,
        nonMuThetaOpt = opt, covMethod = "", calcTables = FALSE,
        handleUninformativeEtas = FALSE, nonMuThetaMaxEval = 40L
      )
    ))
  }

  baseline <- fit("newuoa")
  fisher <- fit("fisher")
  cfg <- attr(fisher$env$saem, "saem.cfg")
  draws <- fisher$env$phiM[, , 2L, cfg$i1 + 1L, drop = FALSE]
  expect_identical(draws, baseline$env$phiM[, , 2L, cfg$i1 + 1L, drop = FALSE])
  expect_identical(fisher$env$saem$phi0OptimizerCounts[["fisherConverged"]], 1L)

  time <- rep(data$TIME, cfg$nmc)
  response <- rep(data$DV, cfg$nmc) - rep(as.numeric(draws), each = 5L)
  slope <- sum(time * response) / sum(time^2)
  variance <- mean((response - slope * time)^2)
  hist <- fisher$env$saem$par_hist
  gain <- 0.5 # first EM step after one burn-in iteration
  expect_equal(unname(fixef(fisher)[["tslope"]]),
               hist[1L, 2L] + gain * (slope - hist[1L, 2L]), tolerance = 5e-4)
  expect_equal(unname(fixef(fisher)[["add.sd"]])^2,
               hist[1L, ncol(hist)]^2 * (1 - gain) + gain * variance, tolerance = 5e-4)
  censored <- data
  censored$CENS <- as.integer(censored$DV < 2.6)
  censored$DV[censored$CENS == 1L] <- 2.6
  censFit <- fit("fisher", censored)
  expect_identical(censFit$env$saem$phi0OptimizerCounts[["fisherConverged"]], 1L)
  censDraws <- censFit$env$phiM[, , 2L, cfg$i1 + 1L, drop = FALSE]
  y <- rep(censored$DV, cfg$nmc)
  isCens <- rep(censored$CENS == 1L, cfg$nmc)
  condNll <- function(par) {
    mean <- rep(as.numeric(censDraws), each = 5L) + par[1L] * time
    sd <- exp(par[2L])
    -sum(dnorm(y[!isCens], mean[!isCens], sd, log = TRUE)) -
      sum(pnorm(y[isCens], mean[isCens], sd, log.p = TRUE))
  }
  censHist <- censFit$env$saem$par_hist
  sdStart <- censHist[1L, ncol(censHist)]
  slopeStart <- censHist[1L, 2L]
  optimum <- optim(c(slopeStart, log(sdStart)), condNll, method = "BFGS",
                   control = list(reltol = 1e-12))$par
  expect_equal(unname(fixef(censFit)[["tslope"]]),
               slopeStart + gain * (optimum[1L] - slopeStart), tolerance = 5e-4)
  expect_equal(unname(fixef(censFit)[["add.sd"]])^2,
               sdStart^2 * (1 - gain) + gain * exp(2 * optimum[2L]), tolerance = 5e-4)

})

test_that("Fisher updates multiple residual endpoints without free structural parameters", {
  mod <- function() {
    ini({
      tbase <- 1.5
      eta.base ~ 0.2
      err1 <- 0.3
      err2 <- 0.5
    })
    model({
      muBase <- tbase + eta.base
      out1 <- muBase
      out2 <- 2 * muBase
      out1 ~ add(err1)
      out2 ~ add(err2)
    })
  }

  data <- expand.grid(TIME = seq_len(5), DVID = 1:2, ID = seq_len(12))
  data <- data[order(data$ID, data$TIME, data$DVID), ]
  data$EVID <- 0L
  data$DV <- with(data, DVID * (1.5 + 0.2 * sin(ID)) + 0.15 * DVID * cos(TIME * ID))
  fit <- suppressWarnings(nlmixr2(
    mod, data, "saem",
    saemControl(
      nBurn = 1L, nEm = 1L, nmc = 2L, seed = 1042L, print = 0,
      nonMuThetaOpt = "fisher", covMethod = "", calcTables = FALSE,
      handleUninformativeEtas = FALSE
    )
  ))
  counts <- fit$env$saem$phi0OptimizerCounts
  expect_identical(counts[["fisherConverged"]], 1L)
  expect_identical(counts[["evaluations"]], 1L)
  expect_true(all(fixef(fit)[c("err1", "err2")] > 0))
})

test_that("joint Fisher scoring uses the draws and responsibilities of SAEM mixtures", {
  mod <- function() {
    ini({
      tka <- log(1.4)
      tcl1 <- log(2)
      tcl2 <- log(5)
      tv <- log(27)
      ecl1 ~ 0.15
      ecl2 ~ 0.15
      p1 <- 0.5
      err <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- mix(exp(tcl1 + ecl1), p1, exp(tcl2 + ecl2))
      v <- exp(tv)
      linCmt() ~ add(err)
    })
  }

  time <- c(0, 0.5, 1, 2, 5, 9, 15, 24)
  data <- do.call(rbind, lapply(seq_len(16L), function(id) {
    k <- c(2, 6)[1L + id %% 2L] * exp(0.15 * sin(id)) / 27
    pred <- 350 / 27 * 1.4 / (1.4 - k) * (exp(-k * time) - exp(-1.4 * time))
    data.frame(ID = id, TIME = time, DV = c(NA, pred[-1L] + 0.2 * cos(id * time[-1L])),
               EVID = c(1L, rep(0L, length(time) - 1L)), AMT = c(350, rep(0, length(time) - 1L)))
  }))
  for (sampler in c("parallel", "msaem")) {
    fit <- suppressWarnings(nlmixr2(
      mod, data, "saem",
      saemControl(
        nBurn = 5L, nEm = 8L, nmc = 2L, seed = 99L, print = 0,
        nonMuThetaOpt = "fisher", mixProbMethod = "annealed", mixSampleMethod = sampler,
        covMethod = "", calcTables = FALSE, handleUninformativeEtas = FALSE
      )
    ))
    raw <- fit$env$saem
    expect_gt(raw$phi0OptimizerCounts[["fisherConverged"]], 0L)
    expect_equal(unname(rowSums(raw$mixWeights)), rep(1, 16L))
    expect_true(all(colSums(raw$mixWeights) > 0))
  }
})
