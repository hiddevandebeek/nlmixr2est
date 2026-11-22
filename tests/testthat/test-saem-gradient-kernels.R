test_that("ODE gradient kernels support analytic and numerical scores", {
  mod <- function() {
    ini({
      tka <- fix(log(1.4))
      tcl <- log(4)
      tv <- fix(log(25))
      eta.cl ~ 0.2
      err <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(err)
    })
  }

  data <- subset(nlmixr2data::theo_sd, ID %in% unique(ID)[1:4])
  fit <- function(method, numerical = FALSE) {
    nu <- c(prior = 0L, rw = 0L, coord = 0L, mala = 0L, hmc = 0L)
    nu[[method]] <- 1L
    penalty <- if (numerical) list(nlmixr2LikPenaltySpec("etaL2", par = c(0, 0))) else NULL
    suppressWarnings(nlmixr2(
      mod, data, "saem",
      saemControl(
        nBurn = 3L, nEm = 4L, nmc = 2L, nu = nu, hmcSteps = 1L,
        seed = 101L, print = 0, covMethod = "", calcTables = FALSE,
        handleUninformativeEtas = FALSE
      ),
      likPenalty = penalty
    ))
  }
  for (numerical in c(FALSE, TRUE)) {
    mala <- fit("mala", numerical)
    hmc <- fit("hmc", numerical)
    score <- if (numerical) "finiteDifference" else "sensitivity"
    expect_identical(mala$env$saem$malaScore, score)
    expect_identical(hmc$env$saem$hmcScore, score)
    expect_gt(mala$env$saem$malaCounts[["accepted"]], 0)
    expect_identical(mala$env$saem$malaCounts, hmc$env$saem$hmcCounts)
    expect_identical(mala$env$phiM, hmc$env$phiM)
  }
})

test_that("SAEM uses event-jump sensitivities for an ETA-dependent dose lag", {
  mod <- function() {
    ini({
      tka <- fix(log(1.4))
      tcl <- log(4)
      tv <- fix(log(25))
      tlag <- log(0.15)
      eta.cl ~ 0.15
      eta.lag ~ 0.1
      err <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      alag(depot) <- exp(tlag + eta.lag)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(err)
    })
  }
  data <- subset(nlmixr2data::theo_sd, ID %in% unique(ID)[1:4])
  ctl <- saemControl(
    nBurn = 3L, nEm = 4L, nmc = 2L, nu = c(mala = 1L), malaStep = 0.08,
    seed = 99L, print = 0, covMethod = "", calcTables = FALSE,
    handleUninformativeEtas = FALSE
  )
  analytic <- suppressWarnings(nlmixr2(mod, data, "saem", ctl))
  numerical <- suppressWarnings(nlmixr2(
    mod, data, "saem", ctl,
    likPenalty = list(nlmixr2LikPenaltySpec("etaL2", par = c(0, 0)))
  ))
  expect_identical(analytic$env$saem$malaScore, "sensitivity")
  expect_identical(numerical$env$saem$malaScore, "finiteDifference")
  expect_equal(analytic$env$phiM, numerical$env$phiM, tolerance = 1e-6)
})
