#' @export
#' @import methods
print.sdmTMB <- function(x, ...) {
  r <- suppressWarnings(tryCatch(x$tmb_obj$report(), error = function(e) NA))
  if (all(is.na(r))) {
    stop("It looks like the model was built with a different version of sdmTMB. ",
      "Please fit your model again with the current sdmTMB.", call. = FALSE)
  }
  # need to initialize the new TMB object once:
  sink(tempfile())
  x$tmb_obj$fn(x$tmb_obj$par)
  lp <- x$tmb_obj$env$last.par.best
  r <- x$tmb_obj$report(lp)
  sink()

  spatial_only <- as.logical(x$tmb_data$spatial_only)

  fit_by <- "ML"
  if ("reml" %in% names(x)) { # for backwards compatibility
    if (isTRUE(x$reml)) fit_by <- "REML" else "ML"
  }

  if (isTRUE(spatial_only)) {
    title <- paste0("Spatial model fit by ", fit_by, " ['sdmTMB']\n")
  } else {
    title <- paste0("Spatiotemporal model fit by ", fit_by, " ['sdmTMB']\n")
  }
  if (x$control$map_rf) {
    title <- paste0("Model fit by ", fit_by, " ['sdmTMB']\n")
  }
  formula <- paste0("Formula: ", deparse(x$call$formula), "\n")
  if (deparse(x$call$time) != "NULL") {
    time <- paste0("Time column: ", deparse(x$call$time), "\n")
  } else {
    time <- NULL
  }
  spde <- paste0("Mesh: ", deparse(x$call$mesh), "\n")
  data <- paste0("Data: ", deparse(x$call$data), "\n")
  family <- paste0("Family: ", paste0(x$family$family, "(link = '", x$family$link, "')"), "\n")
  criterion <- paste0(fit_by, " criterion at convergence: ", mround(x$model$objective, 3), "\n")

  # .formula <- check_and_parse_thresh_params(x$formula, x$data)$formula
  .formula <- x$split_formula$fixedFormula
  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  fe_names <- colnames(x$tmb_data$X_ij)
  fe_names <- fe_names[!fe_names == "offset"]

  pars <- x$model$par
  b_j_exact <- unname(pars[grep("b_j", names(pars))])
  b_j <- round(b_j_exact, 2L)

  if ("ln_phi" %in% names(as.list(pars))) {
    phi <- mround(exp(as.list(pars)$ln_phi), 2L)
    phi <- paste0("Dispersion parameter: ", phi, "\n")
  } else {
    phi <- ""
  }

  sr <- x$sd_report
  sr_est <- as.list(sr, "Estimate")
  if (x$family$family == "tweedie") {
    tweedie_p <- paste0("Tweedie p: ", mround(plogis(sr_est$thetaf) + 1, 2L), "\n")
  } else {
    tweedie_p <- ""
  }

  range <- mround(r$range, 2L)

  pre <- "Spatial SD: "
  if (!is.null(r$sigma_O)) {
    sigma_O <- paste0(pre, mround(r$sigma_O, 2L), "\n")
  } else {
    sigma_O <- ""
  }

  pre <- "Spatiotemporal SD: "
  if (x$tmb_data$spatial_only == 0L) {
    if (!isTRUE(is.na(x$tmb_map$b_epsilon))) {
      sigma_E <- paste0(pre, mround(r$sigma_E, 2L), "\n")
    } else {
      sigma_E <- paste0(pre, mround(r$sigma_E[1], 2L), "\n")
    }
  } else {
    sigma_E <- NULL
  }

  pre <- "Spatiotemporal AR1 correlation (rho): "
  if (!is.null(r$rho) && r$rho != 0L) {
    rho <- paste0(pre, mround(r$rho, 2L), "\n")
  } else {
    rho <- ""
  }

  sr_se <- summary(sr)[, "Std. Error"]
  sr_est <- summary(sr)[, "Estimate"]
  b_j_se <- unname(round(sr_se[grep("b_j", names(sr_se))], 2L))
  b_j <- unname(round(sr_est[grep("b_j", names(sr_est))], 2L))

  mm <- cbind(b_j, b_j_se)
  colnames(mm) <- c("coef.est", "coef.se")
  row.names(mm) <- fe_names

  if (x$tmb_data$has_smooths) {
    bs_se <- unname(round(sr_se[grep("bs", names(sr_se))], 2L))
    bs <- unname(round(sr_est[grep("bs", names(sr_est))], 2L))
    sm <- parse_smoothers(formula = x$formula, data = x$data)
    sm_names <- unlist(lapply(sm$Zs, function(x) attr(x, "s.label")))
    sm_names <- gsub("\\)$", "", gsub("s\\(", "", sm_names))
    sm_names_bs <- paste0("s", sm_names)
    xx <- lapply(sm_names_bs, function(.x) { # split out 2D + smooths
      n_sm <- grep(",", .x) + 1
      if (length(n_sm)) {
        .x <- paste(.x, seq_len(n_sm), sep = "_")
        gsub(",", "", .x)
      } else {
        .x
      }
    })
    sm_names_bs <- unlist(xx)
    sm_names_sds <- paste0("sds(", sm_names, ")")
    mm_sm <- cbind(bs, bs_se)
    row.names(mm_sm) <- sm_names_bs
    mm <- rbind(mm, mm_sm)
    smooth_sds <- round(exp(unname(sr_est[grep("ln_smooth_sigma", names(sr_est))])), 2L)
    re_sm_mat <- matrix(NA_real_, nrow = length(smooth_sds), ncol = 1L)
    re_sm_mat[,1] <- smooth_sds
    rownames(re_sm_mat) <- sm_names_sds
    colnames(re_sm_mat) <- "Std. Dev."
  } else {
    re_sm_mat <- NULL
  }

  sr_se <- as.list(sr, "Std. Error")
  sr_est <- as.list(sr, "Estimate")

  if (x$tmb_data$threshold_func > 0) {
    mm_thresh <- cbind(sr_est$b_threshold, sr_se$b_threshold)
    if (x$threshold_function == 1L) {
      row.names(mm_thresh) <- paste0(x$threshold_parameter, c("-slope", "-breakpt"))
    } else {
      row.names(mm_thresh) <- paste0(x$threshold_parameter, c("-s50", "-s95", "-smax"))
    }
    colnames(mm_thresh) <- c("coef.est", "coef.se")
    mm_thresh[,1] <- round(mm_thresh[,1], 2)
    mm_thresh[,2] <- round(mm_thresh[,2], 2)

    mm <- rbind(mm, mm_thresh)
  } else {
    mm_thresh <- NULL
  }

  .tidy <- tidy(x, "ran_pars")
  if ("tau_G" %in% .tidy$term) {
    re_int_names <- barnames(x$split_formula$reTrmFormulas)
    re_int_mat <- matrix(NA_real_, nrow = length(re_int_names), ncol = 1)
    re_int_mat[,1] <- round(.tidy$estimate[.tidy$term == "tau_G"], 2)
    # rownames(re_int_mat) <- paste(re_int_names, "(Intercept)")
    rownames(re_int_mat) <- re_int_names
    colnames(re_int_mat) <- "Std. Dev."
  } else {
    re_int_mat <- NULL
  }

  if (!is.null(x$time_varying)) {
    tv_names <- colnames(model.matrix(x$time_varying, x$data))
    mm_tv <- cbind(round(as.numeric(sr_est$b_rw_t), 2), round(as.numeric(sr_se$b_rw_t), 2))
    colnames(mm_tv) <- c("coef.est", "coef.se")
    time_slices <- sort(unique(x$data[[x$time]]))
    row.names(mm_tv) <- paste(rep(tv_names, each = length(time_slices)), time_slices, sep = "-")
  } else {
    mm_tv <- NULL
  }

  cat(title,
    formula,
    time,
    spde,
    data,
    family,
    sep = ""
  )

  print(mm)

  if (!is.null(re_int_mat)) {
    cat("\nRandom intercepts:\n")
    print(re_int_mat)
  }

  if (!is.null(re_sm_mat)) {
    cat("\nSmooth terms:\n")
    print(re_sm_mat)
  }

  if (!is.null(x$time_varying)) {
    cat("\nTime-varying parameters:\n" )
    print(mm_tv)
  }

  range_text <- if (x$tmb_data$share_range) {
    paste0("Matern range: ", range[1], "\n")
  } else {
    paste0("Matern range (spatial): ", range[1], "\n",
      "Matern range (spatiotemporal): ", range[2], "\n")
  }
  if (x$control$map_rf) {
    range_text <- NULL
  }

  cat("\n",
    phi,
    tweedie_p,
    range_text,
    sigma_O,
    if (!is.null(sigma_E)) sigma_E,
    rho,
    criterion,
    "\nSee ?tidy.sdmTMB to extract these values as a data frame.\n",
    sep = ""
  )

  invisible(list(fe = mm, tv = mm_tv, thresh = mm_thresh))
}

#' @export
summary.sdmTMB <- function(object, ..., digits) {
  print(object, ...)
}

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

#' Extract the number of observations of a sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats nobs
#' @export
#' @noRd
nobs.sdmTMB <- function(object, ...)
  sum(!is.na(object$data[all.vars(object$formula)[1]]))

#' Extract the log likelihood of a sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats logLik
#' @export
#' @noRd
logLik.sdmTMB <- function(object, ...) {
  val <- -object$model$objective

  nobs <- nobs.sdmTMB(object)
  df <- length(object$model$par) # fixed effects only
  structure(val,
    nobs = nobs, nall = nobs, df = df,
    class = "logLik"
  )
}

#' Extract the AIC of a sdmTMB model
#'
#' @param fit The fitted sdmTMB model
#' @param scale The scale (note used)
#' @param k Penalization parameter, defaults to 2
#' @param ... Anything else
#' @noRd
#'
#' @export
extractAIC.sdmTMB <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L, "df")
  return(c(edf, c(-2 * L + k * edf)))
}
