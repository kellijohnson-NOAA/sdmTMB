#' Plot anisotropy
#'
#' @param object An object from [sdmTMB()].
#'
#' @export
#' @rdname plot_anisotropy
#' @examples
#' \donttest{
#' d <- pcod
#' m <- sdmTMB(data = d,
#'   formula = density ~ 0 + as.factor(year),
#'   time = "year", spde = make_mesh(d, c("X", "Y"), n_knots = 80, type = "kmeans"),
#'   family = tweedie(link = "log"), anisotropy = TRUE,
#'   include_spatial = FALSE)
#' plot_anisotropy(m)
#' }
plot_anisotropy <- function(object) {
  stopifnot(identical(class(object), "sdmTMB"))
  report <- object$tmb_obj$report()
  eig <- eigen(report$H)
  dat <- data.frame(
    x0 = c(0, 0),
    y0 = c(0, 0),
    x1 = eig$vectors[1, , drop = TRUE] * eig$values,
    y1 = eig$vectors[2, , drop = TRUE] * eig$values
  )
  plot(0, xlim = range(c(dat$x0, dat$x1)),
    ylim = range(c(dat$y0, dat$y1)), type = "n", asp = 1, xlab = "", ylab = "")
  graphics::arrows(dat$x0, dat$y0, dat$x1, dat$y1)
  invisible(list(eig = eig, dat = dat, H = report$H))
}


#' Plot diagnostics
#'
#' @param object An object from [sdmTMB()]
#' @param newdata If passed in then predictions are made at the locations of newdata, otherwise
#' predictions are made at the original locations of the data
#' @param plot_type A string deciding which plots to show. If 'all' (default) then all plots are shown.
#' Plots can be shown individually by specifying 'qq' (QQ plot of randomized residuals),
#' 'qq_over_time' (QQ plot of randomized residuals faceted by time step), 'spatial_slope' (shows spatial
#' field for zeta_s, if included), 'spatial_intercept' (shows omega_s, the constant spatial field or spatial
#' intercept, if included), 'st_over_time' (spatiotemporal effects shown across time),
#' 'st_over_space' (spatial plots of epsilon, across all time steps), 'st_over_spacetime' (spatial plots
#' of epsilon, faceted by time steps), 'st_variance' (the standard deviation of epsilon deviations over
#' time steps), 'residuals_over_time' (distribution of randomized residuals by time step),
#' 'residuals_over_spacetime' (randomized residuals over space and faceted by time step)
#' @param plot Whether to show plot the results (default = TRUE)
#'
#' @return
#' A list of ggplot plot objects that can be individually printed to the screen or manipulated
#'
#' @import ggplot2
#' @importFrom stats sd as.formula
#' @export
#' @examples
#' \donttest{
# d <- pcod
# m <- sdmTMB(data = d,
#   formula = density ~ 0 + as.factor(year),
#   time = "year", spde = make_mesh(d, c("X", "Y"), n_knots = 80, type = "kmeans"),
#   family = tweedie(link = "log"), anisotropy = TRUE,
#   include_spatial = FALSE)
# plot_diagnostics(m, plot_type="all")
#' }
plot_diagnostics <- function(object, newdata = NULL, plot_type="all", plot = TRUE) {
  # generate predictions
  pred <- predict(object, newdata = newdata, return_tmb_object = FALSE)
  # generate randomized quantile residuals
  pred$resid <- residuals(object)
  # pred$est_normal = object$family$linkinv(pred$est)
  # if(is.null(newdata)) {
  #   pred$resid_normal = object$response - pred$est_normal
  # }

  plots <- list()
  i <- 0
  # qq plots
  if(plot_type %in% c("all", "qq")) {
    i <- i+1
    plots[[i]] <- ggplot(pred, aes(sample=resid)) +
      stat_qq_line() +
      stat_qq(alpha=0.4) +
      theme_bw() +
      xlab("Theoretical") +
      ylab("Sample")
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if(plot_type %in% c("all", "qq_over_time")) {
    i <- i+1
    plots[[i]] <- ggplot(pred, aes(sample=resid)) +
      stat_qq_line() +
      stat_qq(alpha=0.4) +
      theme_bw() +
      facet_wrap(as.formula(paste("~", object$time))) +
      theme(strip.background =element_rect(fill="white")) +
      xlab("Theoretical") +
      ylab("Sample")
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  # spatial slope
  if(plot_type %in% c("all", "spatial_slope") & sd(pred$zeta_s) != 0) {
    # only grab data from first time slice
    i <- i+1
    plots[[i]] <- ggplot(pred[which(pred[[object$time]] == min(pred[[object$time]])),], aes_string(x = names(pred)[2], y = names(pred)[3], col = "zeta_s")) +
      geom_point(alpha = 0.3) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("")
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  # spatial slope
  if(plot_type %in% c("all", "spatial_intercept") & sd(pred$omega_s) != 0) {
    # only grab data from first time slice
    i <- i+1
    plots[[i]] <- ggplot(pred[which(pred[[object$time]] == min(pred[[object$time]])),], aes_string(x = names(pred)[2], y = names(pred)[3], col = "omega_s")) +
      geom_point(alpha = 0.3) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("")
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  # spatiotemporal components
  if(plot_type %in% c("all", "st_over_time") & sd(pred$epsilon_st) != 0) {
    i <- i+1
    plots[[i]] <- ggplot(pred, aes_string(object$time, "epsilon_st")) +
      geom_hline(aes(yintercept=0),col="red",alpha=0.5) +
      geom_point(position=position_dodge2(0.3),alpha=0.3) +
      theme_bw() +
      xlab("") +
      ylab("Spatiotemporal deviations (epsilon)")
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if(plot_type %in% c("all", "st_over_space") & sd(pred$epsilon_st) != 0) {
    i <- i+1
    plots[[i]] <- ggplot(pred, aes_string(x = names(pred)[2], y = names(pred)[3], col = "epsilon_st")) +
      geom_point(alpha = 0.3) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("")
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if(plot_type %in% c("all", "st_over_spacetime") & sd(pred$epsilon_st) != 0) {
    i <- i+1
    plots[[i]] <- ggplot(pred, aes_string(x = names(pred)[2], y = names(pred)[3], col = "epsilon_st")) +
      geom_point(alpha = 0.5) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("") +
      facet_wrap(as.formula(paste("~", object$time))) +
      theme(strip.background =element_rect(fill="white"))
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if(plot_type %in% c("all", "st_variance") & sd(pred$epsilon_st) != 0) {
    i <- i+1

    # do this with a loop to avoid depedency on dplyr etc
    time_steps = unique(pred[[object$time]])
    plot_df = data.frame("time" = time_steps,
                         "time_sd" = NA)
    for(j in 1:length(time_steps)) {
      plot_df$time_sd[j] = sd(pred[which(pred[[object$time]]==time_steps[j]),"epsilon_st"])
    }

    plots[[i]] <-
      ggplot(plot_df, aes(time,time_sd)) +
      geom_point(size=3) +
      xlab("") +
      ylab("Standard deviation of epsilon") +
      geom_smooth() +
      theme_bw()
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  if(plot_type %in% c("all", "residuals_over_time")) {
    i <- i+1
    plots[[i]] <-
      ggplot(pred, aes_string(object$time, "resid")) +
      geom_hline(aes(yintercept=0),col="red",alpha=0.5) +
      geom_point(position=position_dodge2(0.3),alpha=0.3) +
      theme_bw() +
      xlab("") +
      ylab("Randomized residuals")
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if(plot_type %in% c("all", "residuals_over_spacetime")) {
    i <- i+1
    plots[[i]] <-
      ggplot(pred, aes_string(x = names(pred)[2], y = names(pred)[3], col = "resid")) +
      geom_point(alpha = 0.5) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("") +
      facet_wrap(as.formula(paste("~", object$time))) +
      theme(strip.background =element_rect(fill="white"))
    if(plot==TRUE) {
      print(plots[[i]])
      cat ("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  # a couple other plots to consider would be the first-differenced epsilon values in space,
  # time, or the variance of the first differenced fields.

  # return list of ggplot objects
  return(plots)
}
