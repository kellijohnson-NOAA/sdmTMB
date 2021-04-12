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
plot_diagnostics_sdmTMB <- function(object, newdata = NULL, plot_type="all", plot = TRUE) {
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


#' Plot diagnostics generic
#'
#' @param pred Predicted object or dataframe that can be coerced into a data frame
#' @param X Character string labeling the name of the 'X' variable, defaults to 'X'
#' @param Y Character string labeling the name of the 'Y' variable, defaults to 'Y'
#' @param time Character string labeling the name of the 'time' variable, defaults to 'time'
#' @param predicted Character string labeling the name of the 'predicted' variable, defaults to 'predicted'
#' @param se Character string labeling the name of the 'se' variable, defaults to NULL. This is optional, but
#' can be easily generated for many classes of models (lm, glm, glmmTMB, sdmTMB, gam, etc). If se is included
#' then transparancy of plots will be set to (1/se^2) to highlight residual points that may be large in magnitude
#' and precise. If not included, all points will be shown with the same transparency.
#' @param y Character string labeling the name of the 'y' variable, defaults to 'y'
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
#' @import gstat
#' @importFrom stats sd as.formula
#' @export
#' @examples
#' \donttest{
#' d <- pcod
#' m <- sdmTMB(data = d,
#'  formula = density ~ 0 + as.factor(year),
#'  time = "year", spde = make_mesh(d, c("X", "Y"), n_knots = 80, type = "kmeans"),
#'  family = tweedie(link = "log"), anisotropy = TRUE,
#'  include_spatial = FALSE)
#' d$predicted = predict(m)$est
#' d$obs = d$density
#' plot_diagnostics_generic(d, plot_type="all")
#'
#' }
plot_diagnostics_generic <- function(pred, X = "X", Y = "Y", time = "year",
                predicted = "predicted",
                se = NULL,
                y = "obs",
                plot_type = "all",
                plot = TRUE) {

  # check object
  pred <- as.data.frame(pred)

  pred$size = 1
  size_range <- c(0,1)
  if(!is.na(se)) {
    pred$size <- (1/(pred[[se]]^2))
    size_range <- range(pred$size,na.rm=T) / max(pred$size,na.rm=T) * 3.0
  }

  # rename columns
  colnames(pred)[which(colnames(pred) == X)] <- "X"
  colnames(pred)[which(colnames(pred) == Y)] <- "Y"
  colnames(pred)[which(colnames(pred) == time)] <- "time"
  colnames(pred)[which(colnames(pred) == predicted)] <- "predicted"
  colnames(pred)[which(colnames(pred) == y)] <- "y"

  # sort data by time, and then location
  pred <- pred[with(pred, order(time, X, Y)), ]

  # summarise
  n_raw <- nrow(pred)
  time_steps <- unique(pred[["time"]])
  n_time_steps <- length(time_steps)
  locs <- paste(pred[["X"]], pred[["Y"]], sep = ":")
  unique_locs <- unique(locs)
  n_locs <- length(unique_locs)

  # check to see if predictions have been made at same locations for each time step
  if (n_locs * n_time_steps != n_raw) {
    # create new grid based on all combinations of obsered / predicted
    new_grid <- expand.grid(locs = unique_locs, time = time_steps)
    new_grid$X <- as.numeric(unlist(lapply(strsplit(as.character(new_grid$locs), ":"), getElement, 1)))
    new_grid$Y <- as.numeric(unlist(lapply(strsplit(as.character(new_grid$locs), ":"), getElement, 2)))
    new_grid$y <- NA
    new_grid$predicted <- NA
    # if creating the new grid, we don't have any SE informa
    new_grid$size = 1
    size_range <- c(0,1)
    # then we need to do spatial interpolation. Lots of options -- sdmTMB, mgcv, gstat IDW, fastTps, etc.
    for (t in 1:n_time_steps) {
      idx <- which(pred[["time"]] == time_steps[t])
      f <- gstat(formula = y ~ 1, locations = ~ X + Y, nmax = 5, set = list(idp = 0), data = pred[idx, ])
      # predict to new locs
      idx_pred <- which(new_grid[["time"]] == time_steps[t])
      new_grid$predicted[idx_pred] <- try(predict(f, newdata = new_grid[idx_pred, ])$var1.pred, silent = TRUE)
    }

    # last, need to align observations with new df
    new_grid$y[match(paste(pred$time, pred$X, pred$Y), paste(new_grid$time, new_grid$X, new_grid$Y))] <- pred$y

    # reset pred object to be in regular format
    pred <- new_grid
  }

  # calculate residuals
  pred$resid <- pred$y - pred$predicted

  # Calculate mean spatial field for predictions and residuals. The predictions can be
  # centered by time and space, so that
  # pred = year_fixed_effect + mean_spatial_effect + spatio_temporal_effect
  # This is sort of like modeling in sdmTMB -- but these effects have fixed/random effects confounded
  spatial_mean <- pred[which(pred$time == time_steps[1]), ]
  spatial_mean$predicted <- 0
  spatial_mean$resid <- 0
  for (t in 1:n_time_steps) {
    idx <- which(pred$time == time_steps[t])
    # remove year_fixed_effect by subtracting off mean
    spatial_mean$pred <- spatial_mean$predicted + (pred$predicted[idx] - mean(pred$predicted[idx], na.rm = T))
    spatial_mean$resid <- spatial_mean$resid + (pred$resid[idx] - mean(pred$resid[idx], na.rm = T))
  }
  spatial_mean$predicted <- spatial_mean$predicted / n_time_steps
  spatial_mean$resid <- spatial_mean$resid / n_time_steps

  # calculate first differenced predictions and residuals -- this should help
  # diagnose need for AR(1) and MA(1) components
  # 'pred' data frame already sorted by time, loc x, loc y -- so this
  # should work as long as there's not > 1 observation per location-time
  pred$pred_st <- NA
  pred$diff_pred_st <- NA
  pred$diff_resid <- NA

  for (t in 1:n_time_steps) {
    idx_t <- which(pred$time == time_steps[t])
    # subtract off the year_fixed_effect and mean_spatial_effect
    pred$pred_st[idx_t] <- pred$predicted[idx_t] - mean(pred$predicted[idx_t]) - spatial_mean$predicted
    if (t > 1) {
      idx_t_minus_1 <- which(pred$time == time_steps[t - 1])
      pred$diff_pred_st[idx_t] <- pred$pred_st[idx_t] - pred$pred_st[idx_t_minus_1]
      pred$diff_resid[idx_t] <- pred$resid[idx_t] - pred$resid[idx_t_minus_1]
    }
  }

  # Calculate variance of predictions, and residuals, by time step
  sd_df <- data.frame(
    "time" = time_steps,
    "pred_sd" = NA,
    "resid_sd" = NA
  )
  for (t in 1:length(time_steps)) {
    sd_df$pred_sd[t] <- sd(pred$pred_st[which(pred$time == time_steps[t])])
    sd_df$resid_sd[t] <- sd(pred$resid[which(pred$time == time_steps[t])])
  }

  plots <- list()
  i <- 0
  # qq plots
  if (plot_type %in% c("all", "qq")) {
    i <- i + 1
    plots[[i]] <- ggplot(pred, aes(sample = resid)) +
      stat_qq_line() +
      stat_qq(alpha = 0.4) +
      theme_bw() +
      xlab("Theoretical") +
      ylab("Sample")
    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  # QQ plot by time step is useful to look at whether some of the individual
  # time steps are better / worse than others
  if (plot_type %in% c("all", "qq_over_time")) {
    i <- i + 1
    plots[[i]] <- ggplot(pred, aes(sample = resid)) +
      stat_qq_line() +
      stat_qq(alpha = 0.4) +
      theme_bw() +
      facet_wrap(~time, scale = "free_y") +
      theme(strip.background = element_rect(fill = "white")) +
      xlab("Theoretical") +
      ylab("Sample")
    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  # # spatial intercept -- average across all years, minus fixed effects
  if (plot_type %in% c("all", "spatial_intercept")) {
    # take mean across all time steps

    i <- i + 1
    plots[[i]] <- ggplot(spatial_mean, aes(X, Y, col = pred)) +
      geom_point(aes(size=size), alpha = 0.3) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("") +
      labs(color = "Spatial effect") +
      ggtitle("Mean spatial predictions (centered across time steps)") +
      scale_size(range = size_range)

    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if (plot_type %in% c("all", "spatial_residuals")) {
    # take mean across all time steps
    if (!is.na(sd(sd(spatial_mean$resid)))) {
      i <- i + 1
      plots[[i]] <- ggplot(spatial_mean, aes(X, Y, col = resid)) +
        geom_point(alpha = 0.3) +
        scale_color_gradient2() +
        theme_bw() +
        xlab("") +
        ylab("") +
        labs(color = "Residuals") +
        ggtitle("Mean residuals (centered across time steps)")
      if (plot == TRUE) {
        print(plots[[i]])
        cat("Hit <Return> to see next plot:")
        line <- readline()
      }
    }
  }

  # spatiotemporal components
  if (plot_type %in% c("all", "st_over_time") & n_time_steps > 1) {
    i <- i + 1
    plots[[i]] <- ggplot(pred, aes(time, pred_st)) +
      geom_hline(aes(yintercept = 0), col = "red", alpha = 0.5) +
      geom_point(position = position_dodge2(0.3), alpha = 0.3) +
      theme_bw() +
      xlab("") +
      ylab("Effects") +
      ggtitle("Spatiotemporal effects, centered by time and space")
    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  if (plot_type %in% c("all", "st_over_spacetime")) {
    i <- i + 1
    plots[[i]] <- ggplot(
      pred[which(!is.na(pred$pred_st)), ],
      aes(x = X, y = Y, col = pred_st)) +
      geom_point(aes(size=size), alpha = 0.3) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("") +
      facet_wrap(~time) +
      theme(strip.background = element_rect(fill = "white")) +
      labs(color = "Effects") +
      ggtitle("Spatiotemporal effects, centered by time and space") +
      scale_size(range = size_range)
    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if (plot_type %in% c("all", "st_variance")) {
    i <- i + 1

    plots[[i]] <-
      ggplot(sd_df, aes(time, pred_sd)) +
      geom_point(size = 3) +
      xlab("") +
      ylab("Standard deviation of predictions") +
      geom_smooth() +
      theme_bw()
    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }
  }

  if (plot_type %in% c("all", "residuals_over_time")) {
    i <- i + 1
    plots[[i]] <-
      ggplot(pred, aes(time, resid)) +
      geom_hline(aes(yintercept = 0), col = "red", alpha = 0.5) +
      geom_point(position = position_dodge2(0.3), alpha = 0.3) +
      theme_bw() +
      xlab("") +
      ylab("Residuals")
    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }
  }
  if (plot_type %in% c("all", "residual_variance")) {
    if (length(which(!is.na(sd_df$resid_sd)) > 0)) {
      i <- i + 1

      plots[[i]] <-
        ggplot(sd_df, aes(time, resid_sd)) +
        geom_point(size = 3) +
        xlab("") +
        ylab("Standard deviation of residuals") +
        geom_smooth() +
        theme_bw()
      if (plot == TRUE) {
        print(plots[[i]])
        cat("Hit <Return> to see next plot:")
        line <- readline()
      }
    }
  }
  if (plot_type %in% c("all", "residuals_over_spacetime")) {
    i <- i + 1
    plots[[i]] <-
      ggplot(pred[which(!is.na(pred$resid)), ], aes(x = X, y = Y, color = resid)) +
      geom_point(alpha = 0.5) +
      scale_color_gradient2() +
      theme_bw() +
      xlab("") +
      ylab("") +
      facet_wrap(~time) +
      theme(strip.background = element_rect(fill = "white")) +
      ggtitle("Residuals (observed - predicted)")
    if (plot == TRUE) {
      print(plots[[i]])
      cat("Hit <Return> to see next plot:")
      line <- readline()
    }

    # First differenced predictions are useful at diagnosing correlation, but the fixed effects
    # are also confounded
    if (plot_type %in% c("all", "diff_pred_over_time")) {
      i <- i + 1
      plots[[i]] <-
        ggplot(pred, aes(time, diff_pred_st)) +
        geom_hline(aes(yintercept = 0), col = "red", alpha = 0.5) +
        geom_point(position = position_dodge2(0.3), alpha = 0.3) +
        theme_bw() +
        xlab("") +
        ylab("First differenced predictions")
      if (plot == TRUE) {
        print(plots[[i]])
        cat("Hit <Return> to see next plot:")
        line <- readline()
      }
    }
    if (plot_type %in% c("all", "diff_pred_over_spacetime")) {
      i <- i + 1
      plots[[i]] <-
        ggplot(pred[which(!is.na(pred$diff_pred_st)), ], aes(x = X, y = Y, color = diff_pred_st)) +
        geom_point(alpha = 0.5) +
        scale_color_gradient2() +
        theme_bw() +
        xlab("") +
        ylab("") +
        facet_wrap(~time) +
        theme(strip.background = element_rect(fill = "white")) +
        labs(color = "E[t+1] - E[t]") +
        ggtitle("First differenced predictions")
      if (plot == TRUE) {
        print(plots[[i]])
        cat("Hit <Return> to see next plot:")
        line <- readline()
      }
    }

    # 1st differenced residual. This onlyl works if observations are at same
    # stations or on grid
    if (length(which(!is.na(pred$diff_resid))) > 0) {
      if (plot_type %in% c("all", "diff_resid_over_time")) {
        i <- i + 1
        plots[[i]] <-
          ggplot(pred, aes(time, diff_resid)) +
          geom_hline(aes(yintercept = 0), col = "red", alpha = 0.5) +
          geom_point(position = position_dodge2(0.3), alpha = 0.3) +
          theme_bw() +
          xlab("") +
          ylab("1st differenced residuals")
        if (plot == TRUE) {
          print(plots[[i]])
          cat("Hit <Return> to see next plot:")
          line <- readline()
        }
      }

      if (plot_type %in% c("all", "diff_resid_over_spacetime")) {
        i <- i + 1
        plots[[i]] <-
          ggplot(pred[which(!is.na(pred$diff_resid)), ], aes(x = X, y = Y, color = diff_resid)) +
          geom_point(alpha = 0.5) +
          scale_color_gradient2() +
          theme_bw() +
          xlab("") +
          ylab("") +
          facet_wrap(~time) +
          theme(strip.background = element_rect(fill = "white")) +
          ggtitle("1st differenced residuals")
        if (plot == TRUE) {
          print(plots[[i]])
          cat("Hit <Return> to see next plot:")
          line <- readline()
        }
      }
    }
  }

  # return list of ggplot objects
  return(plots)
}

