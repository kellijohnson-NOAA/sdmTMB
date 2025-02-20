% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stan.R
\name{extract_mcmc}
\alias{extract_mcmc}
\title{Extract MCMC samples from a model fit with \code{\link[tmbstan:tmbstan]{tmbstan::tmbstan()}}.}
\usage{
extract_mcmc(object)
}
\arguments{
\item{object}{Output from \code{\link[tmbstan:tmbstan]{tmbstan::tmbstan()}} run on the \code{tmb_obj}
element of an \code{\link[=sdmTMB]{sdmTMB()}} model. E.g., \code{tmbstan(your_model$tmb_obj)}.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Returns a matrix of parameter samples. Rows correspond to the order
of \code{your_model$tmb_obj$env$last.par.best}. Columns correspond to
posterior samples. Is used internally by \code{\link[=predict.sdmTMB]{predict.sdmTMB()}} to make
fully Bayesian predictions. See the \code{tmbstan_model} argument
in \code{\link[=predict.sdmTMB]{predict.sdmTMB()}}.
}
\examples{

\dontrun{
pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
plot(pcod_spde)

# here we will fix the random field parameters at their approximate
# MLEs (maximum likelihood estimates) from a previous fit
# to improve speed of convergence:
m_tmb <- sdmTMB(density ~ 0 + as.factor(year),
  data = pcod, mesh = pcod_spde, family = tweedie(link = "log"), time = "year",
  control = sdmTMBcontrol(start = list(ln_kappa = rep(-1.58, 2),
    ln_tau_E = -0.15, ln_tau_O = -0.65),
    map = list(ln_kappa = rep(factor(NA), 2),
      ln_tau_E = factor(NA), ln_tau_O = factor(NA))))
m_tmb

# will take 3-5 minutes:
library(tmbstan)
m_stan <- tmbstan(m_tmb$tmb_obj, iter = 200, chains = 1)
print(m_stan, pars = c("b_j", "thetaf", "ln_phi", "omega_s[1]", "epsilon_st[1]"))

post <- extract_mcmc(m_stan)
dim(post)

p <- predict(m_tmb, newdata = qcs_grid, tmbstan_model = m_stan)
p_last <- p[qcs_grid$year == max(qcs_grid$year), ] # just plot last year
pred <- qcs_grid[qcs_grid$year == max(qcs_grid$year), ]
pred$est <- apply(exp(p_last), 1, median)
pred$lwr <- apply(exp(p_last), 1, quantile, probs = 0.1)
pred$upr <- apply(exp(p_last), 1, quantile, probs = 0.9)
pred$cv <- apply(exp(p_last), 1, function(x) sd(x) / mean(x))

library(ggplot2)
ggplot(pred, aes(X, Y, fill = est)) + geom_raster() +
  scale_fill_viridis_c(trans = "log")
ggplot(pred, aes(X, Y, fill = cv)) + geom_raster() +
  scale_fill_viridis_c(trans = "log")

index_quantiles <- get_index_sims(p)
ggplot(index_quantiles, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_line() + geom_ribbon(alpha = 0.5)

index_samples <- get_index_sims(p, return_sims = TRUE)
ggplot(index_samples, aes(as.factor(year), .value)) +
  geom_violin()
}
}
