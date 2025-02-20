% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{sdmTMB}
\alias{sdmTMB}
\title{Fit a spatial or spatiotemporal GLMM with TMB}
\usage{
sdmTMB(
  formula,
  data,
  mesh,
  time = NULL,
  family = gaussian(link = "identity"),
  spatial = c("on", "off"),
  spatiotemporal = c("IID", "AR1", "RW", "off"),
  share_range = TRUE,
  time_varying = NULL,
  spatial_varying = NULL,
  weights = NULL,
  extra_time = NULL,
  reml = FALSE,
  silent = TRUE,
  anisotropy = FALSE,
  control = sdmTMBcontrol(),
  priors = sdmTMBpriors(),
  previous_fit = NULL,
  experimental = NULL,
  do_fit = TRUE,
  spatial_only = deprecated(),
  fields = deprecated(),
  include_spatial = deprecated(),
  spde = deprecated(),
  ...
)
}
\arguments{
\item{formula}{Model formula. See the Details section below for how to specify
offsets and threshold parameters. For index standardization, you may wish
to include \code{0 + as.factor(year)} (or whatever the time column is called)
in the formula. IID random intercepts are possible using \pkg{lme4}
syntax, e.g., \code{+ (1 | g)} where \code{g} is a column with factor levels.
Penalized splines are possible via \pkg{mgcv} with \code{s()}. See examples
and details below.}

\item{data}{A data frame.}

\item{mesh}{An object from \code{\link[=make_mesh]{make_mesh()}}.}

\item{time}{An optional time column name (as character). Can be left as
\code{NULL} for a model with only spatial random fields unless you wish to use
one of the index or center of gravity functions over time.}

\item{family}{The family and link. Supports \code{\link[=gaussian]{gaussian()}}, \code{\link[=Gamma]{Gamma()}},
\code{\link[=binomial]{binomial()}}, \code{\link[=poisson]{poisson()}}, \code{\link[sdmTMB:families]{Beta()}},
\code{\link[sdmTMB:families]{nbinom2()}},
\code{\link[sdmTMB:families]{truncated_nbinom2()}},
\code{\link[sdmTMB:families]{nbinom1()}},
\code{\link[sdmTMB:families]{truncated_nbinom1()}},
\code{\link[sdmTMB:families]{censored_poisson()}},
\code{\link[sdmTMB:families]{student()}}, and
\code{\link[sdmTMB:families]{tweedie()}}]. For binomial family options,
see the 'Binomial families' in the Details section below.}

\item{spatial}{Estimate spatial random fields? Options are
\code{'on'} / \code{'off'} or \code{TRUE} / \code{FALSE}.}

\item{spatiotemporal}{Estimate the spatiotemporal random fields as \code{'IID'}
(independent and identically distributed; default), stationary \code{'AR1'}
(first-order autoregressive), as a random walk (\code{'RW'}), or as fixed at 0
\code{'off'}. Will be set to \code{'off'} if \code{time = NULL}. Note that the
spatiotemporal standard deviation represents the marginal steady-state
standard deviation of the process in the case of the AR1. I.e., it is
scaled according to the correlation. See the \href{https://kaskr.github.io/adcomp/classAR1__t.html}{TMB documentation}. If the AR1
correlation coefficient (rho) is estimated close to 1, say > 0.99, then you
may wish to switch to the random walk \code{"RW"}. Capitalization is ignored. \code{TRUE}
gets converted to \code{'iid'} and \code{FALSE} gets converted to \code{off}.}

\item{share_range}{Logical: estimate a shared spatial and spatiotemporal
range parameter (\code{TRUE}) or independent range parameters (\code{FALSE}).}

\item{time_varying}{An optional one-sided formula describing covariates that
should be modelled as a random walk through time. Be careful not to include
covariates (including the intercept) in both the main and time-varying
formula. I.e., at least one should have \code{~ 0} or \code{~ -1}.}

\item{spatial_varying}{An optional one-sided formula \strong{with a single
predictor} of a coefficient that should varying in space as a random
field. Note that you may want to include a fixed effect for the same
variable to improve interpretability. If the (scaled) time column, will
represent a local-time-trend model. See \doi{10.1111/ecog.05176} and the
\href{https://pbs-assess.github.io/sdmTMB/articles/spatial-trend-models.html}{spatial trends vignette}.
Note this predictor should be centered to have mean zero and have a
standard deviation of approximately 1 (scale by the SD).}

\item{weights}{Optional likelihood weights for the conditional model.
Implemented as in \pkg{glmmTMB}. Weights do not have to sum
to one and are not internally modified. Can also be used for trials with
the binomial family; the weights argument needs to be a vector and not a name
of the variable in the dataframe. See the Details section below.}

\item{extra_time}{Optional extra time slices (e.g., years) to include for
interpolation or forecasting with the predict function. See the
Details section below.}

\item{reml}{Logical: use REML (restricted maximum likelihood) estimation
rather than maximum likelihood? Internally, this adds the fixed effects
to the list of random effects to integrate over.}

\item{silent}{Silent or include optimization details?}

\item{anisotropy}{Logical: allow for anisotropy? See \code{\link[=plot_anisotropy]{plot_anisotropy()}}.}

\item{control}{Optimization control options via \code{\link[=sdmTMBcontrol]{sdmTMBcontrol()}}.}

\item{priors}{Optional penalties/priors via \code{\link[=sdmTMBpriors]{sdmTMBpriors()}}.}

\item{previous_fit}{A previously fitted sdmTMB model to initialize the
optimization with. Can greatly speed up fitting. Note that the data and
model must be set up \emph{exactly} the same way. However, the \code{weights} argument
can change, which can be useful for cross-validation.}

\item{experimental}{A named list for esoteric or in-development options.
Here be dragons.
(Experimental) A column name (as character) of a
predictor of a linear trend (in log space) of the spatiotemporal standard
deviation. By default, this is \code{NULL} and fits a model with a constant
spatiotemporal variance. However, this argument can also be a character
name in the original data frame (a covariate that ideally has been
standardized to have mean 0 and standard deviation = 1). Because the
spatiotemporal field varies by time step, the standardization should be
done by time. If the name of a predictor is included, a log-linear model is
fit where the predictor is used to model effects on the standard deviation,
e.g. \code{log(sd(i)) = B0 + B1 * epsilon_predictor(i)}. The 'epsilon_model' argument may also
be specified. This is the name of the model to use to modeling time-varying epsilon. This
can be one of the following: "trend" (default, fits a linear model without random effects),
"re" (fits a model with random effects in epsilon_st, but no trend), and "trend-re" (a model
that includes both the trend and random effects)}

\item{do_fit}{Fit the model (\code{TRUE}) or return the processed data without
fitting (\code{FALSE})?}

\item{spatial_only}{\strong{Depreciated.} Replaced by \code{spatiotemporal = "off"}.}

\item{fields}{\strong{Depreciated.} Replaced by \code{spatiotemporal}.}

\item{include_spatial}{\strong{Depreciated.} Replaced by \code{spatial}.}

\item{spde}{\strong{Depreciated.} Replaced by \code{mesh}.}

\item{...}{Not currently used.}
}
\value{
An object (list) of class \code{sdmTMB}. Useful elements include:
\itemize{
\item \code{sd_report}: output from \code{\link[TMB:sdreport]{TMB::sdreport()}}
\item \code{gradients}: log likelihood gradients with respect to each fixed effect
\item \code{model}: output from \code{\link[stats:nlminb]{stats::nlminb()}}
\item \code{data}: the fitted data
\item \code{mesh}: the object that was supplied to the \code{mesh} argmument
\item \code{family}: the family object, which includes the inverse link function
\item \code{tmb_params}: The parameters list passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}
\item \code{tmb_map}: The 'map' list passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}
\item \code{tmb_data}: The data list passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}
\item \code{tmb_obj}: The TMB object created by \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}
}
}
\description{
Fit a spatial or spatiotemporal Gaussian Markov random field GLMM with TMB.
This can be useful for (dynamic) species distribution models and relative
abundance index standardization among many other uses.
}
\details{
\strong{Model description}

For now, see the
\href{https://pbs-assess.github.io/sdmTMB/articles/model-description.html}{model description}
vignette for a start. There are also descriptions of particular models in
Anderson et al. (2019) and Barnett et al. (2020) (see reference list below).

\strong{Offsets}

In the model formula, an offset can be included by including \code{+ offset} in
the model formula (a reserved word). The offset will be included in any
prediction. \code{offset} must be a column in \code{data}.

\strong{Binomial families}

Following the structure of \code{\link[stats:glm]{stats::glm()}} and \pkg{glmmTMB}, a binomial
family can be specified in one of 4 ways: (1) the response may be a factor
(and the model classifies the first level versus all others), (2) the
response may be binomial (0/1), (3) the response can be a matrix of form
\code{cbind(success, failure)}, and (4) the response may be the observed
proportions, and the 'weights' argument is used to specify the Binomial size
(N) parameter (\verb{prob ~ ..., weights = N}).

\strong{Smooth terms}

Smooth terms can be included following GAMs (generalized additive models)
using \code{+ s(x)}, which implements a smooth from \code{\link[mgcv:s]{mgcv::s()}}. \pkg{sdmTMB} uses
penalized smooths, constructed via \code{\link[mgcv:smooth2random]{mgcv::smooth2random()}}. This is a similar
approach implemented in \pkg{gamm4} and \pkg{brms}, among other packages.
Within these smooths, the same syntax commonly used in \code{\link[mgcv:s]{mgcv::s()}} can be
applied, e.g. 2-dimensional smooths may be constructed with \code{+ s(x, y)};
smooths can be specific to various factor levels, \code{+ s(x, by = group)}; the
basis function dimensions may be specified, e.g. \code{+ s(x, k = 4)}; and various
types of splines may be constructed such as cyclic splines to model
seasonality, \code{+ s(month, bs = "cc", k = 12)}. Prior to version 0.0.18,
\pkg{sdmTMB} implemented unpenalized splines.

\strong{Threshold models}

A linear break-point relationship for a covariate can be included via
\code{+ breakpt(variable)} in the formula, where \code{variable} is a single covariate
corresponding to a column in \code{data}. In this case, the relationship is linear
up to a point and then constant (hockey-stick shaped).

Similarly, a logistic-function threshold model can be included via
\code{+ logistic(variable)}. This option models the relationship as a logistic
function of the 50\% and 95\% values. This is similar to length- or size-based
selectivity in fisheries, and is parameterized by the points at which f(x) =
0.5 or 0.95. See the vignette.

Note that only a single threshold covariate can be included.

See the
\href{https://pbs-assess.github.io/sdmTMB/articles/threshold-models.html}{threshold vignette}.

\strong{Extra time: forecasting or interpolating}

Extra time slices (e.g., years) can be included for interpolation or
forecasting with the predict function via the \code{extra_time} argument. The
predict function requires all time slices to be defined when fitting the
model to ensure the various time indices are set up correctly. Be careful if
including extra time slices that the model remains identifiable. For example,
including \code{+ as.factor(year)} in \code{formula} will render a model with no data
to inform the expected value in a missing year. \code{\link[=sdmTMB]{sdmTMB()}} makes no attempt
to determine if the model makes sense for forecasting or interpolation. The
options \code{time_varying}, \code{spatiotemporal = "RW"}, and \code{spatiotemporal = "AR1"}
provide mechanisms to predict over missing time slices with process error.

\code{extra_time} can also be used to fill in missing time steps for the purposes
of a random walk or AR1 spatiotemporal field if their inclusion makes the gaps
between time steps even.

\strong{Index standardization}

For index standardization, you may wish to include \code{0 + as.factor(year)}
(or whatever the time column is called) in the formula. See a basic
example of index standardization in the relevant
\href{https://pbs-assess.github.io/sdmTMB/articles/model-description.html}{package vignette}.
You will need to specify the \code{time} argument. See \code{\link[=get_index]{get_index()}} and/or
\code{\link[=get_index_sims]{get_index_sims()}}.

\strong{Regularization and priors}

You can achieve regularization via penalties (priors) on the fixed effect
parameters. See \code{\link[=sdmTMBpriors]{sdmTMBpriors()}}. These should not include \code{offset} terms.
You can fit the model once without penalties and inspect the element
\code{head(your_model$tmb_data$X_ij)} if you want to see how the formula is
translated to the fixed effect model matrix.
}
\examples{
if (inla_installed()) {

# A coarse mesh for example speed:
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 25)
plot(mesh)

# Tweedie:
fit <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = pcod_2011, time = "year", mesh = mesh, family = tweedie(link = "log"))
print(fit)
tidy(fit, conf.int = TRUE)
tidy(fit, effects = "ran_par", conf.int = TRUE)

# Bernoulli:
pcod_binom <- pcod_2011
m_bin <- sdmTMB(present ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  data = pcod_binom, time = "year", mesh = mesh,
  family = binomial(link = "logit"))
print(m_bin)

# Fit a spatial-only model (by not specifying `time` or setting
# `spatiotemporal = 'off'`):
fit <- sdmTMB(
  density ~ depth_scaled + depth_scaled2, data = pcod_2011,
  mesh = mesh, family = tweedie(link = "log"))
print(fit)

# Gaussian:
pcod_gaus <- subset(pcod_2011, density > 0 & year >= 2013)
mesh_gaus <- make_mesh(pcod_gaus, c("X", "Y"), cutoff = 20)
m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  data = pcod_gaus, mesh = mesh_gaus)
print(m_pos)

# With penalized smoothers via mgcv:
m_gam <- sdmTMB(log(density) ~ s(depth),
  data = pcod_gaus, mesh = mesh_gaus
)

# Specifying the basis dimension via `k =`:
m_gam <- sdmTMB(log(density) ~ s(depth, k = 5),
  data = pcod_gaus, mesh = mesh_gaus
)

# Separate smoother by year:
m_gam <- sdmTMB(log(density) ~ s(depth, by = as.factor(year)),
  data = pcod_gaus, mesh = mesh_gaus
)

# Turning off all random fields (creating a regular GLM/GAM):
m_gam <- sdmTMB(log(density) ~ s(depth),
  data = pcod_gaus, mesh = mesh_gaus,
  spatial = "off", spatiotemporal = "off"
)

# With IID random intercepts:
# Simulate some data:
set.seed(1)
x <- runif(500, -1, 1)
y <- runif(500, -1, 1)
loc <- data.frame(x = x, y = y)
mesh_sim <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")
s <- sdmTMB_sim(x = x, y = y, betas = 0, time = 1L,
  phi = 0.1, range = 1.4, sigma_O = 0.2, sigma_E = 0, mesh = mesh_sim)
s$g <- gl(50, 10)
iid_re_vals <- rnorm(50, 0, 0.3)
s$observed <- s$observed + iid_re_vals[s$g]

# Fit it:
fit <- sdmTMB(observed ~ 1 + (1 | g), mesh = mesh_sim, data = s)
print(fit)
tidy(fit, "ran_pars", conf.int = TRUE) # see tau_G
theta <- as.list(fit$sd_report, "Estimate")
plot(iid_re_vals, theta$RE)

\donttest{

# Spatially varying coefficient (year) example:
d <- pcod_2011
d$year_scaled <- as.numeric(scale(d$year))
fit <- sdmTMB(density ~ depth_scaled, data = d,
  mesh = mesh, family = tweedie(link = "log"),
  spatial_varying = ~ 0 + year_scaled, time = "year")
tidy(fit, effects = "ran_par")

# Time-varying effects of depth and depth squared:
fit <- sdmTMB(density ~ 0 + as.factor(year),
  time_varying = ~ 0 + depth_scaled + depth_scaled2,
  data = pcod_2011, time = "year", mesh = mesh,
  family = tweedie(link = "log"))
print(fit)

# See the b_rw_t estimates; these are the time-varying (random walk) effects.
# These could be added to tidy.sdmTMB() eventually.
summary(fit$sd_report)[1:19,]

# Linear breakpoint model on depth:
m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) +
    breakpt(depth_scaled) + depth_scaled2, data = pcod_gaus,
  time = "year", mesh = mesh_gaus)
print(m_pos)
}
}
}
\references{
Main reference/report introducing the package. We plan to write a paper
to cite in the near future:

Anderson, S.C., E.A. Keppel, A.M. Edwards, 2019. A reproducible data synopsis
for over 100 species of British Columbia groundfish. DFO Can. Sci. Advis. Sec.
Res. Doc. 2019/041. vii + 321 p.
\url{https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2019/2019_041-eng.html}

Reference for local trends:

Barnett, L.A.K., E.J. Ward, S.C. Anderson. Improving estimates of species
distribution change by incorporating local trends. Ecography. 44(3):427-439.
\doi{10.1111/ecog.05176}.

Further explanation of the model and application to calculating climate
velocities:

English, P., E.J. Ward, C.N. Rooper, R.E. Forrest, L.A. Rogers, K.L. Hunter,
A.M. Edwards, B.M. Connors, S.C. Anderson. 2021. Contrasting climate velocity
impacts in warm and cool locations show that effects of marine warming are
worse in already warmer temperate waters. In press at Fish and Fisheries.
\doi{10.1111/faf.12613}.

Code for implementing the barrier-SPDE written by Olav Nikolai Breivik and
Hans Skaug.

A number of sections of the original TMB model code were adapted from the
VAST R package:

Thorson, J.T., 2019. Guidance for decisions using the Vector Autoregressive
Spatio-Temporal (VAST) package in stock, ecosystem, habitat and climate
assessments. Fish. Res. 210:143–161.
\doi{10.1016/j.fishres.2018.10.013}.

Code for the \code{family} R-to-TMB implementation, selected parameterizations of
the observation distributions, general package structure inspiration, and the
idea behind the TMB prediction approach were adapted from the glmmTMB R
package:

Mollie E. Brooks, Kasper Kristensen, Koen J. van Benthem, Arni Magnusson,
Casper W. Berg, Anders Nielsen, Hans J. Skaug, Martin Maechler and Benjamin
M. Bolker (2017). glmmTMB Balances Speed and Flexibility Among Packages for
Zero-inflated Generalized Linear Mixed Modeling. The R Journal, 9(2):378-400.
\doi{10.32614/rj-2017-066}.
}
