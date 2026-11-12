# Reproduces the vignette's weights-bias simulation and adds the additive
# weight 1/(tau^2 + v). Truth: group means 0.6 and 2.0 (log10), biological SD
# tau = 0.4, measurement SD = se_fun(true value), U-shaped across the curve.
set.seed(7)
n   <- 20000
mu  <- c(low = 0.6, high = 2.0)
tau <- 0.4
grp <- rep(names(mu), each = n)
tru <- rnorm(2 * n, mean = rep(mu, each = n), sd = tau)
se_fun <- function(y) 0.05 + 0.35 * plogis((0.5 - y) / 0.15) +
                              0.35 * plogis((y - 3.0) / 0.15)
obs <- tru + rnorm(2 * n, 0, se_fun(tru))
v   <- se_fun(obs)^2                       # what the assay reports for each sample

wmean <- function(w) tapply(obs * w, grp, sum)[names(mu)] / tapply(w, grp, sum)[names(mu)]
w_mult <- 1 / v
w_add  <- 1 / (tau^2 + v)
res <- rbind(truth = mu,
             unweighted = tapply(obs, grp, mean)[names(mu)],
             `1/v (multiplicative logic)` = wmean(w_mult),
             `1/(tau2 + v) (additive)`    = wmean(w_add))
print(round(res, 3))
cat("weight max/min: 1/v =", round(max(w_mult) / min(w_mult), 1),
    "  additive =", round(max(w_add) / min(w_add), 1), "\n")
