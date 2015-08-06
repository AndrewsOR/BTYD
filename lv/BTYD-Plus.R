library(devtools)
library(gsl)
devtools::install_github("mplatzer/BTYDplus", dependencies=TRUE)
library(BTYDplus)
library(BTYD)
#demo(package="BTYDplus")
#demo("cdnow")
#demo("gg-nbd")
#demo("mbg-cnbd-k")
#demo("nbd")
#demo("pareto-ggg")
#demo("pareto-nbd-abe")
#demo("timing")

cbs<- as.data.frame(cal.cbs)

# NBD
(params.nbd <- nbd.EstimateParameters(cbs))

# Pareto/NBD (from BTYD package)
(params.pnbd <- BTYD::pnbd.EstimateParameters(cbs))

# Gamma/Gompertz/NBD
#(params.ggnbd <- ggnbd.EstimateParameters(cbs, trace=10)) # Would take several minutes
# Dummy values
(params.ggnbd <- c(r=0.552256, alpha=10.568448, b=0.000035, s=0.609494, beta=0.000417))


# BG/NBD (from BTYD package)
(params.bgnbd <- BTYD::bgnbd.EstimateParameters(cbs))

# MBG/NBD
(params.mbgnbd <- mbgnbd.EstimateParameters(cbs))

# MBG/CNBD-k
#(params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs))
# -> MBG/CNBD-k is identical to MBG/NBD, as no regularity is detected, hence k=1
params.mbgcnbd <- params.mbgnbd

rbind("NBD"=nbd.cbs.LL(params.nbd, cbs),
      "Pareto/NBD"=BTYD::pnbd.cbs.LL(params.pnbd, cbs),
      "GG/NBD"=ggnbd.cbs.LL(params.ggnbd, cbs),
      "BG/NBD"=BTYD::bgnbd.cbs.LL(params.bgnbd, cbs),
      "MBG/NBD"=mbgnbd.cbs.LL(params.mbgnbd, cbs),
      "MBG/CNBD-k"=mbgcnbd.cbs.LL(params.mbgcnbd, cbs))

cbs$nbd <- nbd.ConditionalExpectedTransactions(params.nbd, cbs$T.star, cbs$x, cbs$T.cal)
cbs$pnbd <- BTYD::pnbd.ConditionalExpectedTransactions(params.pnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$ggnbd <- ggnbd.ConditionalExpectedTransactions(params.ggnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$bgnbd <- BTYD::bgnbd.ConditionalExpectedTransactions(params.bgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$mbgnbd <- mbgnbd.ConditionalExpectedTransactions(params.mbgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$mbgcnbd <- mbgcnbd.ConditionalExpectedTransactions(params.mbgcnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)


#Estimate P(alive)

cbs$palive.nbd <- 1
cbs$palive.pnbd <- BTYD::pnbd.PAlive(params=params.pnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.ggnbd <- ggnbd.PAlive(params=params.ggnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.bgnbd <- BTYD::bgnbd.PAlive(params=params.bgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.mbgnbd <- mbgnbd.PAlive(params=params.mbgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.mbgcnbd <- mbgcnbd.PAlive(params=params.mbgcnbd, cbs$x, cbs$t.x, cbs$T.cal)



MAPE <- function(a, f) sum(abs(a-f)/sum(a))
RMSE <- function(a, f) sqrt(mean((a-f)^2))
MSLE <- function(a, f) mean(((log(a+1) - log(f+1)))^2)
BIAS <- function(a, f) sum(f)/sum(a)-1
bench <- function(cbs, models) {
  acc <- t(sapply(models, function(model) c(MAPE(cbs$x.star, cbs[, model]),
                                            RMSE(cbs$x.star, cbs[, model]),
                                            MSLE(cbs$x.star, cbs[, model]),
                                            BIAS(cbs$x.star, cbs[, model]))))
  colnames(acc) <- c("MAPE", "RMSE", "MSLE", "BIAS")
  round(acc, 3)
}

bench(cbs, c("nbd", "pnbd", "ggnbd", "bgnbd", "mbgnbd", "mbgcnbd"))

bench(cbs, c("pnbd", "bgnbd"))

#Estimate Pareto/GGG model via MCMC - this will take several minutes

# Note: For keeping the runtime of this demo short, we limit the MCMC to 500 
# iterations for both chains. In fact, the chains should be run longer to 
# ensure convergence, and collect enough samples

set.seed(1)

pggg.draws <- pggg.mcmc.DrawParameters(cbs, mcmc=500, burnin=100, chains=2, thin=10)

plot(pggg.draws$level_2, density=FALSE)
plot(pggg.draws$level_2, trace=FALSE)

coda::gelman.diag(pggg.draws$level_2)
# -> MCMC chains have not converged yet

(summary(pggg.draws$level_2)$quantiles[, "50%"])

pggg.mcmc.plotRegularityRateHeterogeneity(pggg.draws)
# -> very narrow distribution around k=1; 
# CDNow customers purchases indeed seem to follow Poisson process

round(coda::effectiveSize(pggg.draws$level_2))
# -> effective sample size are small for such a short chain


Estimate Future Transactions & P(active) for MCMC models

# draw future transaction
pnbd.xstar <- mcmc.DrawFutureTransactions(cbs, pnbd.draws, T.star=cbs$T.star)
pggg.xstar <- mcmc.DrawFutureTransactions(cbs, pggg.draws, T.star=cbs$T.star)

# calculate mean over future transaction draws for each customer
cbs$pnbd.mcmc <- apply(pnbd.xstar, 2, mean)
cbs$pggg.mcmc <- apply(pggg.xstar, 2, mean)

# forecasting accuracy
bench(cbs, c("pnbd", "pnbd.mcmc", "pggg.mcmc"))

# calculate P(active)
cbs$pactive.pnbd.mcmc <- apply(pnbd.xstar, 2, function(x) mean(x>0))
cbs$pactive.pggg.mcmc <- apply(pggg.xstar, 2, function(x) mean(x>0))

# calculate P(alive)
cbs$palive.pnbd.mcmc <- mcmc.PAlive(cbs, pnbd.draws)
cbs$palive.pggg.mcmc <- mcmc.PAlive(cbs, pggg.draws)
