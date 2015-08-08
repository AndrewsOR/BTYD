source("setup.R")

#Parameter estimation
params <- pnbd.EstimateParameters(cal.cbs)
params

LL <- pnbd.cbs.LL(params, cal.cbs);
LL

p.matrix <- c(params, LL);
for (i in 1:2){
  params <- pnbd.EstimateParameters(cal.cbs, params);
  LL <- pnbd.cbs.LL(params, cal.cbs);
  p.matrix.row <- c(params, LL);
  p.matrix <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("r", "alpha", "s", "beta", "LL");
rownames(p.matrix) <- 1:3;
p.matrix;


pnbd.PlotTransactionRateHeterogeneity(params)
dev.copy(png,'TransactionRateHeterogeneity.png')
dev.off()
pnbd.PlotDropoutRateHeterogeneity(params)
dev.copy(png,'DropoutRateHeterogeneity.png')
dev.off()

# Individual level estimations

pnbd.Expectation(params, t=52)

# Expected behavior from a particular customer
custName <- sample(cal.cbs[,1],1)
custName

cal.cbs[custName,]

x<-cal.cbs[custName,"x"]
t.x <-cal.cbs[custName,"t.x"]
T.cal <- cal.cbs[custName,"T.cal"]

pnbd.ConditionalExpectedTransactions(params, T.star = 52, x, t.x, T.cal)
pnbd.PAlive(params, x, t.x, T.cal)

# To visualize the distribution of P(Alive) across customers:
p.alives <- pnbd.PAlive(params, cal.cbs[,"x"], cal.cbs[,"t.x"], cal.cbs[,"T.cal"])

ggplot(as.data.frame(p.alives),aes(x=p.alives))+
  geom_histogram(colour="grey",fill="orange")+
  ylab("Number of Customers")+
  xlab("Probability Customer is 'Live'")+
  theme_minimal()
dev.copy(png,'PAlive.png')
dev.off()

# Goodness of Fit
censor <- 7
pnbd.PlotFrequencyInCalibration(params, cal.cbs, censor)
dev.copy(png,'FreqInCalibration.png')
dev.off()

# Verify in Holdout Period
x.star <- hold.cbs[,"x.star"]
comp <- pnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star=52, cal.cbs, x.star, censor)
rownames(comp) <- c("act", "exp", "bin")
comp
dev.copy(png,'FreqVsCondExpFreq.png')
dev.off()

pnbd.PlotRecVsConditionalExpectedFrequency(params, cal.cbs, T.star=52, x.star)
dev.copy(png,'RecVsCondExpFreq.png')
dev.off()

# Plot Actual V/s Expected Transactions on a weekly basis
tot.cbt <- dc.CreateFreqCBT(elog)
head(tot.cbt)

# ...Completed Freq CBT
d.track.data <- rep(0, 7 * 105)
origin <- as.Date("2013-01-01")
for (i in colnames(tot.cbt)){
  date.index <- difftime(as.Date(i), origin) + 1;
  d.track.data[date.index] <- sum(tot.cbt[,i]);
}
w.track.data <- rep(0, 105)
for (j in 1:105){
  w.track.data[j] <- sum(d.track.data[(j*7-6):(j*7)])
}

T.cal <- cal.cbs[,"T.cal"]
T.tot <- 105
n.periods.final <- 105
inc.tracking <- pnbd.PlotTrackingInc(params, T.cal,
                                     T.tot, w.track.data,
                                     n.periods.final)
inc.tracking[,20:25]
dev.copy(png,'TrackingInc.png')
dev.off()

cum.tracking.data <- cumsum(w.track.data)
cum.tracking <- pnbd.PlotTrackingCum(params, T.cal,
                                     T.tot, cum.tracking.data,
                                     n.periods.final)
cum.tracking[,20:25]
dev.copy(png,'TrackingCum.png')
dev.off()

##DERT plot (Discounted Expected Residual Plots)

d <- 0.0027  ### 15\% compounded annually has been converted to 0.0027 compounded continously,
# as we are dealing with weekly data and not annual data

pnbd.Plot.DERT(params, x=0:14, t.x=0:77, T.cal=77.86, d, type="persp")
dev.copy(png,'DERT-persp.png')
dev.off()

pnbd.Plot.DERT(params, x=0:14, t.x=0:77, T.cal=77.86, d, type="contour")
dev.copy(png,'DERT-iso.png')
dev.off()

ave.spend <- as.vector(simData$cust.data$m.x)
ave.spend[ave.spend<0]<-0
tot.trans <- cal.cbs[,"x"]

spend.params <-spend.EstimateParameters(ave.spend, tot.trans)
spend.params
# What is the expected value of a customer with average transactions of $25 and 2 tx during calibration period?
spend.expected.value(spend.params, m.x=25, x=1)

spend.plot.average.transaction.value(spend.params, ave.spend, tot.trans)
dev.copy(png,'SpendPlot.png')
dev.off()
                                     