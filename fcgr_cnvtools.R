# using R 2.15.0 in VM Fedora 15 sandbox environment

# read in FCGR2A:FCGR2C CNP z-scores - standardised to mean 0 standard deviation 1.  
# Post-QC.

# define matrix and assign sample ID's based on original data table.
library(CNVtools)
library(ggplot2)
library(foreign)
zdata <- read.table("/Users/michaelm/Desktop/FCGR/CNVtools/BRAGGSS eular_resp zcnp.txt",
                    h=T, stringsAsFactors=F)
z.signal <- zdata$signal

pdf("/Users/michaelm/Desktop/FCGR/znorm_cnp2a.pdf")
hist(z.signal, breaks=50, main="Z-score 2A CNP", cex.lab=1.3)
dev.off()

#show distribution of data
# run model selection for number of components - optimal should be 3 based on CN<2, CN=2 and CN>2.
# set batches based on cohort (GORA control, BRAGGSS1 and BRAGGSS4).
# set sample IDs based on original file/data table
# run model selection with variance constrained to be proportional to copy number with 3 iterations

batches <- factor(zdata$batch)
sample <- factor(zdata$sample_id)
set.seed(0)
z.results <- CNVtest.select.model(signal=z.signal, 
                                  batch=batches, sample=sample, 
                                  n.H0=3, method="BIC", v.ncomp=1:5, 
                                  v.model.component=rep('gaussian',5), 
                                  v.model.mean=rep("~ strata(cn)",5), 
                                  v.model.var=rep("~strata(cn)",5))
ncomp <- 3
pdf("z.model_select.pdf")
plot(-z.results$BIC, xlab="n comp", ylab="-BIC", type="b", lty=2, col="red", pch='+')
dev.off()

# select 3 component model based on lowest BIC
# run model under null hypothesis to check EM algorithm is converging properly, 
# run with 6 iterations.
z.fit <- CNVtest.binary(signal=z.signal,
                     sample=sample, batch=batches, 
                     ncomp=ncomp, n.H0=6, n.H1=0, 
                     model.var='~ strata(cn,batch)')

print(z.fit$status.H0)

# EM converges, have a look at posterior probability distributions, 
# check nothing funny is going on
pdf("z.fit_posteriors.pdf")
cnv.plot(z.fit$posterior.H0, batch="BRAGGSS1", main="BRAGGSS1", breaks=50, col="red")
cnv.plot(z.fit$posterior.H0, batch="BRAGGSS4", main="BRAGGSS4", breaks=50, col="red")
dev.off()

# all cohorts have posteriors fit representing 3 copy number states
# 2 d.f.test of association with no trend assumed, iterations under H0=6, 
# iterations under H1=3
# assign disease status as GORA=0, BRAGGSS1 and BRAGGSS4=1
trait <- zdata$eular_resp
show(trait)
z.test <- CNVtest.binary(signal=z.signal, sample=sample, 
                      batch=batches, disease.status=trait, ncomp=3, 
                      n.H0=6, n.H1=3, model.disease="~as.factor(cn)", 
                      model.var="~ strata(cn)")
print(z.test$status.H0)
print(z.test$status.H1)

# EM converges under both H0 and H1
# print model statistics

print(z.test$model.H0)
print(z.test$model.H1)

# calculate OR for each comparison against WT (CN=2)
p1 <- (z.test$model.H1$pdc[1,1])
p2 <- (z.test$model.H1$pdc[2,1])
p3 <- (z.test$model.H1$pdc[3,1])
OR1 <- (p1*(1-p2))/(p2*(1-p1))
print(OR1)
OR2 <- (p1*(1-p3))/(p3*(1-p1))
print(OR2)

# write out data to .csv file
write.table(z.test$posterior.H1, file="/Users/michaelm/Desktop/FCGR/posterior_z.csv", 
            sep=",", row.names=FALSE)

# test association with LR test
LR.statistic <- -2*(z.test$model.H0$lnL - z.test$model.H1$lnL)
print(LR.statistic)
1 - pchisq(LR.statistic, df=2)


# fit models of WT vs del and WT vs dup separately.
z.deltest <- CNVtest.binary(signal=z.signal, sample=sample, 
                            batch=batches, disease.status=trait, 
                            ncomp=3, n.H0=6, n.H1=3, model.disease="~I(cn%/%2)", 
                            model.var="~ strata(cn)")
print(z.deltest$status.H0)
print(z.deltest$status.H1)

# model statistics
print(z.deltest$model.H0)
print(z.deltest$model.H1)

p1 <- (z.deltest$model.H1$pdc[1,1])
p2 <- (z.deltest$model.H1$pdc[2,1])
OR <- (p1*(1-p2))/(p2*(1-p1))
print(OR)

LR.statdel <- -2*(z.deltest$model.H0$lnL - z.deltest$model.H1$lnL)
print(LR.statdel)
1-pchisq(LR.statdel, df=1)

# WT vs dup
z.duptest=CNVtest.binary(signal=z.signal, sample=sample, batch=batches, disease.status=trait, ncomp=3, n.H0=6, n.H1=3, model.disease="~I(cn%/%3)", model.var="~ strata(cn)")
print(z.deltest$status.H0)
print(z.deltest$status.H1)

# model statistics
print(z.duptest$model.H0)
print(z.duptest$model.H1)

p1 <- (z.duptest$model.H1$pdc[1,1])
p2 <- (z.duptest$model.H1$pdc[3,1])
OR <- (p1*(1-p2))/(p2*(1-p1))
print(OR)

LR.statdup <- -2*(z.duptest$model.H0$lnL - z.duptest$model.H1$lnL)
print(LR.statdup)
1-pchisq(LR.statdup, df=1)

# Use the posterior probability estimates to test CNV on DAS28 CRP
clin.data <- read.dta("/Users/michaelm/Desktop/FCGR/BRAGGSS/BRAGGSS_fcgr_cnv_clinical.dta")
clin.merge <- merge(z.fit$posterior.H0, clin.data, by.x=c("subject"),
                    by.y=c("sample_id"))

# regress deletion posterior on delta DAS28-CRP, adjust for baseline
del.lm <- glm("CHDAS28 ~ DAS28BL + p1_2c", data=clin.merge,
             family=gaussian(link="identity"))

summary(del.lm)

wt.lm <- glm("CHDAS28 ~ DAS28BL + p2_2c", data=clin.merge,
              family=gaussian(link="identity"))

summary(wt.lm)

dup.lm <- glm("CHDAS28 ~ DAS28BL + p3_2c", data=clin.merge,
              family=gaussian(link="identity"))

summary(dup.lm)