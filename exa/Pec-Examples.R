# {{{ testing seed
library(pec)
set.seed(17)
d <- SimSurv(100)
f <- coxph(Surv(time,status)~X2,data=d)
set.seed(13)
a=pec(f,splitMethod="bootcv",B=3,M=63,keep.index=TRUE,verbose=F)
b=a$splitMethod$index
set.seed(13)
c=resolvesplitMethod(splitMethod="bootcv",N=100,M=63,B=3)$index
stopifnot(all.equal(b,c))
# }}}

## library(party)
## h <- cforest(Surv(time,status)~X2,data=d)
## isS4(h)
## f <- cph(Surv(time,status)~X2,data=d,surv=TRUE)
## set.seed(19)
## A <- pec(f,splitMethod="cv5",B=1,M=63,keep.index=TRUE,verbose=F)

# {{{ testing ipcw
set.seed(18)
A=pec(f,B=30,splitMethod="bootcv")
set.seed(18)
A1=pec(f,B=30,ipcw.refit=T,splitMethod="bootcv")
cbind(A$BootCvErr$CoxModel,A1$BootCvErr$CoxModel)
## graphics::plot(A,xlim=c(0,100))
## graphics::plot(A1,add=TRUE,lty=3)
B=pec(f,cens.model="cox")

# }}}

# {{{ testing splitMethods: cvk, loocv

# Comparing the predictive performance of some standard
# survival regression models 
# --------------------------------------------------------------------

library(pec)
library(rms) # thanks to Frank Harrell 
data(GBSG2)
GBSG2$status <- GBSG2$cens
GBSG2 <- GBSG2[order(GBSG2$time,-GBSG2$status),]
GBSG2$grade.bin <- as.factor(as.numeric(GBSG2$tgrade!="I"))
levels(GBSG2$grade.bin) <- c("I","II/III")

GBSG2=GBSG2[sample(1:NROW(GBSG2),size=50),]
# Kaplan-Meier (survival package)
system.time(pbc.fit0 <- survfit(Surv(time,status)~1,data=GBSG2))
# Kaplan-Meier (prodlim package) is faster
system.time(pbc.fit0 <- prodlim(Hist(time,status)~1,data=pbc) )
## Cox model (rms)
Cox=cph(Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2,x=TRUE,y=TRUE,surv=TRUE,se.fit=FALSE)

set.seed(17)
check.code("pec")
loocv <- pec.list(object=list(Cox),formula=Surv(time,status)~1,data=GBSG2,cens.model="marginal",splitMethod="loocv",B=1,keep.matrix=TRUE,verbose=TRUE)
cv5 <- pec.list(object=list(Cox),formula=Surv(time,status)~1,data=GBSG2,cens.model="marginal",splitMethod="cv5",B=2,keep.matrix=TRUE,verbose=TRUE)
# }}}

# {{{ testing splitMethods: boot632, boot632+

## d <- SimSurv(300)
## f2 <- coxph(Surv(time,status)~rcs(X1)+X2,data=d)

library(pec)
library(rms) # thanks to Frank Harrell 
data(GBSG2)
GBSG2$status <- GBSG2$cens
GBSG2 <- GBSG2[order(GBSG2$time,-GBSG2$status),]
GBSG2$grade.bin <- as.factor(as.numeric(GBSG2$tgrade!="I"))
levels(GBSG2$grade.bin) <- c("I","II/III")

# {{{ Parallel computing
library(pec)
library(doParallel)
library(prodlim) 
library(survival)
set.seed(130971) 
dat <- SimSurv(100) 
Models <-list("Cox.X1"=coxph(Surv(time,status)~X1,data=dat,y=TRUE))
## ,"Cox.X2"=coxph(Surv(time,status)~X2,data=dat,y=TRUE),"Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=dat,y=TRUE))
cl <- makeCluster(2) 
registerDoParallel(cl)
PredError.632plus <- pec(object=Models,formula=Surv(time,status)~X1+X2,data=dat,exact=TRUE,cens.model="marginal",splitMethod="bootcv",B=2,verbose=TRUE)
PredError.632plus <- pec(object=Models,formula=Surv(time,status)~X1+X2,data=dat,exact=TRUE,cens.model="marginal",splitMethod="Boot632plus",B=2,verbose=TRUE,cause=1)
stopCluster(cl) #####################################

# }}}
         
