normal.cox <- coxph(Surv(time,status)~edema,data=pbc)

update.cox <- function(object,tstar,data){
  object$call$data <- data[data$time>tstar,]
  update <- eval(object$call)
  class(update) <- "dynamicCox"
  update
}

predictProb.dynamicCox <- function(object,newdata,cutpoints,learn.data,...){
  p <- matrix(1,nrow=NROW(newdata),ncol=length(cutpoints))
  p
}


pec.c <- pec(object=list(c),
             formula=Surv(time,status)~1,
             data=pbc,
             exact=TRUE,
             method="ipcw",
             cens.model="marg",
             B=0,
             verbose=TRUE)

