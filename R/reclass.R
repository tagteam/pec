##' Tabulate grouped risks predicted by two different methods, models, algorithms
##'
##' All risks are multiplied by 100 before 
##' @title Risk reclassification table
##' @param list A list with two elements. Each element should either
##' be a vector with probabilities, or an object for which \code{predictStatusProb}
##' can extract predicted risk based on newdata.
##' @param newdata Passed on to \code{predictStatusProb}
##' @param cuts Risk quantiles to group risk
##' @param digits Number of digits to show for the predicted risks
##' @return reclassification table 
##' @seealso predictStatusProb
##' @examples
##' library(survival)
#' set.seed(40)
#' d <- prodlim::SimSurv(400)
#' nd <- prodlim::SimSurv(400)
#' Models <- list("Cox.X2"=coxph(Surv(time,status)~X2,data=d),
#'                "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=d))
#' rc <- reclass(Models,Surv(time,status)~1,newdata=nd,time=5)
#' print(rc)
#' plot(rc)
#'
#' set.seed(40)
#' library(riskRegression)
#' library(prodlim)
#' d <- prodlim::SimCompRisk(400)
#' nd <- prodlim::SimCompRisk(400)
#' crModels <- list("Cox.X2"=CSC(Hist(time,event)~X2,data=d),
#'                "Cox.X1.X2"=CSC(Hist(time,event)~X1+X2,data=d))
#' rc <- reclass(crModels,Hist(time,event)~1,newdata=nd,time=5)
#' print(rc)
#'
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
reclass <- function(list,formula,newdata,time,cause,cuts=seq(0,100,25),digits=1){
    stopifnot(length(list)==2)
    # {{{ response
    ## histformula <- formula
    ## if (histformula[[2]][[1]]==as.name("Surv")){
    ## histformula <- update(histformula,paste("prodlim::Hist","~."))
    ## histformula[[2]][[1]] <- as.name("prodlim::Hist")
    ## }
    ## print(histformula)
    ## m <- model.frame(histformula,data,na.action=na.fail)
    m <- model.frame(formula,newdata,na.action=na.omit)
    response <- model.response(m)
    if (match("Surv",class(response),nomatch=0)!=0){
        attr(response,"model") <- "survival"
        attr(response,"cens.type") <- "rightCensored"
        model.type <- "survival"
    }
    model.type <- attr(response,"model")
    if (model.type=="competing.risks"){
        predictHandlerFun <- "predictEventProb"
        if (missing(cause))
            cause <- attr(response,"state")[1]
    }
    else{
        predictHandlerFun <- "predictSurvProb"
    }
    # }}}
    ## for competing risks find the cause of interest.
    if (predictHandlerFun=="predictEventProb"){
        event <- prodlim::getEvent(response,mode="character")
        availableCauses <- unique(event[event!="unknown"])
        if (!match(cause, availableCauses,nomatch=FALSE))
            stop("Cause ",cause," is not among the available causes: ",paste(availableCauses,collapse=", "))
        event <- event==cause
    }
    predrisk <- lapply(list,function(x){
                           if (class(x)[[1]]=="numeric")
                               round(100*x,digits)
                           else{
                               if (predictHandlerFun=="predictEventProb"){
                                   do.call(predictHandlerFun,list(x,newdata=newdata,times=time,cause=cause))
                               } else{
                                     do.call(predictHandlerFun,list(x,newdata=newdata,times=time))
                                 }
                           }
                       })
    ## observed reclassification table
    predriskCut <- lapply(predrisk,function(p)cut(round(100*p,digits),cuts,
                                                  include.lowest=TRUE,
                                                  labels=paste(paste(cuts[-length(cuts)],cuts[-1],sep="-"),"%",sep="")))
    retab <- table(predriskCut[[1]],predriskCut[[2]])
    ## expected reclassification frequencies
    ## edat <- data.frame(cbind(response,do.call("cbind",lapply(predriskCut,as.character))))
    edat <- data.frame(cbind(do.call("cbind",predriskCut),response))
    N <- NROW(edat)
    names(edat)[1:2] <- c("P1","P2")
    cells <- split(edat,list(edat$P1,edat$P2))
    km.time <- predict(prodlim::prodlim(formula,data=edat),times=time,type="cuminc")
    all.comb <- apply(expand.grid(1:(length(cuts)-1),1:(length(cuts)-1)),1,paste,collapse=".")
    nn <- names(list)
    ## Apply Bayes' theorem to calculate expected number of persons with an event in each
    ## cell of the table
    ## P(Pred=x|T<t) = P(Pred=x,T<t) /P(T<t)
    ##               = P(T<t|Pred=x) P(Pred=x) /P(T<t)
    if (predictHandlerFun=="predictEventProb"){
        eventtable <- lapply(availableCauses,function(cause){
                                 expectedevents <- lapply(all.comb,function(j){
                                                              x <- cells[[j]]
                                                              if (NROW(x)>0){predict(prodlim::prodlim(formula,data=x),times=time,cause=cause,type="cuminc") * (NROW(x)/N)/ km.time } else 0
                                                          })
                                 e <- matrix(unlist(expectedevents),ncol=length(cuts)-1)*retab
                                 dimnames(e) <- dimnames(retab)
                                 if (!is.null(nn) & length(nn)==2){
                                     names(dimnames(e)) <- nn
                                 }
                                 e
                             })
        names(eventtable) <- paste("cause",availableCauses,sep=":")
        eventfreetable <- retab - Reduce("+",eventtable)
        dimnames(eventfreetable) <- dimnames(retab)
        if (!is.null(nn) & length(nn)==2){
            names(dimnames(retab)) <- names(dimnames(eventfreetable)) <- nn
        }
    }
    else{
        expectedevents <- lapply(all.comb,function(j){
                                     x <- cells[[j]]
                                     if (NROW(x)>0){predict(prodlim::prodlim(formula,data=x),times=time,type="cuminc") * (NROW(x)/N)/ km.time } else 0
                                 })
        eventtable <- matrix(unlist(expectedevents),ncol=length(cuts)-1)*retab
        eventfreetable <- retab-eventtable
        dimnames(eventtable) <- dimnames(eventfreetable) <- dimnames(retab)
        if (!is.null(nn) & length(nn)==2){
            names(dimnames(retab)) <- names(dimnames(eventtable)) <- names(dimnames(eventfreetable)) 
        }
    }
    out <- list(predRisk=predrisk,
                reclassification.table=retab,
                eventtable=eventtable,
                eventfreetable=eventfreetable,
                cuts=cuts)
    class(out) <- "riskReclassification"
    out
}

##' @S3method print riskReclassification
print.riskReclassification <- function(x,...){
    print(x$table)
    print(x$etable)
}

##' @S3method plot riskReclassification
plot.riskReclassification <- function(x,xlim=c(0,100),ylim=c(0,100),xlab,ylab,grid=TRUE,grid.col=gray(0.9),...){
    if (missing(xlab)) xlab <- paste("Risk (%):",names(dimnames(x$table))[[1]])
    if (missing(ylab)) ylab <- paste("Risk (%):",names(dimnames(x$table))[[2]])
    plot(x$predRisk[[1]],x$predRisk[[2]],axes=FALSE,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
    axis(1,at=x$cuts)
    axis(2,at=x$cuts)
    if (grid==TRUE)
        abline(h = x$cuts, v = x$cuts, col = gray(0.9))
}
