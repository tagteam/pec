##' Retrospective table of risks predicted by two different methods, models, algorithms
##'
##' All risks are multiplied by 100 before 
##' @title Retrospective risk reclassification table
##' @param list A list with two elements. Each element should either
##' be a vector with probabilities, or an object for which
##' \code{predictStatusProb} can extract predicted risk based on data.
##' @param formula
##' @param data Used to extract the response from the data and passed
##' on to \code{predictEventProb} to extract predicted event
##' probabilities.
##' @param time Time interest for prediction.
##' @param cause For competing risk models the cause of
##' interest. Defaults to all available causes.
##' @param cuts Risk quantiles to group risks.
##' @param #'formula A survival formula as obtained either #' with
##' \code{prodlim::Hist} or \code{survival::Surv} which defines the
##' response #' in the \code{data}.
##' @param digits Number of digits to show for the predicted risks
##' @return reclassification tables: overall table and one conditional table for each cause and for subjects event free at time interest.
##' @seealso predictStatusProb
##' @examples
##' library(survival)
#' set.seed(40)
#' d <- prodlim::SimSurv(400)
#' nd <- prodlim::SimSurv(400)
#' Models <- list("Cox.X2"=coxph(Surv(time,status)~X2,data=d),
#'                "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=d))
#' rc <- reclass(Models,Surv(time,status)~1,data=nd,time=5)
#' print(rc)
#' plot(rc)
#'
#' set.seed(40)
#' library(riskRegression)
#' library(prodlim)
#' dcr <- prodlim::SimCompRisk(400)
#' ndcr <- prodlim::SimCompRisk(400)
#' crPred5 <- list("X2"=predictEventProb(CSC(Hist(time,event)~X2,data=dcr),newdata=ndcr,times=5),
#'                 "X1+X2"=predictEventProb(CSC(Hist(time,event)~X1+X2,data=dcr),newdata=ndcr,times=5))
#' rc <- reclass(crPred5,Hist(time,event)~1,data=ndcr,time=3)
#' print(rc)
#' 
#' reclass(crPred5,Hist(time,event)~1,data=ndcr,time=5,cuts=100*c(0,0.05,0.1,0.2,1))
#'
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
reclass <- function(list,formula,data,time,cause,cuts=seq(0,100,25)){
    stopifnot(length(list)==2)
    NC <- length(cuts)
    NR <- NC-1 ## dimension of reclassification tables is NR x NR
    # {{{ response
    ## histformula <- formula
    ## if (histformula[[2]][[1]]==as.name("Surv")){
    ## histformula <- update(histformula,paste("prodlim::Hist","~."))
    ## histformula[[2]][[1]] <- as.name("prodlim::Hist")
    ## }
    ## print(histformula)
    ## m <- stats::model.frame(histformula,data,na.action=na.fail)
    m <- stats::model.frame(formula,data,na.action=na.omit)
    response <- model.response(m)
    if (match("Surv",class(response),nomatch=0)!=0){
        attr(response,"model") <- "survival"
        attr(response,"cens.type") <- "rightCensored"
        model.type <- "survival"
    }
    model.type <- attr(response,"model")
    if (model.type=="competing.risks"){
        predictHandlerFun <- "predictEventProb"
        availableCauses <- attr(response,"states")
        if (missing(cause))
            cause <- availableCauses
    }
    else{
        predictHandlerFun <- "predictSurvProb"
    }
    # }}}
    ## for competing risks find the cause of interest.
    if (predictHandlerFun=="predictEventProb"){
        ncauses <- length(availableCauses)
        if (!all(thecauses <- match(cause, availableCauses,nomatch=FALSE)))
            stop("Cause ",cause," is not among the available causes: ",paste(availableCauses,collapse=", "))
    }
    predrisk <- lapply(list,function(x){
                           if (class(x)[[1]]=="numeric")
                               x
                           else{
                               if (predictHandlerFun=="predictEventProb"){
                                   do.call(predictHandlerFun,list(x,newdata=data,times=time,cause=cause))
                               } else{
                                     do.call(predictHandlerFun,list(x,newdata=data,times=time))
                                 }
                           }
                       })
    ## overall reclassification table
    predriskCut <- lapply(predrisk,function(p){
                              stopifnot(min(p)>=min(cuts))
                              stopifnot(max(p)<=max(cuts))
                              cut(100*p,cuts,
                                  include.lowest=TRUE,
                                  labels=paste(paste(cuts[-NC],cuts[-1],sep="-"),"%",sep=""))
                          })
    retab <- table(predriskCut[[1]],predriskCut[[2]])
    ## reclassification frequencies conditional on outcome
    edat <- data.frame(cbind(do.call("cbind",predriskCut),response))
    edat$event[edat$status==0] <- 0
    N <- NROW(edat)
    names(edat)[1:2] <- c("P1","P2")
    cells <- split(edat,list(factor(edat$P1,levels=1:NR),factor(edat$P2,levels=1:NR)))
    all.comb <- apply(expand.grid(1:(NR),1:(NR)),1,paste,collapse=".")
    nn <- names(list)
    if (!is.null(nn) & length(nn)==2){
        names(dimnames(retab)) <- nn
    }
    ## --------------------------------------------------------------------------------------
    ## Apply Bayes' theorem to calculate expected reclassification probabilities
    ##                      conditional on outcome
    ## --------------------------------------------------------------------------------------
    if (predictHandlerFun=="predictEventProb"){
        ## --------------------------------------------------------------------------------------
        ## Competing risk
        ##
        ## P(X=x|T<=t, cause=j) = P(X=x,T<=t,cause=j)        / P(T<=t,cause=j)
        ##                      = P(T<=t,cause=j|X=x) P(X=x) / P(T<=t,cause=j)
        ##                      = cuminc.x H.x / cuminc
        ## 
        ##           P(X=x|T>t) = P(X=x,T>t) / P(T>t)
        ##                      = P(T>t|X=x) P(X=x) / P(T>t)
        ##                      = efreesurv.x H.x / efreesurv
        ## --------------------------------------------------------------------------------------
        eformula <- Hist(time,event)~1
        Hx <- unlist(lapply(cells,NROW))/N
        cuminc.x <- do.call("rbind",
                            lapply(names(cells),
                                   function(cc){
                                       x <- cells[[cc]]
                                       if (NROW(x)>0){
                                           ## warn if too short followup
                                           if (all(x$time<time)) {
                                               warning(call.=FALSE,paste0("pec::reclass: Cell row ",
                                                           sub("\\."," column ",cc),
                                                           " no subject was followed until time ",
                                                           time,
                                                           ". Result is NA (not available)."))
                                               rep(NA,length(thecauses))
                                           } else{
                                                 fit.x <- prodlim::prodlim(eformula,data=x)
                                                 fitted.causes <- attr(fit.x$model.response,"states")
                                                 nstates <- length(fitted.causes)
                                                 sapply(thecauses,function(j){
                                                            ## it may happen that cause j
                                                            ## does not occur in this cell
                                                            if (sum(x$event==j)>0){
                                                                ## check if there is more than one cause
                                                                if (nstates<length(availableCauses)){
                                                                    if (nstates==1){
                                                                        ## only one cause
                                                                        stats::predict(fit.x,times=time,type="cuminc")
                                                                    } else{
                                                                          ## competing causes but less than all causes
                                                                          ## need to change the value of cause
                                                                          xj.cause <- match(j,fitted.causes,nomatch=0)
                                                                          if (xj.cause==0)
                                                                              stop(paste0("Cause ",j,"does not appear in fit. Fitted are causes: ",fitted.causes))
                                                                          else{
                                                                              stats::predict(fit.x,times=time,cause=xj.cause,type="cuminc")
                                                                          }
                                                                      }
                                                                }else{
                                                                     stats::predict(fit.x,times=time,cause=j,type="cuminc")
                                                                 }
                                                            } else {
                                                                  ## warn if no event of type j
                                                                  jstring <- j
                                                                  if (as.character(j)%in%availableCauses[[j]])
                                                                      jstring <- paste0(j," (",availableCauses[[j]],")")
                                                                  warning(call.=FALSE,paste0("pec::reclass: Cell row ",
                                                                              sub("\\."," column ",cc),
                                                                              " no event of type ",
                                                                              jstring,". Result is 0."))
                                                                  return(0)
                                                              }
                                                        })
                                             }
                                       } else{
                                             ## empty cell
                                             rep(0,length(thecauses))
                                         }}))
        fit <- prodlim::prodlim(eformula,data=edat)
        cuminc <- unlist(lapply(thecauses,function(j){stats::predict(fit,times=time,cause=j,type="cuminc")}))
        Px <- apply(cuminc.x * Hx,1,function(p){p/ cuminc})
        ## rownames(Px) <- paste("cause",thecauses,sep=":")
        efreesurv <- 1-sum(cuminc)
        efreesurv.x <- 1-rowSums(cuminc.x,na.rm=TRUE)
        Px <- rbind(Px,"eventfree"=efreesurv.x * Hx / efreesurv)
        event.retab <- lapply(1:NROW(Px),function(xx){
                                        matrix(Px[xx,],ncol=NR)
                                        matrix(Px[xx,],ncol=NR,dimnames=dimnames(retab))
                                    })
        names(event.retab) <- c(paste("Event:",availableCauses),"Event-free")
    } else{
          ## --------------------------------------------------------------------------------------
          ## Survival
          ##
          ## P(X=x|T<=t) = P(X=x,T<=t) /P(T<=t)
          ##             = P(T<=t|X=x) P(X=x) /P(T<=t)
          ##             = cuminc.x * Hx / cuminc
          ##
          ## P(X=x|T>t)  = P(X=x,T>t) /P(T>t)
          ##             = surv.x * Hx / surv
          ## 
          ## --------------------------------------------------------------------------------------
          eformula <- Hist(time,status)~1
          Hx <- unlist(lapply(cells,NROW))/N
          cuminc <- stats::predict(prodlim::prodlim(eformula,data=edat),times=time,type="cuminc")
          cuminc.x <- sapply(cells,function(x){
                                 if (NROW(x)>0){
                                     ## warn if too short followup
                                     if (all(x$time<time)) {
                                         warning(call.=FALSE,paste0("pec::reclass: Cell row ",
                                                     sub("\\."," column ",cc),
                                                     " no subject was followed until time ",
                                                     time,
                                                     ". Result is NA (not available)."))
                                         NA
                                     }else{
                                          stats::predict(prodlim::prodlim(eformula,data=x),times=time)
                                      }
                                 } else {
                                       ## empty cell
                                       0
                                   }
                             })
          surv <- 1-cuminc
          surv.x <- 1-cuminc.x
          event.retab <- list("event"=matrix(cuminc.x*Hx/surv,ncol=NR,dimnames=dimnames(retab)),
                              "eventfree"=matrix(surv.x * Hx / surv,ncol=NR,dimnames=dimnames(retab)))
      }
    out <- list(time=time,
                predRisk=predrisk,
                reclassification=retab,
                event.reclassification=event.retab,
                cuts=cuts,
                model=model.type)
    class(out) <- "riskReclassification"
    out
}

##' @S3method print riskReclassification
print.riskReclassification <- function(x,percent=TRUE,digits=ifelse(percent,1,2),...){
    cat("Observed overall re-classification table:\n\n")
    print(x$reclassification)
    cat("\nExpected re-classification probabilities (%) among subjects with event until time ",x$time,"\n\n",sep="")
    fmt <- paste0("%1.", digits[[1]], "f")
    dnames <- dimnames(x$reclassification)
    dim <- dim(x$reclassification)
    if (percent==TRUE){
        rlist <- lapply(x$event.reclassification,function(x){
                            matrix(sprintf(fmt=fmt,100*c(x)),nrow=dim[1],ncol=dim[2],dimnames=dnames)
                        })
    }else{
         rlist <- lapply(x$event.reclassification,function(x){
                             matrix(sprintf(fmt=fmt,c(x)),nrow=dim[1],ncol=dim[2],dimnames=dnames)
                         })
     }
    if (x$model=="competing.risks"){
        print.listof(rlist[-length(rlist)],quote=FALSE)
    } else{
          print.listof(rlist[1],quote=FALSE)
      }
    cat("\nExpected re-classification probabilities (%) among subjects event-free until time ",x$time,"\n\n",sep="")
    print.listof(rlist[length(rlist)],quote=FALSE)
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
