### plot.calibrationPlot.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 28 2015 (17:32) 
## Version: 
## last-updated: Oct  1 2015 (19:45) 
##           By: Thomas Alexander Gerds
##     Update #: 79
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Calibration plots 
##'
##' @title Plot objects obtained with \code{calPlot}
##' @param x Object obtained with \code{calPlot}
##' @param ... Not used. 
##' @return Nothing
##' @seealso \code{calPlot}
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
plot.calibrationPlot <- function(x,...){
    # {{{ plot an empty frame
    plotFrames <- x$plotFrames
    control <- x$control
    NF <- x$NF
    if (x$add==FALSE && !x$bars){
        do.call("plot",control$plot)
    }
    if (x$diag && !x$bars){
        segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
    }
    # }}}
    # {{{ show calibration 
    showBars <- function(){
        pf <- na.omit(plotFrames[[1]])
        Pred <- pf$Pred
        Obs <- pf$Obs
        if(!x$legend){
            control$barplot$legend.text <- NULL
        }else{
             if (is.null(control$barplot$legend.text)){
                 control$barplot$legend.text <- control$legend$legend
             }else{
                  control$barplot$args.legend <- control$legend
              }
         }
        if (is.null(control$barplot$space))
            control$barplot$space <- rep(c(1,0),length(Pred))
        PredObs <- c(rbind(Pred,Obs))
        control$barplot$height <- PredObs
        if (x$hanging)
            control$barplot$offset <- c(rbind(0,Pred-Obs))
        if (length(x$names)>0 && (x$names[[1]]!=FALSE) && is.character(x$names)){
            if (length(x$names)==(length(control$barplot$height)/2)){
                ## colnames(control$barplot$height) <- x$names
                control$barplot$names.arg <- c(x$names,rbind(rep("",10)))
            }else{
                 if (x$names[[1]]==FALSE){
                     control$barplot$names.arg <- NULL
                 }
             }
        }
        ## if (x$legend) print(control$barplot$args.legend)
        coord <- do.call("barplot",control$barplot)
        ## message(paste0("Bars are located at ",paste(coord,collapse=",")))
        if (x$hanging)
            do.call("abline",control$abline)
        if (x$showFrequencies){
            if(x$hanging){
                text(x=coord,
                     y=(as.vector(rbind(Pred,Pred)) + rep(c(0.03,0.05),times=length(as.vector(coord))/2)),
                     paste(round(100*c(rbind(Pred,Obs)),0),"%",sep=""),xpd=NA)
            }else{
                 text(coord,
                      c(rbind(Pred,Obs)+0.03),
                      paste(round(100*c(rbind(Pred,Obs)),0),"%",sep=""),xpd=NA)
             }
        }
    }
    showCal <- function(f){
        if (is.null(x$pseudo.col)){
            ccrgb=as.list(col2rgb(x$col[f],alpha=TRUE))
            names(ccrgb) <- c("red","green","blue","alpha")
            ccrgb$alpha <- x$jack.density
            jack.col <- do.call("rgb",c(ccrgb,list(max=255)))
        }
        else
            jack.col <- x$pseudo.col
        if (is.null(x$pseudo.pch)) x$pseudo.pch <- 1
        if (x$showPseudo) {
            points(x$pred[,f+1],x$pred[,1],col=jack.col,pch=x$pseudo.pch)
        }
        pf <- x$plotFrames[[f]]
        if(NROW(pf)==1){
            plottype <- "p"
        } else{
              if (x$method=="quantile"){
                  plottype <- "b"
              } else{
                    plottype <- "l"
                }
          }
        pf <- na.omit(pf)
        if (x$model.type=="survival" && x$type!="survival"){
            lines(1-pf$Pred,1-pf$Obs,col=x$col[f],lwd=x$lwd[f],lty=x$lty[f],type=plottype)
        }
        else
            lines(pf$Pred,pf$Obs,col=x$col[f],lwd=x$lwd[f],lty=x$lty[f],type=plottype)
    }
    if (x$bars) {
        stopifnot(NF==1)
        showBars()
    }else{
         nix <- lapply(1:NF,function(f)showCal(f))
         do.call("legend",control$legend)
     }
    # }}}
    # {{{ axes
    if (x$axes){
        if (x$percent){
            control$axis2$labels <- paste(100*control$axis2$at,"%")
            control$axis1$labels <- paste(100*control$axis1$at,"%")
        }
        if (!x$bars)
            do.call("axis",control$axis1)
        mgp2 <- control$axis2$mgp
        if (length(mgp2)>0){
            oldmgp <- graphics::par()$mgp
            graphics::par(mgp=mgp2)
            control$axis2 <- control$axis2[-match("mgp",names(control$axis2),nomatch=0)]
            title(ylab=x$ylab)
        }
        do.call("axis",control$axis2)
        if (length(mgp2)>0){
            graphics::par(mgp=oldmgp)
        }
    }
    invisible(NULL)
    # }}}
}


#----------------------------------------------------------------------
### plot.calibrationPlot.R ends here
