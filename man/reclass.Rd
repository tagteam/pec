% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reclass.R
\name{reclass}
\alias{reclass}
\title{Retrospective risk reclassification table}
\usage{
reclass(
  object,
  reference,
  formula,
  data,
  time,
  cause,
  cuts = seq(0, 100, 25),
  digits = 2
)
}
\arguments{
\item{object}{Either a 
list with two elements. Each element should either
be a vector with probabilities, or an object for which
\code{predictSurvProb} or \code{predictEventProb} can extract predicted risk based on data.}

\item{reference}{Reference prediction model.}

\item{formula}{A survival formula as obtained either with
\code{prodlim::Hist} or \code{survival::Surv} which defines the
response in the \code{data}.}

\item{data}{Used to extract the response from the data and passed
on to \code{predictEventProb} to extract predicted event
probabilities.}

\item{time}{Time interest for prediction.}

\item{cause}{For competing risk models the cause of
interest. Defaults to all available causes.}

\item{cuts}{Risk quantiles to group risks.}

\item{digits}{Number of digits to show for the predicted risks}
}
\value{
reclassification tables: overall table and one conditional table for each cause and for subjects event free at time interest.
}
\description{
Retrospective table of risks predicted by two different methods, models, algorithms
}
\details{
All risks are multiplied by 100 before
}
\examples{
\dontrun{
library(survival)
set.seed(40)
d <- prodlim::SimSurv(400)
nd <- prodlim::SimSurv(400)
Models <- list("Cox.X2"=coxph(Surv(time,status)~X2,data=d,x=TRUE,y=TRUE),
               "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=d,x=TRUE,y=TRUE))
rc <- reclass(Models,formula=Surv(time,status)~1,data=nd,time=5)
print(rc)
plot(rc)

set.seed(40)
library(riskRegression)
library(prodlim)
dcr <- prodlim::SimCompRisk(400)
ndcr <- prodlim::SimCompRisk(400)
crPred5 <- list("X2"=predictEventProb(CSC(Hist(time,event)~X2,data=dcr),newdata=ndcr,times=5),
                "X1+X2"=predictEventProb(CSC(Hist(time,event)~X1+X2,data=dcr),newdata=ndcr,times=5))
rc <- reclass(crPred5,Hist(time,event)~1,data=ndcr,time=3)
print(rc)

reclass(crPred5,Hist(time,event)~1,data=ndcr,time=5,cuts=100*c(0,0.05,0.1,0.2,1))
}
}
\seealso{
predictStatusProb
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
