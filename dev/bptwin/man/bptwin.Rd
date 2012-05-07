\name{bptwin}
\alias{bptwin}
\title{bptwin}
\description{
  bptwin
}
\usage{
bptwin(formula, data, id, zyg, twinnum, DZ, weight = NULL, 
    control = list(trace = 1), type = "ace", eqmean = TRUE, param = 0, 
    robustvar = TRUE, p, debug = FALSE, ...) 
}
\arguments{
  \item{formula}{}
  \item{data}{}
  \item{id}{}
  \item{zyg}{}
  \item{twinnum}{}
  \item{DZ}{}
  \item{weight}{}
  \item{control}{}
  \item{type}{}
  \item{eqmean}{}
  \item{param}{}
  \item{robustvar}{}
  \item{p}{}
  \item{debug}{}
  \item{\dots}{}
}
\value{
  bptwin
}
\details{
  bptwin
}
\seealso{
  The package \code{twinlm}.
}
\examples{
\dontrun{
bptwin(cancer~1,data=prt0,id="tvparnrs",zyg="zygbin",twinnum="tv",control=list(trace=1))
}
}
\author{Klaus K. Holst}
\keyword{regression}

