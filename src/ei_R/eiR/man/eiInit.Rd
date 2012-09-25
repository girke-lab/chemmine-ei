\name{eiInit}
\alias{eiInit}
\title{
   Initialize a compound database
}
\description{
   Takes the raw compound database in whatever format the given
   measure supports and creates a "data" directory.
}
\usage{
   eiInit(compoundDb, measure = atompairMeasure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compoundDb}{
%%     ~~Describe \code{compoundDb} here~~
}
  \item{measure}{
%%     ~~Describe \code{measure} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
   No value is returned. A directory called "data" will have been
   created in the current working directory.
}
\author{
   Kevin Horan
}


\seealso{
   \code{\link{eiMakeDb}}
   \code{\link{eiPerformanceTest}}
   \code{\link{eiQuery}}
}
\examples{
   data(sdfsample)
   eiInit(sdfsample)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
