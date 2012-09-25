\name{eiMakeDb}
\alias{eiMakeDb}
\title{
   Create an embedded database
}
\description{
   Uses the initalized compound data to create an embedded compound
   databbase with \code{r} reference compounds in \code{d} dimensions.
}
\usage{
   eiMakeDb(r, d, measure = atompairMeasure, dir = ".", refIddb = NA, numSamples = cdbSize() * 0.1, cl = makeCluster(1, type = "SOCK"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{measure}{
%%     ~~Describe \code{measure} here~~
}
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{refIddb}{
%%     ~~Describe \code{refIddb} here~~
}
  \item{numSamples}{
%%     ~~Describe \code{numSamples} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
   No value is returned. Creates files in \code{dir} ("run-r-d" by
   default).
}
\author{
   Kevin Horan
}


\seealso{
   \code{\link{eiInit}}
   \code{\link{eiPerformanceTest}}
   \code{\link{eiQuery}}
}
\examples{
   library(snow)

   r<- 50
   d<- 40
   j=1

   #initialize 
   data(sdfsample)
   eiInit(sdfsample)

   #create compound db
   eiMakeDb(r,d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile=""))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
