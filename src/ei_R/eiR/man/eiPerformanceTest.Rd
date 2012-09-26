\name{eiPerformanceTest}
\alias{eiPerformanceTest}
\title{
   Test the performance of LSH search
}
\description{
   Compares LSH search against both exhustive euclidean search and non-indexed
   nearest neighbor. %TODO: CHECK THIS
}
\usage{
   eiPerformanceTest(r, d, measure = atompairMeasure, dir = ".", 
      K = 6, W = 1.39564, M = 19, L = 10, T = 30)
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
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{W}{
%%     ~~Describe \code{W} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
   TBD
}
\author{
   Kevin Horan
}

\seealso{
   \code{\link{eiInit}}
   \code{\link{eiMakeDb}}
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

   eiPerformanceTest(r,d,K=22)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
