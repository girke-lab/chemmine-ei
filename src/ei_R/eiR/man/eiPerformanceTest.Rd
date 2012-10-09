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
      The number of references used to build the database you wish to
      query agains.
   }
  \item{d}{
      The number of dimensions used to build the database you wish to
      query agains.
   }
 \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
  \item{measure}{
      This should be a list containing three functions named "dbBuilder",
      "dbSubset", and "db2dbDistance". See help on "measure" for details.
   }

  \item{K}{
      Number of neighbors
   }
  \item{W}{
     Tunable LSH parameter. See LSHKIT page for details.
     {http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html}
   }
  \item{M}{
     Tunable LSH parameter. See LSHKIT page for details.
     {http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html}
   }
  \item{L}{
     Number of hash tables
   }
  \item{T}{
     Number of probes
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
   dir=file.path(tempdir(),"perf")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)

   #create compound db
   eiMakeDb(r,d,numSamples=20,dir=dir,cl=makeCluster(j,type="SOCK",outfile=""))

   eiPerformanceTest(r,d,dir=dir,K=22)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
