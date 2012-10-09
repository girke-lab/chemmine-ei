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
   eiMakeDb(r, d, measure = atompairMeasure, dir = ".", refIddb = NA, 
               numSamples = cdbSize(dir) * 0.1, cl = makeCluster(1, type = "SOCK"))
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
  \item{refIddb}{
      An Iddb formatted file containing the reference IDs to use when
      creating this database. If not given, the reference IDs will be
      randomly chosen.
  }
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
  \item{measure}{
      This should be a list containing three functions named "dbBuilder",
      "dbSubset", and "db2dbDistance". See help on "measure" for details.
   }

  \item{numSamples}{
      The number of non-reference samples to be chosen now to be used
      later by the eiPerformanceTest function.
   }
  \item{cl}{
     A SNOW cluster can be given here to run this function in
     parrallel.
   }
}
\details{
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
   dir=file.path(tempdir(),"makedb")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)

   #create compound db
   eiMakeDb(r,d,numSamples=20,dir=dir,cl=makeCluster(j,type="SOCK",outfile=""))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
