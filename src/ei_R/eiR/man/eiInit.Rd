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
   eiInit(compoundDb, dir=".", measure = atompairMeasure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compoundDb}{
     The filename of a list of compounds in a format understood by the
     given measure. For the default measure, this will be SDF.
   }
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
  \item{measure}{
      This should be a list containing three functions named "dbBuilder",
      "dbSubset", and "db2dbDistance". See help on "measure" for details.
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
   dir=file.path(tempdir(),"init")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
