\name{eiQuery}
\alias{eiQuery}
\title{
   Perform a query on an embedded database
}
\description{
   Finds similar compounds for each query.
}
\usage{
   eiQuery(r, d, refIddb, queries, dir = ".", measure = atompairMeasure, 
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
  \item{refIddb}{
      An Iddb formatted file containing the reference IDs of the database
      you wish to query against. This should almost always be the file:
      run-\code{r}-\code{d}/<long random string>.cdb.
   }
  \item{queries}{
      This can be either an SDFset, or a file containg 1 or more
      query compounds. The file should be in whatever format is supported
      by the given measure. The default is SDF.
   }
  \item{dir}{
      %%     ~~Describe \code{dir} here~~
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
  \item{measure}{
      %%     ~~Describe \code{measure} here~~
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
   Returns a data frame with columns 'query', 'target', and
   'distance'. 'query' and 'target' are the compound IDs and
   distance is the distance between them, as computed by
   db2dbDistance. Query IDs are repeated for each matching target
   found.
}
\author{
   Kevin Horan
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
   \code{\link{eiInit}}
   \code{\link{eiMakeDb}}
   \code{\link{eiPerformanceTest}}
}
\examples{

   library(snow)
   r<- 50
   d<- 40
   j<- 1

   #initialize
   data(sdfsample)
   dir=file.path(tempdir(),"query")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)

   #create compound db
   eiMakeDb(r,d,numSamples=20,dir=dir,cl=makeCluster(j,type="SOCK",outfile=""))

   runDir<-file.path(dir,paste("run",r,d,sep="-"))
   matches<-dir(runDir,pattern=".cdb$",full.names=TRUE)
   refIddb = matches[1]

   #find compounds similar two each query
   results = eiQuery(r,d,refIddb,sdfsample[1:2],K=15,dir=dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
