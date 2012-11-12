\name{eiAdd}
\alias{eiAdd}
\title{
   Add new compounds
}
\description{
   Add additional compounds to and existing database
}
\usage{
   eiAdd(r,d,refIddb,additions,dir=".",measure=atompairMeasure)
}
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
      you wish to append to. This should almost always be the file:
      run-\code{r}-\code{d}/<long random string>.cdb.
   }
  \item{additions}{
      The compounds to add. This can be either an SDFset, or a file
      containg 1 or more query compounds. The file should be in whatever
      format is supported by the given measure. The default is SDF.
   }
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
  \item{measure}{
      This should be a list containing three functions named "dbBuilder",
      "dbSubset", and "db2dbDistance". See description in 
      \code{\link{eiInit}} 
   }
}
\details{
   New Compounds can be added to an existing database, however, the 
   reference compounds cannot be changed. 
}
\examples{

   library(snow)
   r<- 50
   d<- 40
   j<- 1

   #initialize
   data(sdfsample)
   dir=file.path(tempdir(),"add")
   dir.create(dir)
   eiInit(sdfsample[1:99],dir=dir)

   #create compound db
   eiMakeDb(r,d,numSamples=20,dir=dir,
      cl=makeCluster(j,type="SOCK",outfile=""))

   runDir<-file.path(dir,paste("run",r,d,sep="-"))
   matches<-dir(runDir,pattern=".cdb$",full.names=TRUE)
   refIddb = matches[1]

   #find compounds similar two each query
   eiAdd(r,d,refIddb,sdfsample[100],dir=dir)


}
\author{ Kevin Horan }
\seealso{
   \code{\link{eiMakeDb}}
   \code{\link{eiPerformanceTest}}
   \code{\link{eiQuery}}
}

