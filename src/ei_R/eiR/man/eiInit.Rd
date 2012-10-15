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
      "dbSubset", and "db2dbDistance".
   }
}
\details{
   EiInit can take either an SDFset, or a filename.  SDF is supported
   by default. Other formats can be supported through customization.
   It might complain if your SDF file does not
   follow the SDF specification. If this happens, you can create an
   SDFset with the \code{read.SDFset} command and then use that
   instead of the filename.  EiInit will create  a folder called
   'data'. Commands should always be executed in the folder containing
   this directory (ie, the parent directory of "data"), or else
   specify the location of that directory with the \code{dir} option.

   All \code{eiXXXX} methods accept a \code{measure} parameter.
   A measure object handles dealing with the compound representations
   and also computing distances between compounds.
   It is a list with three fields, named: dbBuilder, dbSubset, and
   db2dbDistance. Each of these is a function.

   As an example, here is the default measure, which works for sdf
   formated files and computes distances using atompair descriptors:

   \preformatted{
atompairMeasure = list(
dbBuilder = function(input,output)
   batch_sdf_parse(input,output),
dbSubset = function(db,iddb,output)
   db_subset(db,iddb,output),
db2dbDistance = function(db,db2=NA,iddb1=NA,iddb2=NA,file=NA)
{
   if(!is.na(file)){
      if(is.na(db2) && ! is.na(iddb1) && ! is.na(iddb2)){
         db2db_distance2file(db,iddb1,iddb2,file)
      }else if(!is.na(db2) &&  is.na(iddb1) &&  is.na(iddb2)){
         db2db_distance2file(db,db2,file)
      }else{
         stop("bad argument list")
      }
   }else{
      if(is.na(db2) && ! is.na(iddb1) && ! is.na(iddb2)){
         return(.Call("db2db_distance_iddb",as.character(db),
            as.character(iddb1),as.character(iddb2)))
      }else if(!is.na(db2) &&  is.na(iddb1) &&  is.na(iddb2)){
         return(.Call("db2db_distance_db",as.character(db),
            as.character(db2)))
      }else{
         stop("bad argument list\n")
      }
   }
}
)
   }

    

   \code{DbBuilder} transforms the raw input file into whatever
   format will be used by dbSubset and db2dbDistance. It should take
   the name of an input file and the name of an output file to create.

   \code{DbSubset} extracts a subset of compounds from one file and
   puts them in another file. The compounds to extract are defined by 
   an IDDB file. It should take the name of a compound
   file, in the format output by dbBuilder, an IDDB file, and the name
   of the output file

   \code{Db2dbDistance} computes distances between compounds. There
   are two ways to use this function. The first is to specify both
   \code{db} and \code{db2}. This should compute all pairwise
   distances between compounds in these two files. The second way is
   to specify \code{db}, \code{iddb1}, and \code{iddb2}. In this
   case distances are computed only between compounds given in the two
   IDDB files. \code{db} should contain all compounds referenced by
   the IDDB files. In both cases, the \code{file} can optionally be
   given. If given, the output should be written to a file. The output
   should be lines of whitespace separated numbers. Each line should
   contain the distances of the first element of the first db against
   each member of the second db. If \code{file} is not given, the
   result should be returned as a matrix.



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
