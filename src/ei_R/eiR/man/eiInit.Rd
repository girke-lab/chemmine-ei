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
	eiInit(compoundDb,dir=".",format="sdf",descriptorType="ap",append=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compoundDb}{
	  Either a filename of and SDF file, or an SDFset.
   }
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
	\item{format}{
		The format of the data in compoundDb. Currenly only "sdf" is 
		supported.
	}
	\item{descriptorType}{
		The format of the descriptor. Currently supported values are "ap" for atom pair, and 
		"fp" for fingerprint.
	}
	\item{append}{
		If true the given compounds will be added to an existing database
		and the <data-dir>/Main.iddb file will be udpated with the new
		compound id numbers. This should not normally be used directly, use
		\code{\link{eiAdd}} instead to add new compounds to a database.
	}

}
\details{


   EiInit can take either an SDFset, or a filename.  SDF is supported
   by default.
   It might complain if your SDF file does not
   follow the SDF specification. If this happens, you can create an
   SDFset with the \code{read.SDFset} command and then use that
   instead of the filename.  
	
	EiInit will create  a folder called
   'data'. Commands should always be executed in the folder containing
   this directory (ie, the parent directory of "data"), or else
   specify the location of that directory with the \code{dir} option.

}
\value{
   A directory called "data" will have been created in the current working directory.
	The generated compound ids of the given compounds will be returned. These can be used to 
	reference a compound or set of compounds in other functions, such as \code{\link{eiQuery}}.
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
