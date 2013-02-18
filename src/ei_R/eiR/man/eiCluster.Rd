\name{eiCluster}
\alias{eiCluster}
\title{
	Cluster compounds
}
\description{
	Uses Jarvis-Patrick clustering to cluster the compound database using 
	the LSH algorithm to quickly find nearest neighbors.
}
\usage{
	eiCluster(r, d, K, minNbrs, dir = ".", cutoff = NA, descriptorType = "ap", 
		distance = apDistance, W = 1.39564, M = 19, L = 10, T = 30, type = "cluster")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
      The number of references used to build the database you wish to
      query against.
	}
  \item{d}{
      The number of dimensions used to build the database you wish to
      query against.
	}
  \item{K}{
	  The number of neighbors to consider for each compound.
	}
  \item{minNbrs}{
	  The minimum number of neighbors that two comopunds must have in common in order to be joined.
	}
	\item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
	\item{descriptorType}{
		The format of the descriptor. Currently supported values are "ap" for atom pair, and 
		"fp" for fingerprint.
	}
	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors.
	}

  \item{cutoff}{
	  Distance cutoff value. Compounds having a distance larger this this value will not
	  be included in the nearest neighbor table. Note that this is a distance value, not a similarity
	  value, as is often used in other ChemmineR functions.
	}

	\item{W}{
     Tunable LSH parameter. See LSHKIT page for details.
     \url{http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html}
   }
  \item{M}{
     Tunable LSH parameter. See LSHKIT page for details.
     \url{http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html}
   }
  \item{L}{
     Number of hash tables
   }
  \item{T}{
     Number of probes
   }
   \item{type}{
		If "cluster", returns a clustering, else, if "matrix", returns a nearest neighbor matrix.
	}
}
\details{

}
\value{
		If \code{type} is "cluster", returns a clustering.
		This will be a vector in which the names are the compound names, and the values are the cluster
		labels.
		Otherwise, if \code{type} is "matrix", returns a nearest neighbor matrix.
		This will be a matrix with a row for each compound. Each row will contain the index value of
		the neibhboring compounds. If there are not \code{K} neibhbors for a 
		compound, that row will be padded with NAs.

}
\author{
	Kevin Horan
}
}
\examples{

}
