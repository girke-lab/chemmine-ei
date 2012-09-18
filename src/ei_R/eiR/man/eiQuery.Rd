\name{eiQuery}
\alias{eiQuery}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eiQuery(r, d, refIddb, queryFile, dir = ".", measure = atompairMeasure, K = defK, W = defW, M = defM, L = defL, T = defT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{refIddb}{
%%     ~~Describe \code{refIddb} here~~
}
  \item{queryFile}{
%%     ~~Describe \code{queryFile} here~~
}
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{measure}{
%%     ~~Describe \code{measure} here~~
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (r, d, refIddb, queryFile, dir = ".", measure = atompairMeasure, 
    K = defK, W = defW, M = defM, L = defL, T = defT) 
{
    tmpDir = tempdir()
    workDir = file.path(dir, paste("run", r, d, sep = "-"))
    queryDb = file.path(tmpDir, "queries.db")
    refDb = refDb(refIddb, measure)
    query2RefDistFile = file.path(tmpDir, "query2refs.dist")
    numQueries = measure$dbBuilder(queryFile, queryDb)
    allNames = readLines(file.path(dir, paste(ChemDb, "names", 
        sep = ".")))
    queryNames = readLines(paste(queryDb, "names", sep = "."))
    embeddedQueries = embedFromRefs(r, d, refIddb, measure, queryDb, 
        db2 = refDb)
    if (debug) 
        print(embeddedQueries)
    matrixFile = file.path(workDir, sprintf("matrix.\%d-\%d", r, 
        d))
    hits = search(embeddedQueries, matrixFile, queryDb, measure, 
        K = K, W = W, M = M, L = L, T = T)
    if (debug) 
        print(hits)
    results = data.frame(query = rep(NA, K * numQueries), target = rep(NA, 
        K * numQueries), distance = rep(NA, K * numQueries))
    i = 1
    lapply(1:numQueries, function(queryIndex) lapply(1:K, function(hitIndex) {
        results[i, "query"] <<- queryNames[queryIndex]
        results[i, "target"] <<- allNames[hits[[queryIndex]][hitIndex, 
            1]]
        results[i, "distance"] <<- hits[[queryIndex]][hitIndex, 
            2]
        i <<- i + 1
    }))
    if (debug) 
        print(results)
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
