\name{eiPerformanceTest}
\alias{eiPerformanceTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eiPerformanceTest(r, d, measure = atompairMeasure, dir = ".", K = defK, W = defW, M = defM, L = defL, T = defT)
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
function (r, d, measure = atompairMeasure, dir = ".", K = defK, 
    W = defW, M = defM, L = defL, T = defT) 
{
    workDir = file.path(dir, paste("run", r, d, sep = "-"))
    genTestQueryResults(measure)
    eucsearch2file(file.path(workDir, sprintf("matrix.\%s-\%s", 
        r, d)), file.path(workDir, sprintf("matrix.query.\%s-\%s", 
        r, d)), 50000, file.path(workDir, sprintf("eucsearch.\%s-\%s", 
        r, d)))
    matrixFile = file.path(workDir, sprintf("matrix.\%d-\%d", r, 
        d))
    coordQueryFile = file.path(workDir, sprintf("coord.query.\%d-\%d", 
        r, d))
    embeddedTestQueries = t(as.matrix(read.table(coordQueryFile)))
    hits = search(embeddedTestQueries, matrixFile, ChemDb, measure, 
        K = K, W = W, M = M, L = L, T = T)
    out = file(file.path(workDir, "indexed"), "w")
    for (x in hits) cat(paste(x[, 1], x[, 2], sep = ":", collapse = " "), 
        "\n", file = out)
    close(out)
    return(hits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
