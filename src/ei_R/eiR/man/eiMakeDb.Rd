\name{eiMakeDb}
\alias{eiMakeDb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eiMakeDb(r, d, measure = atompairMeasure, dir = ".", refIddb = NA, numSamples = cdbSize() * 0.1, cl = makeCluster(1, type = "SOCK"))
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
  \item{refIddb}{
%%     ~~Describe \code{refIddb} here~~
}
  \item{numSamples}{
%%     ~~Describe \code{numSamples} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
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
function (r, d, measure = atompairMeasure, dir = ".", refIddb = NA, 
    numSamples = cdbSize() * 0.1, cl = makeCluster(1, type = "SOCK")) 
{
    workDir = file.path(dir, paste("run", r, d, sep = "-"))
    if (!file.exists(workDir)) 
        dir.create(workDir)
    queryIds = NA
    if (is.na(refIddb)) {
        prefix <- paste(sample(c(0:9, letters), 32, replace = TRUE), 
            collapse = "")
        refIddb = file.path(workDir, paste(prefix, "cdb", sep = "."))
        queryIds = genRefs(r, numSamples, refIddb)
    }
    if (is.na(queryIds[1])) 
        queryIds = readIddb(TestQueries)
    selfDistFile <- paste(refIddb, "distmat", sep = ".")
    coordFile <- paste(selfDistFile, "coord", sep = ".")
    ref2AllDistFile <- paste(refIddb, "distances", sep = ".")
    embeddedFile <- file.path(workDir, sprintf("coord.\%d-\%d", 
        r, d))
    embeddedQueryFile <- file.path(workDir, sprintf("coord.query.\%d-\%d", 
        r, d))
    coords <- if (file.exists(coordFile)) {
        table.read(coordsFile)
    }
    else {
        if (!file.exists(selfDistFile)) {
            print("generateding selfDistFile")
            measure$db2dbDistance(ChemDb, iddb1 = refIddb, iddb2 = refIddb, 
                file = selfDistFile)
        }
        selfDist <- read.table(selfDistFile)
        coords <- cmdscale(selfDist, d)
        write.table(coords, file = coordFile, row.names = F, 
            col.names = F)
        coords
    }
    if (!file.exists(ref2AllDistFile)) 
        measure$db2dbDistance(ChemDb, iddb1 = Main, iddb2 = refIddb, 
            file = ref2AllDistFile)
    solver <- getSolver(r, d, coords)
    distConn <- file(ref2AllDistFile, "r")
    numJobs = length(cl)
    jobSize = as.integer(cdbSize()/numJobs + 1)
    dataBlocks = Map(function(x) strsplit(readLines(distConn, 
        jobSize), "\\s+"), 1:numJobs)
    clusterApply(cl, 1:numJobs, function(i) {
        solver <- getSolver(r, d, coords)
        data = sapply(dataBlocks[[i]], function(x) embedCoord(solver, 
            d, as.numeric(x)))
        write.table(t(data), file = file.path(workDir, paste(r, 
            d, i, sep = "-")), row.names = F, col.names = F)
        selected = queryIds[queryIds \%in\% ((i - 1) * jobSize + 
            1):(i * jobSize)] - ((i - 1) * jobSize)
        qd = if (length(selected) == 1) 
            t(data[, selected])
        else t(data)[selected, ]
        write.table(qd, file = file.path(workDir, paste("q", 
            r, d, i, sep = "-")), row.names = F, col.names = F)
    })
    close(distConn)
    system(paste("cat", paste(Map(function(x) file.path(workDir, 
        paste(r, d, x, sep = "-")), 1:numJobs), collapse = " "), 
        ">", embeddedFile))
    system(paste("cat", paste(Map(function(x) file.path(workDir, 
        paste("q", r, d, x, sep = "-")), 1:numJobs), collapse = " "), 
        ">", embeddedQueryFile))
    Map(function(x) unlink(file.path(workDir, paste(r, d, x, 
        sep = "-"))), 1:numJobs)
    Map(function(x) unlink(file.path(workDir, paste("q", r, d, 
        x, sep = "-"))), 1:numJobs)
    binaryCoord(embeddedFile, file.path(workDir, sprintf("matrix.\%d-\%d", 
        r, d)), d)
    binaryCoord(embeddedQueryFile, file.path(workDir, sprintf("matrix.query.\%d-\%d", 
        r, d)), d)
    return(file.path(workDir, sprintf("matrix.\%d-\%d", r, d)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
