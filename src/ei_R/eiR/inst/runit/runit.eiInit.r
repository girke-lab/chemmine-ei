
library(eiR)
library(snow)

options(warn=2)
r<- 50
d<- 40
N<- 122
j=1
runDir<-paste("run",r,d,sep="-")

test.aa.eiInit <- function() {
   data(example_compounds)

   cat(paste(paste(example_compounds,collapse="\n"),"\n",sep=""),file="example_compounds.sdf")
   eiInit("example_compounds.sdf")
   checkTrue(file.exists(file.path("data","chem.db")))
   checkTrue(file.exists(file.path("data","chem.db.names")))
   checkTrue(file.exists(file.path("data","main.iddb")))
   i <- readLines(file.path("data","main.iddb"))
   checkEquals(length(i),N)
}

test.ca.eiQuery <- function(){

   data(example_compounds)
   cat(
      paste(
         paste(
            example_compounds[1:which(example_compounds=="$$$$")[2]],
            collapse="\n"),
         "\n",sep=""),
      file="example_queries.sdf")
   matches<-dir(runDir,pattern=".cdb$",full.names=T)
   checkEquals(length(matches),1)
   refIddb = matches[1]
   results = eiQuery(r,d,refIddb,"example_queries.sdf",K=15)
   checkTrue(length(results$distance) != 0)
   checkTrue(all(results$distance <= 1))
   checkEquals(results$distance[1],0)
}
test.ba.eiMakeDb <- function() {

     runChecks = function(){
      checkMatrix(".cdb$",r,1)
      checkMatrix(".cdb.distmat$",r,r)
      checkMatrix(".cdb.distmat.coord$",r,d)
      checkMatrix(".cdb.distances$",N,r)
      checkMatrix(sprintf("coord.%d-%d",r,d),N,d)
      checkMatrix(sprintf("coord.query.%d-%d",r,d),20,d)
      checkTrue(file.info(file.path(runDir,sprintf("matrix.%d-%d",r,d)))$size>0)
      checkTrue(file.info(file.path(runDir,sprintf("matrix.query.%d-%d",r,d)))$size>0)
      Map(function(x)
         checkTrue(!file.exists(file.path(runDir,paste(r,d,x,sep="-")))),1:j)
      Map(function(x)
         checkTrue(!file.exists(file.path(runDir,paste("q",r,d,x,sep="-")))),1:j)
   }
   t=system.time(eiMakeDb(r,d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile="")))
   runChecks()
   print(t)
   #system.time(eiMakeDb(r,d,cl=makeCluster(j,type="SOCK",outfile="")))
   #runChecks()
}
test.da.eiPerformanceTest <- function() {
   r = eiPerformanceTest(r,d,K=22)
   checkMatrix("eucsearch.50-40",20,N)
   checkMatrix("indexed",20,22)
}

test.aaaaa.cleanup<- function(){
   junk <- c("data","example_compounds.sdf","example_queries.sdf","run-50-40")
   unlink(junk,recursive=T)
}

checkMatrix <- function(pattern,x,y){
   matches<-dir(runDir,pattern=pattern,full.names=T)
   checkEquals(length(matches),1)
   file <- matches[1]
   checkTrue(file.info(file)$size>0)
   checkEquals(dim(read.table(file)),c(x,y))
}

#test.snow = function() {
#   options(warn=2)
#   options(error=traceback)
#   j=4
#
#   f=function(i) {
#      library(eiR)
#			#.Call("embedCoord",4,as.integer(3),as.double(1:5))
#			embedCoord(5,2,1:5)
#			t(sapply(1:d,function(x) x*x))
#         3
#		}
#   serFile=file("fun","w")
#   serialize(f,serFile)
#   close(serFile)
#   serFile=file("fun","r")
#   f2=unserialize(serFile)
#   checkEquals(f,f2)
#
#
#   #cl=makeCluster(j,type="SOCK",outfile="")
#   #clusterApply(cl,1:length(cl),function(x) x*2)
#   #clusterApply(cl,1:length(cl),f)
#   #cat(paste(result,collapse=" "))
#}


