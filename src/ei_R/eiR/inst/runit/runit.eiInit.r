
library(eiR)

test.eiInit <- function() {
   data(example_compounds)

   cat(paste(paste(example_compounds,collapse="\n"),"\n",sep=""),file="example_compounds.sdf")
   eiInit("example_compounds.sdf")
   checkTrue(file.exists(file.path("data","chem.db")))
   checkTrue(file.exists(file.path("data","chem.db.names")))
   checkTrue(file.exists(file.path("data","main.iddb")))
   i <- readLines(file.path("data","main.iddb"))
   checkEquals(length(i),122)
}


test.eiMakeDb <- function() {
   options(warn=2)
   r<- 50
   d<- 40
   runDir<-paste("run",r,d,sep="-")
   eiMakeDb(r,d);

   checkMatrix <- function(pattern,x,y){
      matches<-dir(runDir,pattern=pattern,full.names=T)
      checkEquals(length(matches),1)
      file <- matches[1]
      checkTrue(file.info(file)$size>0)
      checkEquals(dim(read.table(file)),c(x,y))
   }
   checkMatrix(".cdb$",r,1)
   checkMatrix(".cdb.distmat$",r,r)
   checkMatrix(".cdb.distmat.coord$",r,d)
   checkMatrix(".cdb.distances$",122,r)
}

test.aaaaa.cleanup<- function(){
   junk <- c("data","example_compounds.sdf","run-50-40")
   unlink(junk,recursive=T)
}
