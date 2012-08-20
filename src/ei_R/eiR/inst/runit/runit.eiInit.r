
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


