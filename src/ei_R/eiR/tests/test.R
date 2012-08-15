
cat("testing stuff 123")
library(eiR)

data(example_compounds)
ls()

cat(paste(example_compounds,collapse="\n"),file="example_compounds.sdf")
eiInit("example_compounds.sdf")


cat("done testing")
