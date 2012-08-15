
db_builder.atompair <- function(input,output)
	batch_sdf_parse(input,output)
eiInit <- function(compoundDb,measure=NA,db_builder = db_builder.atompair)
{
	cat("eiInit")
	if(!file.exists("main"))
		dir.create("data")

	#if(!is.na(measure)){
		##write in config file
	#}
	chemdb <- file.path("data","chem.db")
	if(!file.exists(chemdb)){
		numCompounds = db_builder(compoundDb,chemdb)
		mainIddb = file.path("data","main.iddb")
		write.table(1:numCompounds,mainIddb)

	}

}
