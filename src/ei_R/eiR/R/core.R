
library(tools)

TestQueries = "test_query.iddb"
ChemDb = "chem.db"
Main = "main.iddb"
#numSamples=1000
numSamples=20

cdbCachedSize=NA
cdbSize <- function() {
	if(is.na(cdbCachedSize))
		cdbCachedSize = length(readLines(file.path("data",Main)))
	cdbCachedSize
}

db_builder.atompair <- function(input,output)
	batch_sdf_parse(input,output)

db2db_distance.atompair <- function(db,db2=NA,iddb1=NA,iddb2=NA,file)
{
	if(is.na(db2) && ! is.na(iddb1) && ! is.na(iddb2)){
		cat(" iddb files\n")
		db2db_distance2file(db,iddb1,iddb2,file)
	}else if(!is.na(db2) &&  is.na(iddb1) &&  is.na(iddb2)){
		cat("2 real dbs\n")
		db2db_distance2file(db,db2,file)
	}else{
		cat("bad argument list\n")
	}
}


eiInit <- function(compoundDb,measure=NA,db_builder = db_builder.atompair)
{
	cat("eiInit")
	if(!file.exists("data"))
		dir.create("data")

	#if(!is.na(measure)){
		##write in config file
	#}
	chemdb <- file.path("data",ChemDb)
	if(!file.exists(chemdb)){
		numCompounds = db_builder(compoundDb,chemdb)
		mainIddb = file.path("data",Main)
		writeIddb(1:numCompounds,mainIddb)

	}

}

eiMakeDb <- function(r,d,db2dbDistance=db2db_distance.atompair,
									dir=".",refIddb=NA)
{
	workDir=file.path(dir,paste("run",r,d,sep="-"))
	if(!file.exists(workDir))
		dir.create(workDir)

	#get reference compounds
	if(is.na(refIddb)){
		prefix<- paste(sample(c(0:9,letters),32,replace=TRUE),collapse="")
		refIddb=file.path(workDir,paste(prefix,"cdb",sep="."))
		genRefs(r,numSamples,refIddb)
	}
	
	selfDistFile <- paste(refIddb,"distmat",sep=".")
	coordFile <- paste(selfDistFile,"coord",sep=".")
	ref2AllDistFile <- paste(refIddb,"distances",sep=".")

	#embed references in d dimensional space 
	coords <- if(file.exists(coordFile)){
		table.read(coords)
	}else{
		#compute pairwise distances for all references
		if(! file.exists(selfDistFile)){
			cat("generateding selfDistFile")
			db2dbDistance(file.path("data",ChemDb),iddb1=refIddb,iddb2=refIddb,file=selfDistFile)
		}
		selfDist<-read.table(selfDistFile)

		#use MDS and pairwise distances to embed references
		coords <- cmdscale(selfDist,d)
		write.table(coords,file=coordFile,row.names=F,col.names=F)
		coords
	}

	#compute dist between refs and all compounds
	if(!file.exists(ref2AllDistFile))
		db2dbDistance(file.path("data",ChemDb),iddb1=file.path("data",Main),iddb2=refIddb,file=ref2AllDistFile)
	
	#each job needs: R, D, coords, a chunk of distance data

	#output: coord.r-d, coord.query.r-d: subset wth only testQueries
	#binarize coords

}
writeIddb <- function(data, file)
		write.table(data,file,quote=FALSE,col.names=FALSE,row.names=FALSE)
readIddb <- function(file) as.numeric(readLines(file))


# randomly select n reference compounds. Also sample and stash away
# numSamples query compounds that are not references for later
# testing
genRefs <- function(n,numSamples,refFile)
{
	testQueryFile <-file.path("data",TestQueries)
	mainIds <- readIddb(file.path("data",Main))
	queryIds <- if(file.exists(testQueryFile)) {
			readIddb(testQueryFile)
		}else{
			ids <- sort(sample(mainIds,numSamples))
			writeIddb(ids,testQueryFile)
			ids
		}
	refIds = sort(sample(setdiff(mainIds,queryIds),n))
	writeIddb(refIds,refFile)
}


