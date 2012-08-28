
library(tools)
library(snow)
library(snowfall)

TestQueries = "test_query.iddb"
DataDir = "data"
ChemDb = file.path(DataDir,"chem.db")
Main = file.path(DataDir,"main.iddb")
#numSamples=1000
numSamples=20

cdbCachedSize=NA
cdbSize <- function() {
	if(is.na(cdbCachedSize))
		cdbCachedSize = length(readLines(Main))
	cdbCachedSize
}
embedCoord <- function(s,len,coords) {
	.Call("embedCoord",s,as.integer(len),as.double(coords))
}
embedCoordTest <- function(r,d,refCoords,coords) {
	.Call("embedCoordTest",as.integer(r),as.integer(d),as.double(refCoords),as.double(coords))
}

db_builder.atompair <- function(input,output)
	batch_sdf_parse(input,output)
db_subset.atompair = function(db,iddb,output)
	db_subset(db,iddb,output)

db2db_distance.atompair <- function(db,db2=NA,iddb1=NA,iddb2=NA,file=NA)
{
	if(!is.na(file)){
		if(is.na(db2) && ! is.na(iddb1) && ! is.na(iddb2)){
			cat(" iddb files to file\n")
			db2db_distance2file(db,iddb1,iddb2,file)
		}else if(!is.na(db2) &&  is.na(iddb1) &&  is.na(iddb2)){
			cat("2 real dbs to file\n")
			db2db_distance2file(db,db2,file)
		}else{
			stop("bad argument list\n")
		}
	}else{
		if(is.na(db2) && ! is.na(iddb1) && ! is.na(iddb2)){
			cat(" iddb files\n")
			return(.Call("db2db_distance_iddb",as.character(db),as.character(iddb1),as.character(iddb2)))
		}else if(!is.na(db2) &&  is.na(iddb1) &&  is.na(iddb2)){
			cat("2 real dbs\n")
			return(.Call("db2db_distance_db",as.character(db),as.character(db2)))
		}else{
			stop("bad argument list\n")
		}
	}
}


eiInit <- function(compoundDb,measure=NA,dbBuilder = db_builder.atompair)
{
	cat("eiInit")
	if(!file.exists("data"))
		dir.create("data")

	#if(!is.na(measure)){
		##write in config file
	#}
	if(!file.exists(ChemDb)){
		numCompounds = dbBuilder(compoundDb,ChemDb)
		writeIddb(1:numCompounds,Main)

	}
}

eiMakeDb <- function(r,d,db2dbDistance=db2db_distance.atompair,
									dir=".",refIddb=NA,cl=makeCluster(1,type="SOCK"))
{
	workDir=file.path(dir,paste("run",r,d,sep="-"))
	if(!file.exists(workDir))
		dir.create(workDir)

	queryIds=NA
	#get reference compounds
	if(is.na(refIddb)){
		prefix<- paste(sample(c(0:9,letters),32,replace=TRUE),collapse="")
		refIddb=file.path(workDir,paste(prefix,"cdb",sep="."))
		queryIds=genRefs(r,numSamples,refIddb)
	}
	if(is.na(queryIds[1]))
		queryIds=readIddb(file.path("data",TestQueries))
	
	selfDistFile <- paste(refIddb,"distmat",sep=".")
	coordFile <- paste(selfDistFile,"coord",sep=".")
	ref2AllDistFile <- paste(refIddb,"distances",sep=".")
	embeddedFile <- file.path(workDir,sprintf("coord.%d-%d",r,d))
	embeddedQueryFile <- file.path(workDir,sprintf("coord.query.%d-%d",r,d))

	#embed references in d dimensional space 
	coords <- if(file.exists(coordFile)){
		table.read(coordsFile)
	}else{
		#compute pairwise distances for all references
		if(! file.exists(selfDistFile)){
			cat("generateding selfDistFile\n")
			db2dbDistance(ChemDb,iddb1=refIddb,iddb2=refIddb,file=selfDistFile)
		}
		selfDist<-read.table(selfDistFile)

		#use MDS and pairwise distances to embed references
		coords <- cmdscale(selfDist,d)
		write.table(coords,file=coordFile,row.names=F,col.names=F)
		coords
	}

	#compute dist between refs and all compounds
	if(!file.exists(ref2AllDistFile))
		db2dbDistance(ChemDb,iddb1=Main,iddb2=refIddb,file=ref2AllDistFile)
	
	#each job needs: R, D, coords, a chunk of distance data
	solver <- getSolver(r,d,coords)	
	distConn <- file(ref2AllDistFile,"r")

	numJobs=length(cl)
	jobSize = as.integer(cdbSize() / numJobs + 1) #make last job short

	dataBlocks = Map(function(x)
		strsplit(readLines(distConn,jobSize),"\\s+"),1:numJobs)

	clusterApply(cl,1:numJobs, 
		function(i) { # job i has indicies [(i-1)*jobSize+1, i*jobSize]
			solver <- getSolver(r,d,coords)	
			data = sapply(dataBlocks[[i]],function(x) 
								embedCoord(solver,d,as.numeric(x)))

			write.table(t(data),
				file=file.path(workDir,paste(r,d,i,sep="-")),
				row.names=F,col.names=F)

			#list indexes for this job, see which of them are queries, then shift indexes back to
			# this jobs range before selecting from data.
			selected = queryIds[queryIds %in% ((i-1)*jobSize+1):(i*jobSize)]-((i-1)*jobSize)
			# R magically changes the data type depending on the size, yay!
			qd = if(length(selected)==1)  t(data[,selected]) else t(data)[selected, ]
			write.table(qd ,
				file=file.path(workDir,paste("q",r,d,i,sep="-")),
				row.names=F,col.names=F)
		})
	close(distConn)

	system(paste("cat",
					 paste(Map(function(x) file.path(workDir,paste(r,d,x,sep="-")),1:numJobs),collapse=" "),
					 ">",embeddedFile))
	system(paste("cat",
					 paste(Map(function(x) file.path(workDir,paste("q",r,d,x,sep="-")),1:numJobs),collapse=" "),
					 ">",embeddedQueryFile))

	Map(function(x) unlink(file.path(workDir,paste(r,d,x,sep="-"))),1:numJobs)
	Map(function(x) unlink(file.path(workDir,paste("q",r,d,x,sep="-"))),1:numJobs)

	binaryCoord(embeddedFile,file.path(workDir,sprintf("matrix.%d-%d",r,d)),d)
	binaryCoord(embeddedQueryFile,
		file.path(workDir,sprintf("matrix.query.%d-%d",r,d)),d)

	return(file.path(workDir,sprintf("matrix.%d-%d",r,d)))
}
eiQuery <- function(r,d,refIddb,queryFile,
	dir=".",db2dbDistance=db2db_distance.atompair,
	dbBuilder=db_builder.atompair,dbSubset=db_subset.atompair)
{
		tmpDir=tempdir()
		queryDb = file.path(tmpDir,"query.db")
		refDb = refDb(refIddb,dbSubset)
		query2RefDistFile = file.path(tmpDir,"query2refs.dist")

		numQueries = dbBuilder(queryFile,queryDb)

		allNames = readLines(file.path(dir,paste(ChemDb,"names",sep=".")))
		queryNames = readLines(paste(queryDb,"names",sep="."))

# embedding
		query2RefDists = db2dbDistance(queryDb,db2=refDb)
		coordFile=paste(refIddb,"distmat","coord",sep=".")
		coords = as.matrix(read.table(coordFile))
		solver = getSolver(r,d,coords)
		embedded = t(apply(query2RefDists,c(1),
			function(x) embedCoord(solver,d,x)))
# end embedding

		print(embedded)


		unlink(tmpDir,recursive=T)
}
writeIddb <- function(data, file)
		write.table(data,file,quote=FALSE,col.names=FALSE,row.names=FALSE)
readIddb <- function(file) as.numeric(readLines(file))
readNames <- function(file) as.numeric(readLines(file))


# randomly select n reference compounds. Also sample and stash away
# numSamples query compounds that are not references for later
# testing
genRefs <- function(n,numSamples,refFile)
{
	testQueryFile <-file.path("data",TestQueries)
	mainIds <- readIddb(Main)
	queryIds <- if(file.exists(testQueryFile)) {
			readIddb(testQueryFile)
		}else{
			ids <- sort(sample(mainIds,numSamples))
			writeIddb(ids,testQueryFile)
			ids
		}
	refIds = sort(sample(setdiff(mainIds,queryIds),n))
	writeIddb(refIds,refFile)
	queryIds
}
refDb <- function(refIddb,dbSubset,refDb=paste(refIddb,"db",sep="."))
{
	if(!file.exists(refDb))
		dbSubset(ChemDb,refIddb,refDb)
	return(refDb)
}


