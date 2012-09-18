
library(tools)
library(snow)
library(snowfall)
library(ChemmineR)

DataDir = "data"
TestQueries = file.path(DataDir,"test_query.iddb")
TestQueryResults=file.path(DataDir,"chemical-search.results")
ChemDb = file.path(DataDir,"chem.db")
Main = file.path(DataDir,"main.iddb")

#default lsh parameters
defK=6
defW = 1.39564
defM=19
defL=10 
defT=30 

debug=T

cdbCachedSize=NA
cdbSize <- function() {
	if(is.na(cdbCachedSize))
		cdbCachedSize = length(readLines(Main))
	cdbCachedSize
}
embedCoord <- function(s,len,coords) 
	.Call("embedCoord",s,as.integer(len),as.double(coords))

embedCoordTest <- function(r,d,refCoords,coords) 
	.Call("embedCoordTest",as.integer(r),as.integer(d),as.double(refCoords),as.double(coords))

# requires one query per column, not per row
lshsearch <- function(queries,matrixFile,
	W=NA,H=NA,M=NA,L=NA,K=NA,T=NA,R=NA) 
{
	print(c(W,H,M,L,K,T,R))
	.Call("lshsearch",queries,as.character(matrixFile),
		as.double(W),as.integer(H),as.integer(M),as.integer(L),
		as.integer(K),as.integer(T), as.double(R))
}
atompairMeasure = list(
	dbBuilder = function(input,output)
		batch_sdf_parse(input,output),
	dbSubset = function(db,iddb,output)
		db_subset(db,iddb,output),
	db2dbDistance = function(db,db2=NA,iddb1=NA,iddb2=NA,file=NA)
	{
		if(!is.na(file)){
			if(is.na(db2) && ! is.na(iddb1) && ! is.na(iddb2)){
				#if(debug) print(" iddb files to file")
				db2db_distance2file(db,iddb1,iddb2,file)
			}else if(!is.na(db2) &&  is.na(iddb1) &&  is.na(iddb2)){
				#if(debug) print("2 real dbs to file")
				db2db_distance2file(db,db2,file)
			}else{
				stop("bad argument list")
			}
		}else{
			if(is.na(db2) && ! is.na(iddb1) && ! is.na(iddb2)){
				#if(debug) print(" iddb files")
				return(.Call("db2db_distance_iddb",as.character(db),as.character(iddb1),as.character(iddb2)))
			}else if(!is.na(db2) &&  is.na(iddb1) &&  is.na(iddb2)){
				#if(debug) print("2 real dbs")
				return(.Call("db2db_distance_db",as.character(db),as.character(db2)))
			}else{
				stop("bad argument list\n")
			}
		}
	}
)

eiInit <- function(compoundDb,measure=atompairMeasure)
{
	if(debug) print("eiInit")
	if(!file.exists(DataDir))
		dir.create(DataDir)

	if(!file.exists(ChemDb)){
		numCompounds = measure$dbBuilder(compoundDb,ChemDb)
		writeIddb(1:numCompounds,Main)

	}
}

eiMakeDb <- function(r,d,measure=atompairMeasure,
				dir=".",refIddb=NA,numSamples=cdbSize()*0.1,
				cl=makeCluster(1,type="SOCK"))
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
		queryIds=readIddb(TestQueries)
	
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
			print("generateding selfDistFile")
			measure$db2dbDistance(ChemDb,iddb1=refIddb,iddb2=refIddb,file=selfDistFile)
		}
		selfDist<-read.table(selfDistFile)

		#use MDS and pairwise distances to embed references
		coords <- cmdscale(selfDist,d)
		write.table(coords,file=coordFile,row.names=F,col.names=F)
		coords
	}

	#compute dist between refs and all compounds
	if(!file.exists(ref2AllDistFile))
		measure$db2dbDistance(ChemDb,iddb1=Main,iddb2=refIddb,file=ref2AllDistFile)
	
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
		dir=".",measure=atompairMeasure,
		K=defK, W = defW, M=defM,L=defL ,T=defT )
{
		tmpDir=tempdir()
		workDir=file.path(dir,paste("run",r,d,sep="-"))
		queryDb = file.path(tmpDir,"queries.db")
		refDb = refDb(refIddb,measure)
		query2RefDistFile = file.path(tmpDir,"query2refs.dist")

		#reformat query file
		numQueries = measure$dbBuilder(queryFile,queryDb)

		#read names
		allNames = readLines(file.path(dir,paste(ChemDb,"names",sep=".")))
		queryNames = readLines(paste(queryDb,"names",sep="."))

		#embed queries in search space
		embeddedQueries = embedFromRefs(r,d,refIddb,measure,queryDb,db2=refDb)

		#search for nearby compounds
		if(debug) print(embeddedQueries)
		matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
		hits = search(embeddedQueries,matrixFile,
							queryDb,measure,K=K,W=W,M=M,L=L,T=T)
		if(debug) print(hits)

		#fetch names for queries and hits and put in a data frame
		results = data.frame(query=rep(NA,K*numQueries),
								  target = rep(NA,K*numQueries),
								  distance=rep(NA,K*numQueries))
		i=1
		lapply(1:numQueries,function(queryIndex)
			lapply(1:K,function(hitIndex){
				results[i,"query"]<<-queryNames[queryIndex]
				results[i,"target"]<<-allNames[hits[[queryIndex]][hitIndex,1]]
				results[i,"distance"]<<- hits[[queryIndex]][hitIndex,2]
				i<<-i+1
			}))
		if(debug) print(results)
		return(results)
}
#expects one query per column
search <- function(queries,matrixFile,queryDb,measure,K,...)
{
		neighbors = lshsearch(queries,matrixFile,K=K,...)

		tmpDir=tempdir()
		#compute distance between each query and its candidates	
		Map(function(x) refine(neighbors[x,,],queryDb,x,K,measure,tmpDir),
			1:(dim(queries)[2]))
}
#fetch coords from refIddb.distmat.coord and call embed
embedFromRefs <- function(r,d,refIddb,measure,...)
{
		coordFile=paste(refIddb,"distmat","coord",sep=".")
		coords = as.matrix(read.table(coordFile))
		embed(r,d,coords,measure,...)
}
#take referenace coords, and args to compute distance matrix
#return the emedding of each row in the matrix
embed <- function(r,d,coords, measure,...)
{
		query2RefDists = measure$db2dbDistance(...)
		solver = getSolver(r,d,coords)
		embeddedQueries = apply(query2RefDists,c(1),
			function(x) embedCoord(solver,d,x))
#
}
refine <- function(lshNeighbors,queriesCdb,queryIndex,limit,measure,tmpDir=tempdir())
{
	queryIddb=file.path(tmpDir,"query.iddb")
	queryDb=file.path(tmpDir,"query.db")
	candidatesIddb=file.path(tmpDir,"candidates.iddb")
	candidatesDb=file.path(tmpDir,"candidates.db")

	queryIddb=file.path(tmpDir,"query.iddb")
	writeIddb(c(queryIndex),queryIddb)
	measure$dbSubset(queriesCdb,queryIddb,queryDb)

	queryIddb=file.path(tmpDir,"query.iddb")
	writeIddb(lshNeighbors[,1],candidatesIddb)
	measure$dbSubset(ChemDb,candidatesIddb,candidatesDb)

	queryIddb=file.path(tmpDir,"query.iddb")
	d=measure$db2dbDistance(queryDb,db2=candidatesDb)
	#if(debug) print(str(d))
	lshNeighbors[,2]=d #measure$db2dbDistance(queryDb,db2=candidatesDb)
	limit = min(limit,length(lshNeighbors[,2]))
	queryIddb=file.path(tmpDir,"query.iddb")
	lshNeighbors[order(lshNeighbors[,2])[1:limit],]
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
	testQueryFile <-TestQueries
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
refDb <- function(refIddb,measure,refDb=paste(refIddb,"db",sep="."))
{
	if(!file.exists(refDb))
		measure$dbSubset(ChemDb,refIddb,refDb)
	return(refDb)
}
genTestQueryResults <- function(measure)
{
	if(file.exists(TestQueryResults))
		return()

	out=file(TestQueryResults,"w")
	d=measure$db2dbDistance(ChemDb,iddb1=TestQueries,iddb2=Main)
	for(i in dim(d)[1])
		cat(paste(
				paste(1:dim(d)[2],d[i,],sep=":")[order(d[i,])[1:50000]],
				collapse=" "),"\n",file=out)	
	close(out)
}
eiPerformanceTest <- function(r,d,measure=atompairMeasure,
	dir=".",K=defK, W = defW, M=defM,L=defL ,T=defT )
{
	workDir=file.path(dir,paste("run",r,d,sep="-"))
	genTestQueryResults(measure)
	eucsearch2file(file.path(workDir,sprintf("matrix.%s-%s",r,d)),
				 file.path(workDir,sprintf("matrix.query.%s-%s",r,d)),
				 50000,
				 file.path(workDir,sprintf("eucsearch.%s-%s",r,d)))

	#evaluator TestQueryResuts eucsearch-r-d recall

	matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
	coordQueryFile =file.path(workDir,sprintf("coord.query.%d-%d",r,d))

	embeddedTestQueries = t(as.matrix(read.table(coordQueryFile)))
	hits = search(embeddedTestQueries,matrixFile,
						ChemDb,measure,K=K,W=W,M=M,L=L,T=T)
	out=file(file.path(workDir,"indexed"),"w")
	#if(debug) print(hits)
	for(x in hits)
		cat(paste(x[,1],x[,2],sep=":",collapse=" "),"\n",file=out)
	close(out)

	#indexed_evalutator TestQueryResults indexed indexed.performance

	return(hits)
	
}


