
library(tools)
library(snow)
library(snowfall)
library(ChemmineR)

DataDir = "data"
TestQueries = file.path(DataDir,"test_query.iddb")
TestQueryResults=file.path(DataDir,"chemical-search.results")
ChemPrefix="chem"
ChemDb = file.path(DataDir,paste(ChemPrefix,".db",sep=""))
ChemIndex = file.path(DataDir,paste(ChemPrefix,".index",sep=""))
Main = file.path(DataDir,"main.iddb")

#default lsh parameters
defK=6
defW = 1.39564
defM=19
defL=10 
defT=30 

#debug=TRUE
debug=FALSE

cdbCachedSize=NA
cdbSize <- function() {
	if(is.na(cdbCachedSize))
		cdbCachedSize = getSegmentSize(ChemDb)
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

addToChemIndex <- function(filename,size)
	cat(paste(filename,"\t",size,"\n",sep=""),file=ChemIndex,append=TRUE)
getSegmentSize <- function(filename)
	read.table(ChemIndex,sep="\t",row.names=1)[filename,]

applyOverIndex <- function(indexValues,f)
{
	index=read.table(ChemIndex,sep="\t",row.names=1)
	names=rownames(index)
	sums=sapply(1:length(index[[1]]),function(x) c(x,sum(index[[1]][1:x]))) 
	owners=names[sapply(indexValues,
					function(x) sums[,sums[2,]>=x][1])]

	indexesUsed=unique(owners)

	results = rep(NA,length(indexValues))
	for(i in 1:length(indexesUsed)){
		if(debug) print(i)
		ownerIndex=which(names==indexesUsed[i])
		offset=if(ownerIndex==1) 0 else sums[2,ownerIndex-1]
		#if(debug) print(paste("offset",offset))

		indexSet=which(owners == indexesUsed[i])
		results[indexSet]=f(indexesUsed[i],indexValues[indexSet]-offset)
	}
	#if(debug) print("results:")
	#if(debug) print(results)
	return(results)
}

getIndexOwners <- function(indexValues)
{
	index=read.table(ChemIndex,sep="\t",row.names=1)
	sums=sapply(1:length(index[[1]]),function(x) c(x,sum(index[[1]][1:x]))) 
	list(names=rownames(index),sums=sums,
		  owners=rownames(index)[sapply(indexValues,function(x) sums[,sums[2,]>=x][1])])
}
getNames <- function(indexValues)
{
	index=read.table(ChemIndex,sep="\t",row.names=1)
	sums=sapply(1:length(index[[1]]),function(x) c(x,sum(index[[1]][1:x]))) 
}
eiInit <- function(compoundDb,measure=atompairMeasure)
{
	if(!file.exists(DataDir))
		dir.create(DataDir)

	if(!file.exists(ChemDb)){
		numCompounds = measure$dbBuilder(toSdfFile(compoundDb),ChemDb)
		writeIddb(1:numCompounds,Main)
		addToChemIndex(ChemDb,numCompounds)
	}
}

eiMakeDb <- function(r,d,measure=atompairMeasure,
				dir=".",refIddb=NA,numSamples=cdbSize()*0.1,
				cl=makeCluster(1,type="SOCK"))
{
	workDir=file.path(dir,paste("run",r,d,sep="-"))
	if(!file.exists(workDir))
		dir.create(workDir)
	matrixFile = file.path(workDir,sprintf("matrix.%d-%d",r,d))

	if(file.exists(matrixFile))
		stop(paste("found existing",matrixFile),"stopping")

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
		message(paste("re-using coordfile",coordFile))
		as.matrix(read.table(coordFile))
	}else{
		#compute pairwise distances for all references
		if(! file.exists(selfDistFile)){
			message("generating selfDistFile")
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

	binaryCoord(embeddedFile,matrixFile,d)
	binaryCoord(embeddedQueryFile,
		file.path(workDir,sprintf("matrix.query.%d-%d",r,d)),d)

	return(file.path(workDir,sprintf("matrix.%d-%d",r,d)))
}
eiQuery <- function(r,d,refIddb,queries,
		dir=".",measure=atompairMeasure,
		K=defK, W = defW, M=defM,L=defL ,T=defT )
{
		tmpDir=tempdir()
		workDir=file.path(dir,paste("run",r,d,sep="-"))
		queryDb = file.path(tmpDir,"queries.db")
		refDb = refDb(refIddb,measure)

		queryFile=toSdfFile(queries)
#		querySdfSet=toSdfSet(queries)
		#reformat query file
		numQueries = measure$dbBuilder(queryFile,queryDb)

		#read names
		queryNames = readLines(paste(queryDb,"names",sep="."))

		#embed queries in search space
		embeddedQueries = embedFromRefs(r,d,refIddb,measure,queryDb,db2=refDb)

		#search for nearby compounds
		if(debug) print(embeddedQueries)
		matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
		hits = search(embeddedQueries,matrixFile,
							queryDb,measure,K=K,W=W,M=M,L=L,T=T)
		if(debug) print(hits)

		targetIds=unlist(lapply(1:length(hits),function(x) hits[[x]][,1]))
		targetIds=targetIds[targetIds!=-1]
		targetNames=as.matrix(getNames(targetIds))
		rownames(targetNames)=targetIds
		print("targetIds")
		print(targetIds)
		print(targetNames)


#		chemSdfSet = toSdfSet(ChemDb)
#		results =lapply(1:numQueries,function(queryIndex){
#				chemSdfSet[hits[[queryIndex]][,1]]
#			}
#		)

		numHits=sum(sapply(hits,function(x) sum(x[,1]==-1)))
		#fetch names for queries and hits and put in a data frame
		results = data.frame(query=rep(NA,numHits),
								  target = rep(NA,numHits),
								  distance=rep(NA,numHits))
		i=1
		lapply(1:numQueries,function(queryIndex)
			lapply(1:K,function(hitIndex){
				if(hits[[queryIndex]][hitIndex,1]!=-1){
					results[i,"query"]<<-queryNames[queryIndex]
					results[i,"target"]<<-targetNames[
							as.character(hits[[queryIndex]][hitIndex,1]),1]
					results[i,"distance"]<<- hits[[queryIndex]][hitIndex,2]
					i<<-i+1
				}
			}))
		if(debug) print(results)
		return(results)
}

eiAdd <- function(r,d,refIddb,additions,dir=".",
		measure=atompairMeasure)
{
		tmpDir=tempdir()
		workDir=file.path(dir,paste("run",r,d,sep="-"))
		additionsDb = nextChemDb(dir)

		#reformat query file
		numAdditions=measure$dbBuilder(toSdfFile(additions),additionsDb)
		#additionNames = readLines(paste(additionsDb,"names",sep="."))
		addToChemIndex(additionsDb,numAdditions)

		#embed queries in search space
		embeddedAdditions= embedFromRefs(r,d,refIddb,
									measure,additionsDb,db2=refDb(refIddb,measure))
		print(dim(embeddedAdditions))
		print(embeddedAdditions)
		embeddedFile <- file.path(workDir,sprintf("coord.%d-%d",r,d))

		#add additions to existing coord and names files
		write.table(t(embeddedAdditions),
				file=embeddedFile, append=TRUE,row.names=F,col.names=F)
		binaryCoord(embeddedFile,file.path(workDir,sprintf("matrix.%d-%d",r,d)),d)
		#write.table(additionNames,append=TRUE,row.names=F,col.names=F,
				#quote=F,file=file.path(dir,paste(ChemDb,"names",sep=".")))
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

	#name additions as "chem-x.db"
	#create index file, "chem.index"
	#each line has name of db file and number of entries in it
	#	order of files is important, must match order of lines in 
	#	coord/matrix files.

#	writeIddb(lshNeighbors[,1],candidatesIddb)
#	measure$dbSubset(ChemDb,candidatesIddb,candidatesDb)
#	d=measure$db2dbDistance(queryDb,db2=candidatesDb)

	print(paste("query index:",queryIndex))
	d=multiFileDistance(queryDb,lshNeighbors[,1],measure)

	print("result distance: ")
	if(debug) print(str(d))
	lshNeighbors[,2]=d #measure$db2dbDistance(queryDb,db2=candidatesDb)
	limit = min(limit,length(lshNeighbors[,2]))
	lshNeighbors[order(lshNeighbors[,2])[1:limit],]
}
getNames <- function(indexes)
	applyOverIndex(indexes,function(filename,indexSet)
		readLines(paste(filename,".names",sep=""))[indexSet]
	)

multiFileDistance <- function(db1,indexes,measure)
{
	applyOverIndex(indexes,function(filename,indexSet){
		tempDir=tempdir()
		writeIddb(indexSet, file.path(tempDir,"temp.iddb"))
		measure$dbSubset(filename,file.path(tempDir,"temp.iddb"),
			file.path(tempDir,"tempdb"))
		measure$db2dbDistance(db1,db2=file.path(tempDir,"tempdb"))
	})
}
#multiFileDistance <- function(db1,indexes,measure)
#{
#	#find out which file each index is from
#	#create a subset from that file
#	#compute dist(db1,subset)
#	#append all distance results together
#	if(debug) print(paste("indexes",paste(indexes,collapse=" ")))
#	indexOwners=getIndexOwners(indexes)
#	print("owners")
#	if(debug) print(indexOwners)
#	indexesUsed=unique(indexOwners$owners)
#	print("used")
#	if(debug) print(indexesUsed)
#	tempDir=tempdir()
#
#	allDists = rep(0,length(indexes))
#	#allDists=unlist(sapply(1:length(indexesUsed),function(i){
#	for(i in 1:length(indexesUsed)){
#		if(debug) print(i)
#		ownerIndex=which(indexOwners$names==indexesUsed[i])
#		offset=if(ownerIndex==1) 0 else indexOwners$sums[2,ownerIndex-1]
#		print(paste("offset",offset))
#
#		indexSet=which(indexOwners$owners == indexesUsed[i])
#		writeIddb(indexes[indexSet]-offset, file.path(tempDir,"temp.iddb"))
#		print(paste("using",indexesUsed[i]))
#		print(indexes[indexSet]-offset)
#
#		measure$dbSubset(indexesUsed[i],file.path(tempDir,"temp.iddb"),
#			file.path(tempDir,"tempdb"))
#		d=measure$db2dbDistance(db1,db2=file.path(tempDir,"tempdb"))
#		if(debug) print(d)
#		allDists[indexSet]=d
#	}
#	#}))
#	#if(length(allDists) != length(indexes))
#	#	stop("failed to compute distance for each candidate")
#	allDists
#
#}
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

nextChemDb <- function(dir=".")
{
	regex=paste(ChemPrefix,"-\\d+\\.db$",sep="")
	sub.regex=paste(ChemPrefix,"-(\\d+)\\.db$",sep="")
	values = sub(paste(".*",sub.regex,sep=""),"\\1",
					grep(regex,
						  dir(file.path(dir,DataDir),pattern=regex,full.names=TRUE),
						  value=TRUE))
	nextIndex = if(length(values)==0) 1 else max(as.numeric(values))+1

	return(file.path(dir,DataDir,paste(ChemPrefix,"-",nextIndex,".db",sep="")))
}
#if input is a file, do nothing.
#if input is an SDFset, write it to a file and return the name
#else cause an error
toSdfFile <- function(input)
{
	tempName=input
	if(class(input)[1] == "SDFset"){
		tempName = tempfile()
		origData=datablock(input)
		datablock(input) = rep("",length(input))
		write.SDF(input,tempName)
		datablock(input)=origData
	}
	if(!is.character(tempName)){
		stop("unknown type of input while trying to get an sdf file")
	}

	return(tempName)
}

#like toSdfFile, but return an SDFset
toSdfSet <- function(input)
{
	sdfset = input
	if(is.character(input)){
		sdfset=read.SDFset(input)
	}
	if(class(sdfset)!="SDFset")
		stop("unknown type of input while trying to get an SDFset")
	return(sdfset)
}

