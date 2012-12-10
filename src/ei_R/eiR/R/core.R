
library(snow)

DataDir = "data"
TestQueries = file.path(DataDir,"test_query.iddb")
TestQueryResults=file.path(DataDir,"chemical-search.results")
ChemPrefix="chem"
ChemDb = file.path(DataDir,paste(ChemPrefix,".db",sep=""))
ChemIndex = file.path(DataDir,paste(ChemPrefix,".index",sep=""))
Main = file.path(DataDir,"main.iddb")


debug=TRUE
#debug=FALSE

# Notes
#  Need function to produce descriptors from sdf or smile
#  Need function to compute distances between descriptors

cdbSize <- function(dir=".") {
#	getSegmentSize(file.path(dir,ChemDb),dir)
	#TODO: make this more efficient
	length(readIddb(file.path(dir,Main)))
}
embedCoord <- function(s,len,coords) 
	.Call("embedCoord",s,as.integer(len),as.double(coords))

embedCoordTest <- function(r,d,refCoords,coords) 
	.Call("embedCoordTest",as.integer(r),as.integer(d),as.double(refCoords),as.double(coords))

# requires one query per column, not per row
lshsearch <- function(queries,matrixFile,
	W=NA,H=NA,M=NA,L=NA,K=NA,T=NA,R=NA) 
{
	if(!file.exists(matrixFile)) stop(paste("could not find matrix file:",matrixFile))
	.Call("lshsearch",queries,as.character(matrixFile),
		as.double(W),as.integer(H),as.integer(M),as.integer(L),
		as.integer(K),as.integer(T), as.double(R))
}


# functions needed for sql backend:
# distance

# descriptorStr  raw format (sdf,smile) -> descriptor object -> string
# str2Descriptor  string -> descriptor object
# also need descrptor type, ie, "ap" "fpap", etc.

# and compound format, ei, "SDF", "SMILE", etc.
# X optionally: compound -> string and string -> compound

eiInit <- function(compoundDb,dir=".",format="sdf",descriptorType="ap",append=FALSE)
{
	if(!file.exists(file.path(dir,DataDir)))
		dir.create(file.path(dir,DataDir))

	descriptorFunction = function(set)
		data.frame(descriptor=getTransform(buildType(format,descriptorType))$toString(set),
					  descriptor_type=descriptorType)
	

	conn = initDb(file.path(dir,ChemDb))
	if(tolower(format) == "sdf"){
		compoundIds = loadSdf(conn,compoundDb, descriptors=descriptorFunction)
	}else if(tolower(format) == "smile"){
		compoundIds = loadSmile(conn,compoundDb,descriptors=descriptorFunction)
	}else{
		stop(paste("unknown input format:",format," supported formats: SDF, SMILE"))
	}
	print(paste(length(compoundIds)," loaded by eiInit"))

	writeIddb(compoundIds,file.path(dir,Main),append=append)
	compoundIds
}

eiMakeDb <- function(refs,d,distance=apDistance,
				dir=".",numSamples=cdbSize(dir)*0.1,
				cl=makeCluster(1,type="SOCK"))
{
	workDir=NA
	createWorkDir <- function(r){
		workDir<<-file.path(dir,paste("run",r,d,sep="-"))
		if(!file.exists(workDir))
			dir.create(workDir)
	}

	if(is.character(refs)){ #assume its a filename
		refIds=readIddb(refs)
		r=length(refIds)
		createWorkDir(r)
		refIddb=file.path(workDir,basename(refs))
		file.copy(refs,workDir,overwrite=TRUE)
	}else if(is.numeric(refs)){
		if(length(refs)==0){ #assume its the number of refs to use
			stop(paste("variable refs must be posative, found ",refs))
		}else if(length(refs)==1){ #assume its the number of refs to use
			r=refs
			createWorkDir(r)
			refIddb=genRefName(workDir)
			refIds=genRefs(r,refIddb,dir)
		}else{ #refs is a vector of compound indexes to use a referances
			refIds=refs
			r=length(refIds)
			createWorkDir(r)
			refIddb=genRefName(workDir)
			writeIddb(refIds,refIddb)
		}
	}else{
		stop(paste("don't know how to handle refs:",str(refs)))
	}

	matrixFile = file.path(workDir,sprintf("matrix.%d-%d",r,d))

	if(file.exists(matrixFile))
		stop(paste("found existing",matrixFile),"stopping")


	queryIds=genTestQueryIds(numSamples,dir,refIds)


	selfDistFile <- paste(refIddb,"distmat",sep=".")
	coordFile <- paste(selfDistFile,"coord",sep=".")
	ref2AllDistFile <- paste(refIddb,"distances",sep=".")
	embeddedFile <- file.path(workDir,sprintf("coord.%d-%d",r,d))
	embeddedQueryFile <- file.path(workDir,sprintf("coord.query.%d-%d",r,d))

allIds = readIddb(file.path(dir,Main))
IddbVsIddbDist(file.path(dir,ChemDb),allIds,allIds,distance,file="all2allR.dist")

	#embed references in d dimensional space 
	coords <- if(file.exists(coordFile)){
		message(paste("re-using coordfile",coordFile))
		as.matrix(read.table(coordFile))
	}else{
		#compute pairwise distances for all references
		if(! file.exists(selfDistFile)){
			message("generating selfDistFile")
			#sql_db2dbDistance(file.path(dir,ChemDb),iddb1=refIds,iddb2=refIds,file=selfDistFile)
			IddbVsIddbDist(file.path(dir,ChemDb),refIds,refIds,distance,file=selfDistFile)
		}
		selfDist<-read.table(selfDistFile)

		#use MDS and pairwise distances to embed references
		coords <- cmdscale(selfDist,d)
		write.table(coords,file=coordFile,row.names=F,col.names=F)
		coords
	}
	#compute dist between refs and all compounds
	if(!file.exists(ref2AllDistFile))
		#sql_db2dbDistance(file.path(dir,ChemDb),iddb1=readIddb(file.path(dir,Main)),iddb2=refIds,file=ref2AllDistFile)
		IddbVsIddbDist(file.path(dir,ChemDb),readIddb(file.path(dir,Main)),refIds,distance,file=ref2AllDistFile)
	
	#each job needs: R, D, coords, a chunk of distance data
	solver <- getSolver(r,d,coords)	
	distConn <- file(ref2AllDistFile,"r")

	numJobs=length(cl)
	jobSize = as.integer(cdbSize(dir) / numJobs + 1) #make last job short

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
eiQuery <- function(r,d,refIddb,queries,format="sdf",
		dir=".",descriptorType="ap",distance=apDistance,
		K=200, W = 1.39564, M=19,L=10,T=30)
{
		tmpDir=tempdir()
		workDir=file.path(dir,paste("run",r,d,sep="-"))
		descriptorInfo = getTransform(buildType(format,descriptorType)
												)$toObject(queries)
		queryDescriptors = descriptorInfo$descriptors
		#print("query descriptors")
		#print(ap(queryDescriptors))

		refIds = readIddb(refIddb)
		#queryDb = file.path(tmpDir,"queries.db")
		#refDb = refDb(refIddb,measure,dir=dir)

		#queryFile=toSdfFile(queries)
		#reformat query file
		#numQueries = measure$dbBuilder(queryFile,queryDb)
		numQueries = length(queryDescriptors)

		#read names
		#queryNames = readLines(paste(queryDb,"names",sep="."))
		queryNames = descriptorInfo$names
		#print("queryNames:")
		#print(queryNames)

		#embed queries in search space
		qd=IddbVsGivenDist(file.path(dir,ChemDb),
															refIds,queryDescriptors,distance)
		print("distances")
		print(t(qd))
		#embeddedQueries = embedFromRefs(r,d,refIddb,measure,queryDb,db2=refDb)
		embeddedQueries = embedFromRefs(r,d,refIddb, t(qd))

		#search for nearby compounds
		if(debug) print(embeddedQueries)
		matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
		hits = search(embeddedQueries,matrixFile,
							queryDescriptors,distance,dir,K=K,W=W,M=M,L=L,T=T)
		if(debug) print("hits")
		if(debug) print(hits)

		targetIds=unlist(lapply(1:length(hits),function(x) hits[[x]][,1]))
		targetIds=targetIds[targetIds!=-1]
		targetNames=as.matrix(getNames(targetIds,dir))
		rownames(targetNames)=targetIds
		#print(paste(targetIds,targetNames))



		numHits=sum(sapply(hits,function(x) sum(x[,1]!=-1)))
		#print(paste("numHits:",numHits))
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
		if(debug) print("results:")
		if(debug) print(results)
		return(results)
}

eiAdd <- function(r,d,refIddb,additions,dir=".",format="SDF",
		descriptorType="ap",distance=apDistance)
{
		tmpDir=tempdir()
		workDir=file.path(dir,paste("run",r,d,sep="-"))
		#additionsDb = nextChemDb(dir)

		# add additions to database
		conn = initDb(file.path(dir,ChemDb))
		print("adding")
		compoundIds = eiInit(additions,dir,format,descriptorType,append=TRUE)
		print("getting descriptors")
		additionDescriptors=getDescriptors(conn,descriptorType,compoundIds)
		numAdditions = length(compoundIds)
		print(numAdditions)

		#reformat query file
		#numAdditions=measure$dbBuilder(toSdfFile(additions),additionsDb)
		#addToChemIndex(additionsDb,numAdditions,dir)

		refIds = readIddb(refIddb)
		#embed queries in search space

		embeddedAdditions= embedFromRefs(r,d,refIddb,
									t(IddbVsGivenDist(file.path(dir,ChemDb),refIds,additionDescriptors,distance)))
		if(debug) print(dim(embeddedAdditions))
		if(debug) print(embeddedAdditions)
		embeddedFile <- file.path(workDir,sprintf("coord.%d-%d",r,d))

		#add additions to existing coord and names files
		write.table(t(embeddedAdditions),
				file=embeddedFile, append=TRUE,row.names=F,col.names=F)
		binaryCoord(embeddedFile,file.path(workDir,sprintf("matrix.%d-%d",r,d)),d)
}
#expects one query per column
search <- function(queries,matrixFile,queryDescriptors,distance,K,dir,...)
{
		neighbors = lshsearch(queries,matrixFile,K=2*K,...)
		#print(paste("got ",paste(dim(neighbors),callapse=","),"neighbors back from lshsearch"))
		print("neighbors:")
		print(neighbors)

		#compute distance between each query and its candidates	
		Map(function(i) refine(neighbors[i,,],queryDescriptors[i],K,distance,dir),
			1:(dim(queries)[2]))
}
#fetch coords from refIddb.distmat.coord and call embed
embedFromRefs <- function(r,d,refIddb,query2RefDists)
{
		coordFile=paste(refIddb,"distmat","coord",sep=".")
		coords = as.matrix(read.table(coordFile))
		embed(r,d,coords,query2RefDists)
}
#take referenace coords, and distance matrix
#return the emedding of each row in the matrix
embed <- function(r,d,coords, query2RefDists)
{
		solver = getSolver(r,d,coords)
		embeddedQueries = apply(query2RefDists,c(1),
			function(x) embedCoord(solver,d,x))
}
refine <- function(lshNeighbors,queryDescriptors,limit,distance,dir)
{
	tmpDir=tempdir()
	#queryIddb=file.path(tmpDir,"query.iddb")
	#queryDb=file.path(tmpDir,"query.db")

	#queryIddb=file.path(tmpDir,"query.iddb")
	#writeIddb(c(queryIndex),queryIddb)
	#measure$dbSubset(queriesCdb,queryIddb,queryDb)

	#name additions as "chem-x.db"
	#create index file, "chem.index"
	#each line has name of db file and number of entries in it
	#	order of files is important, must match order of lines in 
	#	coord/matrix files.


	#if(debug) print(paste("query index:",queryIndex))
	#d=multiFileDistance(queryDb,lshNeighbors[,1],measure,dir)
	d = t(IddbVsGivenDist(file.path(dir,ChemDb),lshNeighbors[,1],queryDescriptors,distance))

	#if(debug) print("result distance: ")
	#if(debug) print(str(d))
	lshNeighbors[,2]=d 
	limit = min(limit,length(lshNeighbors[,2]))
	#print(paste("num dists:",length(lshNeighbors[,2]),
			#"limit:",limit))
	lshNeighbors[order(lshNeighbors[,2])[1:limit],]
}
getNames <- function(indexes,dir)
	getCompoundNames(initDb(file.path(dir,ChemDb)),indexes)

writeIddb <- function(data, file,append=FALSE)
		write.table(sort(data),file,quote=FALSE,append=append,col.names=FALSE,row.names=FALSE)
readIddb <- function(file) as.numeric(readLines(file))
readNames <- function(file) as.numeric(readLines(file))


# randomly select n reference compounds. Also sample and stash away
# numSamples query compounds that are not references for later
# testing
genTestQueryIds <- function(numSamples,dir,refIds=c())
{
	testQueryFile <-file.path(dir,TestQueries)
	mainIds <- readIddb(file.path(dir,Main))
	set=setdiff(mainIds,refIds)
	if(numSamples < 0 || numSamples > length(set)) stop(paste("trying to take more samples than there are compounds available",numSamples,length(set)))
	queryIds <- sort(sample(set,numSamples))
	writeIddb(queryIds,testQueryFile)
	queryIds
}
genRefs <- function(n,refFile,dir,queryIds=c())
{
	mainIds <- readIddb(file.path(dir,Main))
	set=setdiff(mainIds,queryIds)
	if(n < 0 || n > length(set)) stop(paste("found more refereneces than compound candidates",n,length(set)))
	refIds = sort(sample(set,n))
	writeIddb(refIds,refFile)
	refIds
}
genRefName <- function(workDir)
	file.path(workDir,
				 paste(paste(sample(c(0:9,letters),32,replace=TRUE),
								 collapse=""),
						 "cdb",sep="."))
#refDb <- function(refIddb,measure,refDb=paste(refIddb,"db",sep="."),dir)
#{
#	if(!file.exists(refDb))
#		measure$dbSubset(file.path(dir,ChemDb),refIddb,refDb)
#	return(refDb)
#}
genTestQueryResults <- function(distance,dir)
{
	if(file.exists(file.path(dir,TestQueryResults)))
		return()

	out=file(file.path(dir,TestQueryResults),"w")
	#d=sql_db2dbDistance(file.path(dir,ChemDb),
	#	iddb1=readIddb(file.path(dir,TestQueries)),iddb2=readIddb(file.path(dir,Main)))
	d=IddbVsIddbDist(file.path(dir,ChemDb),
		readIddb(file.path(dir,TestQueries)),
		readIddb(file.path(dir,Main)),distance)
	if(debug) print(paste("dim(d): ",dim(d)))
	maxLength=min(dim(d)[2],50000)
	for(i in 1:(dim(d)[1]))
		cat(paste(
				paste(1:dim(d)[2],d[i,],sep=":")[order(d[i,])[1:maxLength]],
				collapse=" "),"\n",file=out)	
	close(out)
}
eiPerformanceTest <- function(r,d,distance=apDistance,descriptorType="ap",
	dir=".",K=200, W = 1.39564, M=19,L=10,T=30)
{
	workDir=file.path(dir,paste("run",r,d,sep="-"))
	eucsearch=file.path(workDir,sprintf("eucsearch.%s-%s",r,d))
	genTestQueryResults(distance,dir)
	eucsearch2file(file.path(workDir,sprintf("matrix.%s-%s",r,d)),
				 file.path(workDir,sprintf("matrix.query.%s-%s",r,d)),
				 50000,eucsearch)

	#evaluator TestQueryResuts eucsearch-r-d recall
	evaluator(file.path(dir,TestQueryResults),eucsearch,
		file.path(workDir,"recall"))

	matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
	coordQueryFile =file.path(workDir,sprintf("coord.query.%d-%d",r,d))

	conn = initDb(file.path(dir,ChemDb))
	testQueryDescriptors=getDescriptors(conn,descriptorType,readIddb(file.path(dir,TestQueries)) )
	embeddedTestQueries = t(as.matrix(read.table(coordQueryFile)))
	hits = search(embeddedTestQueries,matrixFile,
						testQueryDescriptors,distance,dir,K=K,W=W,M=M,L=L,T=T)
	indexed=file.path(workDir,"indexed")
	out=file(indexed,"w")
	#if(debug) print(hits)
	for(x in hits)
		cat(paste(x[,1],x[,2],sep=":",collapse=" "),"\n",file=out)
	close(out)

	#indexed_evalutator TestQueryResults indexed indexed.performance
	write.table(compareSearch(file.path(dir,TestQueryResults),indexed),
			file=file.path(workDir,"indexed.performance"),
			row.names=F,col.names=F,quote=F)
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

#distance functions
# all vs descriptors     db,descriptors
# all vs subset          db,iddb
# subset vs descriptors  db,iddb,descriptors
# subset vs subset       db,iddb1,iddb2

desc2descDist <- function(desc1,desc2,dist)
	as.matrix(sapply(desc2,function(x) sapply(desc1,function(y) dist(x,y))))


#dbVsGivenDist<- function(db,descriptors,dist,descriptorType="ap",file=NA){
#	conn = initDb(db)
#
#	preProcess = getTransform(descriptorType)$toObject
#	descriptors=preProcess(descriptors)
#	process = function(record){
#		bufferResultSet(dbSendQuery(conn,"SELECT descriptor FROM descriptors JOIN descriptor_types 
#												USING(descriptor_type_id) WHERE descriptor_type='ap' ORDER BY compound_id"),
#			function(data){
#				data=preProcess(data[1][[1]])
#				record(desc2descDist(data,descriptors,dist))
#		},closeRS=TRUE)
#	}
#	output(file,cdbSize,length(descriptors),process)
#}
dbVsIddbDist<- function(db,iddb,dist,descriptorType="ap",file=NA){
	dbVsGivenDist(db,getDescriptors(initDb(db),descriptorType,iddb),dist,descriptorType,file)
}
IddbVsGivenDist<- function(db,iddb,descriptors,dist,descriptorType="ap",file=NA){

	conn = initDb(db)
	preProcess = getTransform(descriptorType)$toObject
	descriptors=preProcess(descriptors)

	process = function(record){
		batchByIndex(iddb,function(ids){
			outerDesc=getDescriptors(conn,descriptorType,ids)
			#print("innerDesc")
			#print(ids)
			#print(innerDesc)
			outerDesc = preProcess(outerDesc)
			record(desc2descDist(outerDesc,descriptors,dist))
		})
	}
	output(file,length(iddb),length(descriptors),process)
}
IddbVsIddbDist<- function(db,iddb1,iddb2,dist,descriptorType="ap",file=NA){
	conn = initDb(db)

	#print(paste("iddb2:",paste(iddb2,collapse=",")))
	preProcess = getTransform(descriptorType)$toObject
	descriptors = preProcess(getDescriptors(conn,descriptorType,iddb2))
	#print("descriptors")
	#print(str(descriptors))
	process = function(record){
		batchByIndex(iddb1,function(ids){
			#print(paste("processing",paste(ids,collapse=",")))
			#outerDesc = preProcess(dbGetQuery(conn,selectDescriptors(descriptorType,ids))$descriptor)
			outerDesc = preProcess(getDescriptors(conn,descriptorType,ids))
			record(desc2descDist(outerDesc,descriptors,dist))
		})
	}
	output(file,length(iddb1),length(iddb2),process)
}

#choose whether to output a file or a matrix
output <- function(filename,nrows,ncols,process)
	if(!is.na(filename)){ #return result as matrix
		toFile(filename,process)	
	}else{ #write result to file
		toMatrix(nrows,ncols,process)
	}

#send data produced by body to a file
toFile <- function(filename,body){
	f = file(filename,"w")
	body(function(data) write.table(data,file=f,quote=F,sep="\t",row.names=F,col.names=F))
	close(f)
}

#send data produced by body to a matrix
toMatrix <- function(nrows,ncols,body){

	allDists = matrix(NA,nrows,ncols)
	rowCount = 1

	body(function(data){
#		if(debug) print(paste("recording data. rowCount: ",rowCount,", dim:",paste(dim(data),collapse=",")))
#		if(debug) print(dim(allDists))
		allDists[rowCount:(rowCount+dim(data)[1]-1),] <<- data
		rowCount <<- rowCount + dim(data)[1]
	})

	allDists
}

selectDescriptors <- function(type,ids){
	q=paste("SELECT compound_id, descriptor FROM descriptors JOIN descriptor_types USING(descriptor_type_id) WHERE ",
				" descriptor_type='",type,"' AND compound_id IN (", paste(ids,collapse=","),") ORDER
				BY compound_id",sep="")
	#print(q)
	q
}
getDescriptors <- function(conn,type,idList){
	data = selectInBatches(conn,idList,function(ids) selectDescriptors(type,ids))
	n=data$descriptor
	names(n)=data$compound_id
	ordered=n[as.character(idList)]
	#write.table(n,file="descriptors.out")
	ordered
}
	

