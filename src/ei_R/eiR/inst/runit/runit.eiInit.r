
library(eiR)
library(snow)

options(warn=2)
r<- 50
d<- 40
#N<- 122
N<- 100
j=1
runDir<-paste("run",r,d,sep="-")

test.aa.eiInit <- function() {
	DEACTIVATED("slow")
   data(sdfsample)
   compoundIds = eiInit(sdfsample)
   checkTrue(file.exists(file.path("data","chem.db")))
   checkTrue(file.exists(file.path("data","main.iddb")))
   i <- readLines(file.path("data","main.iddb"))
   checkEquals(length(i),N)
	checkEquals(length(compoundIds),N)
	sdfFromDb = getCompounds(initDb(file.path("data","chem.db")),compoundIds)
	checkEquals(length(sdfFromDb),N)
}

testRefs <- function(){
	200+c(1,2,5,8,9,10,11,17,18,19,20,23,24,25,26,29,31,33,34,36,38,43,45,46,47,48,49,51,53,66,67,70,71,72,73,74,75,77,78,79,80,81,82,83,87,88,89,91,99,100)
}
test.ba.eiMakeDb <- function() {

	#DEACTIVATED("slow")
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

	print("by file name")
	eiR:::writeIddb((1:r)+200,"reference_file.cdb")
   eiMakeDb("reference_file.cdb",d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile=""))
   runChecks()
	unlink(runDir,recursive=TRUE)

	print("by number")
   eiMakeDb(r,d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile=""))
   runChecks()
	unlink(runDir,recursive=TRUE)

	print("by vector")
   eiMakeDb(testRefs(),d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile=""))
   runChecks()

	
}
test.ca.eiQuery <- function(){

	DEACTIVATED("slow")
   data(sdfsample)
   refIddb = findRefIddb(runDir)
   results = eiQuery(r,d,refIddb,sdfsample[1:2],K=15)
   checkTrue(length(results$distance) != 0)
   checkTrue(all(results$distance <= 1))
   checkEquals(results$distance[16],0)


	results=eiQuery(r,d,refIddb,203:204,format="compound_id",K=15)
   checkEquals(results$distance[1],0)

	results=eiQuery(r,d,refIddb,c("650002","650003"), format="name",K=15)
   checkEquals(results$distance[1],0)
   #checkEquals(results$distance[9],0) # not reliable

}

test.da.eiPerformanceTest <- function() {
	DEACTIVATED("slow")
   r = eiPerformanceTest(r,d,K=22)
   checkMatrix("chemical-search.results$",20, N,"data")
   checkMatrix("eucsearch.50-40",20,N)
   checkMatrix("^indexed$",20,22)
   checkMatrix("indexed.performance",20,1)
}
test.ea.eiAdd<- function(){

	DEACTIVATED("slow")
   data(example_compounds)
   cat(paste(paste(example_compounds,collapse="\n"),"\n",sep=""),file="example_compounds.sdf")
   options(warn=-1)
   examples=read.SDFset("example_compounds.sdf")
   options(warn=2)
   eiAdd(r,d,findRefIddb(runDir),examples[1:2])

   results = eiQuery(r,d,findRefIddb(runDir),examples[1:2])
   print(results)
   checkEquals(results$distance[1],0)

   eiAdd(r,d,findRefIddb(runDir),examples[4:8])
   results = eiQuery(r,d,findRefIddb(runDir),examples[4])
   checkEquals(results$distance[1],0)
   print(results)
}
test.fa.eiCluster <- function(){
#	DEACTIVATED("off")
	numNbrs=5
	minNbrs=2


	clustering=eiCluster(r,d,K=numNbrs,minNbrs=minNbrs,cutoff=0.5)
	checkTrue(length(clustering) >= N) #eiAdd will add some stuff

	conn = initDb("data/chem.db")
	compoundIds=names(clustering)
	compoundNames=getCompoundNames(conn,compoundIds)
	names(clustering)=compoundNames
	print(sort(clustering))
	print(clusterSizes(clustering))
}
test.fn.cluster_comparison <- function(){

	#DEACTIVATED("off")
	numNbrs=10
	minNbrs=2
	#numNbrs=20
	#minNbrs=17
	fast=TRUE
	cutoff=0.5

	dir="."
	#dir="/home/khoran/runs/drug_bank_1000"
	#r=300
	#d=100

	#clnnm=eiCluster(r,d,K=numNbrs,minNbrs=minNbrs)
	#clustering = jarvisPatrick(clnnm,j=numNbrs,k=minNbrs,fast=fast)
	#compoundIds=eiR:::readIddb(eiR:::Main)

	clustering=eiCluster(r,d,K=numNbrs,minNbrs=minNbrs,dir=dir,cutoff=cutoff)
								#L = 60, T = 50,  M = 9,  W = 101.019)
								#L = 60,T = 50 ,M = 4,  W = 63.8275)
#	checkTrue(length(clustering) >= N) #eiAdd will add some stuff
#
	conn = initDb(file.path(dir,"data/chem.db"))
	compoundIds=names(clustering)
	compoundNames=getCompoundNames(conn,compoundIds)
	names(clustering)=compoundNames
	print(sort(clustering))
#

	#### non lsh clustering

###### gen cl2
	preProcess = eiR:::getTransform("ap")$toObject
	aps=as(preProcess(eiR:::getDescriptors(conn,"ap",compoundIds)),"APset")
	#cid(aps)=compoundNames
	#cid(aps)=as.character(compoundIds)
	cid(aps)=as.character(1:length(compoundIds))

	cl2nnm = jarvisPatrick(aps,j=numNbrs,k=minNbrs,type="matrix",cutoff=cutoff)
	d=dim(cl2nnm)
	cl2nnm=as.numeric(cl2nnm)
	dim(cl2nnm)=d
	rownames(cl2nnm)=cid(aps)

	#print(sapply(seq(dim(cl2nnm)[1]),function(i) 
	#		 sapply(seq(along=cl2nnm[i,]),function(j)
	#				  if(is.na(cl2nnm[i,j])) -1 else cmp.similarity(aps[i],aps[cl2nnm[i,j]]))))
	print(tail(cl2nnm))

	#cl2 = jarvisPatrick(cl2nnm,j=numNbrs,k=minNbrs,fast=fast)

	cl2 = jarvisPatrick_c(cl2nnm,minNbrs,fast=fast)
	#print(cl2)
	#names(cl2)=cid(aps)
	#names(cl2)=as.character(compoundIds)
	names(cl2)=compoundNames
	#print(cl2)
##############


	write.table(clustering,file="sample_lsh_cutoff.clstr",quote=F,sep="\t",col.names=F)
	write.table(cl2,file="sample_true.cutoff.clstr",quote=F,sep="\t",col.names=F)
	print(sort(clustering))
	print(sort(cl2))


	source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/clusterIndex.R")
	print(clusterSizes(clustering))
	print(clusterSizes(cl2))
	ci <- cindex(clV1=clustering, clV2=cl2, minSZ=0, method="jaccard") 
	rand <- cindex(clV1=clustering, clV2=cl2, minSZ=0,
						method="rand")
	rand=rand$Rand_Index
	arand <- cindex(clV1=clustering, clV2=cl2, minSZ=0,
						 method="arand")
	arand=arand$Adjusted_Rand_Index
	print(paste("cluster similarity:",ci$Jaccard_Index,rand,arand))

}

test.fo.nnm_test  <- function(){
	DEACTIVATED("off")
	numNbrs=10
	minNbrs=2
	fast=TRUE
	

	dir="/home/khoran/runs/drug_bank_1000"
	r=300
	d=100
	#dir="."
	cutoff = 0.5

	compoundIds=eiR:::readIddb(file.path(dir,eiR:::Main))

	#compoundIds = compoundIds[1:100]

	print("computeing lsh results")
	lshnnm=lshNnm(r=r,d=d,K=numNbrs,minNbrs=minNbrs,dir=dir,cutoff=cutoff)

#	lshnnm = queriedNnm(compoundIds,r,d,numNbrs,dir)
	#print(lshnnm)

	print("lsh NNM:")
	print(tail(lshnnm,n=30))
	#print(lshnnm)

	print("computing true results")
	truennm = trueNnm(compoundIds,numNbrs,minNbrs,dir=dir,cutoff=cutoff)

	print("true NNM: ")
	print(tail(truennm,n=30))
	#print(truennm)

	#compute precsion/recall
	results=sapply(1:nrow(lshnnm),function(i){
			intersection =  intersect(truennm[i,],lshnnm[i,])
			t=length( intersection[!is.na(intersection)])
			#numFetched = sum(lshnnm[i,] != -1)
			numFetched = sum(!is.na(lshnnm[i,])) 
			numTrue = sum(!is.na(truennm[i,])) 
			#print(paste(t,numFetched,numTrue))
			p=t/numFetched
			r=t/numTrue
			#print(paste(t,numFetched,numTrue,p,r))
			c(p,r,  2* (p*r)/(p+r)  )
				 })
	print(results)
	print(paste(sum(results[1,])/length(results[1,]), sum(results[2,])/length(results[2,])
					,  sum(results[3,])/length(results[3,])))

}
lshNnm <- function(...){

	nnm=eiCluster(...,type="matrix")
	nnm
}
queriedNnm <- function(compoundIds,r,d,numNbrs,dir){

		refIddb=findRefIddb(file.path(dir,paste("run",r,d,sep="-")))
		results = eiQuery(r,d,refIddb,compoundIds,format="compound_id",K=numNbrs,dir=dir)
		cidToPosition = 1:length(compoundIds)
		names(cidToPosition) = as.character(compoundIds)
		print(results)

		t(sapply(seq(along=compoundIds),function(i){
					#print(paste("i:",i,"cid: ",compoundIds[i],paste(which(results$query==compoundIds[i]),collapse=",")))
					cidToPosition[as.character(results[results$query==compoundIds[i],"target_ids"])]
				 }))


}
trueNnm <- function(compoundIds,numNbrs,minNbrs,dir,cutoff=NA){

	conn = initDb(file.path(dir,"data/chem.db"))
	preProcess = eiR:::getTransform("ap")$toObject
	aps=as(preProcess(eiR:::getDescriptors(conn,"ap",compoundIds)),"APset")
	#cid(aps)=compoundNames
	#cid(aps)=as.character(compoundIds)
	cid(aps)=as.character(1:length(compoundIds))



	nnm = jarvisPatrick(aps,type="matrix",cutoff,j=numNbrs,k=minNbrs)
	d=dim(nnm)
	nnm=as.numeric(nnm)
	dim(nnm)=d
	rownames(nnm)=cid(aps)

	nnm
}

clusterSizes <- function(clustering) {
	sizes=Reduce(rbind,lapply(unique(clustering),function(cid)
							  cbind(cid=cid,size=sum(clustering==cid))))
	#sizes[order(sizes[,1]),]
	sizes[sizes[,2]>1,]
}



test.aaaaa.cleanup<- function(){
   #junk <- c("data","example_compounds.sdf","example_queries.sdf","run-50-40")
   junk <- c("example_compounds.sdf","example_queries.sdf","run-50-40")
   unlink(junk,recursive=T)
}
findRefIddb <- function(runDir){
   matches<-dir(runDir,pattern=".cdb$",full.names=T)
   checkEquals(length(matches),1)
   matches[1]
}
checkMatrix <- function(pattern,x,y,dir=runDir){
#	print(paste("searching for ",pattern))
   matches<-dir(dir,pattern=pattern,full.names=T)
#	print(matches)
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


