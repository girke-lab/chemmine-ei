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

#expects two binary vectors
fingerprintDistance <- function(v1,v2)
{

}

rawFingerprintMeasure = list(
	dbBuilder = function(input,output){
		file.copy(input,output,overwrite=TRUE)
	},
	#we assume one FP per line
	dbSubset = function(db,iddb,output){
		ids=readIddb(iddb)
		f=file(otuput,"w")
		dbh=file(db,"r")
		lastId=0
		for(id in ids){
			line = scan(dbh,skip=id-lastId-1,nlines=1)
			lastId=id
			cat(line,f)
		}
	},
	db2dbDistance = function(db,db2=NA,iddb1=NA,iddb2=NA,file=NA){

	}

)


apDistance <- function(d1,d2){
	1-cmp.similarity(d1,d2)
}

# functions needed for sql backend:
# distance
# descriptorStr  raw format (sdf,smile) -> descriptor object -> string
# str2Descriptor  string -> descriptor object
# also need descrptor type, ie, "ap" "fpap", etc.
# X optionally: compound -> string and string -> compound


Translations = list()

addTransform <- function(name,toString,toObject){
	Translations[[name]] <<- list(toString=toString,toObject=toObject)
}
getTransform <- function(name){
	Translations[[name]]
}

buildType <- function(format,descriptorType) tolower(paste(format,descriptorType,sep="-"))

addTransform(buildType("sdf","ap"),
	# Any sdf source -> ap string
	toString = function(input,dir=".") 
		getTransform("ap")$toString(
			getTransform(buildType("sdf","ap"))$toObject(input)$descriptors),
	# Any sdf source -> APset
	toObject = function(input,dir="."){

		sdfset=if(is.character(input) && file.exists(input)){
			read.SDFset(input)
		}else if(inherits(input,"SDFset")){
			input
		}else{
			stop(paste("unknown type for 'input', or filename does not exist. type found:",class(input)))
		}
		list(names=sdfid(sdfset),descriptors=sdf2ap(sdfset))
	}
)
addTransform("ap",  
   # APset -> string,
	toString = function(apset,dir="."){
		unlist(lapply(ap(apset), function(x) paste(x,collapse=", ")))
	},
   # string or list -> AP set list
	toObject= function(v,dir="."){ 
		if(inherits(v,"list") || length(v)==0)
			return(v)

		as( if(!inherits(v,"APset")){
				names(v)=as.character(1:length(v));  
				read.AP(v,type="ap",isFile=F)
			} else v,
			"list")  
	}
)

lapply(c("ap"),function(descriptorType)
	addTransform(buildType("compound_id",descriptorType),
		# compound_id -> ap string
		toString = function(ids,dir="."){
			getDescriptors(initDb(file.path(dir,ChemDb)),descriptorType,ids)
		},
		# compound_id -> AP list object
		toObject = function(ids,dir="."){
			descInfo = getTransform(buildType("compound_id",descriptorType))$toString(ids,dir)
			list(names=names(descInfo),
				  descriptors=getTransform(descriptorType)$toObject(descInfo))
		}
	)
)
lapply(c("ap"),function(descriptorType)
	addTransform(buildType("name",descriptorType),
		# name -> ap string
		toString = function(names,dir="."){
			conn=initDb(file.path(dir,ChemDb))
			getDescriptors(conn,descriptorType,findCompoundsByName(conn,names,keepOrder=TRUE))
		},
		# name -> AP list object
		toObject = function(names,dir="."){
			descInfo = getTransform(buildType("name",descriptorType))$toString(names,dir)
			list(names=names,
				  descriptors=getTransform(descriptorType)$toObject(descInfo))
		}
	)
)
