

defaultDistances = list()

#	ap= function(d1,d2) 1-cmp.similarity(d1,d2),
#	fp= function(d1,d2) 1-fpSim(d1,d2)
#)

setDefaultDistance <- function(name,distance)
	defaultDistances[[name]] <<- distance

getDefaultDist <- function(descriptorType){
	d=defaultDistances[[descriptorType]]
	if(is.null(d))
		warn("no default distance found for desciptor type ",descriptorType)
	d
}

setDefaultDistance("ap", function(d1,d2) 1-cmp.similarity(d1,d2) )
setDefaultDistance("fp", function(d1,d2) 1-fpSim(d1,d2) )


Translations = list()

addTransform <- function(name,isDescriptorType=TRUE,toString,toObject){


	Translations[[name]] <<- list(toString=toString,toObject=toObject)

	if( isDescriptorType ){
		# add extra handlers for compound_id and name types
		addTransform(buildType("compound_id",name),isDescriptorType=FALSE,
			# compound_id -> ap string
			toString = function(ids,dir="."){
				getDescriptors(initDb(file.path(dir,ChemDb)),name,ids)
			},
			# compound_id -> AP list object
			toObject = function(ids,dir="."){
				descInfo = getTransform(buildType("compound_id",name))$toString(ids,dir)
				list(names=names(descInfo),
					  descriptors=getTransform(name)$toObject(descInfo))
			}
		)

		addTransform(buildType("name",name),isDescriptorType=FALSE,
			# name -> ap string
			toString = function(names,dir="."){
				conn=initDb(file.path(dir,ChemDb))
				getDescriptors(conn,name,findCompoundsByName(conn,names,keepOrder=TRUE))
			},
			# name -> AP list object
			toObject = function(names,dir="."){
				descInfo = getTransform(buildType("name",name))$toString(names,dir)
				list(names=names,
					  descriptors=getTransform(name)$toObject(descInfo))
			}
		)
	}	



}
getTransform <- function(name){
	t=Translations[[name]]
	if(is.null(t))
		stop("transform ",name," not defined")
	t
}

buildType <- function(format,descriptorType) tolower(paste(format,descriptorType,sep="-"))

addTransform(buildType("sdf","ap"),isDescriptorType=FALSE,
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

addTransform(buildType("sdf","fp"),isDescriptorType=FALSE,
	# Any sdf source -> fp string
	toString = function(input,dir=".") {
		getTransform("fp")$toString(
			getTransform(buildType("sdf","fp"))$toObject(input)$descriptors)
	},
	# Any sdf source -> FPset
	toObject = function(input,dir="."){
		apList = getTransform(buildType("sdf","ap"))$toObject(input,dir)
		apList$descriptors = desc2fp(apList$descriptors)
		apList
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
addTransform("fp",  
   # FPset -> string,
	toString = function(fpset,dir="."){
		sapply(1:length(fpset), function(i) as(fpset[i],"character") )
	},
   # string or list -> FP set list
	toObject= function(v,dir="."){ 
		if(inherits(v,"list") || length(v)==0)
			return(v)

		as( if(!inherits(v,"FPset")){
				#names(v)=as.character(1:length(v));  
				read.AP(v,type="fp",isFile=F)
			} else v,
			"FP")  
	}
)

