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
