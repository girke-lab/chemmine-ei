# This is an automatically generated file by the R module for SWIG.

##   Generated via the command line invocation:
##	 swig -c++ -module ChemmineRpp -o ChemmineRpp/src/r_wrap.cc -r ../swig.i


#                         srun.swg                            #
#
# This is the basic code that is needed at run time within R to
# provide and define the relevant classes.  It is included
# automatically in the generated code by copying the contents of
# srun.swg into the newly created binding code.


# This could be provided as a separate run-time library but this
# approach allows the code to to be included directly into the
# generated bindings and so removes the need to have and install an
# additional library.  We may however end up with multiple copies of
# this and some confusion at run-time as to which class to use. This
# is an issue when we use NAMESPACES as we may need to export certain
# classes.

######################################################################

if(length(getClassDef("RSWIGStruct")) == 0) 
  setClass("RSWIGStruct", representation("VIRTUAL"))



if(length(getClassDef("ExternalReference")) == 0) 
# Should be virtual but this means it loses its slots currently
#representation("VIRTUAL")
  setClass("ExternalReference", representation( ref = "externalptr"))



if(length(getClassDef("NativeRoutinePointer")) == 0) 
  setClass("NativeRoutinePointer", 
              representation(parameterTypes = "character",
                             returnType = "character",
                             "VIRTUAL"), 
              contains = "ExternalReference")

if(length(getClassDef("CRoutinePointer")) == 0) 
  setClass("CRoutinePointer", contains = "NativeRoutinePointer")


if(length(getClassDef("EnumerationValue")) == 0) 
  setClass("EnumerationValue", contains = "integer")


if(!isGeneric("copyToR")) 
 setGeneric("copyToR",
            function(value, obj = new(gsub("Ref$", "", class(value)))) 
               standardGeneric("copyToR"
           ))

setGeneric("delete", function(obj) standardGeneric("delete"))


SWIG_createNewRef = 
function(className, ..., append = TRUE)
{
  f = get(paste("new", className, sep = "_"), mode = "function")

  f(...)
}

if(!isGeneric("copyToC")) 
 setGeneric("copyToC", 
             function(value, obj = RSWIG_createNewRef(class(value)))
              standardGeneric("copyToC"
            ))


# 
defineEnumeration =
function(name, .values, where = topenv(parent.frame()), suffix = "Value")
{
   # Mirror the class definitions via the E analogous to .__C__
  defName = paste(".__E__", name, sep = "")
  assign(defName,  .values,  envir = where)

  if(nchar(suffix))
    name = paste(name, suffix, sep = "")

  setClass(name, contains = "EnumerationValue", where = where)
}

enumToInteger <- function(name,type)
{
   if (is.character(name)) {
   ans <- as.integer(get(paste(".__E__", type, sep = ""))[name])
   if (is.na(ans)) {warning("enum not found ", name, " ", type)}
   ans
   } 
}

enumFromInteger =
function(i,type)
{
  itemlist <- get(paste(".__E__", type, sep=""))
  names(itemlist)[match(i, itemlist)]
}

coerceIfNotSubclass =
function(obj, type) 
{
    if(!is(obj, type)) {as(obj, type)} else obj
}


setClass("SWIGArray", representation(dims = "integer"), contains = "ExternalReference")

setMethod("length", "SWIGArray", function(x) x@dims[1])


defineEnumeration("SCopyReferences",
                   .values = c( "FALSE" = 0, "TRUE" = 1, "DEEP" = 2))

assert = 
function(condition, message = "")
{
  if(!condition)
    stop(message)

  TRUE
}


if(FALSE) {
print.SWIGFunction =
function(x, ...)
 {
 }
}


#######################################################################

R_SWIG_getCallbackFunctionStack =
function()
{
    # No PACKAGE argument as we don't know what the DLL is.
  .Call("R_SWIG_debug_getCallbackFunctionData")
}

R_SWIG_addCallbackFunctionStack =
function(fun, userData = NULL)
{
    # No PACKAGE argument as we don't know what the DLL is.
  .Call("R_SWIG_R_pushCallbackFunctionData", fun, userData)
}


#######################################################################


setClass('C++Reference', contains = 'ExternalReference')
setClass('_p_Descriptors', contains = 'C++Reference')
setClass('_p_Database', contains = 'C++Reference')



setMethod('[', "ExternalReference",
function(x,i,j, ..., drop=TRUE) 
if (!is.null(x$"__getitem__")) 
sapply(i, function(n) x$"__getitem__"(i=as.integer(n-1))))

setMethod('[<-' , "ExternalReference",
function(x,i,j, ..., value) 
if (!is.null(x$"__setitem__")) {
sapply(1:length(i), function(n) 
x$"__setitem__"(i=as.integer(i[n]-1), x=value[n]))
x
})

setAs('ExternalReference', 'character',
function(from) {if (!is.null(from$"__str__")) from$"__str__"()})

setMethod('print', 'ExternalReference',
function(x) {print(as(x, "character"))})

# Start of new_Descriptors

`Descriptors` = function()
{
  ans = .Call('R_swig_new_Descriptors', PACKAGE='ChemmineRpp')
  class(ans) <- "_p_Descriptors"
  
  ans
  
}

attr(`Descriptors`, 'returnType') = '_p_Descriptors'
class(`Descriptors`) = c("SWIGFunction", class('Descriptors'))

# Start of Descriptors_parse_sdf

`Descriptors_parse_sdf` = function(self, sdf, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Descriptors") 
  sdf = as(sdf, "character") 
  .Call('R_swig_Descriptors_parse_sdf', self, sdf, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`Descriptors_parse_sdf`, 'returnType') = 'numeric'
attr(`Descriptors_parse_sdf`, "inputTypes") = c('_p_Descriptors', 'character')
class(`Descriptors_parse_sdf`) = c("SWIGFunction", class('Descriptors_parse_sdf'))

# Start of Descriptors_parse_sdfile

`Descriptors_parse_sdfile` = function(self, sdfile, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Descriptors") 
  sdfile = as(sdfile, "character") 
  .Call('R_swig_Descriptors_parse_sdfile', self, sdfile, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`Descriptors_parse_sdfile`, 'returnType') = 'numeric'
attr(`Descriptors_parse_sdfile`, "inputTypes") = c('_p_Descriptors', 'character')
class(`Descriptors_parse_sdfile`) = c("SWIGFunction", class('Descriptors_parse_sdfile'))

# Start of Descriptors_parse_smiles

`Descriptors_parse_smiles` = function(self, smile, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Descriptors") 
  smile = as(smile, "character") 
  .Call('R_swig_Descriptors_parse_smiles', self, smile, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`Descriptors_parse_smiles`, 'returnType') = 'numeric'
attr(`Descriptors_parse_smiles`, "inputTypes") = c('_p_Descriptors', 'character')
class(`Descriptors_parse_smiles`) = c("SWIGFunction", class('Descriptors_parse_smiles'))

# Start of Descriptors_get_descriptor

`Descriptors_get_descriptor` = function(self, i, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Descriptors") 
  i = as.numeric(i) 
  
  assert(length(i) == 1 && i >= 0, "All values must be non-negative")
  
  .Call('R_swig_Descriptors_get_descriptor', self, i, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`Descriptors_get_descriptor`, 'returnType') = 'numeric'
attr(`Descriptors_get_descriptor`, "inputTypes") = c('_p_Descriptors', 'numeric')
class(`Descriptors_get_descriptor`) = c("SWIGFunction", class('Descriptors_get_descriptor'))

# Start of Descriptors_get_len

`Descriptors_get_len` = function(self, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Descriptors") 
  .Call('R_swig_Descriptors_get_len', self, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`Descriptors_get_len`, 'returnType') = 'numeric'
attr(`Descriptors_get_len`, "inputTypes") = c('_p_Descriptors')
class(`Descriptors_get_len`) = c("SWIGFunction", class('Descriptors_get_len'))

# Start of delete_Descriptors

`delete_Descriptors` = function(self)
{
  self = coerceIfNotSubclass(self, "_p_Descriptors") 
  .Call('R_swig_delete_Descriptors', self, PACKAGE='ChemmineRpp')
  
}

attr(`delete_Descriptors`, 'returnType') = 'void'
attr(`delete_Descriptors`, "inputTypes") = c('_p_Descriptors')
class(`delete_Descriptors`) = c("SWIGFunction", class('delete_Descriptors'))

# Start of accessor method for Descriptors
setMethod('$', '_p_Descriptors', function(x, name)

{
  accessorFuns = list('parse_sdf' = Descriptors_parse_sdf, 'parse_sdfile' = Descriptors_parse_sdfile, 'parse_smiles' = Descriptors_parse_smiles, 'get_descriptor' = Descriptors_get_descriptor, 'get_len' = Descriptors_get_len)
  idx = pmatch(name, names(accessorFuns))
  if(is.na(idx)) 
  return(callNextMethod(x, name))
  f = accessorFuns[[idx]]
  formals(f)[[1]] = x
  f
}


)
# end of accessor method for Descriptors
setMethod('delete', '_p_Descriptors', function(obj) {delete_Descriptors(obj)})
# Start of similarity

`similarity` = function(d1, d2, .copy = FALSE)
{
  d1 = coerceIfNotSubclass(d1, "_p_Descriptors") 
  d2 = coerceIfNotSubclass(d2, "_p_Descriptors") 
  .Call('R_swig_similarity', d1, d2, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`similarity`, 'returnType') = 'numeric'
attr(`similarity`, "inputTypes") = c('_p_Descriptors', '_p_Descriptors')
class(`similarity`) = c("SWIGFunction", class('similarity'))

# Start of batch_parse

`batch_parse` = function(sdfile, dbfile, .copy = FALSE)
{
  sdfile = as(sdfile, "character") 
  dbfile = as(dbfile, "character") 
  .Call('R_swig_batch_parse', sdfile, dbfile, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`batch_parse`, 'returnType') = 'numeric'
attr(`batch_parse`, "inputTypes") = c('character', 'character')
class(`batch_parse`) = c("SWIGFunction", class('batch_parse'))

# Start of new_Database

`Database` = function()
{
  ans = .Call('R_swig_new_Database', PACKAGE='ChemmineRpp')
  class(ans) <- "_p_Database"
  
  ans
  
}

attr(`Database`, 'returnType') = '_p_Database'
class(`Database`) = c("SWIGFunction", class('Database'))

# Start of Database_open

`Database_open` = function(self, filename, mode, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Database") 
  filename = as(filename, "character") 
  mode = as(mode, "character");     
  .Call('R_swig_Database_open', self, filename, mode, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`Database_open`, 'returnType') = 'numeric'
attr(`Database_open`, "inputTypes") = c('_p_Database', 'character', 'character')
class(`Database_open`) = c("SWIGFunction", class('Database_open'))

# Start of Database_close

`Database_close` = function(self)
{
  self = coerceIfNotSubclass(self, "_p_Database") 
  .Call('R_swig_Database_close', self, PACKAGE='ChemmineRpp')
  
}

attr(`Database_close`, 'returnType') = 'void'
attr(`Database_close`, "inputTypes") = c('_p_Database')
class(`Database_close`) = c("SWIGFunction", class('Database_close'))

# Start of Database_next

`Database_next` = function(self, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Database") 
  ans = .Call('R_swig_Database_next', self, as.logical(.copy), PACKAGE='ChemmineRpp')
  class(ans) <- "_p_Descriptors"
  
  ans
  
}

attr(`Database_next`, 'returnType') = '_p_Descriptors'
attr(`Database_next`, "inputTypes") = c('_p_Database')
class(`Database_next`) = c("SWIGFunction", class('Database_next'))

# Start of Database_store

`Database_store` = function(self, d, .copy = FALSE)
{
  self = coerceIfNotSubclass(self, "_p_Database") 
  d = coerceIfNotSubclass(d, "_p_Descriptors") 
  .Call('R_swig_Database_store', self, d, as.logical(.copy), PACKAGE='ChemmineRpp')
  
}

attr(`Database_store`, 'returnType') = 'numeric'
attr(`Database_store`, "inputTypes") = c('_p_Database', '_p_Descriptors')
class(`Database_store`) = c("SWIGFunction", class('Database_store'))

# Start of delete_Database

`delete_Database` = function(self)
{
  self = coerceIfNotSubclass(self, "_p_Database") 
  .Call('R_swig_delete_Database', self, PACKAGE='ChemmineRpp')
  
}

attr(`delete_Database`, 'returnType') = 'void'
attr(`delete_Database`, "inputTypes") = c('_p_Database')
class(`delete_Database`) = c("SWIGFunction", class('delete_Database'))

# Start of accessor method for Database
setMethod('$', '_p_Database', function(x, name)

{
  accessorFuns = list('open' = Database_open, 'close' = Database_close, 'next' = Database_next, 'store' = Database_store)
  idx = pmatch(name, names(accessorFuns))
  if(is.na(idx)) 
  return(callNextMethod(x, name))
  f = accessorFuns[[idx]]
  formals(f)[[1]] = x
  f
}


)
# end of accessor method for Database
setMethod('delete', '_p_Database', function(obj) {delete_Database(obj)})

