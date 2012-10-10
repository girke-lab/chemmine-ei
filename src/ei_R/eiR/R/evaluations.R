

evaluator <- function(reference,result,output=NA)
{
   maxQueries=1000
   maxK=1000
   rs = c(1,2,5,10,20,30,40,50)
   ks = c(1,10,20,50,100,200,500,1000)

   counters=matrix(0,length(rs),maxK)


   updateCounter <- function(counterIndex, ratio)
   {
      captured=c()

      for(cut in 1:maxK)
      {
         cut2=(cut)*ratio
         if(cut2 > length(x_ref)) 
            break

         contrib=0
         if(cut==1){
            if(x_ref[cut] %in% names(d_target) && 
                  d_target[as.character(x_ref[cut])] <= cut2)
               contrib = contrib+1
            captured=c(captured,contrib)
         }
         else{
            if(x_ref[cut] %in% names(d_target) && 
                  d_target[as.character(x_ref[cut])] <= cut2)
               contrib = contrib+1
            for(i in 1:ratio){
               if(x_target[cut2-i+1]  %in% names(d_ref) && 
                     d_ref[as.character(x_target[cut2-i+1])] < cut)
                  contrib=contrib+1
            }
            #print(paste("rest:",contrib))
            captured=c(captured,captured[cut-1]+contrib)
         }
      }
      #print("captured:")
      #print(captured)
      #print(paste("counterIndex:",counterIndex))
      for(i in 1:length(captured))
         counters[counterIndex,i]<<-counters[counterIndex,i]+captured[i]
      #print(counters[counterIndex,])
      #print("----------------")

   }

   getRanking <- function(x)
      as.numeric(sapply(unlist(strsplit(x,"\\s+")),
                     function(y) strsplit(y,":",fixed=TRUE)[[1]][1]))

   ###########################################

   refFile=gzfile(reference,"r")
   targetFile=gzfile(result,"r")

   for(lineNum in 1:maxQueries)
   {
      refLine=readLines(refFile,n=1)
      targetLine=readLines(targetFile,n=1)

      if(length(refLine)==0 || length(targetLine)==0)
         break

      x_ref=getRanking(refLine)
      d_ref=1:length(x_ref)
      names(d_ref)=as.character(x_ref)

      x_target=getRanking(targetLine)
      d_target=1:length(x_target)
      names(d_target)=as.character(x_target)

      for(i in 1:length(rs)){
         updateCounter(i,rs[i])
      }

      #print("=============================================")
      
   }
   close(targetFile)
   close(refFile)


   ratios=sapply(1:length(rs), function(i){
      x=counters[i,]/1:maxK/maxQueries
      if(!is.na(output))
         write.table(x,
                  file=paste(output,".ratio-",rs[i],sep=""),
                  col.names=F,row.names=F,quote=F)
      x
   })
   colnames(ratios)=rs

   ratioSummary=sapply(1:length(ks),function(j) counters[,ks[j]]/ks[j]/maxQueries)
   rownames(ratioSummary)=rs
   colnames(ratioSummary)=ks
   if(!is.na(output))
      write.csv(ratioSummary,file=paste(output,"csv",sep="."),quote=FALSE)
   return(list(ratios=ratios,ratioSummary=ratioSummary))
}
