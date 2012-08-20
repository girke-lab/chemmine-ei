
library(RUnit)
library(eiR)


path = if(Sys.getenv("RCMDCHECK")=="FALSE"){ 
         file.path("eiR","inst","runit")
       }else{
         file.path(.path.package(package="eiR"),"runit")
      }

testsuite <- defineTestSuite("eiR tests", dirs = path)
testResult <- runTestSuite(testsuite)
printTextProtocol(testResult)
