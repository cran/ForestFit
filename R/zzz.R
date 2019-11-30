zzz<- function(){ 
msg <- c(paste0(
"Welcome to
     _______ ________ ________  _________ _________ _________    _______ __ _________
    / _____// ____  // ____  / / _______// _______//___  ___/   / _____// //___  ___/
   / /___  / /   / // /___/ / / /______ / /______     / /      / /___  / /    / /    
  / ____/ / /   / // __   _/ / _______//______  /    / /      / ____/ / /    / /     
 / /     / /___/ // /  | |  / /______ _______/ /    / /      / /     / /    / /          
/_/     /_______//_/   |_| /________//________/    /_/      /_/     /_/    /_/      version ",
 
packageVersion("ForestFit")),"\nType 'citation(\"ForestFit\")' for citing this R package in publications.")
 return(msg)
}
#   unlockBinding("ForestFit", asNamespace("ForestFit")) 
.onAttach <- function(libname, pkgname) {
  mess <- zzz()
  if(!interactive())
    mess[1] <- paste("Package 'ForestFit' version", packageVersion("ForestFit"))
    packageStartupMessage(mess) 
  invisible()
  }
#  suppressPackageStartupMessages(library(ars))	
#  library(ForestFit, quietly = TRUE)
#  unloadNamespace(ars)
#suppressPackageStartupMessages(library(ForestFit))