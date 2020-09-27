fileload <- function(project){
  inFile = paste('data/',as.character(project),'.RData',sep = '')
  load(inFile)
  loaddata=scrna
  return(loaddata)
}