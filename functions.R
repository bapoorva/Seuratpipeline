#Function to process  the data
createobj <- function(name,org,files,mincells,mingenes,aggr){
  try(if(length(files)==0) stop("No files"))
  
  if(length(files)==1){
    # Load the dataset
    inputdata <- Read10X_h5(files)
    if(aggr==T){
      colnames(inputdata) <- paste0(colnames(inputdata), name)
    }
    
    # Initialize the Seurat object with the raw (non-normalized data).  
    scrna <- CreateSeuratObject(raw.data = inputdata, min.cells = mincells, min.genes = mingenes,project = name)
    scrna@meta.data$var_sample <- name
  }else{
    #Initialize the first object with the raw (non-normalized data) and add rest of the data 
    inputdata <- Read10X_h5(files[1])
    scrna <- CreateSeuratObject(raw.data = inputdata,  min.cells = mincells, min.genes = mingenes, project = 'Rep1')
    #scrna@meta.data$var_sample <- name
    cat('Rep1', length(scrna@cell.names), "\n")
    for(i in 2:length(files)){
      tmp.data <- Read10X_h5(files[i])
      tmp.scrna <- CreateSeuratObject(raw.data = tmp.data, min.cells = mincells, min.genes = mingenes, project = paste0('Rep',i))
      #tmp.scrna@meta.data$var_sample <- name
      cat('Rep', i, ": ", length(tmp.scrna@cell.names), "\n", sep="")
      scrna <- MergeSeurat(scrna, tmp.scrna, do.normalize = FALSE, min.cells = 0, min.genes = 0, add.cell.id2 = paste0('Rep',i))
    }
    cat("merged: ", length(scrna@cell.names), "\n", sep="")
  }
  # calculate the percent.mito values.
  if(org=='mouse'){
    mito.genes <- grep(pattern = "^mt-", x = rownames(x = scrna@data), value = TRUE)
  }else{
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna@data), value = TRUE)
  }
  percent.mito <- Matrix::colSums(scrna@raw.data[mito.genes, ])/Matrix::colSums(scrna@raw.data)
  
  # AddMetaData adds columns to object@meta.data. metadata is a great place to stash QC stats
  scrna <- AddMetaData(object = scrna, metadata = percent.mito, col.name = "percent.mito")
  return(scrna)
}

#Function to get pval from jackstraw slot
jackstrawpval=function(object,PCs){
  
  
  score.thresh = 1e-5
  pAll <- object@dr$pca@jackstraw@emperical.p.value
  pAll <- pAll[, PCs, drop = FALSE]
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  pAll.l <- reshape2::melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  qq.df <- NULL
  score.df <- NULL
  for (i in PCs) {
    q <- qqplot(x = pAll[, i], y = runif(n = 1000), plot.it = FALSE)
    #pc.score=mean(q$y[which(q$x <=score.thresh)])
    pc.score <- suppressWarnings(prop.test(
      x = c(
        length(x = which(x = pAll[, i] <= score.thresh)),
        floor(x = nrow(x = pAll) * score.thresh)
      ),
      n = c(nrow(pAll), nrow(pAll))
    )$p.val)
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
      pc.score <- 1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
    } else {
      score.df <- rbind(score.df, data.frame(PC = paste0("PC",i), Score = pc.score))
    }
    if (is.null(x = qq.df)) {
      qq.df <- data.frame(x = q$x, y = q$y, PC = paste0("PC", i))
    } else {
      qq.df <- rbind(qq.df, data.frame(x = q$x, y = q$y, PC = paste0("PC", i)))
    }
  }
  return(score.df)
}  



