library(shinydashboard)
library(shiny)
library(shinyBS)
library(shinyjs)
library(dashboardthemes)
library(shinyFiles)
library(Seurat)
library(DT)
library(rhandsontable)
library(dplyr)
library(readr)
source("functions.R")
cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")
mouseorthologfile <- "data/mouse_human.csv"

server <- function(input, output,session) {
 
    #Process the input data
  preprocess = reactive({
    withProgress(session = session, message = 'Preprocessing...',detail = 'Please Wait...',{
      validate(need(input$project,"Please enter project name"))
      files=input$dir$datapath
        scrna=createobj(name=input$project,org=input$organism,files=files,input$mincells,input$mingenes,aggr=F)
    })
  })

  output$filelist = renderTable(input$dir)

  ######################################################################################################
  ######################################################################################################
  ###########################################   QC PLOTS ###############################################
  ######################################################################################################
  ######################################################################################################
  #Move to QC plots tab when pre-process button is clicked
  observeEvent(input$preprocess, {
    updateTabsetPanel(session, "tabset1",
                      selected = "QC Plots")
  })
  
  #Create subtabs for the QC plots
  output$plotsubtab <- renderUI({
    tabsetPanel(id = "subTabPanel1",
                tabPanel("Violin Plot",plotOutput("violinplot")),
                tabPanel("Gene Plot",plotOutput("geneplot"))
    )
   
  })
  
  #Create Violin plot
  output$violinplot = renderPlot({
    withProgress(session = session, message = 'Preprocessing...',detail = 'Please Wait...',{
      scrna=preprocess()
      VlnPlot(object = scrna, features.plot = c("nGene", "nUMI", "percent.mito"),do.return = T)
    })
  })
  
  #Create geneplot
  output$geneplot = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      scrna=preprocess()
      par(mfrow = c(1, 2))
      GenePlot(object = scrna, gene1 = "nUMI", gene2 = "percent.mito")
      GenePlot(object = scrna, gene1 = "nUMI", gene2 = "nGene")
    })
  })
  
  ######################################################################################################
  ######################################################################################################
  ################################# FILTERATION/NORMALIZATION/SCALING  ################################
  ######################################################################################################
  ######################################################################################################
  #Move to Filteration tab when button is clicked
  observeEvent(input$proceedtofilt, {
    updateTabsetPanel(session, "tabset1",
                      selected = "Filtration and Normalization")
  })
  
  #Dropdown for variables to filter by
  output$filtvars = renderUI({
    scrna=preprocess()
    withProgress(session = session, message = 'Generating gene names...',detail = 'Please Wait...',{
    genes=rownames(scrna@data)
    selectInput('varstofilt', 'Options',c("nGene","percent.mito",genes), multiple=TRUE, selectize=TRUE, selected = c("nGene","percent.mito"))
    })
  })
  #Filtercells, then normalize and scale
  filtercells= reactive({
  scrna=preprocess()
  validate(need(input$varstofilt,"Select variables to filter"))
  v=vector()
  for(i in 1:length(input$varstofilt)){v=c(v,input$varstofilt[i]) }
  lowthres=input$lowthres
  lt=as.character(strsplit(lowthres,",")[[1]])
  highthres=input$highthres
  ht=as.character(strsplit(highthres,",")[[1]])
  withProgress(session = session, message = 'Filtering Cells...',detail = 'Please Wait...',{
  scrna <- FilterCells(object = scrna, subset.names = v,low.thresholds = as.numeric(lt), high.thresholds =as.numeric(ht))
    })
  #normalize data
  withProgress(session = session, message = 'Normalizing Data...',detail = 'Please Wait...',{
  scrna <- NormalizeData(object = scrna, normalization.method = "LogNormalize",scale.factor = 10000)
  })
  #detection of variable genes
  #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
  withProgress(session = session, message = 'Detecting Variable Genes...',detail = 'Please Wait...',{
    pdf(NULL)
  scrna <- FindVariableGenes(object = scrna, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  dev.off()
  length(x = scrna@var.genes)
  })
  withProgress(session = session, message = 'Scaling Data...',detail = 'Please Wait...',{
  if(input$ccscale=='yes'){
    # if(input$organism=='human'){
    #   #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
    #   scrna <- CellCycleScoring(object = scrna, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
    # }else 
      if(input$organism=='mouse'){
      m2h <- read_csv(mouseorthologfile)
      cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
      cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
    }
    #Scaling the data and removing unwanted sources of variation
    scrna <- CellCycleScoring(object = scrna, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
    scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
  }else{
    scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito"))
  }
  return(scrna)
    })
  })
  
  #Run filtration/normalization and scaling
  observeEvent(input$runnorm, {
    filtercells()
  })
  
  #Move to plots tab when proceed button is clicked
  observeEvent(input$proceedtopca, {
    updateTabsetPanel(session, "tabset1",
                      selected = "PCA")
  })
  
  ######################################################################################################
  ######################################################################################################
  ############################################ PCA ####################################################
  ######################################################################################################
  ######################################################################################################
  
  runpca = reactive({
    withProgress(session = session, message = 'Running PCA...',detail = 'Please Wait...',{
    scrna=filtercells()
    scrna <- RunPCA(object = scrna, pc.genes = scrna@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5,pcs.compute=input$pcslider)
    })
  })
  
  #Generate drop down menu to populate the max number of dimensions used in the scRNA analysis
  output$ndim = renderUI({
    scrna=runpca()
    maxdim="NA"
    maxdim=scrna@calc.params$RunPCA$pcs.compute
    validate(need(is.na(maxdim)==F,"PCA dimensional Reduction has not been computed"))
    var=1:maxdim
    selectInput("ndim","Choose number of dimensions",var,selected = 1)
  })
  
  #Create subtabs for the PCA plots
  output$pcasubtab <- renderUI({
    tabsetPanel(id = "pcasubtab",
                tabPanel("Viz Plot",
                         fluidRow( column(6,uiOutput("ndim")),
                         column(6,sliderInput("ngenes", "Number of genes:",min = 10, max = 50, value = 10,step=5))),
                        plotOutput("vizplot",height = 600)),
                tabPanel("PCA Plot",plotOutput("pcaplot")),
                tabPanel("PCA Heatmap",plotOutput("pcaheatmap",height = 3800))
    )
    
  })
  
  #Create Viz plot
  output$vizplot = renderPlot({
    withProgress(session = session, message = 'Generating Plot...',detail = 'Please Wait...',{
      scrna=runpca()
      dim=input$ndim
      validate(need(dim,"PCA dimensional Reduction has not been computed"))
      VizPCA(object = scrna, pcs.use = dim:dim,nCol=1,font.size = 1,num.genes = input$ngenes)
    })
  })
  
  #Create pcaplot
 output$pcaplot = renderPlot({
    withProgress(session = session, message = 'Generating Plot...',detail = 'Please Wait...',{
      scrna=runpca()
      PCAPlot(object = scrna, dim.1 = 1, dim.2 = 2)
    })
  })
 
 #Create PCA Heatmap
 output$pcaheatmap = renderPlot({
   withProgress(session = session, message = 'Generating Plot...',detail = 'Please Wait...',{
     scrna=runpca()
     PCHeatmap(object = scrna, pc.use = 1:input$pcslider, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
   })
 })
 
 #Move to jackstraw tab when proceed button is clicked
 observeEvent(input$proceedtojs, {
   updateTabsetPanel(session, "tabset1",
                     selected = "JackStraw")
 })
 ######################################################################################################
 ######################################################################################################
 ############################################ JACKSTRAW ###############################################
 ######################################################################################################
 ######################################################################################################
 #Run jackstraw
 runjackstraw= reactive({
   withProgress(session = session, message = 'Running jackstraw...',detail = 'Please Wait...',{
   scrna=runpca()
   scrna <- JackStraw(object = scrna, num.replicate = 100,num.pc=input$pcslider,do.par=T)
 })
 })
 
 #Create subtabs for the PCA plots
 output$jssubtab <- renderUI({
   tabsetPanel(id = "jssubtab",
               tabPanel("JackStraw Plot",
                        plotOutput("jsplot",height = 2000),
                        DT::dataTableOutput("pcpval")),
               tabPanel("Elbow Plot",plotOutput("elbowplot"))
   )
 })
 
 #Create JackStraw Plot
 output$jsplot = renderPlot({
   withProgress(session = session, message = 'Generating Plot...',detail = 'Please Wait...',{
     scrna=runjackstraw()
     JackStrawPlot(object = scrna, PCs = 1:input$pcslider)
   })
 })
 
 #Get a table of PC values
 pcpval= reactive({
   scrna=runjackstraw()
   df=jackstrawpval(scrna,PCs=1:input$pcslider)
 })
 
 #Render table
 output$pcpval = DT::renderDataTable({
   withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     DT::datatable(pcpval(),
                   extensions = c('Buttons','Scroller'),
                   options = list(dom = 'Bfrtip',
                                  searchHighlight = TRUE,
                                  pageLength = 10,
                                  lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                  scrollX = TRUE,
                                  buttons = c('copy', 'print')
                   ),rownames=FALSE,caption= "P-Values of PC's",selection = list(mode = 'single', selected =1),escape = F)
   })
 })
 
 
 #Create Elbow Plot
 output$elbowplot = renderPlot({
   withProgress(session = session, message = 'Generating Plot...',detail = 'Please Wait...',{
     scrna=runjackstraw()
     PCElbowPlot(object = scrna,num.pc =input$pcslider)
   })
 })
 
 #Move to Clustering tab when proceed button is clicked
 observeEvent(input$proceedtocluster, {
   updateTabsetPanel(session, "tabset1",
                     selected = "Clustering")
 })
 ######################################################################################################
 ######################################################################################################
 ###################################### CLUSTERING ####################################################
 ######################################################################################################
 ######################################################################################################
 #Find maximum dims used
 output$finddim = renderUI({
   sliderInput("maxdim", label = "Dimension", min = 5,max =input$pcslider, value = input$pcslider,step=1)
 })
 
  #Run Clustering
 runclust= reactive({
   withProgress(session = session, message = 'Finding Clusters...',detail = 'Please Wait...',{
   scrna=runjackstraw()
   scrna <- RunTSNE(object = scrna, dims.use = 1:input$maxdim, do.fast = TRUE)
   scrna <- RunUMAP(object = scrna, dims.use = 1:input$maxdim, min_dist=0.5,n_neighbors = 15,metric = 'correlation')
   scrna <- FindClusters(scrna, reduction.type = "pca",dims.use = 1:input$maxdim, resolution = input$maxres)
 })
 })
 
 #Make tsne and umap plots
 #Create tSNE Plot
 output$tSNEplot = renderPlot({
   withProgress(session = session, message = 'Generating tSNE Plot...',detail = 'Please Wait...',{
     scrna=runclust()
     TSNEPlot(scrna, do.label = T,do.return=T,colors.use = cpallette,pt.size = input$pointsize)
   })
 })
 
 #Create UMAP Plot
 output$uMAPplot = renderPlot({
   withProgress(session = session, message = 'Generating uMAP...',detail = 'Please Wait...',{
     scrna=runclust()
     DimPlot(scrna, reduction.use = "umap", do.label = T,do.return=T,cols.use = cpallette,pt.size = input$pointsize)
   })
 })

 #Move to Markergenes tab when proceed button is clicked
 observeEvent(input$proceedtomarker, {
   updateTabsetPanel(session, "tabset1",
                     selected = "Marker Genes")
 })
 
 #Save data
 output$saverdata <- downloadHandler(
   filename = function() {
     paste0(input$project,"_",input$maxdim,"_",input$maxres,".RData",sep="")
   },
   content = function(file){
     saveRDS(runclust(),file=file)
   })
 
 ######################################################################################################
 ######################################################################################################
 ####################################### MARKER GENES ################################################
 ######################################################################################################
 ######################################################################################################
 #Create dropdown for marker genes
 output$markergenes <- renderUI({
   withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
   scrna=runclust()
   genes= rownames(scrna@data)
   selectizeInput('markergenes', 'Enter Marker genes (comma separated)', choices = genes, multiple = TRUE,options = list(maxItems = 12))
 })
 })
 
 #Create subtabs for the PCA plots
 output$plotclusttab <- renderUI({
   tabsetPanel(id = "plotclusttab",
               tabPanel("Marker Genes - tSNE",plotOutput("marktSNEplot",height = 900)),
               tabPanel("Marker Genes - UMAP",plotOutput("markuMAPplot",height = 900))
   )
   
 })
 
 #Make tsne and umap plots for marker genes
 #Create tSNE Plot
 output$marktSNEplot = renderPlot({
   withProgress(session = session, message = 'Generating tSNE Plot...',detail = 'Please Wait...',{
     validate(need(input$markergenes,"Enter at least one gene name"))
     scrna=runclust()
     a=FeaturePlot(object = scrna, features.plot = input$markergenes, cols.use = c("lightgrey", "blue"),do.return = T,pt.size = input$ptsize)
     plot_grid(plotlist=a)
   })
 })
 
 #Create UMAP Plot
 output$markuMAPplot = renderPlot({
   withProgress(session = session, message = 'Generating uMAP...',detail = 'Please Wait...',{
     scrna=runclust()
     a=FeaturePlot(object = scrna, features.plot = input$markergenes, cols.use = c("lightgrey", "blue"),do.return = T,reduction.use = "umap",pt.size = input$ptsize)
     plot_grid(plotlist=a)
   })
 })
 
 #Move to Assign cell types tab when proceed button is clicked
 observeEvent(input$proceedtoct, {
   updateTabsetPanel(session, "tabset1",
                     selected = "Assign Cell types")
 })
 ######################################################################################################
 ######################################################################################################
 ##################################### ASSIGN CELL TYPES ##############################################
 ######################################################################################################
 ######################################################################################################
 #Create table to assign cell-types
 celltype = reactive({
   withProgress(session = session, message = 'Generating table...',detail = 'Please Wait...',{
   scrna=runclust()
   a= levels(scrna@ident)
   b=vector(mode="character",length(a))
   df=data.frame(a,b,stringsAsFactors = F)
   colnames(df)=c("cluster","celltype")
   return(df)
 })
 })
 
 #Create datatable to enter celltypes and process it to reset table
 output$celltype = renderRHandsontable({
   table=celltype()
   rhandsontable(as.data.frame(table)) %>% hot_col("cluster",readOnly=TRUE)
   })
 
 #Switch ident based on user-defined cell-types
 changecelltype = reactive({
   scrna=runclust()
   scrna@meta.data= scrna@meta.data %>% rename("var_cluster" := paste("res.",input$maxres,sep=""))
   if(input$ct=="yes"){
     validate(need(values$data,"Enter values for celltype"))
     #validate(need(length(values$data$celltype)==length(values$data$cluster),"Please enter celltype for"))
     celltype=as.data.frame(values$data)
     colnames(celltype)=c("clust","var_celltype")
     meta=scrna@meta.data
     scrna@meta.data=left_join(meta,celltype,by=c("var_cluster"="clust"))
     rownames(scrna@meta.data)=rownames(meta)
     scrna <- SetAllIdent(object = scrna, id = "var_celltype")
   }else{
     #Else set the cluster of selected resolution as the identity/main group of comparison
     scrna <- AddMetaData(scrna, scrna@ident,col.name= "var_cluster")
   }
   return(scrna)
 })
 
 
 values <- reactiveValues()
 observeEvent(input$runButton, {
   withProgress(session = session, message = 'Generating table...',detail = 'Please Wait...',{
    # validate(need(length(input$celltype$data)==length(values$data$cluster),"Please enter celltype for"))
   values$data <-  hot_to_r(input$celltype)
   #changecelltype()
   })
 })

 output$vargenes <- renderText(length(values$data$celltype[is.na(values$data$celltype)==F]))

 ######################################################################################################
 ######################################################################################################
 ###################################   RUN FIND MARKERS ###############################################
 ######################################################################################################
 ######################################################################################################
 
 #Create dropdown for setting category
 output$setcat <- renderUI({
   withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
     scrna=changecelltype()
     metadata=as.data.frame(scrna@meta.data)
     metadata=metadata %>% dplyr::select(starts_with("var_"))
     var=c(colnames(metadata))
     selectInput("setcategory","Choose category to run Find Markers with",var,"pick one")
   })
 })
 
 #run find markers
runfindmarkers <- reactive({
  withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=changecelltype()
    scrna <- SetAllIdent(object = scrna, id = input$setcategory)
    scrna@misc=NA
    scrna@misc <-  vector(mode="list", length=length(levels(scrna@ident)))
    names(scrna@misc)=levels(scrna@ident)
    for(c in levels(scrna@ident)){
      scrna@misc[[c]] <- FindMarkers(scrna,ident.1 = c) %>% tibble::rownames_to_column('gene_name')
      rownames(scrna@misc[[c]])=scrna@misc[[c]]$gene_name
    } 
    return(scrna)
  })
 })

observeEvent(input$findmar, {
  withProgress(session = session, message = 'Generating table...',detail = 'Please Wait...',{
    runfindmarkers()
  })
})

#Save data
output$saverdata2 <- downloadHandler(
  filename = function() {
    paste0(input$project,"_dim",input$maxdim,"_res",input$maxres,".RData",sep="")
  },
  content = function(file){
    saveRDS(runfindmarkers(),file=file)
  })
 ######################################################################################################
 ######################################################################################################
 #################################### AGGREGATE INPUTS  ###############################################
 ######################################################################################################
 ######################################################################################################
# Enter as many fileinput boxes as it takes to aggregate data
output$aggrinput <- renderUI({
  num <- as.integer(input$aggrfiles)  
  lapply(1:(num), function(i) {
    fluidRow(
    column(6,fileInput(paste0("aggrfile",i,sep=""), label = paste0("Choose input file ",i,sep=""),multiple = F)),
    column(6,textInput(paste0("txt", i),label = paste0("Name",i,sep=""),value=""))
    )
  }) #end of lapply
}) #

#Preprocess aggregate data
preprocess2 = reactive({
  withProgress(session = session, message = 'Preprocessing...',detail = 'Please Wait...',{
    num <- as.integer(input$aggrfiles)
    scrna <-  vector(mode="list", length=num)
    lapply(1:(num), function(i) {
    validate(
      need(input[[as.character(paste0("aggrfile",i,sep=""))]],paste0("Please Select input",i,"file",sep="")),
      need(input[[as.character(paste0("txt",i,sep=""))]],paste0("Please Select input",i,"name",sep=""))
      )
      ipfile=input[[as.character(paste0("aggrfile",i,sep=""))]]$datapath
      scrna[[i]]=assign(input[[as.character(paste0("txt",i,sep=""))]],createobj(name=input[[as.character(paste0("txt",i,sep=""))]],org=input$organism2,files=ipfile,mincells=input$mincells2,mingenes=input$mingenes2,aggr = T))
      #scrna=createobj(name=input$txt1,org=input$organism2,files=input$aggrfile1$datapath,mincells=input$mincells2,mingenes=input$mingenes2)
    })
  })
})

#Move to QC plots tab when pre-process button is clicked
observeEvent(input$preprocess2, {
  updateTabsetPanel(session, "tabset2",
                    selected = "QC Plots")
})


#Create subtabs for the QC plots
output$qcsubtab <- renderUI({
  num=as.numeric(input$aggrfiles)
  tabsetPanel(id = "qcsubTabPanel",
              tabPanel("Violin Plot",plotOutput("violinplot2",height = num*400)),
              tabPanel("Gene Plot",plotOutput("geneplot2",height = num*400))
  )
  
})

#Create Violin plot
output$violinplot2 = renderPlot({
  withProgress(session = session, message = 'Preprocessing...',detail = 'Please Wait...',{
    scrna=preprocess2()
    num <- as.integer(input$aggrfiles)
    vln <-  vector(mode="list", length=num)
    for (i in 1:num){
      vln[[i]]=VlnPlot(object = scrna[[i]], features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,do.return=TRUE)
    }
    plot_grid(plotlist = vln,align = 'v',axis="tb",ncol=1)
  })
})

#Create geneplot
output$geneplot2 = renderPlot({
  withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=preprocess2()
    num <- as.integer(input$aggrfiles)
    par(mfrow = c(num, 2))
    for (i in 1:num){
    GenePlot(object = scrna[[i]], gene1 = "nUMI", gene2 = "percent.mito")
    GenePlot(object = scrna[[i]], gene1 = "nUMI", gene2 = "nGene")
    }
  })
})

######################################################################################################
######################################################################################################
################################# FILTERATION/NORMALIZATION/SCALING AGGR  ############################
######################################################################################################
######################################################################################################
#Move to Filteration tab when button is clicked
observeEvent(input$proceedtofilt2, {
  updateTabsetPanel(session, "tabset2",
                    selected = "Filtration and Normalization")
})

#Dropdown for variables to filter by
output$filtvars2 = renderUI({
  scrna=preprocess2()
  withProgress(session = session, message = 'Generating gene names...',detail = 'Please Wait...',{
    num <- as.integer(input$aggrfiles)
    genes=vector(mode="list", length=num)
    for (i in 1:num){
      genes[[i]]=rownames(scrna[[i]]@data)
    }
    genes=Reduce(intersect, genes)
    selectInput('varstofilt2', 'Options',c("nGene","percent.mito",genes), multiple=TRUE, selectize=TRUE, selected = c("nGene","percent.mito"))
  })
})


#Filtercells, then normalize and scale
filtercells2= reactive({
  scrna=preprocess2()
  num=input$aggrfiles
  validate(need(input$varstofilt2,"Select variables to filter"))
  v=vector()
  for(i in 1:length(input$varstofilt2)){v=c(v,input$varstofilt2[i]) }
  lowthres=input$lowthres2
  lt=as.character(strsplit(lowthres,",")[[1]])
  highthres=input$highthres2
  ht=as.character(strsplit(highthres,",")[[1]])
  for(i in 1:num){
  withProgress(session = session, message = 'Filtering Cells...',detail = 'Please Wait...',{
    scrna[[i]] <- FilterCells(object = scrna[[i]], subset.names = v,low.thresholds = as.numeric(lt), high.thresholds =as.numeric(ht))
  })
  #normalize data
  withProgress(session = session, message = 'Normalizing Data...',detail = 'Please Wait...',{
    scrna[[i]] <- NormalizeData(object = scrna[[i]], normalization.method = "LogNormalize",scale.factor = 10000)
  })
  #detection of variable genes
  #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
  withProgress(session = session, message = 'Detecting Variable Genes...',detail = 'Please Wait...',{
    scrna[[i]] <- FindVariableGenes(object = scrna[[i]], do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    length(x = scrna[[i]]@var.genes)
  })
  withProgress(session = session, message = 'Scaling Data...',detail = 'Please Wait...',{
    if(input$ccscale2=='yes'){
        if(input$organism2=='mouse'){
        m2h <- read_csv(mouseorthologfile)
        cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
        cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
      }
      #Scaling the data and removing unwanted sources of variation
      scrna[[i]] <- CellCycleScoring(object = scrna[[i]], s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
      scrna[[i]] <- ScaleData(object = scrna[[i]], vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
    }else{
      scrna[[i]] <- ScaleData(object = scrna[[i]], vars.to.regress = c("nUMI", "percent.mito"))
    }
  })
}

  return(scrna)
})
#Run filtration/normalization and scaling
observeEvent(input$runnorm2, {
  filtercells2()
})

#Move to multicca when proceed button is clicked
observeEvent(input$proceedtocca, {
  updateTabsetPanel(session, "tabset2",
                    selected = "MultiCCA")
})

######################################################################################################
######################################################################################################
############################################ MultiCCA ################################################
######################################################################################################
######################################################################################################

#Run multicca
multicca = reactive({
  if(input$radiofileup=='cont'){
  ob.list <- filtercells2()}
  else if(input$radiofileup=="upload"){
    file=input$ccafileupload
    load(file$datapath)
  }
  genes.use <- c()
  for (i in 1:length(ob.list)) {
    genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
  }
  genes.use <- names(which(table(genes.use) > 1))
  for (i in 1:length(ob.list)) {
    genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
  }
  withProgress(session = session, message = 'Running MultiCCA...',detail = 'Please Wait...',{
  combined <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = input$numccs)
  })
  return(combined)
})


#Create subtabs for the CCA plots
output$ccasubtab <- renderUI({
  tabsetPanel(id = "ccasubTabPanel",
              tabPanel("Dim Plot",plotOutput("dimplot")),
              tabPanel("MetageneBicorPlot,",plotOutput("metagenebicorplot"))
              
  )
  
})

#Create Dim plot and Violin plot for cca
output$dimplot = renderPlot({
  withProgress(session = session, message = 'Generating plot...',detail = 'Please Wait...',{
    combined=multicca()
    p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "var_sample", 
                  pt.size = 0.5, do.return = TRUE)
    p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "var_sample", 
                  do.return = TRUE)
    plot_grid(p1, p2)
  })
})

#Create metagenebicorplot
output$metagenebicorplot = renderPlot({
  withProgress(session = session, message = 'Generating plot...',detail = 'Please Wait...',{
    combined=multicca()
    MetageneBicorPlot(combined, grouping.var = "var_sample", dims.eval = 1:input$numccs)
  })
})

#Move to Clustering when proceed button is clicked
observeEvent(input$proceedtoclust2, {
  updateTabsetPanel(session, "tabset2",
                    selected = "Clustering")
})
######################################################################################################
######################################################################################################
###################################### CLUSTERING AGGR ###############################################
######################################################################################################
######################################################################################################
#Find maximum dims used
output$maxdim2 = renderUI({
  sliderInput("maxdim2", label = "Dimension", min = 5,max =input$numccs, value = input$numccs,step=1)
})

#Run Clustering
runclustcca= reactive({
  withProgress(session = session, message = 'Finding Clusters...',detail = 'Please Wait...',{
    combined=multicca()
    combined <- CalcVarExpRatio(object = combined, reduction.type = "pca",grouping.var = "var_sample", dims.use = 2:input$maxdim2)
    scrna <- SubsetData(combined, subset.name = "var.ratio.pca",accept.low = 0.5)
    scrna <- AlignSubspace(scrna,reduction.type = "cca",grouping.var = "var_sample",dims.align = 1:input$maxdim2)
    
    # t-SNE and Clustering
    scrna <- RunTSNE(scrna,reduction.use = "cca.aligned",dims.use = 1:input$maxdim2,check_duplicates = FALSE,nthreads = 10,max_iter = 2000)
    
    ### run UMAP
    scrna <- RunUMAP(scrna, reduction.use = "cca.aligned", dims.use = 1:input$maxdim2)
    #Find clusters
    scrna <- FindClusters(scrna, reduction.type = "cca.aligned",dims.use = 1:input$maxdim2, save.SNN = T, resolution = input$maxres2)
  })
})

#Make tsne and umap plots
#Create tSNE Plot
output$tSNEplot2 = renderPlot({
  withProgress(session = session, message = 'Generating tSNE Plot...',detail = 'Please Wait...',{
    scrna=runclustcca()
    TSNEPlot(scrna, do.label = T,do.return=T,colors.use = cpallette,pt.size = input$pointsize2)
  })
})

#Create UMAP Plot
output$uMAPplot2 = renderPlot({
  withProgress(session = session, message = 'Generating uMAP...',detail = 'Please Wait...',{
    scrna=runclustcca()
    DimPlot(scrna, reduction.use = "umap", do.label = T,do.return=T,cols.use = cpallette,pt.size = input$pointsize2)
  })
})

#Create tSNE Plot with samples
output$tSNEplotsamp = renderPlot({
  withProgress(session = session, message = 'Generating tSNE Plot...',detail = 'Please Wait...',{
    scrna=runclustcca()
    DimPlot(scrna,reduction.use = "tsne",group.by='var_sample', do.label = T,do.return=T,colors.use = cpallette,pt.size = input$pointsize2)
  })
})

#Create UMAP Plot with samples
output$uMAPplotsamp = renderPlot({
  withProgress(session = session, message = 'Generating uMAP...',detail = 'Please Wait...',{
    scrna=runclustcca()
    DimPlot(scrna, reduction.use = "umap",group.by='var_sample', do.label = T,do.return=T,cols.use = cpallette,pt.size = input$pointsize2)
  })
})
#Move to Markergenes tab when proceed button is clicked
observeEvent(input$proceedtomarker2, {
  updateTabsetPanel(session, "tabset2",
                    selected = "Marker Genes")
})

#Save data
output$saverdataaggr <- downloadHandler(
  filename = function() {
    paste0(input$project,"_",input$maxdim,"_",input$maxres,"_cca.RData",sep="")
  },
  content = function(file){
    saveRDS(runclustcca(),file=file)
  })

######################################################################################################
######################################################################################################
#################################### MARKER GENES AGGR ###############################################
######################################################################################################
######################################################################################################

#Create dropdown for marker genes
output$markergenes2 <- renderUI({
  withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=runclustcca()
    genes= rownames(scrna@data)
    selectizeInput('markergenes2', 'Enter Marker genes', choices = genes, multiple = TRUE,options = list(maxItems = 12))
  })
})

#Create subtabs for the PCA plots
output$plotclusttab2 <- renderUI({
  tabsetPanel(id = "plotclusttabaggr",
              tabPanel("Marker Genes - tSNE",plotOutput("marktSNEplot2",height = 900)),
              tabPanel("Marker Genes - UMAP",plotOutput("markuMAPplot2",height = 900))
  )
  
})

#Make tsne and umap plots for marker genes
#Create tSNE Plot
output$marktSNEplot2 = renderPlot({
  withProgress(session = session, message = 'Generating tSNE Plot...',detail = 'Please Wait...',{
    validate(need(input$markergenes2,"Enter at least one gene name"))
    scrna=runclustcca()
    a=FeaturePlot(object = scrna, features.plot = input$markergenes2, cols.use = c("lightgrey", "blue"),do.return = T,pt.size = input$ptsize2)
    plot_grid(plotlist=a)
  })
})

#Create UMAP Plot
output$markuMAPplot2 = renderPlot({
  withProgress(session = session, message = 'Generating uMAP...',detail = 'Please Wait...',{
    validate(need(input$markergenes2,"Enter at least one gene name"))
    scrna=runclustcca()
    a=FeaturePlot(object = scrna, features.plot = input$markergenes2, cols.use = c("lightgrey", "blue"),do.return = T,reduction.use = "umap",pt.size = input$ptsize2)
    plot_grid(plotlist=a)
  })
})

#Move to Assign cell types tab when proceed button is clicked
observeEvent(input$proceedtoct2, {
  updateTabsetPanel(session, "tabset2",
                    selected = "Assign Cell types")
})

######################################################################################################
######################################################################################################
################################ ASSIGN CELL TYPES AGGR ##############################################
######################################################################################################
######################################################################################################
#Create table to assign cell-types
celltype2 = reactive({
  withProgress(session = session, message = 'Generating table...',detail = 'Please Wait...',{
    scrna=runclustcca()
    a= levels(scrna@ident)
    b=vector(mode="character",length(a))
    df=data.frame(a,b,stringsAsFactors = F)
    colnames(df)=c("cluster","celltype")
    return(df)
  })
})

#Create datatable to enter celltypes and process it to reset table
output$celltype2 = renderRHandsontable({
  table=celltype2()
  rhandsontable(as.data.frame(table)) %>% hot_col("cluster",readOnly=TRUE)
})

values2 <- reactiveValues()
observeEvent(input$runButton2, {
  withProgress(session = session, message = 'Generating table...',detail = 'Please Wait...',{
    values2$data <-  hot_to_r(input$celltype2)
    #changecelltype()
  })
})

#Switch ident based on user-defined cell-types
changecelltype2 = reactive({
  scrna=runclustcca()
  scrna@meta.data= scrna@meta.data %>% rename("var_cluster" := paste("res.",input$maxres2,sep=""))
  if(input$ct2=="yes"){
    validate(need(values2$data,"Enter values for celltype"))
    celltype=as.data.frame(values2$data)
    colnames(celltype)=c("clust","var_celltype")
    meta=scrna@meta.data
    scrna@meta.data=left_join(meta,celltype,by=c("var_cluster"="clust"))
    rownames(scrna@meta.data)=rownames(meta)
    scrna <- SetAllIdent(object = scrna, id = "var_celltype")
  }else{
    #Else set the cluster of selected resolution as the identity/main group of comparison
    scrna <- AddMetaData(scrna, scrna@ident,col.name= "var_cluster")
  }
  return(scrna)
})

######################################################################################################
######################################################################################################
###################################   RUN FIND MARKERS ###############################################
######################################################################################################
######################################################################################################

#Create dropdown for setting category
output$setcat2 <- renderUI({
  withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=changecelltype2()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% dplyr::select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setcategory2","Choose category to run Find Markers with",var,"pick one")
  })
})

#run find markers
runfindmarkers2 <- reactive({
  withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=changecelltype2()
    scrna <- SetAllIdent(object = scrna, id = input$setcategory2)
    scrna@misc=NA
    scrna@misc <-  vector(mode="list", length=length(levels(scrna@ident)))
    names(scrna@misc)=levels(scrna@ident)
    for(c in levels(scrna@ident)){
      scrna@misc[[c]] <- FindMarkers(scrna,ident.1 = c) %>% tibble::rownames_to_column('gene_name')
      rownames(scrna@misc[[c]])=scrna@misc[[c]]$gene_name
    }
    return(scrna)
  })
})

observeEvent(input$findmar2, {
  withProgress(session = session, message = 'Generating table...',detail = 'Please Wait...',{
    runfindmarkers2()
  })
})

#Save data
output$saverdataggr2 <- downloadHandler(
  filename = function() {
    paste0(input$project,"_dim",input$maxdim,"_res",input$maxres,"_cca.RData",sep="")
  },
  content = function(file){
    saveRDS(runfindmarkers2(),file=file)
  })


}