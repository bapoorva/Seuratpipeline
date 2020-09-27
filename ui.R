library(shinydashboard)
library(shiny)
library(shinyBS)
library(shinyjs)
library(visNetwork)
library(dashboardthemes)
library(shinyFiles)
library(rhandsontable)
options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize=50*1024^2) 
ui <- dashboardPage(skin = "red",
  dashboardHeader(title = "SeUrAt PiPeLiNe"),
  dashboardSidebar(div(style="overflow-y: scroll"),
                   sidebarMenu(
                     menuItem("Single", tabName = "single", icon = icon("user")),
                     menuItem("Aggregate", tabName = "aggregate", icon = icon("users"))
                     )#End sidebar menu
  ),#end dashboardSidebar

  ######################################################################################################################################
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    tags$head(tags$style(HTML(".sidebar { height: 420vh; overflow-y: auto; }
                                             .shiny-notification{position: fixed;top: 33%;left: 45%;right: 30%;}" )
                         )),
    useShinyjs(),
    tabItems(
   tabItem(tabName = "single",
    tabBox(id = "tabset1",width=12,
      tabPanel("Input", 
               fileInput("dir", label = "Choose input file",multiple = T),
               textInput("project", label = "Enter Project Name"),
                 #tableOutput("filelist"),
               selectInput("organism", label = "Choose an organism",choices = list("Mouse" = 'mouse', "Human" = 'human'),selected = 1),
               sliderInput("mincells", label = "Minimum number of cells that should express a gene", min = 1, max = 500, value = 3),
               sliderInput("mingenes", label = "Mininum number of genes expressed in a cell", min = 100, max = 500, value = 200),
               actionButton("preprocess", label = "Preprocess Data")
      ),
      tabPanel("QC Plots", uiOutput("plotsubtab"),
               actionButton("proceedtofilt", label = "Proceed to Filtration")),
      tabPanel("Filtration and Normalization", h4("Select variables to filter by"),
               uiOutput("filtvars"),
               textInput("lowthres", label = "Low threshold", value = "500,-Inf"),
               textInput("highthres", label = "High Threshold",value = "Inf,0.05"),
               radioButtons("ccscale", label = "Regress out cell cycle",inline=T,choices = list("Yes" = 'yes', "No" = 'no'),selected = 'yes'),
               actionButton("runnorm", label = "Run"),br(),br(),
               conditionalPanel(
                 condition = "input.runnorm ==true",actionButton("proceedtopca", label = "Proceed to PCA")
               )
               ),
      tabPanel("PCA",
               sliderInput("pcslider", label = "Total number of PCS's to compute", min = 5,max = 100, value = 80),
               uiOutput("pcasubtab"),
               actionButton("proceedtojs", label = "Proceed to JackStraw")),
      tabPanel("JackStraw",
               uiOutput("jssubtab"),
               actionButton("proceedtocluster", label = "Proceed to Clustering")
        ),
      tabPanel("Clustering",
               fluidRow(
               column(6,uiOutput('finddim')),
               column(6,sliderInput("maxres", label = "Resolution", min = 0.4,max = 1.2, value = 0.6,step=0.2))
               ),
               fluidRow(
               column(12,sliderInput("pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25))),
               fluidRow(
                 column(6,plotOutput("tSNEplot")),
                 column(6,plotOutput("uMAPplot"))
               ),
               actionButton("proceedtomarker", label = "Proceed to Marker Genes"),
               downloadButton("saverdata", label = "Save Dataset")
        ),
    tabPanel("Marker Genes",
             uiOutput('markergenes'),
             sliderInput("ptsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
             uiOutput("plotclusttab"),
             actionButton("proceedtoct", label = "Proceed to Assign Celltypes")
    ),
    tabPanel("Assign Cell types",
    radioButtons("ct", label = "Do you want to assign celltypes?",inline=T,choices = list("Yes" = 'yes', "No" = 'no'),selected = 'no'),
    conditionalPanel(
      condition = "input.ct == 'yes'",
      rHandsontableOutput('celltype'),br(),
      actionButton("runButton","Assign celltypes"),br(),br(),
      uiOutput('setcat'),
      actionButton("findmar","Run Find Markers"),br(),br()
    ),
    #verbatimTextOutput('vargenes'),
    #DT::dataTableOutput('vargenes'),
    downloadButton("saverdata2", label = "Save Dataset")
    )
    )#End of tabbox
   ),#End of tab single
   ######################################################################################################################################
   tabItem(tabName = "aggregate",
           tabBox(id = "tabset2",width=12,
                  tabPanel("Input", 
                     selectInput("aggrfiles", label = "Enter how many input files to aggregate",choices = 1:10,selected = 1),
                     uiOutput("aggrinput"),
                     selectInput("organism2", label = "Choose an organism",choices = list("Mouse" = 'mouse', "Human" = 'human'),selected = 1),
                     sliderInput("mincells2", label = "Minimum number of cells that should express a gene", min = 1, max = 500, value = 3),
                     sliderInput("mingenes2", label = "Mininum number of genes expressed in a cell", min = 100, max = 500, value = 200),
                     actionButton("preprocess2", label = "Preprocess Data")) ,
                  tabPanel("QC Plots", uiOutput("qcsubtab"),
                              actionButton("proceedtofilt2", label = "Proceed to Filtration")),
                  tabPanel("Filtration and Normalization", h4("Select variables to filter by"),
                           uiOutput("filtvars2"),
                           textInput("lowthres2", label = "Low threshold", value = "500,-Inf"),
                           textInput("highthres2", label = "High Threshold",value = "Inf,0.05"),
                           radioButtons("ccscale2", label = "Regress out cell cycle",inline=T,choices = list("Yes" = 'yes', "No" = 'no'),selected = 'yes'),
                           actionButton("runnorm2", label = "Run"),br(),br(),
                           conditionalPanel(
                             condition = "input.runnorm2 ==true",actionButton("proceedtocca", label = "Combine and Proceed to MultiCCA")
                           ),
                           verbatimTextOutput("console")),
                tabPanel("MultiCCA", 
                         radioButtons("radiofileup", label = "File input type",choices = list("Continue from previous tab " = "cont", "Upload new Rdata" = "upload"),selected = "cont"),
                         conditionalPanel(
                           condition = "input.radiofileup =='upload'",fileInput('ccafileupload', 'Upload RData File')),
                         sliderInput("numccs", label = "Number of CC's to run", min =5,max = 60, value =5,step=1),
                         uiOutput("ccasubtab"),
                        actionButton("proceedtoclust2", label = "Proceed to Clustering")),
                tabPanel("Clustering",
                         fluidRow(
                           column(6,uiOutput('maxdim2')),
                           column(6,sliderInput("maxres2", label = "Resolution", min = 0.4,max = 1.2, value = 0.6,step=0.2))
                         ),
                         fluidRow(
                           column(12,sliderInput("pointsize2", "Point Size:",min = 0, max = 5, value = 1,step=.25))),
                         fluidRow(
                           column(6,plotOutput("tSNEplot2")),
                           column(6,plotOutput("uMAPplot2"))
                         ),
                         fluidRow(
                           column(6,plotOutput("tSNEplotsamp")),
                           column(6,plotOutput("uMAPplotsamp"))
                         ),
                         actionButton("proceedtomarker2", label = "Proceed to Marker Genes"),
                         downloadButton("saverdataaggr", label = "Save Dataset")
                ),
                tabPanel("Marker Genes",
                         uiOutput('markergenes2'),
                         sliderInput("ptsize2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                         uiOutput("plotclusttab2"),
                         actionButton("proceedtoct2", label = "Proceed to Assign Celltypes")
                ),
                tabPanel("Assign Cell types",
                         radioButtons("ct2", label = "Do you want to assign celltypes?",inline=T,choices = list("Yes" = 'yes', "No" = 'no'),selected = 'no'),
                         conditionalPanel(
                           condition = "input.ct2 == 'yes'",
                           rHandsontableOutput('celltype2'),br(),
                           actionButton("runButton2","Assign celltypes"),br(),br(),
                           uiOutput('setcat2'),
                           actionButton("findmar2","Run Find Markers"),br(),br()
                         ),
                         downloadButton("saverdataggr2", label = "Save Dataset")
                )
           )#End of tabbox
)#End of tab aggregate
    )#End if tab items 
  )#end of dashboard body
  )#end of dashboard page
