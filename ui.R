library(shiny)
library(markdown)


shinyUI(pageWithSidebar(
  
  # application title
  headerPanel("Tox21 Activity Profiler"),
  
  sidebarPanel(
    
    tags$head(tags$meta(`http-equiv`="pragma", content="no-cache"), 
              tags$meta(`http-equiv`="Cache-control", content="no-cache, no-store")),
    
    # profiling options
    h4('Profiling'),
    
    radioButtons("proftype", "Profile type:", 
                 list("signal", "activity"), selected = "activity"),
    tags$hr(),
    
    conditionalPanel(
      condition = "input.proftype == 'signal'",
      radioButtons("sigtype", "Signal type:",
                   list("wAUC" = "wauc.logit"))
    ),
    
    tags$br(),
    
    conditionalPanel(
      condition = "input.proftype == 'activity'",
      radioButtons("acttype", "Activity type:",
                   list("wAUC" = "nwauc.logit",  "POD" = "npod", "EC50" = "nec50"), selected="npod")
    ),
    
    tags$br(),
    
    conditionalPanel(
      h4('Activity filtering (not applicable if loading the data matrix)'),
      
      condition = "(input.acttype == 'npod' || input.acttype == 'nec50' || input.acttype == 'nwauc.logit'  ) && input.proftype == 'activity'",
      sliderInput("nwauc_thres", 
                  "wAUC threshold", min=0, max=1, value=0, step=0.05),
      tags$br(),
      sliderInput("ncmax_thres", 
                  "Emax threshold", min=0, max=100, value=0, step=5),
      tags$br(),
      sliderInput("npod_thres", 
                  "POD threshold",  min=3, max=10, value=3, step=0.5),
      tags$br(),
      sliderInput("nec50_thres", 
                  "EC50 threshold", min=3, max=10, value=3, step=0.5),
      tags$br(),
      sliderInput("pod_diff_thres", 
                  "log10(ratio of signal to cytotoxicity)(inhibition-type assays only)", min=0, max=3, value=0, step=0.2),
      tags$br(),
      checkboxInput("cytofilter", "apply pipeline's cytotoxicity filter", TRUE),
      
      tags$br(),
      checkboxInput("nocyto", "no cytotoxicity observed in tested concentration range (inhibition-type assays only)", FALSE),
      
      tags$br(),
      checkboxInput("isgoodcc2", "only curve class 1.1, 1.2, 2.1", FALSE),
      
      tags$br(),
      checkboxInput("nohighcv", "exclude high activity variation between sources", TRUE),
      
      tags$br(),
      checkboxInput("noinconlab", "exclude inconclusive label", TRUE)
    ),
    
    tags$hr(),
    
    # arrange the compounds
    h4('Sort compounds by ...'),
    radioButtons("sort_method", "Method:",
                 list('structure similarity' = 'chemclust',
                      'activity similarity' = 'actclust',
                      'toxicity score (only for activity)' = 'toxscore')),
    
    tags$hr(),
    
    # filter the assays
    h4('Filter assays'), 
    textInput('reg_sel', 'names (regular expression)', 'cytotoxicity'),
    checkboxInput("inv_sel", "invert your selection", TRUE),
    
    tags$hr(),
    
    # input control
    h4('Input'),
    tags$textarea(id="cmpds", rows=3, cols=1, ""),
    helpText("Note: copy & paste from excel file with two columns: CAS|GSID & Cluster"),
    
    h6('or'),
    fileInput('file1', 'Import the data matrix', multiple=FALSE),
    #    helpText("Note: a tab-delimited file with two columns: CAS & Cluster|userClust"),
    #    h6('or'),
    #    selectInput("dataset", "Choose a pre-defined set:", 
    #                choices = c("no selection", "polycyclic aromatic hydrocarbons (PAHs)", "flame retardants (FRs)")),
    tags$hr(),
    
    # miscellaneous functions
    h4('Others'),

    checkboxInput("showdendro", "show compound similarity dendrogram", FALSE),
    checkboxInput("keepsize", "keep heatmap size one page", FALSE),
    
    tags$br(),
    
    # fontsize
    sliderInput("fontsize", 
                "fontsize", min = 2,max = 28, value = 16, step=2),
    helpText("tip: lower the fontsize when saving the plot"),
  
    # output functions
    br(),
    downloadButton('downloadData', 'Download Activities'),
    downloadButton('downloadPlot', 'Save Plot')
    
    
  ),
  mainPanel(
    tabsetPanel(
     
      tabPanel( 'Input chemicals', dataTableOutput('contents')),
      tabPanel( "Profile", plotOutput("profiling", height=1000, width="500%")), # i think the height don't affect
      tabPanel( "Potency boxplot", plotOutput("box",  height=1000, width="500%")),
      tabPanel( 'Activity data', dataTableOutput('dd')),
      tabPanel( 'Assays', dataTableOutput('assay_info')),
      tabPanel('About', includeHTML("README.html"))
    )
  )
))

