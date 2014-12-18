library(shiny)
library(markdown)
# todo: need to be replaced by markdown page
#tabPanelAbout <- source(paste(getwd(), "/source/about.R", sep=""))$value

shinyUI(pageWithSidebar(
  
  # application title
  headerPanel("Compound signal/activity profiling (by Tox21 qHTS assays)"),
  
  sidebarPanel(
    
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
    
    # profiling options
    h4('Profiling'),
    
    radioButtons("proftype", "Profile type:", 
                 list("signal", "activity")),
    tags$hr(),
    
    # todo: wauc matrix
    conditionalPanel(
      condition = "input.proftype == 'signal'",
      radioButtons("sigtype", "Signal type:",
                   list("wAUC" = "wauc.logit"))
    ),
    
    tags$br(),
    
    # todo: nwauc, pod, ac50 matrix
    conditionalPanel(
      condition = "input.proftype == 'activity'",
      radioButtons("acttype", "Activity type:",
                   list("wAUC" = "nwauc.logit",  "POD" = "npod", "AC50" = "nac50"))
    ),
    
    tags$br(),
    
    # todo: nwauc, emax, pod, ac50, pod_med_diff, call matrix
    conditionalPanel(
      h4('Activity filtering (not applicable if loading the data matrix'),
      
      condition = "(input.acttype == 'npod' || input.acttype == 'nac50' || input.acttype == 'nwauc.logit'  ) && input.proftype == 'activity'",
      sliderInput("nwauc_thres", 
                  "wAUC threshold", min=0, max=1, value=0.05, step=0.05),
      tags$br(),
      sliderInput("nemax_thres", 
                  "Emax threshold", min=0, max=100, value=0, step=5),
      tags$br(),
      sliderInput("npod_thres", 
                  "POD threshold",  min=3, max=10, value=3, step=0.5),
      tags$br(),
      sliderInput("nac50_thres", 
                  "AC50 threshold", min=3, max=10, value=3, step=0.5),
      tags$br(),
      sliderInput("pod_diff_thres", 
                  "log10(ratio of signal to cytotoxicity)(inhibition-type assays only)", min=0, max=3, value=0, step=0.2),
      checkboxInput("nocyto", "not cytotoxic (inhibition-type assays only)", FALSE),
      #tags$br(),
      #checkboxInput("isstrong", "wauc above the median of least potent (10 uM)", FALSE),
      
      tags$br(),
      checkboxInput("isgoodcc2", "only curve class 1.1, 1.2, 2.1", FALSE),
      
      tags$br(),
      checkboxInput("nohighcv", "exclude high activity variation between sources", TRUE)
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
#    helpText("e.g. to view only assays used Tox21 Library version#1: cytotoxicity|pparg_antagonism|ppard|are|hse plus inverse"),
#    helpText("to select exact assays: \bmitotox\b|er_antagonism|ar_antagonism(mdakb2)|\baromatase\b|\bare\b"),
    
    tags$hr(),
    
    # miscellaneous functions
    h4('Others'),

    # todo: adjust the variable
    # to show the dendrogram
    checkboxInput("showdendro", "show compound similarity dendrogram ", FALSE),
    
    tags$br(),
    
    # fontsize
    sliderInput("fontsize", 
                "fontsize", min = 2,max = 28, value = 18, step=2),
    helpText("tip: lower the fontsize when saving the plot"),
  
    # output functions
    br(),
    downloadButton('downloadData', 'Download'),
    downloadButton('downloadPlot', 'Save Plot')
    
    
  ),
  mainPanel(
    tabsetPanel(
     
      tabPanel( 'Input chemicals', dataTableOutput('contents')),
      #tabPanel( 'Input chemicals', htmlOutput('contents')),
      tabPanel( "Profile", plotOutput("profiling", height=1000, width="500%")), # i think the height don't affect
      tabPanel( "Potency boxplot", plotOutput("box",  height=1000, width="500%")),
      tabPanel( 'Data', dataTableOutput('assay_des')),
      #tabPanelAbout()
      tabPanel('About', includeHTML("test.html"))
    )
  )
))

