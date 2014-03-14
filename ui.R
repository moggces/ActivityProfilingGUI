
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

#setwd("~/ShinyApps/profiling/")
tabPanelAbout <- source(paste(getwd(), "/source/about.R", sep=""))$value

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Tox21 qHTS assays: signal/activity profiling"),
  
  sidebarPanel(
    h4('Input'),
    fileInput('file1', 'Import a list of chemcials', multiple=FALSE),
    helpText("Note: a tab-delimited file with two columns: CAS & Cluster"),
    h6('or'),
    selectInput("dataset", "Choose a pre-defined set:", 
                choices = c("no selection", "polycyclic aromatic hydrocarbons (PAHs)", "flame retardants (FRs)")),
    tags$hr(),
    
    h4('Profiling'),
    
    radioButtons("proftype", "Profile type:", 
                 list("signal", "activity")),
    tags$hr(),
    
    conditionalPanel(
      condition = "input.proftype == 'signal'",
      radioButtons("sigtype", "Signal type:",
                   list("wAUC" = "signal_wauc"))
    ),
      
    conditionalPanel(
      condition = "input.proftype == 'activity'",
      radioButtons("acttype", "Activity type:",
                   list("wAUC" = "nwauc",  "POD" = "pod", "AC50" = "ac50")),
      radioButtons("actstrict", "Activity stringency:",
                   list("loose" = "loose",  "medium" = "medium", "tight" = "tight"), selected = "medium")
    ),
    
    
    tags$hr(),
    
    h4('Sort columns by ...'),
    radioButtons("sort_method", "Method:",
                 list('structure similarity' = 'chemclust',
                      'activity similiarty' = 'actclust',
                      'toxicity score (not for signal)' = 'toxscore')),
                 
    #checkboxInput("chemclust", "chemical similarity", TRUE),
    
    tags$hr(),
    
    h4('Options'),
    checkboxInput("recyto", "remove cytotoxicity assays", TRUE),
    checkboxInput("showheat", "show the heatmap", TRUE),
    
    tags$br(),
    
    sliderInput("fontsize", 
                "fontsize", min = 2,max = 28, value = 18, step=2),
    
    conditionalPanel(
      condition = "input.acttype == 'pod' || input.acttype == 'ac50'" ,
      sliderInput("nwauc_thres", 
                  "wAUC threshold for plotting", min=0, max=1, value=0.1, step=0.05)
    ),
  
    
    br(),
    downloadButton('downloadData', 'Download')
    #downloadButton('downloadPlot', 'Save Plot')
    
    
  ),
  mainPanel(
    tabsetPanel(
     
      tabPanel( 'Input chemicals', dataTableOutput('contents')),
      #tabPanel( 'Input chemicals', htmlOutput('contents')),
      tabPanel( "Profile", plotOutput("profiling", height=1000, width="500%")),
      tabPanelAbout()
    )
  )
))

