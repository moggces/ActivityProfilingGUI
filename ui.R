library(shiny)
library(markdown)

shinyUI(
  bootstrapPage(
    tags$head(
      tags$meta(`http-equiv`="pragma", content="no-cache"),
      tags$meta(`http-equiv`="Cache-control", content="no-cache, no-store")
    ),
    tags$head(tags$script(src="extra.js")),
    tags$head(tags$link(rel="stylesheet", type="text/css", href="headerfooter.css")),
    tags$head(tags$link(rel="stylesheet", type="text/css", href="custom.css")),
    includeHTML(path="./www/header.html"),
    tabsetPanel(
      tabPanel('About', tags$div(includeMarkdown("README.Rmd"), class='container')),
      tabPanel('Application',
        pageWithSidebar(
          headerPanel('', windowTitle='Tox21 Activity Profiler'),
          sidebarPanel(
            
            # input control
            h3('Required inputs'),
            wellPanel(
              h4('CAS & Cluster (optional)'),
              tags$textarea(id="cmpds", class='col-xs-12', rows=5, ""),
              helpText("Note: copy & paste from excel file with two columns: "),
              
              h6('or'),
              fileInput('file1', 'Import a data matrix', multiple=FALSE),
              tags$hr(),
              
              h4('Activity type'),
              conditionalPanel(
                condition = "input.proftype == 'activity'",
                radioButtons("acttype", "select one of the metrics",
                             list("wAUC" = "nwauc.logit",  "POD" = "npod", "EC50" = "nec50"), selected="npod")
              ),
              
              radioButtons("proftype", "select the profiling type", list("signal", "activity"), selected = "activity"),
              conditionalPanel(
                condition = "input.proftype == 'signal'",
                radioButtons("sigtype", "Signal type:",
                             list("wAUC" = "wauc.logit"))
              )
              
            ),
            
            h3('Other settings'),
            wellPanel(
              h4('Limit activities'),
              tags$br(),
              
              conditionalPanel(
                condition = "(input.acttype == 'npod' || input.acttype == 'nec50' || input.acttype == 'nwauc.logit'  ) && input.proftype == 'activity'",
                sliderInput("nwauc_thres",
                            "wAUC threshold", min=0, max=1, value=0, step=0.05),
                tags$br(),
                sliderInput("ncmax_thres",
                            "Emax threshold", min=0, max=100, value=0, step=5),
                tags$br(),
                #sliderInput("npod_thres",
                #            "POD threshold",  min=3, max=10, value=3, step=0.5),
                numericInput("npod_thres", label = "POD threshold (uM)", value = NA),
                tags$br(),
                #sliderInput("nec50_thres",
                #           "EC50 threshold", min=3, max=10, value=3, step=0.5),
                numericInput("nec50_thres", label = "EC50 threshold (uM)", value = NA),
                
                tags$br(),
                sliderInput("pod_diff_thres",
                            "log10(ratio of signal to cytotoxicity)(inhibition-type assays only)", min=0, max=3, value=0, step=0.2),
                tags$hr(),
              
                
                checkboxInput("cytofilter", "exclude activity due to cytotoxicity (loose)", TRUE),
                
                tags$br(),
                checkboxInput("nocyto", "exclude activity due to cytotoxicity (strict)", FALSE),
                
                tags$br(),
                checkboxInput("isgoodcc2", "exclude activity with suboptimal NCATS fits", FALSE),
                
                tags$br(),
                checkboxInput("nohighcv", "exclude activity varied between sources", TRUE)
          
              )
              
            ),
            wellPanel(
              h4('Filter assays'),
              tags$br(),
              textInput('reg_sel', 'names (regular expression)', 'cytotoxicity'),
              checkboxInput("inv_sel", "invert your selection", TRUE)
              
            ),
            wellPanel(
              h4('Sort chemicals'),
              tags$br(),
              radioButtons("sort_method", "Method:",
                           list('structure similarity' = 'chemclust',
                                'activity similarity' = 'actclust',
                                'toxicity score (sum of activities)' = 'toxscore'), selected="actclust")
            ),
            h3('Formatting & outputs'),
            wellPanel(
              
              checkboxInput("noinconlab", "make inconclusive as inactive", TRUE),
              checkboxInput("showdendro", "show compound similarity dendrogram", FALSE),
              checkboxInput("keepsize", "keep heatmap size one page", FALSE),
              
              tags$br(),
              # fontsize
              sliderInput("fontsize",
                          "fontsize", min = 2,max = 28, value = 16, step=2),
              helpText("tip: lower the fontsize when saving the plot"),
              
              # output functions
              tags$hr(),
              downloadButton('downloadData', 'Download Activities'),
              downloadButton('downloadPlot', 'Save Plot'),
              downloadButton('downloadEnrich', 'Download Enrichment analysis')
              
            )    
            ),
          
          mainPanel(
            tabsetPanel(
              tabPanel('Input chemicals', dataTableOutput('contents')),
              tabPanel("Profile", plotOutput("profiling", height=1000, width="500%")), # i think the height don't affect
              tabPanel("Potency boxplot", plotOutput("box",  height=1000, width="500%")),
              tabPanel('Activity data', dataTableOutput('dd')),
              tabPanel('Enrichment analysis', dataTableOutput('enrich')),
              tabPanel('Assays', dataTableOutput('assay_info'))
            )
          )
        )
      )
    ),
    includeHTML(path="./www/footer.html")
  )
)

