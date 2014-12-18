
# shiny note: 1) it can't discriminate R vs. r in the r script file 
# 2) the renderTable has trouble with encoding issue (cannot recognize ppp file iconv -t UTF-8 -f ISO-8859-1)
# looks like the new read.table can automatically select best encoding
# 3) in shiny server, once you delete a file but replace a file with same name. somwhow don't know how to refresh its
# but if you update .R file, you can refresh to get new functions


library(shiny)
library(plyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)
options(stringsAsFactors = FALSE)


#Sys.setlocale(locale="C")
#setwd("~/ShinyApps/profiling/")
source(paste(getwd(), "/source/customized.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/pheatmap_display_number.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/get.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/load.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/mis.R", sep=""), local=TRUE)
environment(pheatmap_new_label) <- environment(pheatmap)

# the example files
pah_file <- './data/pah_60_lead_clust_v3.txt'
fr_file <- './data/flame_retardant_clusters_selected_v4.txt'

# load assay related parameters
logit_para_file <- './data/tox21_assay_collection.txt'  #logit_para_file call_match_pathway_name
assay_names <- load_profile(logit_para_file) # global, dataframe output

# load chemical information (will include purity later)
profile_file <- './data/tox21_compound_id_v5a.txt'
master <- load_profile(profile_file) # global, dataframe output

# ??
removeAbnormalDirection <- FALSE

# load the activities (all data) and the structure fp matrix
activities_rdata <- './data/activities.RData'
struct_mat_rdata <- './data/struct_mat.RData'
load(activities_rdata) # global, matrix output, activities
load(struct_mat_rdata) # global, matrix output, struct_mat


# heatmap settings
wauc_breaks <- c( -1, -0.75, -0.5, -0.25, -0.1, -0.02, 0, 0.0001, 0.1, 0.25, 0.5, 0.75, 1) # upper is filled , lower is empty 
wauc_colors <-  c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE", "#D1E5F0", "#F7F7F7", "gray", "#FDDBC7" ,"#F4A582" ,"#D6604D" ,"#B2182B", "#67001F"  ) #RdBu
wauc_leg_breaks <- c(-1, -0.75, -0.5, -0.25,  0,   0.25, 0.5, 0.75, 1 )
wauc_leg_labels <- c("-1", "-0.75", "-0.5", "-0.25",  "0", "0.25", "0.5", "0.75", "1")
potency_breaks <- c(-10, -9, -7.5, -5, -4.5, -4, -0.02, 0, 0.0001, 4, 4.5, 5, 7.5, 9, 10)
potency_colors <- c("darkorange","#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "gray", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30", "chartreuse") #BrBG
potency_leg_breaks <- c(-10, -9, -7.5, -5, -4.5, -4,  0,  4, 4.5, 5, 7.5, 9,10 )
potency_leg_labels <- c("-10", "-9", "-7.5", "-5", "-4.5", "-4",  "0",  "4", "4.5", "5", "7.5", "9", "10")


shinyServer(function(input, output) {
   
  chemical_loader <- reactive({
    
    result <- NULL
    inFile <- input$file1
    path <- NULL
#     path <- switch(input$dataset,
#                      "no selection" = NULL, 
#                      "polycyclic aromatic hydrocarbons (PAHs)" = pah_file,
#                      "flame retardants (FRs)" = fr_file)
    textdata <- input$cmpds
    if (! is.null(inFile)) { path <- inFile$datapath; filen <- inFile$name }
    if (textdata != '' ) result <- load_text_2_df(textdata)
    if (! is.null(path)) result <- load_data_matrix(path, filen) # as long as path or file has something it will override
   
    return(result)

  })
  
  matrix_subsetter <- reactive({
    partial <- NULL
    reg_sel <- input$reg_sel # select the assays
    inv_sel <- input$inv_sel # inverse the selection
    rename_assay <- TRUE
    
    # get all chemical information
    id_info <- chemical_loader()
    
    chem_id_df <- get_lookup_list(id_info[['id']], master)
    #ip <- subset(chem_id_df, ! is.na(StructureID), select=c(CAS, Cluster))
    
    # the basic identifies , GSID + Cluster
    ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
    
    # collect all the matrices
    
    full <- list()
    full <- activities 
    
    # if it is a data matrix
    if (length(id_info) > 1)
    {
      full <- id_info[! grepl('id', names(id_info))]
      chemical_name_ref <- conversion(master, inp='CAS', out='GSID')
      rownames(full[[1]]) <- chemical_name_ref[as.character(rownames(full[[1]]))]
      rename_assay <- FALSE
    }
    full[['struct']] <- struct_mat
    
    # subset the matrices by chemicals
    partial <- get_input_chemical_mat(ip, full)
    
    # rename the assays
    partial <- rename_mat_col_row(partial,  master, assay_names, rename_assay=rename_assay)
    
    # subset the matrices by assay names
    partial <- get_assay_mat(partial, reg_sel, invSel=inv_sel)
    
    # sort the matrix
    partial <- sort_matrix(partial)
    
    return(partial)
  })
  
  activity_filter <- reactive({
    profile_type <- input$proftype
    activity_type <- input$acttype
    nwauc_thres <- input$nwauc_thres
    nemax_thres <- input$nemax_thres
    npod_thres <- input$npod_thres
    nac50_thres <- input$nac50_thres
    pod_diff_thres <- input$pod_diff_thres
    #isstrong <- input$isstrong
    nocyto <- input$nocyto
    isgoodcc2 <- input$isgoodcc2
    nohighcv <- input$nohighcv
    
    partial <- matrix_subsetter()
    
    # if it is data matrix
    if (length(partial) == 2) return(partial)
   
    act_mat_names <- c('npod', 'nac50', 'nwauc.logit')
    # reverse direction of mitotox could be meaningful
    partial <- fix_mitotox_reverse(partial,act_mat_names=act_mat_names )
    
#    # sort the matrix
#    partial <- sort_matrix(partial)
    
    # filtering
    partial <- filter_activity_by_type(partial, 'nwauc.logit', nwauc_thres, act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'nemax', nemax_thres,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'npod', npod_thres,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'nac50', nac50_thres,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'pod_med_diff', pod_diff_thres,act_mat_names=act_mat_names)
    #partial <- filter_activity_by_type(partial, 'hitcall', thres=NULL, decision=isstrong,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'pod_med_diff', thres=NULL, decision=nocyto,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'cc2', thres=NULL, decision=isgoodcc2,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'cv.wauc', thres=NULL, decision=nohighcv,act_mat_names=act_mat_names)
    
    return(partial)
  })
  
  matrix_editor <- reactive({
    
    partial <- activity_filter()
    # if it is data matrix
    if (length(partial) == 2) return(partial)
    
    # create CV marks
    cv_mark <- get_cv_mark_mat(partial[['cv.wauc']], partial[['nwauc.logit']])
    partial[['cv_mark']] <- cv_mark
    
    # make activities matrix as 0.0001
    act_mat_names <- c('npod', 'nac50', 'nwauc.logit')
    partial <- assign_reverse_na_number(partial, act_mat_names=act_mat_names)
    acts <- partial[c( act_mat_names, 'wauc.logit', 'struct', 'cv_mark')]
    return(acts)
  })
  
  
  heatmap_para_generator <- reactive({
    sort_meth <- input$sort_method
    profile_type <- input$proftype
    activity_type <- ''
    
    # get all chemical information
    chem_id_df <- get_lookup_list(chemical_loader()[['id']], master)
    #ip <- subset(chem_id_df, ! is.na(StructureID), select=c(CAS, Cluster))
    
    # the basic identifies , GSID + Cluster
    ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
    
    # the cleaned matrices
    dt <- matrix_editor() # c('npod', 'nac50', 'nwauc.logit','wauc.logit', 'cv_mark', 'struct') 
    
    # if the input is data matrix
    if (length(dt) == 2 ) 
    {
      activity_type <- names(dt)[1]
      act <- dt[[1]]
      cv <-  matrix("", nrow(act), ncol(act), dimnames=dimnames(act))
    } else
    {
      if (profile_type == 'activity')
      {
        activity_type <- input$acttype
        act <- dt[[activity_type]]
        
      } else
      {
        act <- dt[['wauc.logit']]
      }
      
      cv <- dt[['cv_mark']]
      
    }
    
    # struct matrix
    struct <- dt[['struct']]
    
    # first, cluster the chemicals
    dcols <- dist(struct, method = "binary") ## chemicals
    
    annotation <- get_heatmap_annotation(dcols, ip, master, dmat=dt, actType=activity_type) #data.frame output
    annt_colors <- get_heatmap_annotation_color(annotation,  actType=activity_type)
    
    # cluster compounds by various methods
    if (sort_meth == 'actclust')
    {
      dcols <- dist(act, method = "euclidean") ## chemicals by assays
      
    } else if (sort_meth == 'toxscore' )
    {
      tox_order <- rownames(annotation)[order(annotation$toxScore)]
      act <- act[tox_order, ]
      cv <- cv[tox_order, ]
    }
    # cluster assays by similarity 
    drows <- dist(t(act) , method = "euclidean") ## assays
    return(list(dcols=dcols, drows=drows, annotation=annotation, annt_colors=annt_colors, act=act, struct=struct, cv=cv))
    
  })

  
    select_plot <- reactive({
    showDendrogram <- input$showdendro
    profile_type <- input$proftype
    sort_meth <- input$sort_method
    fsize <- input$fontsize
      
    color <- wauc_colors
    breaks <- wauc_breaks
    leg_labels <- wauc_leg_labels
    leg_breaks <- wauc_leg_breaks
    
    if (profile_type == 'activity')
    {
      activity_type <- input$acttype
      if (activity_type != 'nwauc.logit')
      {
        color <- potency_colors
        breaks <- potency_breaks
        leg_labels <- potency_leg_labels
        leg_breaks <- potency_leg_breaks
      }
    }
    
    if (! is.null(chemical_loader()) )
    {
      # note pheatmap input has to have the same order!!!
      paras <- heatmap_para_generator()
      act <- paras[['act']]
      cv <- paras[['cv']]
      dcols <- paras[['dcols']]
      drows <- paras[['drows']]
      annotation <- paras[['annotation']]
      annt_colors <- paras[['annt_colors']]
      
      if (! showDendrogram)
      {
        if (profile_type == 'signal')
        {
          p <- pheatmap(t(act), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, clustering_distance_rows = drows, clustering_distance_cols = dcols, clustering_method = "average")
        } else if (sort_meth != 'toxscore')
        {
          p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
        } else
        {
          p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
        }
      } else if (sort_meth != 'toxscore' )
      {
        p <- plot(hclust(dcols, method="average"), hang=-1)
      }
    }
    return(p)
  })

  output$contents <- renderDataTable({
    if ( ! is.null(chemical_loader()) ) get_lookup_list(chemical_loader()[['id']], master)
  })

  output$assay_des <- renderDataTable({
    
#     paras <- heatmap_para_generator() #heatmap_para_generator
#     act <- paras[['act']]
#     annotation <- paras[['annotation']]
#     result <- get_output_df(act, annotation)
#     return(result)
    
    # for testing
     paras <- heatmap_para_generator()
     return(data.frame(paras))
  })
  

  getVarWidth <- reactive({
    if ( ! is.null(chemical_loader()) ) df <- get_lookup_list(chemical_loader()[['id']], master)
    ncmpd <- sum(rowSums(apply(df, 2, function(x) x == '' | is.na(x))) == 0)
    if (ncmpd < 40)
    {
      return(1200)
    } else
    {
      return(ncmpd*30)
    }
  })
  
  output$profiling <- renderPlot({
    select_plot()
 
  }, width=getVarWidth)
  
  output$box <- renderPlot({
    profile_type <- input$proftype
    fsize <- input$fontsize
    sort_meth <- input$sort_method
    
    p <- NULL
    if (profile_type == 'activity')
    {
      activity_type <- input$acttype
      if (activity_type == 'npod' | activity_type == 'nac50')
      {
        paras <- heatmap_para_generator()
        act <- paras[['act']]
        annotation <- paras[['annotation']]
        dcols <- paras[['dcols']]
        
        result <- get_output_df(act, annotation)
        p <- get_pod_boxplot(result, fontsize=fsize, sortby=sort_meth, dcols=dcols, global_para=assay_names)
      }
    }
    if (! is.null(p)) print(p)
    
  },  width=getVarWidth)

  output$downloadData <-  downloadHandler(
    filename = function() {
      if (input$proftype == 'profile')
      {
        paste(input$proftype, '_', input$sigtype, '.txt', sep='')
      } else
      {
        paste(input$proftype, '_', input$acttype, '.txt', sep='')
      }
       },
    content = function(file) {
      paras <- heatmap_para_generator()
      act <- paras[['act']]
      annotation <- paras[['annotation']]
      result <- get_output_df(act, annotation)
      write.table(result, file, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=FALSE)
    }
  )

output$downloadPlot <- downloadHandler(
       filename = function() { 
         if (input$proftype == 'profile')
        {
           paste(input$proftype, '_', input$sigtype, '.pdf', sep='')
         } else
         {
           paste(input$proftype, '_', input$acttype, '.pdf', sep='')
         } 
         },
       content = function(file) {
         #png(file, width=9, height=6.5, units="in", res=600)
         pdf(file, width=9, height=6.5)
         select_plot2()
         dev.off()
       }
)

select_plot2 <- function () {
  showDendrogram <- input$showdendro
  profile_type <- input$proftype
  sort_meth <- input$sort_method
  fsize <- input$fontsize
  color <- wauc_colors
  breaks <- wauc_breaks
  leg_labels <- wauc_leg_labels
  leg_breaks <- wauc_leg_breaks
  
  if (profile_type == 'activity')
  {
    activity_type <- input$acttype
    if (activity_type != 'nwauc.logit')
    {
      color <- potency_colors
      breaks <- potency_breaks
      leg_labels <- potency_leg_labels
      leg_breaks <- potency_leg_breaks
    }
  }
  
  if (! is.null(chemical_loader()) )
  {
    # note pheatmap input has to have the same order!!!
    paras <- heatmap_para_generator()
    act <- paras[['act']]
    cv <- paras[['cv']]
    dcols <- paras[['dcols']]
    drows <- paras[['drows']]
    annotation <- paras[['annotation']]
    annt_colors <- paras[['annt_colors']]
    
    if (! showDendrogram)
    {
      if (profile_type == 'signal')
      {
        p <- pheatmap(t(act), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, clustering_distance_rows = drows, clustering_distance_cols = dcols, clustering_method = "average")
      } else if (sort_meth != 'toxscore')
      {
        p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
      } else
      {
        p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
      }
    } else if (sort_meth != 'toxscore' )
    {
      p <- plot(hclust(dcols, method="average"), hang=-1)
    }
  }
  return(p)
}
  

    
  
  
  
})
