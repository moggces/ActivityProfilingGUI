
# shiny note: 1) it can't discriminate R vs. r in the r script file 
# 2) the renderTable has trouble with encoding issue (cannot recognize ppp file iconv -t UTF-8 -f ISO-8859-1)
# looks like the new read.table can automatically select best encoding
# 3) in shiny server, once you delete a file but replace a file with same name. somwhow don't know how to refresh its
# but if you update .R file, you can refresh to get new functions

# chemical_loader() out: list(id, ?(nwauc.logit or npod or nec50 or unknown))
# matrix_subsetter() out: list(activities or ?(nwauc.logit or npod or nec50 or unknown), struct)
# activity_filter() out: same as above 
# matrix_editor() out: list(nwauc.logit, npod, nec50,wauc.logit, struct, cv_mark, label)
# heatmap_para_generator() out: list(dcols, drows, annotation, annt_colors, act=act, struct, cv, label)
# tox21_data_generator()
# cas_data_generator()
# select_plot()

# todo:
# 1. download potency plot
# 2. broaden the "unknown" color scheme
# 6. filter by call meta

library(shiny)
library(plyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)
library(Cairo)
options(stringsAsFactors = FALSE)


#Sys.setlocale(locale="C")
#setwd("~/ShinyApps/profiling/")
source(paste(getwd(), "/source/customized.R", sep=""), local=TRUE)
#source(paste(getwd(), "/source/pheatmap_display_number.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/get.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/load.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/mis.R", sep=""), local=TRUE)
#environment(pheatmap_new_label) <- environment(pheatmap) pheatmap v. < 1.0

# load assay related parameters
logit_para_file <- './data/tox21_call_descriptions_v2.txt' #tox21_assay_collection.txt
assay_names <- load_profile(logit_para_file) # global, dataframe output

# load chemical information (will include purity later)
profile_file <- './data/tox21_compound_id_v5a7.txt' #colunm name has to be GSID # v5a3
master <- load_profile(profile_file) # global, dataframe output

# load the activities (all data) and the structure fp matrix
struct_mat_rdata <- './data/struct_mat.RData'
load(struct_mat_rdata, verbose=TRUE) # global, matrix output, struct_mat
activities <- readRDS('./data/activities_combined_170306.rds')

# remove the structures with low purity
#struct_mat <- struct_mat[rownames(struct_mat) %in% rownames(activities[[1]]),]
# very weird!! this line causes no error frozen on shiny public server

# heatmap settings
# the negative direction breaks won't capture wauc with very small values
wauc_breaks <- c( -1, -0.75, -0.5, -0.25, -0.1, -0.02, 0, 0.0001, 0.1, 0.25, 0.5, 0.75, 1) # upper is filled , lower is empty 
wauc_colors <-  c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE", "#D1E5F0", "#F7F7F7", "gray", "#FDDBC7" ,"#F4A582" ,"#D6604D" ,"#B2182B", "#67001F"  ) #RdBu
wauc_leg_breaks <- c(-1, -0.75, -0.5, -0.25,  0,   0.25, 0.5, 0.75, 1 )
wauc_leg_labels <- c("-1", "-0.75", "-0.5", "-0.25",  "0", "0.25", "0.5", "0.75", "1")

potency_breaks <- c(-0.02, 0, 0.0001, 4, 4.5, 5, 7.5, 9, 10)
potency_colors <- c("#F5F5F5", "gray", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30", "chartreuse") #BrBG
potency_leg_breaks <- c( 0,  4, 4.5, 5, 7.5, 9,10 )
potency_leg_labels <- c( "inactive",  "100uM", "30uM", "10uM", "0.3uM", "1nM", "0.1nM")

# potency_breaks <- c(-10, -9, -7.5, -5, -4.5, -4, -0.02, 0, 0.0001, 4, 4.5, 5, 7.5, 9, 10)
# potency_colors <- c("darkorange","#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "gray", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30", "chartreuse") #BrBG
# potency_leg_breaks <- c(-10, -9, -7.5, -5, -4.5, -4,  0,  4, 4.5, 5, 7.5, 9,10 )
# potency_leg_labels <- c("-10", "-9", "-7.5", "-5", "-4.5", "-4",  "0",  "4", "4.5", "5", "7.5", "9", "10")


shinyServer(function(input, output) {

# chemical_loader()
  chemical_loader <- reactive({
    
    result <- NULL
    path <- NULL
    
    # input file
    inFile <- input$file1
    # input textarea
    textdata <- input$cmpds
    
    if (! is.null(inFile)) { path <- inFile$datapath; filen <- inFile$name }
    if (textdata != '' ) result <- load_text_2_df(textdata)
    if (! is.null(path)) result <- load_data_matrix(path, filen) # as long as path or file has something it will override
   
    return(result)

  })

# matrix_subsetter()
  matrix_subsetter <- reactive({
    partial <- NULL
    reg_sel <- input$reg_sel # select the assays
    inv_sel <- input$inv_sel # inverse the selection
    nolowQC <- input$nolowQC # remove low QC 
    rename_assay <- FALSE # use the assay_names df
    
    # get all chemical information
    id_info <- chemical_loader()
    
    chem_id_df <- get_lookup_list(id_info[['id']], master)
    #ip <- subset(chem_id_df, ! is.na(StructureID), select=c(CAS, Cluster))
    
    # the basic identifies , GSID + Cluster
    ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
    
    # collect all the matrices and store in full (list)
    
    full <- list()
    full <- activities$cas_qc 
    if(! nolowQC) full <- activities$cas
    
    # if it is a data matrix input, only CAS ID is allowd
    input_chemical_name <- NULL
    if (length(id_info) > 1) # for loading the data matrix function
    {
      full <- id_info[! grepl('id', names(id_info))]
      chemical_name_ref <- conversion(master, inp='CAS', out='GSID')
      #rownames(full[[1]]) <- chemical_name_ref[as.character(rownames(full[[1]]))]
      avail_name <- chemical_name_ref[as.character(rownames(full[[1]]))]
      full[[1]] <- full[[1]][! is.na(avail_name), ]
      rownames(full[[1]]) <- avail_name[!is.na(avail_name)]
      
      if (! is.null(id_info[['id']]$input_Chemical.Name)) {
        input_chemical_name <- conversion(join(id_info[['id']], master), inp='GSID', out='input_Chemical.Name')
      }
      rename_assay <- FALSE
    }
    
    # the structure fingerprint matrix
    full[['struct']] <- struct_mat
    
    # subset the matrices by chemicals
    partial <- get_input_chemical_mat(ip, full)
    
    # rename the assays & chemicals
    partial <- rename_mat_col_row(partial,  master, assay_names, input_chemical_name, rename_chemical=TRUE, rename_assay=rename_assay)
    
    # subset the matrices by assay names
    partial <- get_assay_mat(partial, reg_sel, invSel=inv_sel)
    #print(partial[['npod']])
    # sort the matrix
    #partial <- sort_matrix(partial)
    
    return(partial)
  })

# activity_filter()  
  activity_filter <- reactive({
    
    # load all the activity filter parameters
    profile_type <- input$proftype
    activity_type <- input$acttype
    nwauc_thres <- input$nwauc_thres
    ncmax_thres <- input$ncmax_thres
    npod_thres <- ifelse(is.na(input$npod_thres), 3, log10(input$npod_thres/1000000)*-1)
    nec50_thres <- ifelse(is.na(input$nec50_thres), 3, log10(input$nec50_thres/1000000)*-1)
    #pod_diff_thres <- input$pod_diff_thres
    wauc_fold_thres <- input$wauc_fold_thres
    #isstrong <- input$isstrong
    nocyto <- input$nocyto
    isgoodcc2 <- input$isgoodcc2
    nohighcv <- input$nohighcv
    cytofilter <- input$cytofilter
    noauto <- input$noauto
    noch2issue <- input$noch2issue
    
    partial <- matrix_subsetter()
    
    # if it is data matrix input, don't change 
    if (length(partial) == 2) return(partial)
   
    act_mat_names <- c('npod', 'nec50', 'nwauc.logit')
    # reverse direction of mitotox could be meaningful
    #partial <- fix_mitotox_reverse(partial,act_mat_names=act_mat_names )
    
    # filtering
    partial <- filter_activity_by_type(partial, 'nwauc.logit', nwauc_thres, act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'ncmax', ncmax_thres,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'npod', npod_thres,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'nec50', nec50_thres,act_mat_names=act_mat_names)
      #partial <- filter_activity_by_type(partial, 'pod_med_diff', pod_diff_thres,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'label_cyto', thres=NULL, decision=cytofilter,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'wauc_fold_change', wauc_fold_thres,act_mat_names=act_mat_names)
      #partial <- filter_activity_by_type(partial, 'hitcall', thres=NULL, decision=isstrong,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'wauc_fold_change', thres=1, decision=nocyto,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'cc2', thres=NULL, decision=isgoodcc2,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'label_autof', thres=NULL, decision=noauto,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'label_ch2', thres=NULL, decision=noch2issue,act_mat_names=act_mat_names)
    
    # it has to be the end
    partial <- filter_activity_by_type(partial, 'cv.wauc', thres=NULL, decision=nohighcv,act_mat_names=act_mat_names)
    #print(partial[['npod']])
    
    return(partial)
  })

# matrix_editor()
  matrix_editor <- reactive({
    
    noincon_label <- input$noinconlab #inconclusive label
    act_mat_names <- c('npod', 'nec50', 'nwauc.logit')
    
    partial <- activity_filter()
    #print(partial[['npod']])
    # if it is data matrix input, skip 
    if (length(partial) == 2) return(partial)
    
    # create CV marks
    cv_mark <- get_cv_mark_mat(partial[['cv.wauc']], partial[['nwauc.logit']])
    partial[['cv_mark']] <- cv_mark
    
    # make activities matrix (< 0 and NA) as 0.0001
    partial <- assign_reverse_na_number(partial, act_mat_names=act_mat_names)
    #print(partial[['npod']])
    
    # remove inconclusive label (0.0001 as 0 ) (but keep the untested ones = 0.0001)
    if (noincon_label) partial <- remove_inconclusive_label(partial, act_mat_names=act_mat_names)

    acts <- partial[c( act_mat_names, 'wauc.logit', 'struct', 'cv_mark', 'label')]
    #print(partial[['npod']])
    return(acts)
  })
  
#heatmap_para_generator()  
  heatmap_para_generator <- reactive({
    sort_meth <- input$sort_method
    profile_type <- input$proftype
    activity_type <- ''
    
    # get all chemical information
    input_chemical_name <- NULL
    chem_id_df <- get_lookup_list(chemical_loader()[['id']], master)
    if (! is.null(chem_id_df$input_Chemical.Name)) {
      input_chemical_name <- conversion(chem_id_df, inp='Chemical.Name', out='input_Chemical.Name')
    }
    # the basic identifies , GSID + Cluster
    # can add the Chemical.Name here
    ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster,Chemical.Name))
    
    # the cleaned matrices
    dt <- matrix_editor() 
    if (is.null(dt)) return(NULL)
    
    # if the input is data matrix, creat a blank CV matrix
    if (length(dt) == 2 ) 
    {
      activity_type <- names(dt)[1]
      act <- dt[[1]]
      cv <-  matrix("", nrow(act), ncol(act), dimnames=dimnames(act))
      label <- matrix("", nrow(act), ncol(act), dimnames=dimnames(act))
    } else
    {
      # it has to be here to add more lines for the duplicates
      dt <- duplicate_chemical_row(dt, ip)
      
      if (profile_type == 'activity')
      {
        activity_type <- input$acttype
        act <- dt[[activity_type]]
        
      } else
      {
        act <- dt[['wauc.logit']]
      }
      
      cv <- dt[['cv_mark']]
      label <- dt[['label']]
      
    }
    
    # struct matrix
    struct <- dt[['struct']]
    # first, cluster the chemicals
    #print(str_c("line271", rownames(struct)))
    dcols <- dist(struct, method = "binary") ## chemicals
    
    # very, very cumbersome functions. better to split, merge dt + activity_type
    annotation <- get_heatmap_annotation(dcols, ip, master, input_chemical_name=input_chemical_name, dmat=dt, actType=activity_type) #data.frame output
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
      label <- label[tox_order, ]
    }
    # cluster assays by similarity 
    drows <- dist(t(act) , method = "euclidean") ## assays
    return(list(dcols=dcols, drows=drows, annotation=annotation, annt_colors=annt_colors, act=act, struct=struct, cv=cv, label=label))
    
  })

  chemical_enricher <- reactive({
    paras <- heatmap_para_generator()
    if (is.null(paras)) return(NULL)
    
    # chemical information
    chem_id_df <- get_lookup_list(chemical_loader()[['id']], master)
    ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster,Chemical.Name))
    
    # parameters
    reg_sel <- input$reg_sel # select the assays
    inv_sel <- input$inv_sel # inverse the selection
    nolowQC <- input$nolowQC # remove the low QC
    rename_assay <- FALSE # use the assay_names df
    profile_type <- input$proftype
    activity_type <- input$acttype
    act_mat_names <- activity_type
    if (profile_type != 'activity') return(NULL)
    
    # get the partial matrix
    partial <- activity_filter()
    # if it is data matrix input, skip 
    if (length(partial) == 2) return(NULL)
    #filtered activies < 0, active >0, inactive =0 or inconclusive in the beginning, NA non tested
    partial[[act_mat_names]][ (is.na(partial[[act_mat_names]]) | partial[[act_mat_names]] == 0.0001)   & ! is.na(partial[['cc2']]) ] <- 0
    
    # add duplicate rows due to duplicate cluster information
    partial <- duplicate_chemical_row(partial, ip)
    #print(str_c("line324", rownames(partial[[act_mat_names]])))
    
    
    # load all the activity filter parameters
    nwauc_thres <- input$nwauc_thres
    ncmax_thres <- input$ncmax_thres
    npod_thres <- ifelse(is.na(input$npod_thres), 3, log10(input$npod_thres/1000000)*-1)
    nec50_thres <- ifelse(is.na(input$nec50_thres), 3, log10(input$nec50_thres/1000000)*-1)
    #pod_diff_thres <- input$pod_diff_thres
    #isstrong <- input$isstrong
    nocyto <- input$nocyto
    isgoodcc2 <- input$isgoodcc2
    nohighcv <- input$nohighcv
    cytofilter <- input$cytofilter
    wauc_fold_thres <- input$wauc_fold_thres
    noauto <- input$noauto
    noch2issue <- input$noch2issue
    
    
    
    full <- activities$cas_qc
    if (! nolowQC) full <- activities$cas
    # subset the matrices by assay names
    
    # rename the assays & chemicals
    full <- rename_mat_col_row(full,  master, assay_names, input_chemical_name=NULL, rename_chemical=FALSE, rename_assay=rename_assay)
    # subset the matrices by assay names
    full <- get_assay_mat(full, reg_sel, invSel=inv_sel)
    
    # filtering
    full <- filter_activity_by_type(full, 'nwauc.logit', nwauc_thres, act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'ncmax', ncmax_thres,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'npod', npod_thres,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'nec50', nec50_thres,act_mat_names=act_mat_names)
    #full <- filter_activity_by_type(full, 'pod_med_diff', pod_diff_thres,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'label_cyto', thres=NULL, decision=cytofilter,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'wauc_fold_change', wauc_fold_thres, act_mat_names=act_mat_names)
    #full <- filter_activity_by_type(full, 'hitcall', thres=NULL, decision=isstrong,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'wauc_fold_change', thres=1, decision=nocyto,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'cc2', thres=NULL, decision=isgoodcc2,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'label_autof', thres=NULL, decision=noauto,act_mat_names=act_mat_names)
    full <- filter_activity_by_type(full, 'label_ch2', thres=NULL, decision=noch2issue,act_mat_names=act_mat_names)
    
    
    # it has to be the end
    full <- filter_activity_by_type(full, 'cv.wauc', thres=NULL, decision=nohighcv,act_mat_names=act_mat_names)
    
    #filtered activies < 0, active >0, inactive =0 or inconclusive in the beginning, NA non tested
    full[[act_mat_names]][ (is.na(full[[act_mat_names]]) | full[[act_mat_names]] == 0.0001)   & ! is.na(full[['cc2']]) ] <- 0
    
    #print(paras[['annotation']])
    #print(rownames(paras[['annotation']]))
    result <- get_clust_assay_enrichment(partial[[act_mat_names]], full[[act_mat_names]], paras[['annotation']], calZscore=FALSE)
    
    return(result)
    
  })
  
    select_plot <- reactive({
    showDendrogram <- input$showdendro
    keepsize <- input$keepsize
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
          #pheatmap v. < 1.0
          #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
          p <- pheatmap(t(act),  fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
        } else
        {
          #pheatmap v. < 1.0
          #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
          p <- pheatmap(t(act),  fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
        }
      } else if (sort_meth != 'toxscore' )
      {
        p <- plot(hclust(dcols, method="average"), hang=-1)
      }
    }
    return(p)
  })

  tox21id_data_generator <- reactive({
    paras <- heatmap_para_generator() #heatmap_para_generator
    actf <- paras[['act']]
    id_info <- chemical_loader()
    id_data <- master
    isUpload <- FALSE
    if(length(id_info) > 1) {
      id_data <- id_info[['id']]
      isUpload <- TRUE
    }
    if (! isUpload)
    {
      result <- get_source_data_long(source_acts=activities$tox21agencyid, chem_id_master=master, filtered_act=actf)
    } else {result <- NULL}
    return(result)
  })
  
  cas_data_generator <- reactive({
    actwithflag <- input$actwithflag
    paras <- heatmap_para_generator() #heatmap_para_generator
    
    id_info <- chemical_loader()
    id_data <- master
    isUpload <- FALSE
    if(length(id_info) > 1) {
      id_data <- id_info[['id']]
      isUpload <- TRUE
    }
    result <- get_output_df(paras, id_data, isUpload=isUpload, actwithflag=actwithflag)
    return(result)
  })
    
  output$contents <- renderDataTable({
    if ( ! is.null(chemical_loader()) ) get_lookup_list(chemical_loader()[['id']], master)
  })

  output$casdata <- renderDataTable({
    
     #return(matrix_subsetter()[['nwauc.logit']])
    return(cas_data_generator())
    # for testing
#       paras <- heatmap_para_generator()
#       return(data.frame(rownames(paras[['act']])))
    
  })
  
  output$tox21iddata <- renderDataTable({

    return(tox21id_data_generator())
  })
  
  output$enrich <- renderDataTable({
    #return(as.data.frame(chemical_enricher()[['modl_acc']]))
    return(chemical_enricher())
  })
  
  output$assay_info <- renderDataTable({
  
    #col_n <- c('common_name','technology','cell_type','species','abbreviation', 'PubChem AID')
    #result <- assay_names[, colnames(assay_names) %in% col_n]
    partial <- matrix_subsetter()
    not_want <- c('_for_FDA_A_name', '_target_type_gene_go.biological.process',	
                  '_target_type_gene_ctd.disease', '_technology_long.description',
                  '_technology_short.description','protocol_call_db.name_parent',
                  'protocol_call_db.name_readout_primary','protocol_CEBS.batch',
                  'protocol_call_db.name_readout_secondary',
                  'protocol_db.name','protocol_time_release',
                  'protocol_slp','protocol_description')
    result <- assay_names[, ! colnames(assay_names) %in% not_want]
    result <- result %>%
        filter(protocol_call_db.name != '')  %>% #the ones with call definition
        filter(protocol_call_db.name %in% colnames(partial[['npod']])) %>%
        #select(noquote(order(colnames(.)))) #reorder the columns alphabetically
        select(protocol_call_db.name, protocol_call_db.name_display.name, 
               starts_with("target"), starts_with("technology"), starts_with("format"),
               starts_with("provider"), starts_with("protocol"))
               
    return(result)
    
  })

  getVarWidth <- reactive({
    ncmpd <- 0
    keepsize <- input$keepsize
    if ( ! is.null(chemical_loader()) & ! keepsize) 
    {
      chem_id_df <- get_lookup_list(chemical_loader()[['id']], master)
      ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
      ncmpd <- nrow(ip)
    }
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
      if (activity_type == 'npod' | activity_type == 'nec50')
      {
        paras <- heatmap_para_generator()
        act <- paras[['act']]
        annotation <- paras[['annotation']]
        dcols <- paras[['dcols']]
        
        id_info <- chemical_loader()
        id_data <- master
        isUpload <- FALSE
        if(length(id_info) > 1) {
          id_data <- id_info[['id']]
          isUpload <- TRUE
        }
        result <- get_output_df(paras, id_data, isUpload,actwithflag=FALSE)
        result <- select(result, -Chemical.Name_original) # remove the new added column after get_output_df
        p <- get_pod_boxplot(result, fontsize=fsize, sortby=sort_meth, dcols=dcols, global_para=assay_names)
      }
    }
    if (! is.null(p)) print(p)
    
  },  width=getVarWidth)

  output$downloadCASData <-  downloadHandler(
    
    filename = function() {
      if (input$proftype == 'signal')
      {
        paste(input$proftype, '_', input$sigtype, '.txt', sep='')
      } else
      {
        paste(input$proftype, '_', input$acttype, '.txt', sep='')
      }
       },
    content = function(file) {
      result <- cas_data_generator()
      #result <- get_published_data_only_commonname(result, assay_names)  # to remove unpublished data
      write.table(result, file, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=FALSE)
    }
  )
  
  output$downloadTox21IDData <-  downloadHandler(
    
    filename = function() {
      paste(as.numeric(as.POSIXct(Sys.time())), ".txt", sep="")
    },
    content = function(file) {
      result <- tox21id_data_generator()
      write.table(result, file, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=FALSE)
    }
  )
  
  output$downloadEnrich <-  downloadHandler(
    filename = function() {
    
      paste(input$proftype, '_', input$acttype, '_enrichment.txt', sep='')
    },
    content = function(file) {
      result <- chemical_enricher()
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
  keepsize <- input$keepsize
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
        #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
        p <- pheatmap(t(act),  fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
      } else
      {
        #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
        p <- pheatmap(t(act), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
      }
    } else if (sort_meth != 'toxscore' )
    {
      p <- plot(hclust(dcols, method="average"), hang=-1)
    }
  }
  return(p)
}
  

    
  
  
  
})
