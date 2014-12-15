#output: list
get_mat <- function (master_df, para, select=c('nwauc', 'call', 'pod', 'ac50'), removeAbnormalDirection=TRUE, isSignal=FALSE)
{
  nwauc_mat <- NULL
  call_mat <- NULL
  pod_mat <- NULL
  ac_mat <- NULL
  
  if ('nwauc' %in% select )
  {
    if (isSignal) 
    {
      select_opt <- 'wauc'
    } else
    {
      select_opt <- 'nwauc'
    }
    
    temp_mat <- master_df[, c('CAS', grep(paste("_", select_opt, sep=""), colnames(master_df), value=TRUE))]
    if (removeAbnormalDirection)
    {
      temp_mat[,-c(1)][temp_mat[,-c(1)] < 0 & ! is.na(temp_mat[,-c(1)])] <- NA
    }
    temp_mat <- melt(temp_mat,  id.vars=c('CAS'), variable.name = 'pathway', value.name = 'nwauc_4toxpi')
    temp_mat$inv_logit <- NA
    temp_mat$inv_logit <- unlist(lapply(1:length(temp_mat$nwauc_4toxpi), function(x) convert_wauc_2logit(select=select_opt, temp_mat$nwauc_4toxpi[x], temp_mat$pathway[x], para))	)
    
    #### sanity check: plotting
    #temp_mat_2 <- merge(temp_mat, para, by.x="pathway", by.y="nwauc")
    #temp_mat_2 <- subset(temp_mat_2, ! is.na(nwauc_4toxpi) )
    #temp_mat_2$activity <- 'inactive'
    #temp_mat_2[abs(temp_mat_2$nwauc_4toxpi) > temp_mat_2$low  , ]$activity <- 'marginal active'
    #temp_mat_2[abs(temp_mat_2$nwauc_4toxpi) > temp_mat_2$high, ]$activity <- 'active'
    
    #dev.new()
    #p <- ggplot(temp_mat_2, aes(x=nwauc_4toxpi, y=inv_logit, color=activity))
    #p + geom_point(alpha=0.5)  + facet_wrap(~ pathway_name, scales = "free_x")  + scale_x_continuous("raw wAUC", limits = c(0, 100)) + scale_y_continuous("normalized wAUC")
    
    # check the distribution
    temp_mat <- dcast(temp_mat, CAS ~ pathway, value.var='inv_logit')
    nwauc_mat <- temp_mat[, -c(1)]
    #mat[is.na(temp_mat) ] <- -0.1
    rownames(nwauc_mat) <- temp_mat$CAS
    nwauc_mat <- nwauc_mat[order(rownames(nwauc_mat)),]
    nwauc_mat <- nwauc_mat[,order(colnames(nwauc_mat))]
  }
  
  if ( 'call' %in%  select )
  {
    
  }
  
  if ( 'pod' %in% select )
  {
    temp_mat <- master_df[, c('CAS', grep(paste("_", 'pod_4toxpi', sep=""), colnames(master_df), value=TRUE))]
    if (removeAbnormalDirection)
    {
      temp_mat[,-c(1)][temp_mat[,-c(1)] < 0 & ! is.na(temp_mat[,-c(1)])] <- NA
    }
    pod_mat <- temp_mat[, -c(1)]
    rownames(pod_mat) <- temp_mat$CAS
    pod_mat <- pod_mat[order(rownames(pod_mat)),]
    pod_mat <- pod_mat[,order(colnames(pod_mat))]
  }
  
  if (  'ac50' %in% select)
  {
    temp_mat <- master_df[, c('CAS', grep(paste("_", 'ac50_4toxpi', sep=""), colnames(master_df), value=TRUE))]
    if (removeAbnormalDirection)
    {
      temp_mat[,-c(1)][temp_mat[,-c(1)] < 0 & ! is.na(temp_mat[,-c(1)])] <- NA
    }
    ac_mat <- temp_mat[, -c(1)]
    rownames(ac_mat) <- temp_mat$CAS
    ac_mat <- ac_mat[order(rownames(ac_mat)),]
    ac_mat <- ac_mat[,order(colnames(ac_mat))]
  }
  
  return(list(act=nwauc_mat, call=call_mat, pod=pod_mat, ac50=ac_mat))
}

# filter the matrix by chemical input 
get_input_chemical_mat <- function (input, full)
{
  partial <- list()
  for (name in names(full))
  {
    partial[[name]] <- full[[name]][as.character(rownames(full[[name]])) %in% as.character(input[['GSID']]),] # CAS here
    
  }
  
  return(partial)
}

# filter the matrix by assays regular expression
get_assay_mat <- function (partial, regSel, invSel=FALSE)
{
  for (name in names(partial))
  {
    if (name != 'struct' ) 
    {
      partial[[name]] <- partial[[name]][,grep(regSel, colnames(partial[[name]]), value = TRUE, invert = invSel)]
    }
  }
  return(partial)
}



get_lookup_list <- function (input, master)
{
  #result <- subset(master, select=c(CAS, Chemical.Name, StructureID))
  #result <- merge(input,result, by='CAS', all.x=TRUE)
  result <- join(input, master)
  return(result)
}


get_heatmap_annotation <- function (d, input, master, cutoff=0.7, method="average", dmat, actType='')
{
  # chemical structure clustering
  hc <- hclust(d, method=method)
  group <- cutree(hc, h=cutoff)
  group_cp <- group
  group_t <- sort(table(group), decreasing=TRUE)
  
  for (i in 1:length(group_t))
  {
    if (group_t[i] == 1)
    {
      group_cp[group == names(group_t)[i]] <- 0
    } else
    {
      group_cp[group == names(group_t)[i]] <- i
    }
  }
  
  # create annotations: chemClust
  annotation <- data.frame(chemClust = as.factor(group_cp))
  rownames(annotation) <- names(group_cp)
  
  # create annoations: userClust
  annotation2 <- data.frame(userClust = as.factor(input[['Cluster']]))
  
  if (nrow(annotation2) > 0)
  {
    rownames(annotation2) <- as.character(input[['GSID']])
    chemical_name_ref <- conversion(master, inp='GSID', out='Chemical.Name')
    rownames(annotation2) <- chemical_name_ref[as.character(rownames(annotation2))]
    
    annotation <- merge(annotation, annotation2, by="row.names")
    rownames(annotation) <- annotation$Row.names
    annotation <- annotation[,-which(colnames(annotation) == 'Row.names')]
  }
  
  # create annotations: toxScore
  annotation3 <- data.frame()
  
  if (actType == 'nwauc.logit' )
  {
    annotation3 <- data.frame(toxScore = rowSums(abs(dmat[[actType]]) ))                           
  } else if (actType == 'npod' | actType == 'nac50' )
  {
    annotation3 <- data.frame(toxScore = unlist(lapply(1:nrow(dmat[[actType]]), function (x) sum(abs(dmat[[actType]][x,])*dmat[['nwauc']][x,]) )))
  }
  
  if (nrow(annotation3) > 0)
  {
    rownames(annotation3) <- rownames(dmat[['nwauc.logit']])
    annotation <- merge(annotation, annotation3, by="row.names")
    rownames(annotation) <- annotation$Row.names
    annotation <- annotation[,-which(colnames(annotation) == 'Row.names')]
  }
  
  return(annotation)
}

get_heatmap_annotation_color <- function(annotation, actType=NULL)
{
  user <- rainbow(length(unique(annotation[['userClust']])))
  names(user) <- levels(annotation[['userClust']])
  chem  <- rainbow(length(unique(annotation[['chemClust']])))
  names(chem) <- levels(annotation[['chemClust']])
  
  #if (actType != '')
  if (! is.null(actType))
  {
    tox <-  c("#F7F4F9", "#E7E1EF", "#D4B9DA", "#C994C7", "#DF65B0", "#E7298A", "#CE1256", "#980043", "#67001F") #PuRd
    return(list(userClust=user, chemClust=chem, toxScore=tox))
   
  } else
  {
    return(list(userClust=user, chemClust=chem))
  }
  
}

get_output_df <- function (act, annotation)
{
  
  act$Chemical.Name <- rownames(act)
  annotation$Chemical.Name <- rownames(annotation)
  result <- join(annotation, act)
  result <- join(subset(master, select=c(CAS, Chemical.Name)), result, type = "inner")
  return(result)
}

get_pod_boxplot <- function (pod, fontsize, sortby, dcols, global_para)
{
  # order the chemical.name 
  h <- hclust(dcols, method='average')
  pod[, 'Chemical.Name'] <- ordered(pod$Chemical.Name, levels=h$label[h$order])
 
  if (sortby == 'toxscore') pod[, 'Chemical.Name'] <- ordered(pod$Chemical.Name, levels=pod$Chemical.Name[order(pod$toxScore)])
  
  # melt the data and exclude the all inactives
  pod_m <-  melt(pod, id.vars = c( 'CAS', 'Chemical.Name', 'chemClust', 'userClust', 'toxScore'), value.name = "pod_value", variable.name = 'pathway')
  pod_m <- subset(pod_m, pod_value > 1)  # Chemical.Name is a factor. So if completely inactve. it won't be removed
  
  mat <- pod_m
  
  #create conversion
  let <- conversion(global_para, inp='pathway_name', out='pathway_abb')
  let2 <- paste(let, names(let), sep="=") # color legend
  names(let2) <- names(let)

  #add a new column
  mat[, 'path_abb'] <- let[as.character(mat$pathway)]

  p <- ggplot(data=mat, aes(x=Chemical.Name, y=pod_value*-1+6)) + 
    geom_boxplot(outlier.size = 0) +
    geom_text(aes(label=path_abb, color=pathway), size=7, alpha=0.7, position="jitter") + 
    scale_color_discrete("",labels=let2) + 
    scale_x_discrete("", drop=FALSE) + # keep the no activity ones
     theme(text=element_text(size=fontsize), 
           axis.text.x = element_text( angle=90, color="black")) + 
    scale_y_continuous('uM', breaks=seq(-10+6, -3+6, by=1), limits=c(-10+6, -3+6), labels = math_format(10^.x)) + 
    annotation_logticks(sides = "l") 
  return(p)
}

get_property_name <- function (master)
{
  col_list <- strsplit(colnames(master), '.', fixed=TRUE)
  #unique(unlist(lapply(col_list, function (x) x[[length(x)]]))) # get the unique assay name
  names <- lapply(col_list, function (x)
  {
    
    if (length(x) == 3)
    {
      return(paste(x[[1]], '.', x[[2]], sep=""))
    } else {return(x[[1]])}
    
  }
  )
  names <- unique(unlist(names))
  return(names)
}

get_cv_mark_mat <- function(cv)
{
  cv_mark <- cv
  cv_mark <- ""
  cv_mark[cv > 1.4] <- "#"
  return(cv_mark)
}
