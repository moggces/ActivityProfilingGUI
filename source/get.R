# specific to activity file generated from KNIME
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

# a wrapper for join, it can detect, CAS, GSID automatically
get_lookup_list <- function (input, master)
{
  #result <- subset(master, select=c(CAS, Chemical.Name, StructureID))
  #result <- merge(input,result, by='CAS', all.x=TRUE)
  result <- join(input, master)
  return(result)
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

# its linked with nwauc.logit matrix results. if active and high CV -> mark
get_cv_mark_mat <- function(cv, nwauc)
{
  cv_mark <- cv
  cv_mark[is.na(cv_mark) ] <- -0.1
  cv_mark[cv_mark > 1.4  & nwauc > 0.0001 ] <- "#"
  cv_mark[cv_mark != "#"] <- ''
  return(cv_mark)
  
}

# dependent on conversion
# d: is the distance matrix
# input: chemical identification (GSID + Cluster)
# master: all mapping info
# dmat

get_heatmap_annotation <- function (d, input, master, input_chemical_name=NULL, cutoff=0.7, method="average", dmat, actType='')
{
  chemical_name_ref <- input_chemical_name
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
    if (is.null(chemical_name_ref)) chemical_name_ref <- conversion(master, inp='GSID', out='Chemical.Name')
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
  } else if (actType == 'npod' | actType == 'nec50' )
  {
    if ( ! is.null(dmat[['nwauc.logit']]))
    {
      annotation3 <- data.frame(toxScore = unlist(lapply(1:nrow(dmat[[actType]]), function (x) sum(abs(dmat[[actType]][x,])*dmat[['nwauc.logit']][x,]) )))
    } else
    {
      annotation3 <- data.frame(toxScore = rowSums(abs(dmat[[actType]]) ))
    }
    
  }
  
  if (nrow(annotation3) > 0)
  {
    rownames(annotation3) <- rownames(dmat[[1]])
    annotation <- merge(annotation, annotation3, by="row.names")
    rownames(annotation) <- annotation$Row.names
    annotation <- annotation[,-which(colnames(annotation) == 'Row.names')]
  }
  
  return(annotation)
}

# rainbow color to generate unique colors
# toxScore is a continuous color
get_heatmap_annotation_color <- function(annotation, actType='')
{
  user <- rainbow(length(unique(annotation[['userClust']])))
  names(user) <- levels(annotation[['userClust']])
  chem  <- rainbow(length(unique(annotation[['chemClust']])))
  names(chem) <- levels(annotation[['chemClust']])
  
  if (actType != '')
    #if (! is.null(actType))
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
  let <- conversion(global_para, inp='common_name', out='letters')
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
    #theme_bw(base_size = fontsize) + 
    annotation_logticks(sides = "l") 
  return(p)
}




