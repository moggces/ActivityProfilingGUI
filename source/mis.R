# to convert mapping from inp to out
conversion <- function (master, inp, out)
{
  result <- master[, out]
  names(result) <- master[, inp]
  return(result)
}

# split the master file into a list of matrices
split_master_2_matrix <- function(master, props, id='CAS')
{
  id_data <- master[, id]
  result <- lapply(as.list(props), function (x)
  {
    col_ids <- grepl(paste('^',x, sep=""), colnames(master))
    mat <- master[,col_ids]
    rownames(mat) <- id_data
    colnames(mat) <- sub(paste(x, '.', sep=""), "", colnames(mat))
    return(mat)
  }
  )
  names(result) <- props
  return(result)
}

# required: conversion, rename the GSID to Chemical.Name (row), rename assay to common_name (col, if TRUE)
rename_mat_col_row <- function (partial, master, assay_names, input_chemical_name=NULL, rename_assay=TRUE)
{
  chemical_name_ref <- input_chemical_name
  if (is.null(chemical_name_ref)) {
    chemical_name_ref <- conversion(master, inp='GSID', out='Chemical.Name')
  } 
  
  for (name in names(partial))
  {
    rownames(partial[[name]]) <-  chemical_name_ref[as.character(rownames(partial[[name]])) ]
    if (  name != 'struct' & rename_assay )
    {
      
      pathway_ref <- conversion(assay_names, inp='assay', out='common_name')
      colnames(partial[[name]]) <-  pathway_ref[as.character(colnames(partial[[name]])) ]
    }
  }
  return(partial)
}

# sort the matrices in a list based on rownames and columnames
sort_matrix <- function (partial)
{
  for (name2 in names(partial))
  {
    if (name2 == 'struct' )
    {
      partial[[name2]] <- partial[[name2]][order(rownames(partial[[name2]])),]
      
    } else
    {
      partial[[name2]] <- partial[[name2]][order(rownames(partial[[name2]])),]
      partial[[name2]] <- partial[[name2]][,order(colnames(partial[[name2]]))]
    }
  }
  return(partial)
}

# make increasing resp in mitotox as active,  dependent on mitotox name: inhibition_MMP
fix_mitotox_reverse <- function(partial, act_mat_names=c('npod', 'nec50', 'nwauc.logit'))
{
  for (name in act_mat_names)
  {
    mitotox_id <- grepl('inhibition_MMP', colnames(partial[[name]]))
    if (sum(mitotox_id) > 0)
    {
      rev_ids <- partial[[name]][, mitotox_id] < 0 &  ! is.na(partial[[name]][, mitotox_id])
      partial[[name]][rev_ids, mitotox_id] <- (partial[[name]][rev_ids, mitotox_id])*-1
    }
  }
  return(partial)
}

# use the information of assay_names (antagonism|inhibition)
filter_activity_by_type <- function(partial, type, thres=NULL, decision=FALSE, act_mat_names=c('npod', 'nec50', 'nwauc.logit'))
{
  if (type == 'nwauc.logit') 
  {
    # for codel 104
    if (thres == 0) thres <- 0.0001
  }
  for (name in act_mat_names)
  {
    ids <- matrix(FALSE, nrow(partial[[name]]), ncol(partial[[name]]))
    if (type == 'pod_med_diff')
    {
      ant_ids <- grepl('antagonism|inhibition', colnames(partial[[name]]))
      if (sum(ant_ids) > 0) 
      {
        ids <- matrix(FALSE, nrow(partial[[name]][, ant_ids]), ncol(partial[[name]][, ant_ids]))
        if (! is.null(thres) ) ids <- (partial[[type]][, ant_ids]*-1) < thres & ! is.na(partial[[type]][, ant_ids]) & ! is.na(partial[[name]][, ant_ids]) & partial[[name]][, ant_ids] > 0.0001
        if ( decision ) ids <- ! is.na(partial[[type]][, ant_ids]) & ! is.na(partial[[name]][, ant_ids]) & partial[[name]][, ant_ids] > 0.0001
        partial[[name]][, ant_ids][ids] <- (partial[[name]][, ant_ids][ids])*-1
      }
      
    } else
    {
      if (type == 'nwauc.logit' | type == 'npod' | type == 'nec50') ids <- partial[[type]] < thres & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'ncmax') ids <- abs(partial[[type]]) < thres & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'hitcall' & isTRUE(decision)) ids <- partial[[type]] != 1 & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'cc2' & isTRUE(decision)) ids <- ( abs(partial[[type]]) != 1.1 & abs(partial[[type]]) != 1.2 & abs(partial[[type]]) != 2.1) & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'cv.wauc' & isTRUE(decision)) ids <- partial[[type]] > 1.4 & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      
      partial[[name]][ids] <- (partial[[name]][ids])*-1
    }
  }
  return(partial)
}

# make < 0 as inconclusive (0.0001), make NA in potency as 0 if inactive, (dependent on nwauc.logit matrix)
assign_reverse_na_number <- function (partial, act_mat_names=c('npod', 'nec50', 'nwauc.logit'))
{
  result <- partial
  for (name in names(partial))
  {
    if (name == 'npod' | name == 'nec50') {
      result[[name]][ partial[['nwauc.logit']] == 0 & ! is.na(partial[['nwauc.logit']]) ] <- 0
      result[[name]][ result[[name]] < 0 |  is.na(result[[name]])   ] <- 0.0001
    }
    if (name == 'nwauc.logit' ) result[[name]][ partial[[name]] < 0 |  is.na(partial[[name]]) ] <- 0.0001
    if (name == 'wauc.logit') result[[name]][  is.na(partial[[name]]) ] <- 0.0001
  }
  
  return(result)
}

# remove inconclusive label (0.0001) but keep the untested as 0.0001
remove_inconclusive_label <- function (partial, act_mat_names=c('npod', 'nec50', 'nwauc.logit'))
{
  result <- partial
  for (name in names(partial))
  {
    if (name == 'npod' | name == 'nec50' | name == 'nwauc.logit' | name == 'wauc.logit') 
    {result[[name]][ partial[[name]] == 0.0001 & ! is.na(partial[['cc2']]) ] <- 0}
  }
  return(result)
}
