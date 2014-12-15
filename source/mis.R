conversion <- function (master, inp, out)
{
  result <- master[, out]
  names(result) <- master[, inp]
  return(result)
}

convert_wauc_2logit <- function (select=c('nwauc', 'wauc'), nwauc, pathway, para)
{
  
  geom_mean <- conversion(para, inp=select, out="geom_mean")
  inactive_inv_logit <- conversion(para, inp=select, out="inactive_inv_logit")
  
  result <- NA
  med_point <- geom_mean[as.character(pathway)]
  if (nwauc >= 0 & ! is.na(nwauc))
  {
    result <-  ( (exp( (nwauc - med_point)/med_point) )/( med_point + (exp( ( nwauc - med_point )/med_point)) ) )  - inactive_inv_logit[as.character(pathway)]
  } else if (nwauc < 0 & ! is.na(nwauc))
  {
    nwauc <- nwauc*-1
    result <- ( (exp( (nwauc - med_point)/med_point) )/( med_point + (exp( ( nwauc - med_point )/med_point)) ) ) - inactive_inv_logit[as.character(pathway)]
    result <- result*-1
  }
  return (result)
}

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


filter_activity_by_type <- function(partial, type, thres=NULL, decision=FALSE, act_mat_names=c('npod', 'nac50', 'nwauc.logit'))
{
  if (type == 'nwauc.logit') 
  {
    if (thres == 0) thres <- 0.0001
  }
  for (name in act_mat_names)
  {
    ids <- matrix(FALSE, nrow(partial[[name]]), ncol(partial[[name]]))
    if (type == 'pod_med_diff')
    {
      ant_ids <- grepl('antagonism|inhibition', colnames(partial[[name]]))
      ids <- (partial[[type]][, ant_ids]*-1) < thres & ! is.na(partial[[type]][, ant_ids]) & ! is.na(partial[[name]][, ant_ids]) & partial[[name]][, ant_ids] > 0.0001
      partial[[name]][, ant_ids][ids] <- (partial[[name]][, ant_ids][ids])*-1
      
    } else
    {
      if (type == 'nwauc.logit' | type == 'npod' | type == 'nac50') ids <- partial[[type]] < thres & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'nemax') ids <- abs(partial[[type]]) < thres & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'hitcall' & isTRUE(decision)) ids <- partial[[type]] == 1 & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'cc2' & isTRUE(decision)) ids <- ( abs(partial[[type]] == 1.1) | abs(partial[[type]] == 1.2) | abs(partial[[type]] == 2.1)) & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      partial[[name]][ids] <- (partial[[name]][ids])*-1
    }
  }
  return(partial)
}

fix_mitotox_reverse <- function(partial, act_mat_names=c('npod', 'nac50', 'nwauc.logit'))
{
  for (name in act_mat_names)
  {
    mitotox_id <- grepl('mitotox', colnames(partial[[name]]))
    if (sum(mitotox_id) > 0)
    {
      rev_ids <- partial[[name]][, mitotox_id] < 0 &  ! is.na(partial[[name]][, mitotox_id])
      partial[[name]][rev_ids, mitotox_id] <- (partial[[name]][rev_ids, mitotox_id])*-1
    }
  }
  return(partial)
}

assign_reverse_na_number <- function (partial, act_mat_names=c('npod', 'nac50', 'nwauc.logit'))
{
  for (name in act_mat_names)
  {
    partial[[name]][partial[[name]] < 0 | is.na(partial[[name]]) ] <- 0.0001
  }
  return(partial)
}

rename_mat_col_row <- function (partial, master, assay_names)
{
  chemical_name_ref <- conversion(master, inp='GSID', out='Chemical.Name')
  
  for (name in names(partial))
  {
    rownames(partial[[name]]) <-  chemical_name_ref[as.character(rownames(partial[[name]])) ]
    if (  name != 'struct')
    {
      
      pathway_ref <- conversion(assay_names, inp='assay', out='common_name')
      colnames(partial[[name]]) <-  pathway_ref[as.character(colnames(partial[[name]])) ]
    }
  }
  return(partial)
}

edit_mat_manual <- function (partial, nwaucThres=0.0001, actType='', regSel='', invSel=FALSE)
{
  
  partial[['nwauc']][is.na(partial[['nwauc']]) ] <- 0.0001  # for plotting
  
  # cv: na -> -0.1 , nwauc < 0.1 (change to 0.05) -> -0.1 , cv > 1.4 -> # , clean CV
  # filter the pathways
  if ( ! is.null( partial[['cv']] ))
  {
    partial[['cv']][is.na(partial[['cv']])] <- -0.1 # just to remove the na 
    partial[['cv']][partial[['nwauc']] < 0.1  ] <- -0.1 # if two small just let it be
    partial[['cv']][partial[['cv']] > 1.4] <- "#"
    partial[['cv']][partial[['cv']] != "*"] <- ""
    ################################ temp for paper
    #partial[['nwauc']] <- partial[['nwauc']][,grep("are|hse|fxr_|pparg_antagonism|ppard_", colnames(partial[['nwauc']]), value = TRUE, invert = TRUE)]
    #partial[['cv']] <- partial[['cv']][,grep("are|hse|fxr_|pparg_antagonism|ppard_", colnames(partial[['cv']]), value = TRUE, invert = TRUE)]
    ###########################################
    #if (removeCytotoxic) partial[['cv']] <- partial[['cv']][,grep("cytotoxicity", colnames(partial[['cv']]), value = TRUE, invert = TRUE)]
    partial[['cv']] <- partial[['cv']][,grep(regSel, colnames(partial[['cv']]), value = TRUE, invert = invSel)]
  }
  
  # nwauc and mitotox 
  for (name in names(partial))
  {
    if (name != 'cv' & name != 'struct' ) 
    {
      partial[[name]] <- partial[[name]][,grep(regSel, colnames(partial[[name]]), value = TRUE, invert = invSel)]
      partial[[name]][is.na(partial[[name]]) ] <- 0.0001  # for plotting
      if (actType == 'nwauc' & name == 'nwauc' )
      {
        if (sum(colnames(partial[['nwauc']]) ==  'mitotox') > 0)
        {
          partial[['nwauc']][, -which(colnames(partial[['nwauc']]) ==  'mitotox')][partial[['nwauc']][, -which(colnames(partial[['nwauc']]) ==  'mitotox')] < 0 ] <- 0.0001
        }
      } else if (actType != '')
      {
        partial[[name]][partial[[name]] < 0 ] <- 0.0001
      }
    }
  }
  
  # the dt40 assays issue
  if (actType != '')
  {
    for (name in c('pod', 'ac50'))
    {
      x <- TRUE
      if (sum(colnames(partial[['nwauc']]) == 'dna_damage(dsb)') > 0)
      {
        partial[[name]][, 'dna_damage(dsb)'] <- partial[['nwauc']][, 'dna_damage(dsb)']
        partial[[name]][, 'dna_damage(dsb)'][partial[[name]][, 'dna_damage(dsb)'] > 0 ] <- 5
      }
      #if (grepl('dna_damage(srf)', colnames(partial[[name]]), fixed=TRUE))
      if (sum(colnames(partial[['nwauc']]) ==  'dna_damage(srf)') > 0)
      {
        partial[[name]][, 'dna_damage(srf)'] <- partial[['nwauc']][, 'dna_damage(srf)']
        partial[[name]][, 'dna_damage(srf)'][partial[[name]][, 'dna_damage(srf)'] > 0 ] <- 5
      }
      
    }
    
    #sort
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
    
    # for nwauc_threshold , activity filter
    for (name in c('pod', 'ac50'))
    {
      partial[[name]][partial[['nwauc']] < nwaucThres & partial[[name]]  != 0.0001 ] <- (partial[[name]][partial[['nwauc']] < nwaucThres & partial[[name]]  != 0.0001 ])*-1
    }
    
  }
  
  
  return(partial)
}

# by using the conversion function to rename the rowname of the matrix
rename_mat <- function (partial, master, para, actType='')
{
  chemical_name_ref <- conversion(master, inp='CAS', out='Chemical.Name')
  
  for (name in names(partial))
  {
    rownames(partial[[name]]) <-  chemical_name_ref[as.character(rownames(partial[[name]])) ]
    if (  name != 'struct')
    {
      temp_name <- name
      if (name == 'nwauc' & actType == '')
      {
        temp_name <- 'wauc'
      }
      if (name == 'cv') temp_name <- 'nwauc'
      pathway_ref <- conversion(para, inp=temp_name, out='pathway_name')
      colnames(partial[[name]]) <-  pathway_ref[as.character(colnames(partial[[name]])) ]
    }
  }
  return(partial)
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

