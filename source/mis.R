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


edit_mat_manual <- function (partial, nwaucThres=0.0001, actType='', regSel='', invSel=FALSE)
{
  
  partial[['nwauc']][is.na(partial[['nwauc']]) ] <- 0.0001  # for plotting
  
  
  if ( ! is.null( partial[['cv']] ))
  {
    partial[['cv']][is.na(partial[['cv']])] <- -0.1 # just to remove the na 
    partial[['cv']][partial[['nwauc']] < 0.1  ] <- -0.1 # if two small just let it be
    partial[['cv']][partial[['cv']] > 1.4] <- "*"
    partial[['cv']][partial[['cv']] != "*"] <- ""
    ################################ temp for paper
    #partial[['nwauc']] <- partial[['nwauc']][,grep("are|hse|fxr_|pparg_antagonism|ppard_", colnames(partial[['nwauc']]), value = TRUE, invert = TRUE)]
    #partial[['cv']] <- partial[['cv']][,grep("are|hse|fxr_|pparg_antagonism|ppard_", colnames(partial[['cv']]), value = TRUE, invert = TRUE)]
    ###########################################
    #if (removeCytotoxic) partial[['cv']] <- partial[['cv']][,grep("cytotoxicity", colnames(partial[['cv']]), value = TRUE, invert = TRUE)]
    partial[['cv']] <- partial[['cv']][,grep(regSel, colnames(partial[['cv']]), value = TRUE, invert = invSel)]
  }
  
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
    
    for (name in c('pod', 'ac50'))
    {
      partial[[name]][partial[['nwauc']] < nwaucThres & partial[[name]]  != 0.0001 ] <- (partial[[name]][partial[['nwauc']] < nwaucThres & partial[[name]]  != 0.0001 ])*-1
    }
    
  }
  
  
  return(partial)
}

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

