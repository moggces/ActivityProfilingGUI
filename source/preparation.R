#### for generating RData in data folder
#### one-time only

##### pre-calculate CV ###

sub_file <- '/tox21_sub_all_pathway_library2_loose_110213.txt'
inconclusiveAsInactive <- FALSE

lib2 <- cal_activity_cv(sub_file, inconclusiveAsInactive=inconclusiveAsInactive)

sub_file <- '/tox21_sub_all_pathway_library1_loose_110213.txt'
inconclusiveAsInactive <- FALSE

lib1 <- cal_activity_cv(sub_file, inconclusiveAsInactive=inconclusiveAsInactive)
lib_cv <- join(lib1, lib2, type='full')

cv_mat <- lib_cv[,-1]
rownames(cv_mat) <- lib_cv$CAS
cv_mat <- cv_mat[order(rownames(cv_mat)),]
cv_mat <- cv_mat[,order(colnames(cv_mat))]
save(cv_mat, file='cv_mat.RData')

cal_activity_cv <- function (sub_file, inconclusiveAsInactive=FALSE)
{
  lib1 <- read.table(paste(getwd(),  sub_file, sep=""), header = TRUE, sep = "\t", quote = "", comment.char = "")
  lib1$unique_id <- paste(lib1$Sample.ID, lib1$Cmpd_Library, sep="@")
  lib1 <- lib1[, c('unique_id', 'CAS',  grep("nwauc_4toxpi$", colnames(lib1), value=TRUE))]
  lib1 <- melt(lib1,  id.vars = c("unique_id","CAS"), na.rm=FALSE, variable.name = "pathway", value.name='nwauc')
  if (inconclusiveAsInactive) lib1[is.na(lib1$nwauc),]$nwauc <- 0 # how should we treat the inconclusive
  lib1_sum <- ddply(lib1, .(CAS, pathway), summarise, mean=mean(nwauc, na.rm=TRUE), std=sd(nwauc, na.rm=TRUE), freq=sum(!is.na(nwauc)), blank=sum(is.na(nwauc)))
  lib1_sum$cv <- lib1_sum$std/lib1_sum$mean
  lib1_sum$cv <- abs(lib1_sum$cv)
  
  lib1 <- dcast(lib1_sum, CAS ~ pathway, value.var='cv')
  return (lib1)
}

##### pre-calculate logit parameters
logit_para_file <- '/logit_para_file.txt'  #logit_para_file call_match_pathway_name
removeAbnormalDirection <- FALSE
isSignal <- FALSE

profile_file <- '/tox21_cas_library1-library2_loose_110213.txt'
select <- c('nwauc', 'pod', 'ac50')
loose <- cal_activity_logit(profile_file, logit_para_file, select=select, isSignal=isSignal, removeAbnormalDirection=removeAbnormalDirection)
save(loose, file='loose.RData')

profile_file <- '/tox21_cas_library1-library2_medium_110213.txt'
select <- c('nwauc', 'pod', 'ac50')
medium <- cal_activity_logit(profile_file, logit_para_file, select=select,isSignal=isSignal, removeAbnormalDirection=removeAbnormalDirection)
save(medium, file='medium.RData')

profile_file <- '/tox21_cas_library1-library2_tight_110213.txt'
select <-  c('nwauc', 'pod', 'ac50')
tight <- cal_activity_logit(profile_file, logit_para_file, select=select, isSignal=isSignal, removeAbnormalDirection=removeAbnormalDirection)
save(tight, file='tight.RData')

profile_file <- '/tox21_cas_library1-library2_loose_110213.txt'
select <- 'wauc'
isSignal <- TRUE
signal_wauc <- cal_activity_logit(profile_file, logit_para_file, select=select, isSignal=isSignal, removeAbnormalDirection=removeAbnormalDirection)
save(signal_wauc, file='signal_wauc.RData')

cal_activity_logit <- function (profile_file, logit_para_file,  select=c('nwauc', 'call', 'pod', 'ac50'), isSignal=FALSE, removeAbnormalDirection=FALSE)
{
  master <- load_profile(profile_file) # global, dataframe output
  para <- load_logit_para(logit_para_file) # global, dataframe output
  # cas as ID
  mat_list <- get_mat(master, para, select=select, isSignal=isSignal, removeAbnormalDirection=removeAbnormalDirection) # list output
  mat_list[['act']] <- mat_list[['act']][,colSums(is.na(mat_list[['act']]))< nrow(mat_list[['act']])]
  return (mat_list)
}

##### pre-load structure fingerprints
profile_file <- '/tox21_cas_library1-library2_loose_110213.txt' # temp ;  bit redundent, just to get the mapping
master <- load_profile(profile_file) # global, dataframe output
structure_fp_base <- '/tox21_8598_fp' # tox21_8598_fp tox21_8306_fp
struct_mat <- load_struc_fp_file(structure_fp_base,master) #global, mat output
save(struct_mat, file='struct_mat.RData')

# split the text file into pieces
profile_file <- 'U:/Projects/TOX21/report/tox21_qHTS_all_111014.txt'
master <- load_profile(profile_file) # the function is modified
props <- get_property_name(master)
activities <- split_master_2_matrix(master, 
      props=grep('label|cc2|cv\\.wauc|pod_med_diff|a_normal|hitcall|npod|nemax|nac50|n?wauc\\.logit|wauc\\.logit)', props, value=TRUE), id='Chemical.ID.GSID')
activities[['label']][, 'label'] <- NULL
lapply(activities, ncol)
save(activities, file='activities.RData')

# load the structure fingerprint
structure_fp_base <- 'U:/Projects/TOX21/Chemical_Curation_from_Ann/20140722/tox21_v5a_leadscope_fp' # tox21_8598_fp tox21_8306_fp
struct_mat <- load_struc_fp_file(structure_fp_base,NULL) #global, mat output
read.excel <- function(header=TRUE,...) {
  read.table("clipboard",sep="\t",header=header,...)
}
dd <- read.excel()
v <- conversion(dd, inp='DSSTox_CID', out='DSSTox_GSID')
v[as.character(rownames(struct_mat))]
save(struct_mat, file='struct_mat.RData')
