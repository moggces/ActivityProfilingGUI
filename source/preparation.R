#split the activity text file into pieces
profile_file <- 'U:/Projects/TOX21/report/tox21_qHTS_all_111415.txt'
master <- load_profile(profile_file)
props <- get_property_name(master)
# 
# # to reduce the size, only used parameters are included
activities <- split_master_2_matrix(master, 
      props=grep('nemax|label|cc2|cv\\.wauc|pod_med_diff|wauc_fold_change|a_normal|hitcall|npod|ncmax|nec50|n?wauc\\.logit|wauc\\.logit|wauc|pod|ec50|cmax|nac50|nemax)', 
                 props, value=TRUE), id='GSID') # change ID (GSID)
activities[['label']][, 'label'] <- NULL # a label column that is accidentally included
lapply(activities, ncol) #the number has to match up with the number in tox21_assay_collection$assay
save(activities, file='activities.RData')
# 
# # load the structure fingerprint
# source(paste(getwd(), "/source/customized.R", sep=""), local=TRUE)
# structure_fp_base <- 'U:/Projects/TOX21/Chemical_Curation_from_Ann/20140722/tox21_v5a_leadscope_fp_extend' # tox21_8598_fp tox21_8306_fp
# struct_mat <- load_struc_fp_file(structure_fp_base,NULL) 
# save(struct_mat, file='struct_mat.RData')

# get the hitcall_simplified
#load("U:/Projects/CurveP/ActivityProfilingGUI/data/activities.RData")
#hitcall <- activities[['hitcall']]
#hitcall[activities[['cv.wauc']] > 1.4 & ! is.na(activities[['cv.wauc']])] <- NA
#hitcall[hitcall < 0 & ! is.na(hitcall)] <- NA
#tox21_assay_collection <- read.delim("U:/Projects/CurveP/ActivityProfilingGUI/data/tox21_assay_collection.txt", stringsAsFactors=FALSE)
#name_vec <- tox21_assay_collection$common_name
#names(name_vec) <- tox21_assay_collection$assay
#name_vec[colnames(hitcall)]
#colnames(hitcall) <- name_vec[colnames(hitcall)]
#write.table(hitcall, file=paste(getwd(), '/data/hitcall_simplified.txt' , sep=""),  row.names = TRUE, col.names = TRUE, sep="\t", quote=FALSE, append=FALSE)