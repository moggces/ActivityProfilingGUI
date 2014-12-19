# split the activity text file into pieces
profile_file <- 'U:/Projects/TOX21/report/tox21_qHTS_all_111014.txt'
master <- load_profile(profile_file)
props <- get_property_name(master)
activities <- split_master_2_matrix(master, 
      props=grep('label|cc2|cv\\.wauc|pod_med_diff|a_normal|hitcall|npod|nemax|nac50|n?wauc\\.logit|wauc\\.logit)', 
                 props, value=TRUE), id='Chemical.ID.GSID')
activities[['label']][, 'label'] <- NULL # a label column that is accidentally included
lapply(activities, ncol)
save(activities, file='activities.RData')

# load the structure fingerprint
source(paste(getwd(), "/source/customized.R", sep=""), local=TRUE)
structure_fp_base <- 'U:/Projects/TOX21/Chemical_Curation_from_Ann/20140722/tox21_v5a_leadscope_fp_extend' # tox21_8598_fp tox21_8306_fp
struct_mat <- load_struc_fp_file(structure_fp_base,NULL) 
save(struct_mat, file='struct_mat.RData')
