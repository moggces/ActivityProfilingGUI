load_input_file <- function (input_file)
{
  if (file.exists(input_file))
  { 
    input <- read.table(input_file, header = TRUE, sep = "\t", quote = '', comment.char = "")
  } else
  {
    input <- read.table(paste(getwd(),  input_file, sep=""), header = TRUE, sep = "\t", quote = '', comment.char = "")
  }
  return(input)
}

load_data_matrix <- function (input_file, file_name)
{
  result <- list()
  
  dm_id <- 'unknown'
  dm_id  <- sub(".*(nwauc\\.logit|npod|nac50).*", "\\1", file_name )
  input <- read.table( input_file, header = TRUE, sep = "\t", quote = '', check.names=FALSE, comment.char = "") 
  if (is.null(input$Cluster)) 
  {
    input[, "Cluster"] <- 'unassigned'
    if ( ! is.null(input$userClust)) input[, "Cluster"] <- input$userClust 
  }
  pot_id_cols <- c('CAS','GSID', 'Chemical.Name','chemClust', 'userClust',	'toxScore', 'Cluster')
  id_df <- input[, colnames(input) %in% pot_id_cols]
  data_df <- input[,! colnames(input) %in% pot_id_cols]
  rownames(data_df) <- id_df$CAS
  result[['id']] <- subset(id_df, select=c(CAS, Cluster))
  result[[dm_id]] <- data_df
  return(result)
}

load_profile <- function (profile_file)
{
  #return ( read.table(paste(getwd(),  profile_file, sep=""), header = TRUE, sep = "\t", quote = '"', comment.char = "") )
  return ( read.table( profile_file, header = TRUE, sep = "\t", quote = '', check.names=FALSE, comment.char = "") )
}

#### deprecated
load_logit_para <- function (logit_para_file)
{
  name_match <- read.table(paste(getwd(),  logit_para_file, sep=""), header = TRUE, sep = "\t", quote = "", comment.char = "")
  name_match$geom_mean <- sqrt(name_match$low*name_match$high)
  name_match$inactive_inv_logit <- (exp( (0-name_match$geom_mean)/name_match$geom_mean ))/(name_match$geom_mean+exp((0-name_match$geom_mean)/name_match$geom_mean))
  
  #geom_mean <- name_match$geom_mean 
  #names(geom_mean) <- name_match$nwauc
  #inactive_inv_logit <- name_match$inactive_inv_logit
  #names(inactive_inv_logit) <- name_match$nwauc
  #return(list(geom_mean=geom_mean, inactive_inv_logit=inactive_inv_logit))
  return(name_match)
}

load_struc_fp_file <- function (structure_fp_base, master=NULL)
{
  d <- readXAfile(structure_fp_base, sep="\t")
  mat <- d$x
  if (! is.null(master))
  {
    cas_ref <- conversion(master, inp='StructureID', out='CAS')
    rownames(mat) <- cas_ref[as.character(rownames(mat))]
  }
  mat <- mat[order(rownames(mat)),]
  return(mat)
}

load_text_2_df <- function (textarea)
{
  rows <- lapply(unlist(strsplit(textarea, "\n")), function (x) unlist(strsplit(x, "\t")))
  mat <- do.call("rbind", rows)
  result <- data.frame(mat[-1,], stringsAsFactors=FALSE)
  colnames(result) <- mat[1,]
  if (is.null(result$Cluster)) result[, "Cluster"] <- 'unassigned'
  return(list(id=result))
}
