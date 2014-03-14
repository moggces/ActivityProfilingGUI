# R function to read descriptor x, a, and xa file
# Sept 22, 2011 moggces@gmail.com

readXAfile<-function(base, newV=FALSE, sep=""){
	if (newV == FALSE)
	{
		xfile <- paste(base, '.x', sep="")
		ds<-read.table(xfile, fill=TRUE, row.names=NULL, header=FALSE,comment.char="",quote="", sep=sep)
		Ncpd<-as.numeric(as.character(ds[1, "V1"]))
		Ndes<-as.numeric(as.character(ds[1, "V2"]))
		dsheader<-as.vector(unlist(ds[2, 1:Ndes]))
		name <- as.vector(unlist(ds[3:(3+Ncpd-1), 2]))
		coff <- NULL
		dsdata<-as.numeric(as.matrix(ds[3:(3+Ncpd-1), -c(1,2)])) 
		dim(dsdata) <- c(length(name), length(dsheader))
		colnames(dsdata)<-dsheader
		rownames(dsdata)<-name
		
		if (nrow(ds) == Ncpd + 4) # this x file has coefficients
		{
			coff <- ds[ c(nrow(ds)-1, nrow(ds)), 1:Ndes]
			names(coff) = dsheader
		}
		afile <- paste(base, '.a', sep="")
		if (file.exists(afile))
		{
			dsa<-read.table(afile,row.names=NULL,col.names=c("name","activity"),
			header=F,comment.char="",quote="", sep=sep)	
			list(x=as.data.frame(dsdata), desName=dsheader, a=dsa$activity,cpdName=as.vector(dsa$name), coff=coff)
		} else 
		{
			list(x=as.data.frame(dsdata), desName=dsheader, cpdName=name)
		}
	} else 
	{
		xafile <- paste(base, '.xa', sep="")
		ds<-read.table(xafile, fill=TRUE, row.names=NULL, header=FALSE,comment.char="",quote="", sep=sep)
		Ncpd<-as.numeric(as.character(ds[1, "V1"]))
		Ndes<-as.numeric(as.character(ds[1, "V2"]))
		dsheader<-as.vector(unlist(ds[2, 1:Ndes]))
		name <- as.vector(unlist(ds[3:(3+Ncpd-1), 2]))
		
		dsdata<-as.numeric(as.matrix(ds[3:(3+Ncpd-1), -c(1,2)])) 
		dim(dsdata) <- c(length(name), length(dsheader))
		colnames(dsdata)<-dsheader
		rownames(dsdata)<-name
		
		activity <- as.character(unlist(ds[3:(3+Ncpd-1), 3]))
		
		if (nrow(ds) == Ncpd + 4) # this xa file has coefficients
		{
			coff <- ds[ c(nrow(ds)-1, nrow(ds)), 1:Ndes]
			names(coff) = dsheader
		}
		list(x=as.data.frame(dsdata), desName=dsheader, a=activity,cpdName=name, coff=coff)
	}
		
}