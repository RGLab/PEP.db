###
# Read alignment from fasta file to get scales and sequences
###
readAlign<-function(filename=NULL)
{
	## EASY MODE:
	## Assume first line is >HXB2, second is hxb2 seq, third is >newSeqID,  fourth is newSeq seq
	if(is.null(filename)){filename<-system.file("extdata/alignment.fasta", package="PEP.db")}
	alFile<-file(filename,open="r")
	lineList<-list()
	while(length(oneLine <- readLines(alFile, n = 2, warn = FALSE))) #n=4 means read four lines (negative value for whole file)
	{
		lineList<-c(lineList,oneLine)
	}
	close(alFile) #
	refSeq<-lineList[[2]]
	
	len<-nchar(refSeq)
	gapCnt<-0
	refScale<-numeric(len)
	for(i in 1:len)
	{
		if(substr(refSeq,i,i)=="-")
		{
			gapCnt<-gapCnt+1
		}
		refScale[i]<-i-gapCnt
	}
	refObj<-list()
	refObj[[1]]<-refScale #the scale
	refObj[[2]]<-refSeq   #the sequence with gaps
	return(refObj)
}


###
# convertPep
#  Changes the position, aligned,, trimmed and peptide columns
###
convertPep<-function(rd=PEP.db::pep_hxb2,filename=NULL,refScale=NULL)
{
	## Read the file
	if(is.null(filename)){filename<-system.file("extdata/hxb2_7subtypes.fasta", package="PEP.db")}
	alFile<-file(filename,open="r")
	lineList<-readLines(alFile, n=-1L)#16) #16 lines, i.e ref+7 subtypes
	close(alFile) #
	
	refSeq<-lineList[2]
	len<-nchar(refSeq)
	gapCnt<-0
	
	##Get the refScale if needed
	if(is.null(refScale))
	{
		refScale<-numeric(len)
		for(i in 1:len)
		{
			if(substr(refSeq,i,i)=="-")
			{
				gapCnt<-gapCnt+1
			}
			refScale[i]<-i-gapCnt
		}
	}
	
	##Create the 7 subtype lists
	lCnt<-3 #no need to loop on the reference
	sTypeIDList<-character()
	sTypeSeq<-numeric()
	while(lCnt < 16) #length(lineList))
	{
		sType<-gsub("> ","",lineList[lCnt])
		sType<-gsub(">","",sType)
		sTypeIDList<-c(sTypeIDList,sType)
		sTypeSeq<-c(sTypeSeq,lineList[lCnt+1])
		lCnt<-lCnt+2
	}
	
	##Convert the positions of the rd
	nrd<-coord2ext(rd,refScale)
	
	#t1<-system.time({
	
	##Get the aligned sequence
	aligned<-reference<-character(length(rownames(nrd)))
	for(ID in sTypeIDList)
	{
		TFvec<-nrd[[ID]]
		idx<-which(TFvec==TRUE)
		newAligns<-sapply(idx, function(x){
					substr(sTypeSeq[which(sTypeIDList==ID)],start(nrd)[x],end(nrd)[x])
				})
		aligned[idx]<-newAligns
		
		newRef<-sapply(idx, function(x){
					substr(refSeq,start(nrd)[x],end(nrd)[x])
				})
		reference[idx]<-newRef
	}
	
	##Get the trimmed sequence
	trimmed<-c()
	for(seqCnt in 1:length(aligned))
	{
		trimmedS<-c()
		for(AACnt in 1:nchar(aligned[seqCnt]))
		{
			if(substr(reference[seqCnt],AACnt,AACnt)!="-")
			{
				trimmedS<-c(trimmedS,substr(aligned[seqCnt],AACnt,AACnt))
			}
		}
		trimmed<-c(trimmed,paste(trimmedS,collapse=""))
	}
	
	#})#end t1
	
	##Set the values for the new rd
	nrd$aligned<-aligned
	nrd$trimmed<-trimmed
	#print(t1)
	return(nrd)
}
