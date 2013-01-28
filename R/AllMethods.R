# Input: A DataFrame like pep_hxb2 or pep_mac239
#        The only requirements are petides as rownames and a `clade` column
# Output: matrix of logical with peptides as rownames and clades as colnames
setMethod("clade",
		signature=signature(object="RangedData"),
		definition=function(object)
		{
			cladeList<-unique(unlist(strsplit(levels(as.factor(object$clade)),","))) #List of all possible clades
			len<-nrow(object)
			retMatrix<-matrix(FALSE, nrow=len, ncol=length(cladeList))
			pepClades<-strsplit(object$clade, split=",") #clades for each peptide
			
			for(pepIdx in 1:len)
			{
				tmpList<-cladeList %in% pepClades[[pepIdx]]
				retMatrix[pepIdx,]<-tmpList
			}
			
			rownames(retMatrix)<-rownames(object)
			colnames(retMatrix)<-cladeList
			
			return(retMatrix)
		})

