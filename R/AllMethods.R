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


#
# Convert the coordinates of an object into the extended coordinate given a scale
# Input: An object and a reference scale
# Output: an object of the same type with coordinates in extended system
setMethod("coord2ext", signature=(obj="numeric"), function(obj, refScale)
		{ 
			extVec<-sapply(obj, function(x){
						min(
								if(x<0)
											x
										else if(x==0)
											which(refScale==1)
										else if(x>refScale[length(refScale)])
											which(refScale==refScale[length(refScale)])
										else if(length(which(refScale==x)))
											which(refScale==x)
										else
											NaN
						)})#sapply#function#min
			extVec<-extVec[!is.na(extVec)]
			return(extVec)
		})

setMethod("coord2ext", signature=(obj="RangedData"), function(obj, refScale)
		{
			if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values
			
			extStart<-coord2ext(start(obj),refScale)
			extEnd<-coord2ext(end(obj),refScale)
			#assign new start coordinates after end to avoid width<0 issues
			end(obj)<-extEnd	
			start(obj)<-extStart
			return(obj)
		})

