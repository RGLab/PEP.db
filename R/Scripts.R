#### The following functions have been used to create some of the 
#### data for the package and will require heavy modifications in 
#### order to be used on other datasets.
#### They are not exported and you should not be reading this.


####
## SIV_2  #It is used to create pep_mac239
## pep_mac239<-.SIV_2("seqFile.fasta", "SIVPeptides-ToOrder.txt")
####
.SIV_2<-function(fastaFile,orderFile)
{ 
	
	
	#Parse fasta file
	IN<-file(fastaFile,open="r")
	lineList<-list()
	while(length(oneLine <- readLines(IN, n = 4, warn = FALSE))) #n=4 means read four lines (negative value for whole file)
	{
		lineList<-c(lineList,oneLine)
	}
	close(IN)
	#Store sequences
	seq1<-lineList[2] #239
	seq2<-lineList[4] #660
	
	#Parse order file
	IN<-file(orderFile,open="r")
	lineList<-list()
	while(length(oneLine <- readLines(IN, n = -1, warn = FALSE)))
	{
		lineList<-c(lineList,oneLine)
	}  
	close(IN)
	lineList<-lineList[which(lineList!="" & lineList!=" ")]
	
	#Lists to create RangedData
	peptideList<-c()
	startList<-c()
	alignedList<-c()
	cladeList<-peptideNb<-c()
	z1sumList<-z2sumList<-z3sumList<-z4sumList<-z5sumList<-c()
	
	idx<-1
	while(idx<=length(lineList))
	{
		if(substr(strsplit(lineList[[idx]],"\\.")[[1]][2],1,1)=="1") #if the peptide is specific to its sequence
		{
			splitLine1<-strsplit(lineList[[idx]]," ")[[1]]
			splitLine1<-splitLine1[which(splitLine1!="")]
			splitLine2<-strsplit(lineList[[idx+1]]," ")[[1]]
			splitLine2<-splitLine2[which(splitLine2!="")]
			
			inc<-2
			#if peptide in E660 without match in mac239
			if(strsplit(splitLine1[3], "_")[[1]][1]!=strsplit(splitLine2[3], "_")[[1]][1])
			{
				inc<-1
			}
			#Get the sequences
			#Couples are given with E660 before mac239
			pep1<-splitLine1[2]
			pep2<-splitLine2[2]
			
			peptideList<-c(peptideList,tail(c(pep2,pep1),inc))
			alignedList<-c(alignedList,tail(c(pep2,pep2),inc)) #mac239 is the ref
			#Get the positions
			if(inc==2)
			{
				start<-regexpr(pep2, seq1)[[1]]
				if(start==-1)
				{
					start<-tail(startList,1)+3
				}
			}
			else
			{
				start<-tail(startList,1)+3
			}
			startList<-c(startList, rep(start, inc)) #the pos are given relative to mac239 coordinates
			#clade
			peptideNb<-c(peptideNb, rep(as.numeric(strsplit(splitLine1[3], "_")[[1]][1]), inc))
			cladeList<-c(cladeList, head(c("E660", "mac239"), inc))
			#zScores
			z1sumList<-c(z1sumList,head(c(.zSum(pep2,"z1"),.zSum(pep1,"z1")),inc))  
			z2sumList<-c(z2sumList,head(c(.zSum(pep2,"z2"),.zSum(pep1,"z2")),inc))
			z3sumList<-c(z3sumList,head(c(.zSum(pep2,"z3"),.zSum(pep1,"z3")),inc))
			z4sumList<-c(z4sumList,head(c(.zSum(pep2,"z4"),.zSum(pep1,"z4")),inc))
			z5sumList<-c(z5sumList,head(c(.zSum(pep2,"z5"),.zSum(pep1,"z5")),inc))
			
			idx<-idx+inc
		}
		else
		{
			splitLine<-strsplit(lineList[[idx]]," ")[[1]]
			#Sequence
			pep<-splitLine[2]
			align<-pep
			peptideList<-c(peptideList,pep)
			alignedList<-c(alignedList,pep)
			#Position
			start<-regexpr(pep, seq1)[[1]]
			if(start==-1)
			{
				print(idx)
			}
			startList<-c(startList,start)
			#clade
			peptideNb<-c(peptideNb, strsplit(splitLine[3], "_")[[1]][1])
			cladeList<-c(cladeList, "E660,mac239")
			#zSum
			z1sumList<-c(z1sumList,.zSum(pep,"z1"))  
			z2sumList<-c(z2sumList,.zSum(pep,"z2"))
			z3sumList<-c(z3sumList,.zSum(pep,"z3"))
			z4sumList<-c(z4sumList,.zSum(pep,"z4"))
			z5sumList<-c(z5sumList,.zSum(pep,"z5"))
			
			idx<-idx+1	
		}
		
		
		
	}
	#Make peptideNb column
	#Make clade column
	
	nrd<-RangedData(ranges=IRanges(start=startList,width=15), aligned=alignedList,
			peptideNb=peptideNb, clade=cladeList,
			z1sum=z1sumList,z2sum=z2sumList,z3sum=z3sumList,z4sum=z4sumList,z5sum=z5sumList)
	rownames(nrd)<-peptideList
#  nrd<-nrd[order(start(nrd)),]
	return(nrd)
}

#
# Returns the zZsum for the AAString
#   ref:Classification of G-protein coupled receptors by 
#       alignment-independent extraction of principal chemical
#       properties of primary amino acid sequences
#
.zSum<-function(AAString, Z=NULL)
{
	#if no Z is specified, the function returns a matrice with the 5zScales
	if(is.null(Z))
		Z<-c("z1","z2","z3","z4","z5")
	rownames<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
	colnames<-c("z1","z2","z3","z4","z5")
	zTable<-c(
			0.24,	-2.32,	0.6,	-0.14,	1.3,
			3.52,	2.5,	-3.5,	1.99,	-1.7,
			3.05,	1.62,	1.04,	-1.15,	1.61,
			3.98,	0.93,	1.93,	-2.46,	0.75,
			0.84,	-1.67,	3.71,	0.18,	-2.65,
			1.75,	0.5,	-1.44,	-1.34,	0.66,
			3.11,	0.26,	-0.11,	-3.04,	-0.25,
			2.05,	-4.06,	0.36,	-0.82,	-0.38,
			2.47,	1.95,	0.26,	3.9,	0.09,
			-3.89,	-1.73,	-1.71,	-0.84,	0.26,
			-4.28,	-1.3,	-1.49,	-0.72,	0.84,
			2.29,	0.89,	-2.49,	1.49,	0.31,
			-2.85,	-0.22,	0.47,	1.94,	-0.98,
			-4.22,	1.94,	1.06,	0.54,	-0.62,
			-1.66,	0.27,	1.84,	0.7,	2,
			2.39,	-1.07,	1.15,	-1.39,	0.67,
			0.75,	-2.18,	-1.12,	-1.46,	-0.4,
			-4.36,	3.94,	0.59,	3.44,	-1.59,
			-2.54,	2.44,	0.43,	0.04,	-1.47,
			-2.59,	-2.64,	-1.54,	-0.85,	-0.02)
	dim(zTable)<-c(20,5)
	dimnames(zTable)<-list(rownames,colnames)
	
	total<-0
	for(AA in 1:nchar(AAString))
	{
		total<-total+zTable[substr(AAString,AA,AA),Z]
	}
	
	return(total)
	
}


### Construction of the new pep_hxb2 with the new clade system
#nPep<-pep_hxb2
#cladCol<-unlist(lapply(rownames(nPep), function(x)
#		{
#			paste(names(which(unlist(nPep@values[["gp160"]][x,3:9])==TRUE)), collapse=",")
#		}))
#
#ndf<-as.data.frame(nPep[,c(1:2,11:16)])
#final<-cbind(ndf, clade=cladCol, stringsAsFactors=FALSE)
#cladeHXB2<-as(final, "RangedData")
#rownames(cladeHXB2)<-rownames(nPep)

