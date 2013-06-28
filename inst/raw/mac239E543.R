library(pepStat)

### Extract the peptides from the gpr files.

ps <- makePeptideSet(path="./data/gprFiles/", mapping.file="./data/mapping.csv", log=TRUE, norm.empty=FALSE, rm.control.list=c("empty"), check.row.order=TRUE)

dat <- values(ps)[[1]]
m239peptides <- unique(dat[dat$featureID%in%c(">SIVmac239 Env",">SIVmac239 Env SIVsmE543 Env gp140"),"peptide"])
write.table(data.frame(peptide=peptides), file="m236PepList", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="")
E543peptides <- unique(dat[dat$featureID%in%c(">SIVsmE543 Env gp140",">SIVmac239 Env SIVsmE543 Env gp140"),"peptide"]) 
write.table(data.frame(peptide=peptides), file="E543PepList", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="")

### Query LANL HIV sequence database
#http://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html
#Download result files: table.txt
tab <- read.table("m239pos.tsv", header=TRUE)
tab <- tab[,c("start", "end", "query_sequence", "polyprotein")]
#mannualy correct missaligned peptides
fix <- data.frame(start=131, end=145, query_sequence = "ASTTSTTASAKVDMV", polyprotein="gp160") #matched with LTR5
tab <- rbind(tab, fix)
tab <- tab[order(tab$start),]

#Get zScores
tab <- cbind(tab, makeZpepMatrix(as.character(tab$query_sequence)))

rd <- RangedData(ranges=IRanges(start= tab$start, end= tab$end, names= tab$query_sequence), space=tab$polyprotein, names=tab$query_sequence, clade="mac239", z1 = tab$z1, z2 = tab$z2, z3 = tab$z3, z4 =  tab$z4, z5 = tab$z5)


tab2 <- read.table("smE543pos.tsv", header=TRUE)
tab2 <- tab2[,c("start", "end", "query_sequence", "polyprotein")]
tab2 <- tab2[order(tab2$start),]
tab2 <- cbind(tab2, makeZpepMatrix(as.character(tab2$query_sequence)))
rd2 <- RangedData(ranges=IRanges(start= tab2$start, end= tab2$end, names= tab2$query_sequence), space=tab2$polyprotein, names=tab2$query_sequence, clade="smE543", z1 = tab2$z1, z2 = tab2$z2, z3 = tab2$z3, z4 =  tab2$z4, z5 = tab2$z5) 

#change clade info
rdo <- rd[rownames(rd)%in%rownames(rd2),]
rdo$clade <- "mac239,smE543"
rd2 <- rd2[!(rownames(rd2) %in% rownames(rdo)),]
rd <- rd[!(rownames(rd) %in% rownames(rdo)),]

pep_m239smE543 <- rbind(rd,rd2,rdo)
save(pep_m239smE543, file="pep_m239smE543.rda")


##mannualy correct missaligned peptides
#fixQS <- c("DQNSNRWKQQKKPEQ", "TTAKSTTSTTTTTVT", "VEDRDQNSNRWKQQK", "AETTTTAKSTTSTTT") #Matched with Vif Tat Vif 
#fixStart <-c( , 131, , 127)
#fixEnd <- c( , 145, , 141)
#fix <- data.frame(start=131, end=145, query_sequence = "ASTTSTTASAKVDMV", polyprotein="gp160")

#STRAT 2
#1. Reconstruct sequences from the peptides
#2. Align the sequence to LANL mac239
#3. locate the peptides on the custom seq
#4. translate into LANL mac239 coords

clade <- "mac239"
#clade <- "smE543"
if(clade=="mac239"){
  peps <- as.character(tab$query_sequence)
  seq <- "MGCLGNQLLIAILLL"
} else{
  peps <- as.character(tab2$query_sequence)
  seq <- "MGCLGNQLLIALLLV"
}

#Check that peptides do overlap
len <- length(peps)
for(i in 2:len){
  res <- grep(substr(peps[i-1], nchar(peps[i-1])-2, nchar(peps[i])), peps[i])
  pairwiseAlignment(peps[i-1], peps[i], type="overlap")
  if(res!=1)
    stop("no overlap")
}

for(i in 2:len){
  res <- pairwiseAlignment(substr(seq, nchar(seq)-14, nchar(seq)), peps[i], type="overlap")
  ovlp <- as.character(res@subject)
  toAdd <- gsub(ovlp, "", peps[i])
  seq <- paste0(seq, toAdd)
  if(i%%10==0)
    print(i)
}

#Get positions
starts <- c()
for(i in 1:len){
  starts <- c(starts, regexpr(peps[i], seq)[[1]])
}
ends <- starts+14

#Get extended positions
if(clade=="mac239"){
  #refScale <- Pviz::readAlign("inst/extdata/alignments/MuscleLANL_239_E543.fasta", seqLine=4)[[1]]
  refScale <- Pviz::readAlign("~/workspace/BioC/PEP.db/inst/extdata/alignments/MuscleLANL_239_E543.fasta", seqLine=4)[[1]]
else{
  #refScale <- Pviz::readAlign("inst/extdata/alignments/MuscleLANL_239_E543.fasta", seqLine=6)[[1]]
  refScale <- Pviz::readAlign("~/workspace/BioC/PEP.db/inst/extdata/alignments/MuscleLANL_239_E543.fasta", seqLine=6)[[1]]
}
extStarts <- Pviz::coord2ext(starts, refScale)
extEnds <- Pviz::coord2ext(ends, refScale)

#Get LANLmac239 positions
LANLscale <- Pviz::readAlign("inst/extdata/alignments/MuscleLANL_239_E543.fasta")[[1]]
LANLstarts <- LANLscale[extStarts]
LANLends <- LANLscale[extEnds]

#Get alignment for each peptide
refSeq <- Pviz::readAlign("inst/extdata/alignments/MuscleLANL_239_E543.fasta")[[2]]
aligned <- c()
for(i in 1:len){
  aligned <- c( aligned, substr(refSeq, extStarts[i], extEnds[i]))
}
zMat <- makeZpepMatrix(peps)
if(clade=="mac239"){
  rdm239 <- RangedData(ranges=IRanges(start=extStarts, end=extEnds, names=peps),space="gp160", aligned=aligned, clade=clade, zMat)
} else{
  rdE543 <- RangedData(ranges=IRanges(start=extStarts, end=extEnds, names=peps),space="gp160", aligned=aligned, clade=clade, zMat)
}

rdo <- rdm239[rownames(rdm239) %in% rownames(rdE543),]
rdo$clade <- "mac239,E543"
rd1 <- rdm239[!rownames(rdm239)%in%rownames(rdo),]
rd2 <- rdE543[!rownames(rdE543)%in%rownames(rdo),]
pep_m239smE543 <- rbind(rd1, rdo, rd2)
save(pep_m239smE543, file="pep_m239smE543.rda")




