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
