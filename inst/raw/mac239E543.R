library(pepStat)

### Extract the peptides from the gpr files.

ps <- makePeptideSet(path="./data/gprFiles/", mapping.file="./data/mapping.csv", log=TRUE, norm.empty=FALSE, rm.control.list=c("empty"), check.row.order=TRUE)

dat <- values(ps)[[1]]
m239peptides <- unique(dat[dat$featureID%in%c(">SIVmac239 Env",">SIVmac239 Env SIVsmE543 Env gp140"),"peptide"])
write.table(data.frame(peptide=peptides), file="m236PepList", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="")

### Query LANL HIV sequence database
#http://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html
#Download result files: table.txt
tab <- read.table("m239pos.tsv", header=TRUE)
tab <- tab[,c("start", "end", "query_sequence", "polyprotein")]
#mannualy correct missaligned peptides
fix <- data.frame(start=131, end=145, query_sequence = "ASTTSTTASAKVDMV", polyprotein="gp160")
tab <- rbind(tab, fix)
tab <- tab[order(tab$start),]

#Get zScores
tab <- cbind(tab, makeZpepMatrix(as.character(tab$query_sequence)))

rd <- RangedData(ranges=IRanges(start= tab$start, end= tab$end, names= tab$query_sequence), space=tab$polyprotein, names=tab$query_sequence, z1 = tab$z1, z2 = tab$z2, z3 = tab$z3, z4 =  tab$z4, z5 = tab$z5)

pep_m239.2 <- rd
save(pep_m239.2, file="pep_m239.2.rda")

