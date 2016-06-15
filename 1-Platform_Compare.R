## Figure 1 for HLA Manuscript ##
## Tally Samples by Platform ##
## Compare HLA Results Across Platforms ##
## March 1, 2016 ##
## Kristopher Standish ##

## Overall Goal:
 # Compare Results of HLA typing from multiple platforms
 	 # SOP = Reads -> SOAP-HLA
   # CHP = SNP Chip -> SNP2HLA
   # SEQ = HaplotypeCaller SNPs -> SNP2HLA
   # LAB = Lab Typing (Gold Standard)

#############################################################
## GAME PLAN ################################################
#############################################################
## (Which questions do I want to answer w/ this figure?)

## A - How many samples were typed by each method?
 # How many samples overlap one another?
 # SOP, CHP, SEQ, LAB
 # Binary Heatmap of Platform vs Sample (typed or not)
## B - Precision of Platforms
 # 2, 4+
 # 5 Genes
 # How to represent diploidy?
   # % of Haplotypes? (e.g., 2*N_Samps)
 # 3-D Table:
   # Counts by Precision x Gene x Platform
   # 2 x 5 x 4 Table
## C - Concordance b/n SNP2HLA approaches
 # 

#############################################################
## GLOBALS ##################################################
#############################################################

## Specify Phenos/Covs
 # Response Phenotypes
PHENOS <- c("DEL_MNe_MN","DEL_lCRP_MNe_MN","DEL_rSJC_MNe_MN","DEL_rTJC_MNe_MN",
	"DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA","RF" )
 # Covariates (ACPA + disease duration)
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)",
	"","")
 # Covariates (BMI + RF + ACPA + disease duration)
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"BMI+ACPA+RF+log(DIS_DUR)","BMI+ACPA+RF+log(DIS_DUR)","BMI+ACPA+RF+log(DIS_DUR)","BMI+ACPA+RF+log(DIS_DUR)",
	"","")
 # Covariates (RF + ACPA + disease duration)
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA+RF+log(DIS_DUR)","ACPA+RF+log(DIS_DUR)","ACPA+RF+log(DIS_DUR)","ACPA+RF+log(DIS_DUR)",
	"","")
PH.COV.1 <- data.frame(PHENOS,COVS)
 # Look at Fewer Phenotypes
PH.COV <- PH.COV.1[c(1,5,9),]
# PH.COV <- PH.COV.1[c(1,2,5,6,9,10),]
PHENOS <- as.character(PH.COV[,1])
COVS <- as.character(PH.COV[,2])

## COLORS: Set Color Scheme
 # FCT: Blend 2 Colors
BLEND <- function( colors ) { colorRampPalette(colors)(3)[2] }
 # Set Color Palette
COLS.list <- c("firebrick2","chocolate2","gold1","chartreuse2","cadetblue2","dodgerblue2","slateblue2","magenta2")
COLS.list.2 <- c("aquamarine3","deeppink2",BLEND(c("gold1","chartreuse2")),BLEND(c("steelblue3","slateblue3")),"tomato2","deepskyblue1",BLEND(c("sienna2","goldenrod1")) )
## COLORS: Pick Specific Colors
 # Multiple Hypothis Correction Color
COLS.cor <- COLS.list[1]
 # Platform Colors
PLATS <- c("SOP","CHP","SEQ","LAB")
COLS.plat <- COLS.list.2[1:length(PLATS)]
names(COLS.plat) <- PLATS
 # Haplotype/AA Frequency Colors
COLS.fr <- COLS.list.2[7]
COLS.AA <- c(colorRampPalette(COLS.list)(26),"black","grey90","grey50") ; names(COLS.AA) <- c(LETTERS,"*",".","?")
 # Phenotype Colors
COLS.ph <- COLS.list.2[c(6,3,2)]
names(COLS.ph) <- PHENOS
 # Null LMM Model Colors
COLS.null <- COLS.list.2[c(1,4,5)]
COLS.beta <- COLS.list.2[1]

#############################################################
## LOAD DATA ################################################
#############################################################
library(xlsx)
library(stringr)
library(gplots)
library(xtable)
## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data Sets & Save Locations (Mac)
PathToSOAP <- "/Users/kstandis/Data/Janssen/Data/HLA/SOAP_HLA_Types/20151211_HLA_Types.Rdata"
PathToS2H <- "/Users/kstandis/Data/Janssen/Data/HLA/SNP2HLA_Types/20160119_Tables/"
PathToS2H <- "/Users/kstandis/Data/Janssen/Data/HLA/SNP2HLA_Types/20160323_Tables/"
PathToLAB <- "/Users/kstandis/Data/Janssen/Data/HLA/Lab_Types/LabCorp HLA typing data_ART3001.xls"
PathToComp <- "/Users/kstandis/Data/Janssen/Data/HLA/Compare/20160323_Tables/1-CompiledTypes.Rdata"
PathToPheno <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20150520_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_ManuHLA_Fig1/",sep="")
PathToRefs <- "/Users/kstandis/Data/Genetics/HLA/Alignments_Rel_3190/"
if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }

## Load Compiled HLA Types
load( PathToComp )
TYPES <- TO_SAVE

## Load Individual HLA Types
 # SOAP-HLA
load( PathToSOAP )
SOP.l <- COMPILE
 # SNP2HLA
CHP.2.l <- read.table( paste(PathToS2H,"CHP_2.txt",sep=""),sep="\t",header=T,colClasses="character" )
CHP.4.l <- read.table( paste(PathToS2H,"CHP_4.txt",sep=""),sep="\t",header=T,colClasses="character" )
SEQ.2.l <- read.table( paste(PathToS2H,"SNP_2.txt",sep=""),sep="\t",header=T,colClasses="character" )
SEQ.4.l <- read.table( paste(PathToS2H,"SNP_4.txt",sep=""),sep="\t",header=T,colClasses="character" )
 # Lab Typing
LAB.l <- read.xlsx( PathToLAB, sheetIndex=1, rowIndex=1:101, header=T, colIndex=1:11)

## Load Phenotypes & Sample Lists
FT <- read.table( PathToPheno, sep="\t",header=T )

## Load Amino Acid Data
FILES <- list.files(PathToRefs)
gsub("_prot.txt","",FILES[grep("_prot",FILES)])

#############################################################
## ORGANIZE PATIENT DATA ####################################
#############################################################

#### SAMPLE LISTS ####

## Get Sample List for each Platform
SAMP.all <- union( as.character(FT$ID), colnames(SOP.l$GENES.2.list$DRB1) )
SAMP.sop <- colnames( SOP.l$GENES.2.list$DRB1 )
SAMP.chp <- rownames( CHP.2.l )
SAMP.seq <- rownames( SEQ.2.l )
SAMP.lab <- unique(as.character( LAB.l$Accession[which( LAB.l$Accession %in% SAMP.all )] ))
SAMP.asn <- setdiff( SAMP.seq, SAMP.chp )
SAMP.art3001 <- as.character( FT$ID )

## LAB: Specify Rows by Sample Name
Which_Samp <- which( LAB.l$Accession %in% SAMP.lab )
Subject_Key <- data.frame( ID=LAB.l$Accession[Which_Samp], NUM=LAB.l$Subject[Which_Samp] )
Subject_Key <- Subject_Key[ which(!duplicated(Subject_Key$ID)), ]
LAB.2 <- merge( Subject_Key, LAB.l, by.x="NUM",by.y="Subject" )
 # Problem w/ Subject Names W367072 & W367073
 # Listed as SAME sample in LAB typing, but different for phenotypes
LAB.2 <- LAB.2[ -which(LAB.2$ID=="W367073"), ]
SAMP.lab <- as.character(unique( LAB.2$ID ))

## Get Samples that Intersect ALL Platforms
SAMP.int <- Reduce( intersect, list(SAMP.sop,SAMP.chp,SAMP.seq,SAMP.lab) )

#### GENE LISTS ####

## Get Gene List for each Platform
GENE.sop <- sort( names( SOP.l$GENES.2.list ) )
GENE.chp <- sort( unique( sapply( colnames(CHP.2.l), function(x) substr(x,1,nchar(x)-1) ) ) ) # sort( gsub("HLA","", unique( sapply( colnames(CHP.2.l), function(x) substr(x,1,nchar(x)-1) ) )) )
GENE.seq <- sort( GENE.chp )
GENE.lab <- sort( unique(as.character(LAB.2$TestName)) ) # sort( gsub("HLA-","", unique(as.character(LAB.2$TestName)) ) )
GENE.all <- gsub("[0-9]", "", GENE.sop[2:length(GENE.sop)] ) ; GENE.all[grep("TAP",GENE.all)] <- c("TAP1","TAP2")
GENE.all <- gsub("[0-9]", "", GENE.sop )
GENE.tab <- array(, c(length(GENE.all),4) )
rownames(GENE.tab) <- GENE.all ; colnames(GENE.tab) <- PLATS
 # Fill Table
for ( row in 1:nrow(GENE.tab) ) {
	name <- rownames(GENE.tab)[row]
	# if ( row < 4 ) {
	which_sop <- union( which(GENE.sop==name | GENE.sop==paste("HLA",name,sep="") | GENE.sop==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.sop ) )
	which_chp <- union( which(GENE.chp==name | GENE.chp==paste("HLA",name,sep="") | GENE.chp==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.chp ) )
	which_seq <- union( which(GENE.seq==name | GENE.seq==paste("HLA",name,sep="") | GENE.seq==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.seq ) )
	which_lab <- union( which(GENE.lab==name | GENE.lab==paste("HLA",name,sep="") | GENE.lab==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.lab ) )
	GENE.tab[row,"SOP"] <- ifelse( length(which_sop)>0, GENE.sop[ which_sop ], NA )
	GENE.tab[row,"CHP"] <- ifelse( length(which_chp)>0, GENE.chp[ which_chp ], NA )
	GENE.tab[row,"SEQ"] <- ifelse( length(which_seq)>0, GENE.seq[ which_seq ], NA )
	GENE.tab[row,"LAB"] <- ifelse( length(which_lab)>0, GENE.lab[ which_lab ], NA )
}
GENE.tab["C","LAB"] <- "HLA-CW"

#############################################################
## A - SAMPLE TYPING BY PLATFORM ############################
#############################################################

## Set Sample & Platform Lists/Parameters
SAMPS <- as.character( FT$ID )

## Compile Table w/ all Samples & Platforms (binary)
A.tab <- array( 0 ,c(length(SAMPS),length(PLATS)) )
colnames(A.tab) <- PLATS
rownames(A.tab) <- SAMPS
A.tab[which(SAMPS%in%SAMP.lab),"LAB"] <- 1
A.tab[which(SAMPS%in%SAMP.chp),"CHP"] <- 2
A.tab[which(SAMPS%in%SAMP.seq),"SEQ"] <- 3
A.tab[which(SAMPS%in%SAMP.sop),"SOP"] <- 4
A.tab.sum <- rowSums(A.tab)
 # How many people typed per platform?
A.plat <- apply( A.tab, 2, function(x)length(which(x!=0)) )

## Heatmap of Binary Platform Table
Temp.names <- paste(colnames(A.tab),"\n(n=",A.plat,")",sep="")
colnames(A.tab) <- Temp.names
# A.tab.2 <- t(A.tab[ order(A.tab[,4],A.tab[,3]), 4:1 ])
A.tab.2 <- t(A.tab[ order(A.tab[,3],A.tab.sum), ])
colnames(A.tab.2) <- NULL
COLS.heat <- c("black",rev(COLS.plat[PLATS]))
BRKS <- seq(-.5,4.5,1)
png(paste(PathToPlot,"1A-Platform_Heat.png",sep=""), height=1200,width=2000,pointsize=32 )
heatmap.2( A.tab.2, scale="none",trace="none",
	Colv=NA,Rowv=NA,dendrogram="none",sepwidth=c(.0001,.01),rowsep=1:nrow(A.tab.2),colsep=1:ncol(A.tab.2),
	col=COLS.heat,key=F,lwid=c(1,20),lhei=c(1,5),
	labCol="",cexRow=1.2,
	main="HLA Typing Platform by Patient",xlab="Patient",ylab="" )
dev.off()

#############################################################
## B - PLATFORM PRECISION BY GENE ###########################
#############################################################

## Set Genes
GENES <- c("A","B","C","DRB","DQB")

## FCT: Darken Platform Color
DARK <- function( color, n_out ) { colorRampPalette(c(color,"black"))(n_out+1)[1:n_out] }

## Specify Table w/ Types across platforms at 4-digit precision
TYPE.4 <- TYPES$DIG4[SAMPS,,GENES]
TYPE.2 <- TYPES$DIG2[SAMPS,,GENES]

## How many haplotypes in population are typed at 4-digit precision (by platform)?
TEMP.4 <- lapply(PLATS, function(p)lapply(GENES, function(g) table(apply(sapply(strsplit(TYPE.4[,p,g],"///"),"[",1:2),2,function(s)length(which(nchar(s)>=4)) )) ))
names(TEMP.4) <- PLATS
for ( p in PLATS ) { names(TEMP.4[[p]]) <- GENES }
## How many haplotypes in population are typed (non-NA, by platform)? 
# TEMP.2 <- lapply(PLATS, function(p)lapply(GENES, function(g) table(apply(sapply(strsplit(TYPE.2[,p,g],"///"),"[",1:2),2,function(s)length(which(!is.na(s))) )) ))
TEMP.2 <- lapply(PLATS, function(p)lapply(GENES, function(g) table(apply(sapply(strsplit(TYPE.2[,p,g],"///"),"[",1:2),2,function(s)length(which(s!="NA")) )) ))
names(TEMP.2) <- PLATS
for ( p in PLATS ) { names(TEMP.2[[p]]) <- GENES }

## Compile Precision Sums into Table
PREC.names <- c("d0","d2","d2+","d4+")
B.tab <- array( 0, c(length(GENES),length(PLATS),length(PREC.names)) )
rownames(B.tab) <- GENES
colnames(B.tab) <- PLATS
dimnames(B.tab)[[3]] <- PREC.names
for ( p in PLATS ) {
	for ( g in GENES ) {
		B.tab[g,p,"d4+"] <- sum( TEMP.4[[p]][[g]]*as.numeric(names(TEMP.4[[p]][[g]])) )
		B.tab[g,p,"d2+"] <- sum( TEMP.2[[p]][[g]]*as.numeric(names(TEMP.2[[p]][[g]])) )
		# (2*A.plat[p]) - B.tab[g,p,"d4+"]
	}
}
B.tab[,,"d2"] <- B.tab[,,"d2+"] - B.tab[,,"d4+"]
B.tab[c("A","B","C"),"LAB","d2"] <- 40
B.tab[c("A","B","C"),"LAB","d4+"] <- 0
B.tab[c("DRB","DQB"),"LAB","d2"] <- 0
B.tab[c("DRB","DQB"),"LAB","d4+"] <- 40

B.tab.d4 <- apply(B.tab[,,"d4+"],1,function(x) x/(2*A.plat) )
B.tab.d2 <- apply(B.tab[,,"d2"],1,function(x) x/(2*A.plat) )
B.tab.d0 <- 1 - (B.tab.d4+B.tab.d2)
## Plot it
 # Parameters
YLIM <- c(0,1.3)
MAIN <- "Precision of HLA Typing by Platform"
XLAB <- "Gene"
YLAB <- "% of Haplotypes Typed by Platform"
COLS.bars <- matrix(unlist( lapply(COLS.plat,function(x)DARK(x,2)) ),ncol=2,byrow=T)
rownames(COLS.bars) <- names(COLS.plat)
 # Open File
png(paste(PathToPlot,"1B-Platform_Precision.png",sep=""), height=1200,width=1600,pointsize=32 )
 # Create Plot
BARS <- barplot( B.tab.d4, col=COLS.bars[,1], beside=T,border=NA, yaxt="n",ylim=YLIM, main=MAIN,xlab=XLAB,ylab=YLAB )
axis( 2, at=seq(0,1,.2),label=seq(0,100,20), las=2 )
abline(h=seq(0,1,.2),lty=3,col="grey50",lwd=1)
 # Plot Data @ 4,2,0-digit precision
barplot( B.tab.d4, col=COLS.bars[,1], beside=T,border=NA, yaxt="n",add=T )
rect(xleft=c(BARS)-.5, ybottom=c(B.tab.d4), xright=c(BARS)+.5, ytop=c(B.tab.d2+B.tab.d4), col=COLS.bars[,2],border=NA )
rect(xleft=c(BARS)-.5, ybottom=c(B.tab.d4+B.tab.d2), xright=c(BARS)+.5, ytop=1, col="black",border=NA )
 # Add Legends
legend("topright",fill=COLS.bars[,1],legend=PLATS,title="Platform",ncol=2 )
legend("topleft",fill=c("grey70","grey30","black"),legend=c("4+ digit","2 digit","Not Typed"),title="Precision",ncol=2 )
 # Add "n=" Labels
text( BARS[,1], .1, label=paste("n=",A.plat,sep=""), srt=90 )
dev.off()

#############################################################
## C - CONCORDANCE b/n COMPUTATIONAL TYPES ##################
#############################################################

#######################################
## Calculate Concordance b/n Methods ##

## FCT: Compare One Platform to Several Others
Calc_Conc <- function( TYPE.arr, Platform, Prim_Comp, Samples, Genes ) {
	## Determine Number of Shared Alleles per Person per Gene
	Samp.int <- intersect( Samples, rownames(TYPE.arr) )
	# TYPE.plat <- TYPE.arr[ Samples,,Genes ]
	TYPE.plat <- TYPE.arr[ Samp.int,,Genes ]
	TYPE.plat.conc <- TYPE.plat
	for ( g in 1:length(Genes) ) {
		gene <- Genes[g]
		for ( c in 1:length(Prim_Comp) ) {
			comp <- Prim_Comp[c]
			# Split Alleles
			split.plat <- strsplit( TYPE.plat[,Platform,gene], "///" )
			split.comp <- strsplit( TYPE.plat[,comp,gene], "///" )
			for ( s in 1:nrow(TYPE.plat) ) {
				samp <- rownames(TYPE.plat)[s]
				temp.plat <- split.plat[[samp]] ; temp.plat <- temp.plat[which(temp.plat!="NA")]
				temp.comp <- split.comp[[samp]] ; temp.comp <- temp.comp[which(temp.comp!="NA")]
				if ( length(temp.plat)==0 ) { 
					TYPE.plat.conc[samp,comp,gene] <- NA
					next
				}
				if ( length(temp.comp)>0 ) {
					temp_conc.fwd <- length(which( temp.plat == temp.comp ))
					temp_conc.rev <- length(which( temp.plat == rev(temp.comp) ))
					TYPE.plat.conc[samp,comp,gene] <- max( temp_conc.fwd, temp_conc.rev )
				}else{
					TYPE.plat.conc[samp,comp,gene] <- NA
				}
			}
		}
	}

	## Calculate % Concordance for Each Gene & Platform
	TYPE.plat.c.a <- array( ,c( length(Prim_Comp), dim(TYPE.plat)[3] ))
	colnames(TYPE.plat.c.a) <- dimnames(TYPE.plat)[[3]]
	rownames(TYPE.plat.c.a) <- Prim_Comp
	 # ...and determine # of Patients in Intersection of Platforms
	TYPE.plat.c.b <- TYPE.plat.c.a
	 # Number of Patients w/ 0,1,2 correct alleles
	TYPE.plat.c.c <- array( ,c( 3, dim(TYPE.plat)[3], length(Prim_Comp) ))
	rownames(TYPE.plat.c.c) <- 0:2
	colnames(TYPE.plat.c.c) <- dimnames(TYPE.plat)[[3]]
	dimnames(TYPE.plat.c.c)[[3]] <- Prim_Comp
	 # Loop through Genes & Patients
	for ( g in 1:length(Genes) ) {
		gene <- Genes[g]
		TYPE.plat.c.a[,gene] <- apply( TYPE.plat.conc[,Prim_Comp,gene], 2, function(x) mean(as.numeric(x),na.rm=T)/2 )
		TYPE.plat.c.b[,gene] <- apply( TYPE.plat.conc[,Prim_Comp,gene], 2, function(x) length(which(!is.na(x))) )
		TYPE.plat.c.c[,gene,] <- rbind( apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==0))), apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==1))), apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==2))) )
	}
	## Compile Outputs
	COMPILE <- list( Type=TYPE.plat, Conc=TYPE.plat.conc, A=TYPE.plat.c.a, B=TYPE.plat.c.b, C=TYPE.plat.c.c )
	return(COMPILE)
}

## SNP Chip (SNP2HLA) Concordance
CHPv <- list()
Platform <- "CHP"
Prim_Comp <- c("SEQ","SOP","LAB","LAB.lik")
Samples <- SAMP.chp
 # Genes A, B, C (2-digit)
Genes <- c("A","B","C","DRB","DQB")
CHPv$p2 <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
 # Genes DRB, DQB (4-digit)
Genes <- c("DRB","DQB")
CHPv$p4 <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )

## HaplotypeCaller (SNP2HLA) Concordance
SEQv <- list()
Platform <- "SEQ"
Prim_Comp <- c("CHP","SOP","LAB","LAB.lik")
Samples <- SAMP.seq
 # Genes A, B, C (2-digit)
Genes <- c("A","B","C","DRB","DQB")
SEQv$p2 <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
 # Genes DRB, DQB (4-digit)
Genes <- c("DRB","DQB")
SEQv$p4 <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )

## HaplotypeCaller (SNP2HLA) Concordance
SOPv <- list()
Platform <- "SOP"
Prim_Comp <- c("CHP","SEQ","LAB","LAB.lik")
Samples <- SAMP.art3001
 # Genes A, B, C (2-digit)
Genes <- c("A","B","C","DRB","DQB")
SOPv$p2 <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
 # Genes DRB, DQB (4-digit)
Genes <- c("DRB","DQB")
SOPv$p4 <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )

#######################################
## Compile Concordance Results ########

## 2-digit Precision (All Genes)
VS <- c("SEQ","SOP")
COMP.2.CHPv <- data.frame( PLAT="CHP",VS=VS, N=apply(CHPv$p2$B[VS,],1,max), CHPv$p2$A[VS,])
VS <- c("CHP","SOP")
COMP.2.SEQv <- data.frame( PLAT="SEQ",VS=VS, N=apply(SEQv$p2$B[VS,],1,max), SEQv$p2$A[VS,] )
COMP.2 <- rbind( COMP.2.CHPv, COMP.2.SEQv )
COMP.2 <- COMP.2[-which(COMP.2$PLAT=="SEQ" & COMP.2$VS=="CHP"),]
rownames(COMP.2) <- paste(COMP.2$PLAT,COMP.2$VS,sep="v")

## 4-digit Precision (DRB,DQB)
VS <- c("SEQ","SOP")
COMP.4.CHPv <- data.frame( PLAT="CHP",VS=VS, N=apply(CHPv$p4$B[VS,],1,max), CHPv$p4$A[VS,] )
VS <- c("CHP","SOP")
COMP.4.SEQv <- data.frame( PLAT="SEQ",VS=VS, N=apply(SEQv$p4$B[VS,],1,max), SEQv$p4$A[VS,] )
COMP.4 <- rbind( COMP.4.CHPv, COMP.4.SEQv )
COMP.4 <- COMP.4[-which(COMP.4$PLAT=="SEQ" & COMP.4$VS=="CHP"),]
rownames(COMP.4) <- paste(COMP.4$PLAT,COMP.4$VS,sep="v")

## Plot Concordance Results
 # Plot Parameters
YLIM <- c(0,1.2)
YLAB.2 <- "% Concordance (2-digit Precision)"
YLAB.4 <- "% Concordance (4-digit Precision)"
MAIN.2 <- "Concordance b/n Computational Typing Approaches (2-digit)"
MAIN.4 <- "Concordance (4-digit)"
LEG.lab <- rownames(COMP.2)
LEG.main <- "Platform Comparison"
XLAB <- "Gene"
 # Open File
png(paste(PathToPlot,"1C-Conc.Comp.png",sep=""), height=1200,width=1800,pointsize=32 )
par(layout( matrix(1:2,ncol=2),width=c(4,2) ))
 # 2-Digit
TEMP.2 <- barplot( data.matrix(COMP.2[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.2$PLAT)],border=NA,yaxt="n",
	ylim=YLIM,ylab=YLAB.2,xlab=XLAB,main=MAIN.2 )
abline(h=seq(0,1,.2),lty=3,col="grey50",lwd=1)
axis( 2, at=seq(0,1,.2),label=seq(0,100,20),las=2 )
barplot( data.matrix(COMP.2[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.2$PLAT)],border=NA,yaxt="n",add=T )
barplot( data.matrix(COMP.2[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.2$VS)],border=COLS.plat[as.character(COMP.2$VS)],yaxt="n",add=T,density=25,angle=45 )
legend( "topright",fill=COLS.plat[as.character(COMP.2$PLAT)],legend=LEG.lab,title=LEG.main,bg=NA,ncol=3)
legend( "topright",fill=COLS.plat[as.character(COMP.2$VS)],border=COLS.plat[as.character(COMP.2$VS)],density=25,angle=45,legend=LEG.lab,title=LEG.main,bg=NA,ncol=3)
text( TEMP.2[,1], .1, srt=90, label=paste("n=(",COMP.2$N,")",sep="") )
 # 4-Digit
TEMP.4 <- barplot( data.matrix(COMP.4[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.4$PLAT)],border=NA,
	yaxt="n",ylim=YLIM,ylab=YLAB.4,xlab=XLAB,main=MAIN.4 )
abline(h=seq(0,1,.2),lty=3,col="grey50",lwd=1)
axis( 2, at=seq(0,1,.2),label=seq(0,100,20),las=2 )
barplot( data.matrix(COMP.4[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.4$PLAT)],border=NA,yaxt="n",add=T )
barplot( data.matrix(COMP.4[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.4$VS)],border=COLS.plat[as.character(COMP.4$VS)],yaxt="n",add=T,density=25,angle=45 )
text( TEMP.4[,1], .1, srt=90, label=paste("n=(",COMP.4$N,")",sep="") )
dev.off()

#############################################################
## D - CONCORDANCE b/n LAB & COMP TYPES #####################
#############################################################

#######################################
## Compile Table for 20 Patients ######
TYPE.lab <- TYPE.4[SAMP.lab,,]

## Make into tall table
TAB.lab.1 <- Reduce( rbind, lapply(1:5,function(x)data.frame(Gene=dimnames(TYPE.lab)[[3]][x],TYPE.lab[,,x])) )
TAB.lab.2 <- TAB.lab.1[,-grep("SOP.alt",colnames(TAB.lab.1))]
TAB.lab.3 <- apply( TAB.lab.2, 2, function(x)gsub("NA///NA","-",x) )
TAB.lab.4 <- data.frame( Samp=rownames(TAB.lab.3), TAB.lab.3, stringsAsFactors=F )
TAB.lab.4[which(duplicated(TAB.lab.4[,"Gene"])),"Gene"] <- ""

## Save Table & Make LaTeX Output
write.table( TAB.lab.4, file=paste(PathToPlot,"T1-LabSamp.table.csv",sep=""),sep=",",row.names=F,col.names=T,quote=F )
writeLines( print(xtable(TAB.lab.4, include.rownames=F )), con=paste(PathToPlot,"T1-LabSamp.LaTeX.txt",sep="") )

#######################################
## Calculate Concordance b/n Methods ##

## SNP Chip (SNP2HLA) Concordance
CHPvl <- list()
Platform <- "CHP"
Prim_Comp <- c("LAB","LAB.lik")
Samples <- SAMP.lab
 # Genes A, B, C (2-digit)
Genes <- c("A","B","C","DRB","DQB")
CHPvl$p2 <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
 # Genes DRB, DQB (4-digit)
Genes <- c("DRB","DQB")
CHPvl$p4 <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )

## HaplotypeCaller (SNP2HLA) Concordance
SEQvl <- list()
Platform <- "SEQ"
Prim_Comp <- c("LAB","LAB.lik")
Samples <- SAMP.lab
 # Genes A, B, C (2-digit)
Genes <- c("A","B","C","DRB","DQB")
SEQvl$p2 <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
 # Genes DRB, DQB (4-digit)
Genes <- c("DRB","DQB")
Samples <- SAMP.chp
SEQvl$p4 <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )

## HaplotypeCaller (SNP2HLA) Concordance
SOPvl <- list()
Platform <- "SOP"
Prim_Comp <- c("LAB","LAB.lik")
Samples <- SAMP.lab
 # Genes A, B, C (2-digit)
Genes <- c("A","B","C","DRB","DQB")
SOPvl$p2 <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
 # Genes DRB, DQB (4-digit)
Genes <- c("DRB","DQB")
SOPvl$p4 <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )


#######################################
## Compile Concordance Results ########

## 2-digit Precision (All Genes)
VS <- c("LAB","LAB.lik")
COMP.2l.CHPv <- data.frame( PLAT="CHP",VS=VS, N=apply(CHPvl$p2$B[VS,],1,max), CHPvl$p2$A[VS,] )
COMP.2l.SEQv <- data.frame( PLAT="SEQ",VS=VS, N=apply(SEQvl$p2$B[VS,],1,max), SEQvl$p2$A[VS,] )
COMP.2l.SOPv <- data.frame( PLAT="SOP",VS=VS, N=apply(SOPvl$p2$B[VS,],1,max), SOPvl$p2$A[VS,] )
COMP.2l <- rbind( COMP.2l.CHPv[1,], COMP.2l.SEQv[1,], COMP.2l.SOPv[1,] )
rownames(COMP.2l) <- paste(COMP.2l$VS,COMP.2l$PLAT,sep="v")

## 4-digit Precision (DRB,DQB)
VS <- c("LAB","LAB.lik")
COMP.4l.CHPv <- data.frame( PLAT="CHP",VS=VS, N=apply(CHPvl$p4$B[VS,],1,max), CHPvl$p4$A[VS,] )
COMP.4l.CHPv.2 <- colSums( COMP.4l.CHPv[,-(1:3)],na.rm=T )
COMP.4l.SEQv <- data.frame( PLAT="SEQ",VS=VS, N=apply(SEQvl$p4$B[VS,],1,max), SEQvl$p4$A[VS,] )
COMP.4l.SEQv.2 <- colSums( COMP.4l.SEQv[,-(1:3)],na.rm=T )
COMP.4l.SOPv <- data.frame( PLAT="SOP",VS=VS, N=apply(SOPvl$p4$B[VS,],1,max), SOPvl$p4$A[VS,] )
COMP.4l.SOPv.2 <- colSums( COMP.4l.SOPv[,-(1:3)],na.rm=T )
TEMP.N <- COMP.2l$N
COMP.4l <- rbind( COMP.4l.CHPv.2, COMP.4l.SEQv.2, COMP.4l.SOPv.2 )
COMP.4l <- data.frame( PLAT=c("CHP","SEQ","SOP"),VS="LAB", N=TEMP.N, COMP.4l )
rownames(COMP.4l) <- paste(COMP.4l$VS,COMP.4l$PLAT,sep="v")

## Plot Concordance Results
 # Plot Parameters
YLIM <- c(0,1.2)
YLAB.2 <- "% Concordance (2-digit Precision)"
YLAB.4 <- "% Concordance (4-digit Precision)"
MAIN.2 <- "Concordance w/ Lab Typing (2-digit)"
MAIN.4 <- "Concordance (4-digit)"
LEG.lab <- as.character(COMP.2l$PLAT)
LEG.main <- "vs Lab"
XLAB <- "Gene"
 # Open File
png(paste(PathToPlot,"1D-Conc.Lab.png",sep=""), height=1200,width=1800,pointsize=32 )
par(layout( matrix(1:2,ncol=2),width=c(4,2) ))
 # 2-Digit
TEMP.2 <- barplot( data.matrix(COMP.2l[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.2l$PLAT)],border=NA,yaxt="n",
	ylim=YLIM,ylab=YLAB.2,xlab=XLAB,main=MAIN.2 )
abline(h=seq(0,1,.2),lty=3,col="grey50",lwd=1)
axis( 2, at=seq(0,1,.2),label=seq(0,100,20),las=2 )
barplot( data.matrix(COMP.2l[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.2l$PLAT)],border=NA,yaxt="n",add=T )
barplot( data.matrix(COMP.2l[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.2l$VS)],border=COLS.plat[as.character(COMP.2l$VS)],yaxt="n",add=T,density=25,angle=45 )
legend( "topright",fill=COLS.plat[as.character(COMP.2l$PLAT)],legend=LEG.lab,title=LEG.main,bg=NA,ncol=3)
legend( "topright",fill=COLS.plat[as.character(COMP.2l$VS)],border=COLS.plat[as.character(COMP.2l$VS)],density=25,angle=45,legend=LEG.lab,title=LEG.main,bg=NA,ncol=3)
text( TEMP.2[,1], .1, srt=90, label=paste("n=(",COMP.2l$N,")",sep="") )
 # 4-Digit
TEMP.4 <- barplot( data.matrix(COMP.4l[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.4l$PLAT)],border=NA,
	yaxt="n",ylim=YLIM,ylab=YLAB.4,xlab=XLAB,main=MAIN.4 )
abline(h=seq(0,1,.2),lty=3,col="grey50",lwd=1)
axis( 2, at=seq(0,1,.2),label=seq(0,100,20),las=2 )
barplot( data.matrix(COMP.4l[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.4l$PLAT)],border=NA,yaxt="n",add=T )
barplot( data.matrix(COMP.4l[,-(1:3)]), beside=T, col=COLS.plat[as.character(COMP.4l$VS)],border=COLS.plat[as.character(COMP.4l$VS)],yaxt="n",add=T,density=25,angle=45 )
text( TEMP.4[,1], .1, srt=90, label=paste("n=(",COMP.4l$N,")",sep="") )
dev.off()

# #############################################################
# ## E - HAPLOTYPE COUNTS for SOAP-HLA (by Gene) ##############
# #############################################################

# ## Get Haplotype Counts by Gene
# HAP_CNT.4 <- apply( TYPE.4[,"SOP",], 2, function(x) table( unlist(strsplit(x,"///")) ) )
# HAP_CNT.2 <- apply( TYPE.2[,"SOP",], 2, function(x) table( unlist(strsplit(x,"///")) ) )

# ## FCT: Plot Aggregate Haplotype Counts for each Gene
# Plot_HapCnt <- function( gene ) {
# 	png(paste(PathToPlot,"2B-",gene,".png",sep=""), height=1400,width=2000,pointsize=28 )
# 	par(mfrow=c(2,1))
# 	## 4-digit
# 	YLIM <- c( 0, max(HAP_CNT.4[[gene]]) ) * c(1,1.3)
# 	TEMP.4 <- barplot( HAP_CNT.4[[gene]], col=COLS.plat["SOP"],ylim=YLIM,
# 		main=paste("Aggregate Haplotype Count (SOAP-HLA 4-digit): HLA-",gene,sep=""),ylab="# Haplotypes",xlab="Haplotype",las=2 )
# 	abline( h=seq(0,400,50),lty=3,col="grey50",lwd=1 )
# 	TEMP.4 <- barplot( HAP_CNT.4[[gene]], col=COLS.plat["SOP"],las=2,add=T )
# 	text( TEMP.4, HAP_CNT.4[[gene]]+YLIM[2]/8, label=paste("(n=",HAP_CNT.4[[gene]],")",sep=""),srt=90,cex=.9 )
# 	## 2-digit
# 	YLIM <- c( 0, max(HAP_CNT.2[[gene]]) ) * c(1,1.3)
# 	TEMP.2 <- barplot( HAP_CNT.2[[gene]], col=COLS.plat["SOP"],ylim=YLIM,
# 		main=paste("Aggregate Haplotype Count (SOAP-HLA 2-digit): HLA-",gene,sep=""),ylab="# Haplotypes",xlab="Haplotype",las=2 )
# 	abline( h=seq(0,400,50),lty=3,col="grey50",lwd=1 )
# 	TEMP.2 <- barplot( HAP_CNT.2[[gene]], col=COLS.plat["SOP"],las=2,add=T )
# 	text( TEMP.2, HAP_CNT.2[[gene]]+YLIM[2]/8, label=paste("(n=",HAP_CNT.2[[gene]],")",sep=""),srt=90,cex=.9 )
# 	dev.off()
# }
# # Plot_HapCnt("A")
# lapply( GENES, Plot_HapCnt )



#############################################################
## END OF DOC ###############################################
#############################################################
