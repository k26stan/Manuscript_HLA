## Table 1 for HLA Manuscript ##
## Summary Baseline Stats for Cohort ##
## November 23, 2016 ##
## Kristopher Standish ##

#############################################################
## LOAD DATA ################################################
#############################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to HLA Data Sets
PathToTypes <- "/Users/kstandis/Data/Janssen/Data/HLA/SOAP_HLA_Types/20151211_HLA_Types.Rdata"
PathToAA <- paste("/Users/kstandis/Data/Janssen/Data/HLA/Amino_Acids/20160126_HLA_AA.Rdata")
PathTo1KG <- "/Users/kstandis/Data/Genetics/HLA/1KG/20140702_hla_diversity.txt"
PathTo1KG.2 <- "/Users/kstandis/Data/Genetics/1KG/Panel_Key.txt"
PathToRefs <- "/Users/kstandis/Data/Genetics/HLA/Alignments_Rel_3190/"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
PathToFT.3002 <- "/Users/kstandis/Data/Janssen/Data/Pheno/Raw_Files/20160315_ART3002/ART3002_Clinical.csv"
PathToSave <- paste("/Users/kstandis/Data/Janssen/Data/HLA/Association/",DATE,"_HLA_Assoc",sep="")
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_ManuHLA_Fig2/",sep="")
dir.create( PathToPlot )

## Load Janssen HLA Results
 # Types
load( PathToTypes )
TYPES.l <- COMPILE
 # Amino Acids
load( PathToAA )
AA.l <- COMPILE

## 1000 Genomes HLA Data
HLA.l <- read.table(PathTo1KG,header=T)
HLA.key <- read.table(PathTo1KG.2,header=T)
HLA <- merge( HLA.key, HLA.l, by.y="id",by.x="sample")

## Get Phenotype Info
FT.l <- read.table( PathToFT, sep="\t",header=T )
FT.3002 <- read.csv( PathToFT.3002 )
colnames(FT.3002)[1] <- "ID"

#############################################################
## ORGANIZE DATA ############################################
#############################################################

## Filter Phenotypic Data (based on length of participation in study)
RM.LT8 <- which( FT.l$IN < 8 )
RM.LT8.samps <- as.character( FT.l$ID[RM.LT8] )
FT <- FT.l[ -RM.LT8, ]

## Specify Samples
SAMPS <- list()
SAMPS$ART3001 <- as.character(FT$ID)
SAMPS$ART3002 <- grep( "-",colnames(PAT_DOS[[1]]),value=T ) # setdiff(colnames(PAT_DOS.4[[1]]),ART3001)
SAMPS$ALL <- colnames(PAT_DOS[[1]])

## Pull out Patient Data
PAT_DOS <- TYPES.l$GENES.2.list
PAT_TYP <- TYPES.l$GENES
PAT_AA <- AA.l$PAT_AA
 # Print Example
X <- 1:10
PAT_AA$DRB1[X,X]

## Patients & Genes
PATS <- colnames(PAT_DOS$A)
N.PATS <- length(PATS)
GENE_LIST.all <- names(AA.l$PAT_AA)
N.GENE.all <- length(GENE_LIST.all) # nrow(PAT_TYP)
 # Specify Important Genes
GENE_LIST <- GENE_LIST.all[1:5]
GENE_LIST <- GENE_LIST.all[c(1:5,7)]
N.GENE <- length(GENE_LIST)

## Cut Patient Types to 2 & 4 Digit Precision
 # PAT_TYP
PAT_TYP.4 <- PAT_TYP.2 <- PAT_TYP
PAT_TYP.4 <- apply( PAT_TYP, 1, function(x) apply( sapply( strsplit(x,":"), "[",1:2),2, function(y) paste(y,collapse=":") ) )
PAT_TYP.4 <- apply( PAT_TYP.4, 1, function(x) gsub("[A-z]","",x) )
PAT_TYP.2 <- apply( PAT_TYP, 1, function(x) sapply( strsplit(x,":"), "[",1) )
PAT_TYP.2 <- apply( PAT_TYP.2, 1, function(x) gsub("[A-z]","",x) )
 # PAT_DOS
PAT_DOS.4 <- PAT_DOS.2 <- list()
for ( g in 1:N.GENE ) { gene <- GENE_LIST[g]
	unique_haps <- unique( PAT_TYP.4[gene,] )
	PAT_DOS.4[[gene]] <- array( 0,c(length(unique_haps),N.PATS) )
	colnames(PAT_DOS.4[[gene]]) <- PATS ; rownames(PAT_DOS.4[[gene]]) <- unique_haps
	unique_haps <- na.omit( unique( PAT_TYP.2[gene,] ) )
	PAT_DOS.2[[gene]] <- array( 0,c(length(unique_haps),N.PATS) )
	colnames(PAT_DOS.2[[gene]]) <- PATS ; rownames(PAT_DOS.2[[gene]]) <- unique_haps
	for ( p in 1:N.PATS ) { pat <- PATS[p]
		which_haps <- PAT_TYP.4[ gene,paste(pat,1:2,sep="_") ]
		if ( length(which(duplicated(which_haps)))>0 ) {
			PAT_DOS.4[[gene]][which_haps,pat] <- 2
		}else{
			PAT_DOS.4[[gene]][which_haps,pat] <- 1
		}
		which_haps <- na.omit( PAT_TYP.2[ gene,paste(pat,1:2,sep="_") ] )
		if ( length(which(duplicated(which_haps)))>0 ) {
			PAT_DOS.2[[gene]][which_haps,pat] <- 2
		}else{
			PAT_DOS.2[[gene]][which_haps,pat] <- 1
		}
	}
}

#############################################################
## GET SUMMARY STATS ########################################
#############################################################

## Demographic Info
 # Sex
D.SEX.cnt <- table( FT$SEX )
D.SEX.prc <- 100*D.SEX.cnt / nrow(FT)
 # Age
D.AGE.mn <- mean( FT$AGE )
D.AGE.sd <- sd( FT$AGE )
D.AGE.qrt <- quantile( FT$AGE, c(0.25,0.75) )
 # Height
D.HT.mn <- mean( FT$HT )
D.HT.sd <- sd( FT$HT )
D.HT.qrt <- quantile( FT$HT, c(0.25,0.75) )
 # Weight
D.WT.mn <- mean( FT$WT )
D.WT.sd <- sd( FT$WT )
D.WT.qrt <- quantile( FT$WT, c(0.25,0.75) )
 # BMI
D.BMI.mn <- mean( FT$BMI )
D.BMI.sd <- sd( FT$BMI )
D.BMI.qrt <- quantile( FT$BMI, c(0.25,0.75) )
 # Disease Duration
D.DIS_DUR.mn <- mean( FT$DIS_DUR )
D.DIS_DUR.sd <- sd( FT$DIS_DUR )
D.DIS_DUR.qrt <- quantile( FT$DIS_DUR, c(0.25,0.75) )
 # Race (Percent White/Non-Hisp)
D.EUR.cnt <- table( FT$RACE,FT$ETHN ) # ["WHITE","NOT HISPANIC OR LATINO"]
D.EUR.prc <- 100*D.EUR.cnt / nrow(FT)

## Baseline Disease State
 # DAS
BL.DAS.mn <- mean( FT$DAS_0wk )
BL.DAS.sd <- sd( FT$DAS_0wk )
BL.DAS.qrt <- quantile( FT$DAS_0wk, c(0.25,0.75) )
 # SJC
BL.SJC.mn <- mean( FT$SJC_0wk )
BL.SJC.sd <- sd( FT$SJC_0wk )
BL.SJC.qrt <- quantile( FT$SJC_0wk, c(0.25,0.75) )
 # TJC
BL.TJC.mn <- mean( FT$TJC_0wk )
BL.TJC.sd <- sd( FT$TJC_0wk )
BL.TJC.qrt <- quantile( FT$TJC_0wk, c(0.25,0.75) )
 # CRP
BL.CRP.mn <- mean( FT$CRP_0wk )
BL.CRP.sd <- sd( FT$CRP_0wk )
BL.CRP.qrt <- quantile( FT$CRP_0wk, c(0.25,0.75) )
 # ACPA
BL.ACPA.cnt <- table( FT$ACPA )
BL.ACPA.prc <- 100*BL.ACPA.cnt / nrow(FT)
 # RF
BL.RF.cnt <- table( FT$RF )
BL.RF.prc <- 100*BL.RF.cnt / nrow(FT)

ACPA/RF









