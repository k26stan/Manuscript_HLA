## Figure 2 for HLA Manuscript ##
## Plot Frequencies of Types, AAs, and Collapsed Haplotypes ##
## April 5, 2016 ##
## Kristopher Standish ##

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

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to HLA Data Sets
PathToTypes <- "/Users/kstandis/Data/Janssen/Data/HLA/SOAP_HLA_Types/20151211_HLA_Types.Rdata"
PathToAA <- paste("/Users/kstandis/Data/Janssen/Data/HLA/Amino_Acids/20160126_HLA_AA.Rdata")
PathTo1KG <- "/Users/kstandis/Data/Genetics/HLA/1KG/20140702_hla_diversity.txt"
PathTo1KG.2 <- "/Users/kstandis/Data/Genetics/1KG/Panel_Key.txt"
PathToRefs <- "/Users/kstandis/Data/Genetics/HLA/Alignments_Rel_3190/"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
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
## PLOT HLA TYPE FREQUENCIES ################################
#############################################################

## FCT: Plot Allele Frequency of Each Haplotype
Plot_AF <- function( Samples, tag ) {
	FREQS <- list()
	for ( g in 1:N.GENE ) {
		gene <- GENE_LIST[g]
		FREQ.4 <- rowSums( PAT_DOS.4[[gene]][,Samples] ) ; names(FREQ.4)[which(names(FREQ.4)==":")]<-"NA" ; FREQ.4 <- FREQ.4[order(names(FREQ.4))]
		names(FREQ.4) <- gsub(":","",names(FREQ.4),fixed=T)
		FREQ.4 <- FREQ.4[ which(names(FREQ.4)!="NA") ]
		FREQ.2 <- rowSums( PAT_DOS.2[[gene]][,Samples] ) ; FREQ.2 <- FREQ.2[order(names(FREQ.2))]
		## Plot it
		# 4 digit
		YLIM <- c( 0, max(FREQ.4) ) * c(1,1.25)
		LINES <- which(!duplicated( substr(names(FREQ.4),1,2) ))
		png( paste(PathToPlot,tag,"-2B-TypeFreq.Pr4.",gene,".png",sep=""), height=800,width=2000,pointsize=30 )
		# par(mfrow=c(2,1))
		par(mar=c(4,4,3,2))
		TEMP.4 <- barplot( FREQ.4, col=COLS.fr,border=NA,ylim=YLIM,
			main=paste("Aggregate Haplotype Count (SOAP-HLA 4-digit): HLA-",gene,sep=""),ylab="# Haplotypes",xlab="Haplotype",las=2 )
		abline( h=seq(0,400,50),v=TEMP.4[LINES]-.6,lty=3,col="grey50",lwd=1 )
		TEMP.4 <- barplot( FREQ.4, col=COLS.fr,border=NA,las=2,add=T )
		text( TEMP.4, FREQ.4+YLIM[2]*.1, label=paste("(n=",FREQ.4,")",sep=""),srt=90,cex=.9 )
		dev.off()
		# 2 digit
		YLIM <- c( 0, max(FREQ.2) ) * c(1,1.25)
		png( paste(PathToPlot,tag,"-2B-TypeFreq.Pr2.",gene,".png",sep=""), height=800,width=1600,pointsize=30 )
		par(mar=c(4,4,3,2))
		TEMP.2 <- barplot( FREQ.2, col=COLS.fr,border=NA,ylim=YLIM,
			main=paste("Aggregate Haplotype Count (SOAP-HLA 2-digit): HLA-",gene,sep=""),ylab="# Haplotypes",xlab="Haplotype",las=2 )
		abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
		TEMP.2 <- barplot( FREQ.2, col=COLS.fr,border=NA,las=2,add=T )
		text( TEMP.2, FREQ.2+YLIM[2]*.1, label=paste("(n=",FREQ.2,")",sep=""),srt=90,cex=.9 )
		dev.off()
		FREQS[[gene]] <- list( FREQ.4=FREQ.4, FREQ.2=FREQ.2 )
	}
	return(FREQS)
}
 # Run Function on Various Sample Lists
FREQS <- list()
FREQS <- lapply( names(SAMPS), function(x) Plot_AF(SAMPS[[x]],x) )
names(FREQS) <- names(SAMPS)

#############################################################
## PLOT AMINO ACID FREQUENCIES ##############################
#############################################################

## Plot AA Level Diversity within Cohort
 # Functions used in calculatiing AA Diversity
Find_Div <- function( aa_vector ) {
	TAB <- table( aa_vector )
	To_Ignore <- which(names(TAB) %in% c("?","*",".") )
	if ( length(To_Ignore)>0 ) { TAB <- TAB[-To_Ignore] }
	return(TAB)
}
allele_bar <- function(pat_tab.prc,pat_tab.names,x) {
	xval <- as.numeric(colnames(pat_tab.prc)[x])
	yvals <- cumsum(pat_tab.prc[,x])
	barplot( t(t(pat_tab.prc[,x])),add=T,xaxt="n",beside=F,border=NA,col=COLS.AA[pat_tab.names[,x]],space=c(xval,0),width=1)
}

Plot_AAF <- function( Samples, tag ) {
	Samples.typ <- paste( rep(Samples,each=2),rep(1:2,length(Samples)),sep="_" )
	# AA_FREQS <- list()
	for ( g in 1:N.GENE ) {
		gene <- GENE_LIST[g]
		# Pull out Patient Data for Gene
		pat_typ <- PAT_TYP[gene,Samples.typ]
		pat_dos <- PAT_DOS[[gene]][,Samples]
		pat_aa.r <- pat_aa <- PAT_AA[[gene]][Samples.typ,] # Raw Table
		ISNA <- which(is.na(pat_aa.r))
		if ( length(ISNA)>0 ) { pat_aa[ISNA] <- "?" } # After converting missing values to "?"

		## Get Amino Acid Frequencies
		 # What is the maximum number of AA's at one position?
		max.temp <- max( unlist(lapply( apply( pat_aa, 2, table ), length)) )
		 # Which Positions have >1 Known AA? (Use Raw Table before "?" were included)
		which.temp <- which( unlist(lapply( apply( pat_aa, 2, Find_Div ), length)) > 1 )
		 # Get AA Frequencies/Proportions & Names
		pat_tab <- sapply(apply( pat_aa, 2, table ), "[", 1:max.temp ) ; colnames(pat_tab) <- gsub("Pos_","",colnames(pat_tab))
		pat_tab.names <- sapply(apply( pat_aa, 2, function(x) names(table(x)) ), "[", 1:max.temp )
		pat_tab.prc <- pat_tab / nrow(pat_aa)
		 # Plot
		XLIM <- range(as.numeric(colnames(pat_tab)))
		png( paste(PathToPlot,tag,"-2C-AAfreq_",gene,".png",sep=""), height=1000,width=2000,pointsize=30 )
		layout( matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(7,1) ) # layout( matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(3,2), heights=c(1,1) )
		plot( 0,0,type="n",xlim=XLIM,ylim=c(0,1), main=paste(tag,": Amino Acid Frequency - HLA-",gene,sep=""),xlab="Amino Acid",ylab="Frequency",xaxt="n")
		axis( 1, at=seq(-1000,1000,20),las=2 )	
		# abline( v=seq(-1000,1000,20),lty=3,col="grey50" )
		abline( h=seq(0,1,.2),lty=3,col="grey50" ) ; abline( h=c(0,1) )
		DWAI <- lapply( which.temp, function(x)allele_bar(pat_tab.prc,pat_tab.names,x) )
		## Amino Acid Key
		Which_AA <- which(names(COLS.AA)%in% pat_aa )
		COLS.aa <- COLS.AA[Which_AA]
		barplot( matrix(rep(1,length(COLS.aa)),ncol=1), beside=F,col=COLS.aa,xaxt="n",yaxt="n",ylab="AA")
		axis( 2, at=1:length(COLS.aa)-.5, label=names(COLS.aa),las=2,cex.axis=.8 )
		dev.off()
		# AA_FREQS[[gene]] <- list( )
	}
}
 # Run Function on Various Sample Lists
SCRAP <- lapply( names(SAMPS), function(x) Plot_AAF(SAMPS[[x]],x) )

#############################################################
## COLLAPSED HAPLOTYPE FREQUENCIES ##########################
#############################################################

## Function to do Haplotype Analysis of Specified Amino Acid Positions
HAP_AN <- function( Positions, tag ) {
	## Pull together Haplotypes from Specified Positions
	N.res <- length(Positions)
	HAP <- apply( PAT_AA$DRB1[,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
	HAP <- gsub("NA","-",HAP)
	HAP.uniq <- sort(unique(HAP))
	HAP.plot <- setdiff( HAP.uniq, paste(rep("-",N.res),collapse=""))
	 # Get Haplotype Frequencies
	HAP.freq <- table(HAP)
	HAP.freq <- HAP.freq[ HAP.plot ]
	#  # Convert to Array for Additive Analysis
	# HAP.arr <- array( 0,c(N.PATS,length(HAP.uniq)) )
	# colnames(HAP.arr) <- HAP.uniq
	# rownames(HAP.arr) <- PATS
	# for ( pat in PATS ) {
	# 	HAP.pat <- HAP[ grep(pat,names(HAP)) ]
	# 	if ( HAP.pat[1]==HAP.pat[2] ) { HAP.arr[pat,HAP.pat[1]] <- 2
	# 	}else{ HAP.arr[pat,HAP.pat] <- 1 }
	# }
	# HAP.arr <- HAP.arr[ SAMPS$ART3001, HAP.plot ]

	## Plot Haplotype Frequency
	YLIM <- c( 0, max(HAP.freq) ) * c(1,1.35)
	png( paste(PathToPlot,"/DRB1_",tag,"_2D-HapFreq.png",sep=""), height=800,width=1000,pointsize=30 )
	TEMP <- barplot( HAP.freq,las=2,col=COLS.fr,border=NA,main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype Frequency"),ylab="# Haplotypes",ylim=YLIM)
	abline(h=seq(0,1000,50),lty=3,col="grey50")
	barplot( HAP.freq,las=2,col=COLS.fr,border=NA,add=T )
	text( TEMP, HAP.freq+.13*YLIM[2], label=paste("(n=",HAP.freq,")",sep=""), srt=90 )
	dev.off()
}

OUT <- list()
##########################################
## POS 11, 71, 74 ##
 # Viatte, et al (2015)
Positions <- c(11,71,74)
tag <- paste(c("p",Positions),collapse="")
OUT[[tag]] <- HAP_AN(Positions,tag)

##########################################
## POS 11, 13, 71, 74 ##
Positions <- 70:74
tag <- "pSE"
OUT[[tag]] <- HAP_AN(Positions,tag)






#############################################################
## END OF DOC ###############################################
#############################################################



















# library(xlsx)
# library(stringr)
# library(gplots)
# library(xtable)
# ## Set Date
# DATE <- gsub("-","",Sys.Date())

# ## Set Paths to Data Sets & Save Locations (Mac)
# PathToSOAP <- "/Users/kstandis/Data/Janssen/Data/HLA/SOAP_HLA_Types/20151211_HLA_Types.Rdata"
# PathToS2H <- "/Users/kstandis/Data/Janssen/Data/HLA/SNP2HLA_Types/20160119_Tables/"
# PathToS2H <- "/Users/kstandis/Data/Janssen/Data/HLA/SNP2HLA_Types/20160323_Tables/"
# PathToLAB <- "/Users/kstandis/Data/Janssen/Data/HLA/Lab_Types/LabCorp HLA typing data_ART3001.xls"
# PathToComp <- "/Users/kstandis/Data/Janssen/Data/HLA/Compare/20160323_Tables/1-CompiledTypes.Rdata"
# PathToPheno <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20150520_Full_Table.txt"
# PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_CompareHLA/",sep="")
# PathToRefs <- "/Users/kstandis/Data/Genetics/HLA/Alignments_Rel_3190/"
# if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }

# ## Load Compiled HLA Types
# load( PathToComp )
# TYPES <- TO_SAVE

# ## Load Individual HLA Types
#  # SOAP-HLA
# load( PathToSOAP )
# SOP.l <- COMPILE
#  # SNP2HLA
# CHP.2.l <- read.table( paste(PathToS2H,"CHP_2.txt",sep=""),sep="\t",header=T,colClasses="character" )
# CHP.4.l <- read.table( paste(PathToS2H,"CHP_4.txt",sep=""),sep="\t",header=T,colClasses="character" )
# SEQ.2.l <- read.table( paste(PathToS2H,"SNP_2.txt",sep=""),sep="\t",header=T,colClasses="character" )
# SEQ.4.l <- read.table( paste(PathToS2H,"SNP_4.txt",sep=""),sep="\t",header=T,colClasses="character" )
#  # Lab Typing
# LAB.l <- read.xlsx( PathToLAB, sheetIndex=1, rowIndex=1:101, header=T, colIndex=1:11)

# ## Load Phenotypes & Sample Lists
# FT <- read.table( PathToPheno, sep="\t",header=T )

# ## Load Amino Acid Data
# FILES <- list.files(PathToRefs)
# gsub("_prot.txt","",FILES[grep("_prot",FILES)])
