## Figure 3 for HLA Manuscript ##
## Plot Association Analyses for HLA Types/Amino Acids ##
## March 7, 2016 ##
## Kristopher Standish ##

## Brainstorming
 # Phenotypes
   # Drug Response (delta-DAS,lCRP,rSJC,rTJC)
   # Disease Severity
   # RF/ACPA Status
 # HLA Predictors
   # ANOVA over all Types for each Gene
     # 2- & 4-Digit Type  
   # Additive/Dominant/Recessive for each Type for each Gene
     # 2- & 4-Digit Type
   # Amino-Acid Level Tests for each Gene (4- or Best-Digit Precision)
     # ANOVA over all Amino Acids at each Position
     # Additive/Dominant/Recessive for each Type for each Gene
   # Collapsed Amino-Acid Haplotypes
     # DRB1
       # Pos 11,71,74
       # Pos 70-74 (Shared Epitope)
       # Pos 11,13,71,74
     # B
       # Pos 9
     # DPB1
       # Pos 9
   # Collapsed Sets of Collapsed Amino-Acid Haplotypes
     # (Viatte 2014)

library(gplots)

#############################################################
## GLOBALS ##################################################
#############################################################

## Specify Phenos/Covs
 # Response Phenotypes
PHENOS <- c("DEL_MNe_MN","DEL_lCRP_MNe_MN","DEL_rSJC_MNe_MN","DEL_rTJC_MNe_MN",
	"DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA","RF" )
 # Covariates (BMI + RF + ACPA + disease duration)
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"BMI+ACPA+RF+log(DIS_DUR)","BMI+ACPA+RF+log(DIS_DUR)","BMI+ACPA+RF+log(DIS_DUR)","BMI+ACPA+RF+log(DIS_DUR)",
	"","")
 # Covariates (RF + ACPA + disease duration)
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA+RF+log(DIS_DUR)","ACPA+RF+log(DIS_DUR)","ACPA+RF+log(DIS_DUR)","ACPA+RF+log(DIS_DUR)",
	"","")
 # Covariates (ACPA + disease duration)
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)",
	"","")
 # Covariates (ACPA)
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA","ACPA","ACPA","ACPA",
	"","")

 # Look at Fewer Phenotypes
PH.COV.1 <- data.frame(PHENOS,COVS)
PH.COV <- PH.COV.1[c(1,5,9),]
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
PathToAA <- "/Users/kstandis/Data/Janssen/Data/HLA/Amino_Acids/20160126_HLA_AA.Rdata"
PathTo1KG <- "/Users/kstandis/Data/Genetics/HLA/1KG/20140702_hla_diversity.txt"
PathTo1KG.2 <- "/Users/kstandis/Data/Genetics/1KG/Panel_Key.txt"
PathToRefs <- "/Users/kstandis/Data/Genetics/HLA/Alignments_Rel_3190/"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
PathToAssocP <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160503_HLA_Assoc_0-P_Precise.Rdata"
PathToAssocB <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160503_HLA_Assoc_0-B_Precise.Rdata"
PathToHLA <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/"
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_ManuHLA_Fig3/",sep="")
dir.create( PathToPlot )

## Load HLA Association Results
load(file=PathToAssocP)
summary(P.out)
load(file=PathToAssocB)
summary(B.out)

## Load Janssen HLA Results
 # Types
load( PathToTypes )
TYPES.l <- COMPILE
 # Amino Acids
load( PathToAA )
AA.l <- COMPILE

## Load HLA Types
HLA_AA.l <- read.table(paste(PathToHLA,"20160208_HLA_Assoc_HLA_AA_Table.txt",sep=""),header=T,sep="\t" )
HLA_TYP.l <- read.table(paste(PathToHLA,"20160208_HLA_Assoc_HLA_Types_Table.txt",sep=""),header=T,sep="\t" )
colnames(HLA_AA.l) <- gsub(".","_",colnames(HLA_AA.l),fixed=T)
colnames(HLA_AA.l) <- gsub("-","_",colnames(HLA_AA.l),fixed=T)
colnames(HLA_TYP.l) <- gsub(".","_",colnames(HLA_TYP.l),fixed=T)
colnames(HLA_TYP.l) <- gsub("-","_",colnames(HLA_TYP.l),fixed=T)

## Get Clinical Phenotype Info
FT.l <- read.table( PathToFT, sep="\t",header=T )

#############################################################
## ORGANIZE DATA ############################################
#############################################################

## Filter Phenotypic Data (based on length of participation in study)
RM.LT8 <- which( FT.l$IN < 8 )
RM.LT8.samps <- as.character( FT.l$ID[RM.LT8] )
FT <- FT.l[ -RM.LT8, ]

## Specify ACPA as 1/0
FT$ACPA <- as.numeric(factor(FT$ACPA=="Positive")) - 1

## Specify Samples
SAMPS <- list()
SAMPS$ART3001 <- as.character(FT$ID)
SAMPS$ART3002 <- grep( "-",colnames(TYPES.l$GENES.2.list$A),value=T ) # setdiff(colnames(PAT_DOS.4[[1]]),ART3001)
SAMPS$ALL <- colnames(TYPES.l$GENES.2.list$A)

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

## Merge Clinical w/ HLA Data
MG.aa <- merge( FT, HLA_AA.l, by="ID" )
MG.typ <- merge( FT, HLA_TYP.l, by="ID" )

#############################################################
## PLOT NULL MODELS #########################################
#############################################################

## Calculate Null (Non-Genetic) Models
FORM.null <- list()
MOD.null <- list()
for ( p in 1:2 ) {
	pheno <- PHENOS[p]
	covs <- COVS[p]
	FORM.null[[pheno]] <- paste(pheno,"~",covs)
	MOD.null[[pheno]] <- lm( as.formula(FORM.null[[pheno]]), data=FT )
}

## Pull out Results
B.null <- lapply( MOD.null,function(x)summary(x)$coefficients[,"Estimate"] )
SE.null <- lapply( MOD.null,function(x)summary(x)$coefficients[,"Std. Error"] )
P.null <- lapply( MOD.null,function(x)summary(x)$coefficients[,"Pr(>|t|)"] )
## Plot Results from Null Models
png( paste(PathToPlot,"Null_1-CovarBetas.png",sep=""),height=1000,width=2000,pointsize=30 )
layout( matrix(1:2,ncol=2),widths=c(2.5,1) )
 # Severity
par(mar=c(7,4,4,1))
YLIM <- c( 0, max(B.null$DAS_BL_MN) ) * c(1,1.3)
GAP <- .1*diff(YLIM)
TEMP <- barplot( B.null$DAS_BL_MN, col=COLS.beta,border=NA,las=2,ylim=YLIM,main="Clinical Model: Disease Severity",ylab="Effect Size" )
abline( h=0:10,lty=3,col="grey50",lwd=1 )
barplot( B.null$DAS_BL_MN, col=COLS.beta,border=NA,las=2,add=T )
arrows( TEMP, B.null$DAS_BL_MN-SE.null$DAS_BL_MN, TEMP, B.null$DAS_BL_MN+SE.null$DAS_BL_MN, lwd=4,length=0 )
text( TEMP, B.null$DAS_BL_MN+GAP*sign(B.null$DAS_BL_MN)+GAP/3, label=paste("B =",round(B.null$DAS_BL_MN,3)) )
text( TEMP, B.null$DAS_BL_MN+GAP*sign(B.null$DAS_BL_MN)-GAP/3, label=paste("p =",formatC(P.null$DAS_BL_MN,format="e",digits=2)) )
 # Response
par(mar=c(7,4,4,1))
YLIM <- range(B.null$DEL_MNe_MN) * c(1.7,2.5)
GAP <- .1*diff(YLIM)
TEMP <- barplot( B.null$DEL_MNe_MN, col=COLS.beta,border=NA,las=2,ylim=YLIM,main="Clinical Model: Response",ylab="Effect Size" )
abline( h=seq(-5,5,.25),lty=3,col="grey50",lwd=1 )
barplot( B.null$DEL_MNe_MN, col=COLS.beta,border=NA,las=2,add=T )
arrows( TEMP, B.null$DEL_MNe_MN-SE.null$DEL_MNe_MN, TEMP, B.null$DEL_MNe_MN+SE.null$DEL_MNe_MN, lwd=4,length=0 )
text( TEMP, B.null$DEL_MNe_MN+GAP*sign(B.null$DEL_MNe_MN)+GAP/3, label=paste("B =",round(B.null$DEL_MNe_MN,3)) )
text( TEMP, B.null$DEL_MNe_MN+GAP*sign(B.null$DEL_MNe_MN)-GAP/3, label=paste("p =",formatC(P.null$DEL_MNe_MN,format="e",digits=2)) )
dev.off()

#############################################################
## PLOT ASSOCIATION RESULTS #################################
#############################################################

## FCT: Modified Nyhold method of calculating Effective Hypotheses
NYHOLT <- function( PRED_TAB ) {
    mcor <- cor(PRED_TAB,method="pearson",use='na.or.complete')
    zz.which <- lapply( 1:ncol(mcor), function(x)length(which( mcor[x,-x]==1 )) )
    zz <- length(which(zz.which>0))
    # zz <- length(which(mcor[upper.tri(mcor)]==1))
	PRED_TAB.2 <- PRED_TAB
    iter <- 0
    while ( zz>0 & iter<1e3 ) {
    	PRED_TAB.2 <- PRED_TAB.2[ , -which.max(zz.which) ]
    	mcor <- cor(PRED_TAB.2,method="pearson",use='na.or.complete')
	    zz.which <- lapply( 1:ncol(mcor), function(x)length(which( mcor[x,-x]==1 )) )
	    zz <- length(which(zz.which>0))
	    iter <- iter+1
    }
    eig <- eigen(mcor,T,only.values=T)
    m <- length(eig$values)
    meff <- 1 + (m-1)*(1-var(eig$values)/m)
    return( meff )
}

## FCT: Pull out & Plot Results for a given Gene
 # (Updated Function w/ Single Plot)
PLOT_RESULTS <- function( P.pr, B.pr, gene, pheno_cov_table, tag ) {
	## PULL RESULTS ##################

	## Names/Tags/Parameters
	PHENOS <- as.character(pheno_cov_table[,"PHENOS"])
	N.PHENOS <- length(PHENOS)
	tag.temp <- paste(tag,gene,sep="_")

	## Actual Results
	 # Pull TYP ANOVA Results
	TYP_AOV <- P.pr$TYP[[gene]]
	 # Pull TYP Additive Results
	TYP_DOS <- P.pr$TYP_DOS[[gene]]
	TYP_DOS.B <- B.pr$TYP_DOS[[gene]]
	 # Pull AA ANOVA Results
	AA_AOV <- P.pr$AA_AOV[[gene]]
	 # Pull AA Additive Results
	AA_DOS <- P.pr$AA_DOS[[gene]]
	 # Pull pat_AA Results
	pat_aa <-  PAT_AA[[gene]] # Raw Table

	## Pull HLA Type Data from Table
	HLA_AA.cols <- grep(paste("^",tag.temp,"_Pos",sep=""),colnames(HLA_AA.l))
	HLA_AA.temp <- HLA_AA.l[,HLA_AA.cols]
	colnames(HLA_AA.temp) <- gsub( paste(tag.temp,"_",sep=""),"",colnames(HLA_AA.temp) )
	HLA_TYP.cols <- grep(paste("^",tag.temp,sep=""),colnames(HLA_TYP.l))
	HLA_TYP.temp <- HLA_TYP.l[,HLA_TYP.cols]
	colnames(HLA_TYP.temp) <- paste("T",gsub("_","",gsub( tag.temp,"",colnames(HLA_TYP.temp) )),sep="")

	## Make Table for TYP Plot(s)
	TYP <- cbind( TYP_AOV, TYP_DOS[,order(colnames(TYP_DOS))] ) ; colnames(TYP)[1] <- "ANOVA"
	TYP <- TYP[PHENOS,]
	TYP.2 <- TYP[PHENOS,-which(colnames(TYP)%in%c("ANOVA","T"))]
	TYP.B <- TYP_DOS.B[,order(colnames(TYP_DOS.B))]
	TYP.B <- TYP.B[PHENOS,]
	if ( length(intersect(c("ANOVA","T"),colnames(TYP.B)))>0 ) {
		TYP.B.2 <- TYP.B[PHENOS,-which(colnames(TYP.B)%in%c("ANOVA","T"))]
	}else{ TYP.B.2 <- TYP.B }
	
	## Amino Acid ANOVA Results Plot(s)
	ISNA <- which( apply(AA_AOV,2,function(x) all(is.na(x)) ) | colnames(AA_AOV)=="TYP" )
	AA_AOV <- AA_AOV[ PHENOS,-ISNA ]

	## PLOT RESULTS ##################

	## Global Plotting Parameters
	NYH.ph <- NYHOLT( FT[,PHENOS] )
	PCH.ph <- rep(1,N.PHENOS)
	PCH.ph[grep("BL_MN",PHENOS)] <- 2
	PCH.ph[grep("RF|ACPA",PHENOS)] <- 3
	names(PCH.ph) <- PHENOS

	## HLA-Type Results Plots
	 # Specific Plotting Parameters
	YLIM <- c(0, max(4,-log10(min(TYP/30))) )
	YLIM.B <- extendrange( TYP.B.2 )
	MAIN.1 <- paste("HLA Type ANOVA:",gene)
	MAIN.2 <- paste("HLA Type Dosage Regression:",gene)
	NYH <- NYH.ph*NYHOLT(HLA_TYP.temp[,colnames(TYP_DOS)])
	# BH <- exp*(.05/NYH.ph)
	# points( -log10(exp), -log10(BH), type="l",lty=3,col="firebrick2",lwd=3)

	png( paste(PathToPlot,tag,"_3ABC_TYP_AAanova",gene,".png",sep=""),height=1200,width=1800,pointsize=28 )
	# layout( matrix(c(1,1,2,2, 3,3,3,4, 5,6,6,7), nrow=3, byrow=TRUE), widths=c(.2,1,3,1),height=c(1,1,1) )
	layout( matrix(c(1,2,2,3,3,4), ncol=3, byrow=TRUE), widths=c(1,3,2),height=c(1,1) )
	 # ANOVA
	barplot( -log10(TYP[,"ANOVA"]), beside=T, las=2,col=COLS.ph,border=NA,ylim=YLIM,main=MAIN.1,ylab="-log10(p)")
	abline( h=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH.ph)),lty=2,col=COLS.cor,lwd=4 )
	barplot( -log10(TYP[,"ANOVA"]), beside=T, las=2,col=COLS.ph,border=NA,add=T )
	 # Regression (P-Values)
	barplot( -log10(TYP.2), beside=T, las=2,col=COLS.ph,border=NA,ylim=YLIM,main=MAIN.2,ylab="-log10(p)")
	abline( h=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH)),lty=2,col=COLS.cor,lwd=4 )
	legend( "topleft", legend=rownames(TYP),title="Phenotype",fill=COLS.ph,border=NA, cex=1.2,ncol=nrow(TYP.2) ) # ncol=ceiling(nrow(TYP)/4),
	barplot( -log10(TYP.2), beside=T, las=2,col=COLS.ph,border=NA,add=T )
	#  # Regression (B-Values)
	# barplot( TYP.B.2, beside=T, las=2,col=COLS.ph,ylim=YLIM.B,main=MAIN.2,ylab="Effect Size")
	# abline( h=0:20,lty=3,col="grey50",lwd=1 )
	# abline( h=-log10(.05/(NYH)),lty=3,col="magenta2",lwd=3 )
	# legend( "topright", legend=rownames(TYP),title="Phenotype",fill=COLS.ph, ncol=ceiling(nrow(TYP)/4),cex=.9 )
	# barplot( TYP.B.2, beside=T, las=2,col=COLS.ph,add=T )
	# dev.off()

	## Amino Acid ANOVA
	 # Specific Plotting Parameters
	YLIM <- c(0, max(4,-log10(min(AA_AOV/30))) )
	XVALS <- gsub("Pos_","",colnames(AA_AOV))
	XVALS <- as.numeric( gsub(".","-",XVALS,fixed=T) )
	XLIM <- range(XVALS)
	MAIN <- paste("Amino Acid ANOVA: HLA",gene)
	NYH.bonf <- prod(dim(AA_AOV))
	 # Manhattan Style Plot Across Gene
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main=MAIN,xlab="Amino Acid Position",ylab="-log10(p)",xaxt="n")
	axis( 1, at=seq(-1000,1000,20),las=1 )	
	abline( v=seq(-1000,1000,10),lty=3,col="grey50")
	abline( h=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH.bonf)),lty=2,col=COLS.cor,lwd=4 )
	# legend( "topright", legend=rownames(AA_AOV),title="Phenotype",col=COLS.ph,pch=PCH.ph, ncol=3,cex=.8,pt.lwd=2 )
	SCRAP <- lapply( PHENOS, function(p) points( XVALS,-log10(AA_AOV[p,]), col=COLS.ph[p],pch=PCH.ph[p],lwd=2 ) )
	 # QQ Plot
	plot( 0,0,type="n",xlim=YLIM,ylim=YLIM, main="QQ P-Vals",xlab="-log10(Exp)",ylab="-log10(Exp)")
	abline( 0,1, lwd=2, col="black" )
	abline( h=0:20,v=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH.bonf)),lty=2,col=COLS.cor,lwd=4 )
	legend( "bottomright", legend=rownames(AA_AOV),title="Phenotype",col=COLS.ph,pch=PCH.ph, cex=1,pt.lwd=2 )
	exp <- 1:ncol(AA_AOV) / ncol(AA_AOV)
	SCRAP <- lapply( PHENOS, function(p) points( -log10(exp),-log10(sort(AA_AOV[p,])), col=COLS.ph[p],pch=PCH.ph[p],lwd=2 ) )
	dev.off()

	## Amino Acid Additive Model
	for ( p in 1:N.PHENOS ) {
		pheno <- PHENOS[p]
		ISNA <- which( unlist(lapply(AA_DOS[[pheno]],function(x) all(x=="NA") )) )
		AAP <- AA_DOS[[pheno]][-ISNA]
		AAP <- AAP[which(names(AAP)!="TYP")]
		AAP <- lapply( AAP, function(x) x[which(x!=0)] )
		HLA_AA.names <- intersect( colnames(HLA_AA.temp), gsub( ".","_",names(unlist(AAP)), fixed=T ) )

		## Plotting Parameters & Data
		 # Data
		exp <- 1:length(unlist(AAP)) / length(unlist(AAP))
		obs <- unlist(AAP)
		 # Parameters
		XVALS <- gsub("Pos_","",names(AAP))
		XVALS <- as.numeric( gsub(".","-",XVALS,fixed=T) )
		XLIM <- range(XVALS)
		YLIM <- c(0, max(4,-log10(as.numeric(Reduce(min,AAP))/30)) )
		NYH <- NYHOLT(HLA_AA.temp[,HLA_AA.names])
		MAIN <- paste("Amino Acid Additive Model: HLA",gene,"-",pheno)
		HEIGHT <- 800
		WIDTH <- 2400 # 2000+5*ncol(AA_AOV) + 200
		PLOT_RATIO <- c( (WIDTH-HEIGHT)/WIDTH, HEIGHT/WIDTH, 200/WIDTH )
		png( paste(PathToPlot,tag,"_3D_AAdose",gene,"_",pheno,".png",sep=""),height=HEIGHT,width=WIDTH,pointsize=34 )
		layout( matrix(c(1,2,3), 1, 3, byrow = TRUE), widths=PLOT_RATIO ) # layout( matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(3,2), heights=c(1,1) )
		 # Manhattan Style Plot
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main=MAIN,xlab="Amino Acid Position",ylab="-log10(p)",xaxt="n")
		axis( 1, at=seq(-1000,1000,20),las=1 )	
		abline( v=seq(-1000,1000,10),lty=3,col="grey50")
		abline( h=0:20,lty=3,col="grey50",lwd=1 )
		abline( h=-log10(.05/(NYH)),lty=2,col=COLS.cor,lwd=4 )
		 # Data
		SCRAP <- lapply( 1:length(AAP), function(x) points(rep(XVALS[x],length(AAP[[x]])),-log10(AAP[[x]]), pch=names(AAP[[x]]),col=COLS.AA[names(AAP[[x]])] ) )
		## QQ Plot
		plot( 0,0,type="n",xlim=YLIM,ylim=YLIM, main="QQ P-Vals",xlab="-log10(Exp)",ylab="-log10(Exp)")
		abline( 0,1, lwd=2, col="black" )
		abline( v=seq(-1000,1000,10),lty=3,col="grey50")
		abline( h=0:20,v=0:20,lty=3,col="grey50",lwd=1 )
		 # Multiple Hypotheses
		abline( h=-log10(.05/(NYH)),lty=2,col=COLS.cor,lwd=4 )
		 # Data
		points( -log10(exp),-log10(sort(obs)), col=COLS.AA[unlist(lapply(AAP,names))][order(obs)],pch=unlist(lapply(AAP,names))[order(obs)],lwd=2 )
		## Amino Acid Key
		Which_AA <- which(names(COLS.AA)%in% pat_aa )
		COLS.aa <- COLS.AA[Which_AA]
		barplot( matrix(rep(1,length(COLS.aa)),ncol=1), beside=F,col=COLS.aa,xaxt="n",yaxt="n",ylab="AA")
		axis( 2, at=1:length(COLS.aa)-.5, label=names(COLS.aa),las=2,cex.axis=.8 )
		dev.off()
	}
}
# PLOT_RESULTS( P.out$Dig_4, "DRB1", "Pr4" )

## Run Through Genes
# WHICH_GENES <- c("A","B","C","DQB1","DRB1","DPA1","DPB1","DQA1")
WHICH_GENES <- c("DQB1","DRB1")
for ( gene in WHICH_GENES ) {
	PLOT_RESULTS( P.out$Dig_4, B.out$Dig_4, gene, PH.COV, "Pr4" )
}

## Run Through Genes
WHICH_GENES <- c("A","B","C","DQB1","DRB1")
for ( gene in WHICH_GENES ) {
	PLOT_RESULTS( P.out$Dig_2, B.out$Dig_2, gene, PH.COV, "Pr2" )
}

#############################################################
## HI-LITE SPECIFIC RESULTS #################################
#############################################################

HEAT_TAB <- function( MG, column ) {
	TEMP <- table( MG$ACPA, MG[,column] )
	COLS <- colorRampPalette(c("chartreuse2","white","firebrick2"))(100)
	png( paste(PathToPlot,"HiLite-",column,".1.png",sep=""),height=500,width=600,pointsize=16 )
	heatmap.2( prop.table(TEMP,1), Colv=NA,Rowv=NA,dendrogram="none",trace="none",
		col=COLS, scale="none",cellnote=TEMP,notecol="black",notecex=2.4,
		lhei=c(1,4),lwid=c(1,6),key=F,cexRow=1.4,cexCol=1.4,
		main=paste("ACPA Status vs",column),xlab=paste(column,"Alleles"),ylab="ACPA (1=Negative)" )
	dev.off()
	png( paste(PathToPlot,"HiLite-",column,".2.png",sep=""),height=500,width=600,pointsize=16 )
	heatmap.2( prop.table(TEMP,2), Colv=NA,Rowv=NA,dendrogram="none",trace="none",
		col=COLS, scale="none",cellnote=TEMP,notecol="black",notecex=2.4,
		lhei=c(1,4),lwid=c(1,6),key=F,cexRow=1.4,cexCol=1.4,
		main=paste("ACPA Status vs",column),xlab=paste(column,"Alleles"),ylab="ACPA (1=Negative)" )
	dev.off()
	return(TEMP)
}

## ACPA vs HLA-DRB1*03:01
 # NOT Shared Epitope
 # Protective against RA
table( MG.typ$ACPA, MG.typ$Pr4_DRB1_03_01 )
round(prop.table( table( MG.typ$ACPA, MG.typ$Pr4_DRB1_03_01 ), 2 ),3)
fisher.test( table( MG.typ$ACPA, MG.typ$Pr4_DRB1_03_01 ) )
HEAT_TAB( MG.typ, "Pr4_DRB1_03_01" )

## ACPA vs HLA-DRB1*01:01
 # IS Shared Epitope
 # RA Susceptibility Haplotype
table( MG.typ$ACPA, MG.typ$Pr4_DRB1_01_01 )
round(prop.table( table( MG.typ$ACPA, MG.typ$Pr4_DRB1_01_01 ), 2),3)
fisher.test( table( MG.typ$ACPA, MG.typ$Pr4_DRB1_01_01 ) )
HEAT_TAB( MG.typ, "Pr4_DRB1_01_01" )

## ACPA vs HLA-DRB1*01:01
 # IS Shared Epitope
 # RA Susceptibility Haplotype
HEAT_TAB( MG.aa, "Pr4_DRB1_Pos_77_N" )
HEAT_TAB( MG.aa, "Pr4_DRB1_Pos_74_R" )
HEAT_TAB( MG.aa, "Pr4_DRB1_Pos_77_T" )

#############################################################
## END OF DOC ###############################################
#############################################################









# #############################################################
# ## HAPLOTYPE ANALYSIS (DRB1) ################################
# #############################################################

# ## Function to do Haplotype Analysis of Specified Amino Acid Positions
# HAP_AN <- function( Positions, tag ) {
# 	## Pull together Haplotypes from Specified Positions
# 	HAP <- apply( PAT_AA$DRB1[,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
# 	HAP <- gsub("NA","-",HAP)
# 	HAP.uniq <- sort(unique(HAP))
# 	N.HAP <- length(HAP.uniq)
# 	 # Get Haplotype Frequencies
# 	HAP.freq <- table(HAP)
# 	HAP.rare <- names(HAP.freq)[which(HAP.freq < 15 )]
# 	 # Convert to Array for Additive Analysis
# 	HAP.arr <- array( 0,c(N.PATS,length(HAP.uniq)) )
# 	colnames(HAP.arr) <- HAP.uniq
# 	rownames(HAP.arr) <- PATS
# 	for ( pat in PATS ) {
# 		HAP.pat <- HAP[ grep(pat,names(HAP)) ]
# 		if ( HAP.pat[1]==HAP.pat[2] ) { HAP.arr[pat,HAP.pat[1]] <- 2
# 		}else{ HAP.arr[pat,HAP.pat] <- 1 }
# 	}
# 	 # Plot Haplotype Frequency
# 	png( paste(PathToPlot,"/DRB1_",tag,"_1-HapFreq.png",sep=""), height=800,width=1600,pointsize=30 )
# 	barplot( HAP.freq,las=2,col=COLS.fr,main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype Frequency"),ylab="# Haplotypes")
# 	abline(h=seq(0,1000,20),lty=3,col="grey50")
# 	barplot( HAP.freq,las=2,col=COLS.fr,add=T)
# 	dev.off()

# 	## Combine Haplotype & Phenotype Data
# 	MG.HAP <- merge( HAP.arr, FT, by.x="row.names",by.y="ID" )
# 	# MG.HAP[,"ACPA"] <- MG.HAP[,"ACPA"]=="Positive"
# 	MG.HAP[,"ACPA"] <- MG.HAP[,"ACPA"]-1
# 	MG.HAP2 <- merge( data.frame(Samp=sapply(strsplit(names(HAP),"_"),"[",1),HAP), FT, by.x="Samp",by.y="ID" )
# 	# MG.HAP2[,"ACPA"] <- MG.HAP2[,"ACPA"]=="Positive"
# 	MG.HAP2[,"ACPA"] <- MG.HAP2[,"ACPA"]-1

# 	## Association w/ Additive by Haplotype
# 	LM.HAP.AA.F <- LM.HAP.AA <- LM.HAP.ANV <- P.HAP.AA <- list()
# 	P.HAP.ANV <- numeric(length(PHENOS)) ; names(P.HAP.ANV) <- PHENOS
# 	for ( p in 1:length(PHENOS) ) {
# 		pheno <- PHENOS[p]
# 		cov <- COVS[p]
# 		# ANOVA Model
# 		if ( cov=="" ) {
# 			formula <- as.formula(paste( pheno,"~ HAP" ))
# 		}else{
# 			formula <- as.formula(paste( pheno,"~",cov,"+HAP" ))
# 		}
# 		if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
# 			LM.HAP.ANV[[pheno]] <- MOD <- chisq.test( table( MG.HAP2[,pheno],as.character(MG.HAP2[,"HAP"]) ) )
# 			P.HAP.ANV[pheno] <- LM.HAP.ANV[[pheno]]$p.value
# 		}else{
# 			LM.HAP.ANV[[pheno]] <- MOD <- lm( formula, data=MG.HAP2, subset=which(!(HAP%in%HAP.rare)) )
# 			P.HAP.ANV[pheno] <- anova(MOD)["HAP","Pr(>F)"]
# 			# Plot Residuals
# 			if ( cov=="" ) {
# 				formula <- as.formula(paste( pheno,"~ HAP" ))
# 				RESID <- MG.HAP2[,pheno]
# 				names(RESID) <- 1:length(RESID)
# 			}else{
# 				formula <- as.formula(paste( pheno,"~",cov ))
# 				RESID <- resid(lm( formula, data=MG.HAP2))
# 			}
# 			png( paste(PathToPlot,"/DRB1_",tag,"_2-BoxPlot_",pheno,".png",sep=""), height=1400,width=1600,pointsize=30 )
# 			par(mfrow=c(2,1))
# 			barplot( table(HAP),las=2,col="grey25",main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype Frequency"),ylab="# Haplotypes")
# 			abline(h=seq(0,1000,20),lty=3,col="grey50")
# 			barplot( table(HAP),las=2,col="grey25",add=T)
# 			boxplot( RESID ~ MG.HAP2$HAP[as.numeric(names(RESID))], las=2,col=COLS.ph[p],main=paste("HLA-DRB1: Pos11,71,74 Haplotype vs",pheno),ylab=paste("Resid: vs",cov) )
# 			abline(h=seq(-10,10,1),lty=3,col="grey50")
# 			boxplot( RESID ~ MG.HAP2$HAP[as.numeric(names(RESID))],add=T,las=2,col=COLS.ph[p] )
# 			points( RESID ~ MG.HAP2$HAP[as.numeric(names(RESID))],pch="+" )
# 			dev.off()
# 		}
# 		# Additive Model
# 		  # Individual Haplotype Dosage
# 		LM.HAP.AA[[pheno]] <- list()
# 		P.HAP.AA[[pheno]] <- numeric(N.HAP) ; names(P.HAP.AA[[pheno]]) <- HAP.uniq
# 		for ( hap in HAP.uniq ) {
# 			if ( all( strsplit(hap,"")[[1]]=="-" ) ) { next }
# 			formula <- as.formula(paste( pheno,"~",cov,"+",hap ))
# 			if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
# 				formula <- as.formula(paste( pheno,"~",hap ))
# 				LM.HAP.AA[[pheno]][[hap]] <- MOD <- glm( formula, data=MG.HAP, family=binomial(logit) )
# 				P.HAP.AA[[pheno]][hap] <- summary(MOD)$coefficients[hap,"Pr(>|z|)"]
# 			}else{
# 				LM.HAP.AA[[pheno]][[hap]] <- MOD <- lm( formula, data=MG.HAP )
# 				P.HAP.AA[[pheno]][hap] <- anova(MOD)[hap,"Pr(>F)"]
# 			}
# 		} # Close Haplotype Loop
# 		  # Full Haplotype Model
# 		HAP.uniq.2 <- grep("-",HAP.uniq,invert=T,value=T)
# 		formula <- as.formula(paste( pheno,"~",cov,"+",paste(HAP.uniq.2,collapse="+") ))
# 		if ( !(pheno %in% c("RF_ACPA","ACPA","RF")) ) {
# 			LM.HAP.AA.F[[pheno]] <- lm( formula, data=MG.HAP )
# 		}
# 	} # Close Pheno Loop

# 	## Plot Results
# 	P.comp <- cbind( Reduce( rbind, P.HAP.AA ), P.HAP.ANV )
# 	colnames(P.comp)[ncol(P.comp)] <- "ANOVA"
# 	rownames(P.comp) <- names(P.HAP.AA)
# 	P.comp <- P.comp[,-which(colnames(P.comp)==paste(rep("-",length(Positions)),collapse=""))]

# 	YLIM <- c( 0,-log10(min(P.comp)/30) )
# 	png( paste(PathToPlot,"/DRB1_",tag,"_3-HapAssoc.png",sep=""), height=800,width=1600,pointsize=30 )
# 	barplot( -log10(P.comp), beside=T, col=COLS.ph,ylim=YLIM,las=2,main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype ANOVA"),ylab="-log10(p)")
# 	abline(h=seq(0,20,1),lty=3,col="grey50" )
# 	abline( h=-log10(.05/(10^(0:6))),lty=2,col=COLS.cor,lwd=2 )
# 	legend( "topright",legend=rownames(P.comp),fill=COLS.ph,ncol=3,cex=.8)
# 	barplot( -log10(P.comp), beside=T, col=COLS.ph,ylim=YLIM,las=2,add=T )
# 	dev.off()
# 	## Compile Outputs
# 	# OUT <- list( MOD.AA=LM.HAP.AA, MOD.ANV=LM.HAP.ANV, P.AA=P.HAP.AA, P.ANV=P.HAP.ANV )
# 	OUT <- list( MOD.AA=LM.HAP.AA, MOD.AA.F=LM.HAP.AA.F, MOD.ANV=LM.HAP.ANV, P=P.comp, HAP=HAP.arr, SAMPS=MG.HAP$Row.names )
# }

# OUT <- list()

# ##########################################
# ## POS 11, 71, 74 ##
#  # Viatte, et al (2015)
# Positions <- c(11,71,74)
# tag <- paste(c("p",Positions),collapse="")
# OUT[[tag]] <- HAP_AN(Positions,tag)

# ##########################################
# ## POS 11, 13, 71, 74 ##
# Positions <- c(11,13,71,74)
# tag <- paste(c("p",Positions),collapse="")
# OUT[[tag]] <- HAP_AN(Positions,tag)

# ##########################################
# ## POS 11, 13, 71, 74 ##
# Positions <- 70:74
# tag <- "pSE"
# OUT[[tag]] <- HAP_AN(Positions,tag)

# MG.hap.se <- merge( FT, OUT$pSE$HAP, by.x="ID",by.y="row.names" )
# HEAT_TAB( MG.hap.se, "QKRAA" )
# TEMP.tab <- HEAT_TAB( MG.hap.se, "QKRGR" )
# fisher.test(TEMP.tab)

# ##########################################
# ## Compile Results & Plot ##

# ## Pull out Relevant Info from Models
# N.ph <- 2 # Number of Phenotypes used
# BETAS <- SES <- PS <- FREQS <- list()
# BETAS.f <- SES.f <- PS.f <- FREQS.f <- list()
# for ( hap in names(OUT) ) {
# 	# Individual Models
# 	PS[[hap]] <- matrix( unlist(lapply( PHENOS[1:N.ph], function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Pr(>|t|)"] )) )), byrow=F,ncol=N.ph )
# 	BETAS[[hap]] <- matrix( unlist(lapply( PHENOS[1:N.ph], function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Estimate"] )) )), byrow=F,ncol=N.ph )
# 	SES[[hap]] <- matrix( unlist(lapply( PHENOS[1:N.ph], function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Std. Error"] )) )), byrow=F,ncol=N.ph )
# 	colnames(BETAS[[hap]]) <- colnames(SES[[hap]]) <- colnames(PS[[hap]]) <- PHENOS[1:N.ph]
# 	rownames(BETAS[[hap]]) <- rownames(SES[[hap]]) <- rownames(PS[[hap]]) <- names(OUT[[hap]]$MOD.AA[[1]])
# 	# Full Models
# 	PS.f[[hap]] <- matrix( unlist(lapply( PHENOS[1:N.ph], function(x)summary(OUT[[hap]]$MOD.AA.F[[x]])$coefficients[rownames(PS[[hap]]),"Pr(>|t|)"] )), byrow=F,ncol=N.ph )
# 	BETAS.f[[hap]] <- matrix( unlist(lapply( PHENOS[1:N.ph], function(x)summary(OUT[[hap]]$MOD.AA.F[[x]])$coefficients[rownames(PS[[hap]]),"Estimate"] )), byrow=F,ncol=N.ph )
# 	SES.f[[hap]] <- matrix( unlist(lapply( PHENOS[1:N.ph], function(x)summary(OUT[[hap]]$MOD.AA.F[[x]])$coefficients[rownames(PS[[hap]]),"Std. Error"] )), byrow=F,ncol=N.ph )
# 	colnames(BETAS.f[[hap]]) <- colnames(SES.f[[hap]]) <- colnames(PS.f[[hap]]) <- colnames(PS.f[[hap]])
# 	rownames(BETAS.f[[hap]]) <- rownames(SES.f[[hap]]) <- rownames(PS.f[[hap]]) <- rownames(PS.f[[hap]])	
# 	# Haplotype Frequencies
# 	FREQS[[hap]] <- colSums( OUT[[hap]]$HAP[ OUT[[hap]]$SAMPS, ] )
# }

# ## Compare Beta to Allele Frequency & Variance
# HAP_OUT <- list()
# for ( hap in names(OUT) ) {
# 	P.n <- colnames(BETAS[[hap]])
# 	P.num <- length(P.n)
# 	H.n <- rownames(BETAS[[hap]])
# 	H.f <- FREQS[[hap]][H.n]
# 	H.b <- BETAS[[hap]][H.n,]
# 	H.s <- SES[[hap]][H.n,]
# 	H.p <- PS[[hap]][H.n,]

# 	## BETA for Response vs Disease Severity
# 	for ( p in 1 ) {
# 		y_pheno <- PHENOS[p]
# 		x_pheno <- PHENOS[p+1]
# 		png( paste(PathToPlot,"DRB1_",hap,"_4-BETA.1.",x_pheno,".png",sep=""),height=1600,width=1600,pointsize=36 )
# 		XLIM <- extendrange(H.b[,x_pheno],f=.2)
# 		YLIM <- extendrange(H.b[,y_pheno],f=.2)
# 		plot( H.b[,x_pheno],H.b[,y_pheno], pch=21,col=COLS.ph[x_pheno],bg=COLS.ph[y_pheno], xlim=XLIM,ylim=YLIM,main=paste("Beta Estimates of",y_pheno,"vs",x_pheno,"-",hap),xlab=paste("Beta:",x_pheno),ylab=paste("Beta:",y_pheno) )
# 		abline( h=seq(-5,5,.5),v=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
# 		abline( h=0,v=0,lty=1,col="grey50",lwd=1 )
# 		arrows( H.b[,x_pheno]+H.s[,x_pheno],H.b[,y_pheno],H.b[,x_pheno]-H.s[,x_pheno],H.b[,y_pheno], code=3,angle=90,lwd=5,col=COLS.ph[x_pheno] )
# 		arrows( H.b[,x_pheno],H.b[,y_pheno]+H.s[,y_pheno],H.b[,x_pheno],H.b[,y_pheno]-H.s[,y_pheno], code=3,angle=90,lwd=5,col=COLS.ph[y_pheno] )
# 		text( H.b[,x_pheno],H.b[,y_pheno]+.02*diff(YLIM), label=rownames(H.b), col=COLS.ph[x_pheno], pos=4,cex=1.2 )
# 		points( H.b[,x_pheno],H.b[,y_pheno], pch=21,col=COLS.ph[x_pheno],bg=COLS.ph[y_pheno], lwd=5,cex=1.2 )
# 		MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno] )
# 		abline(MOD,lwd=6,lty=2,col=COLS.ph[y_pheno] )
# 		text( quantile(XLIM,.1),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col=COLS.ph[y_pheno],cex=1.2 )
# 		dev.off()
# 	}

# 	## BETA vs SE vs FREQ
# 	png( paste(PathToPlot,"DRB1_",hap,"_5-BETA.2.png",sep=""),height=1600,width=1600,pointsize=36 )
# 	# png( paste(PathToPlot,"DRB1-2_PAIRS.",hap,".png",sep=""),height=800,width=2400,pointsize=36 )
# 	# par(mfrow=c(1,3))
# 	 # BETA vs FREQ
# 	# plot( c(H.b) ~ rep(H.f,P.num), ylab="Beta",xlab="Frequency",main=paste("Effect Size vs Haplotype Frequency:",hap),col=COLS.ph[rep(P.n,each=P.num)],pch=as.numeric(as.factor(rep(H.n,P.num) )) )
# 	 # SE vs FREQ
# 	XLIM <- c(0,max(H.f)) ; YLIM <- c(0,max(H.s))
# 	plot( c(H.s) ~ rep(H.f,P.num), xlim=XLIM,ylim=YLIM,ylab="SE",xlab="Frequency",main=paste("Standard Error vs Haplotype Frequency:",hap),col=COLS.ph[rep(P.n,each=P.num)],pch=as.numeric(as.factor(rep(H.n,P.num) )), cex=3*sqrt(abs(c(H.b))),lwd=4 )
# 	abline(h=seq(-5,5,.1),v=seq(0,1000,50),lty=3,col="grey50",lwd=1)
# 	legend( quantile(XLIM,.65),quantile(YLIM,1), fill=COLS.ph,legend=names(COLS.ph) )
# 	legend( quantile(XLIM,.465),quantile(YLIM,1), pch=as.numeric(as.factor(H.n)),legend=H.n )
# 	legend( quantile(XLIM,.3),quantile(YLIM,1), pch=16,pt.cex=3*sqrt(c(.1,.5,1)),legend=c(.1,.5,1),title="Effect Size" )
# 	 # SE vs BETA
# 	# plot( c(H.b) ~ c(H.s), ylab="Beta",xlab="SE",main=paste("Effect Size vs Standard Error:",hap),col=COLS.ph[rep(P.n,each=P.num)],pch=as.numeric(as.factor(rep(H.n,P.num) )) )
# 	dev.off()

# 	## Weighted Regression of BETAs (Response vs Severity)
# 	 # ...and also plot of Effect Size + Frequency (like Viatte, 2015)
# 	for ( p in 1 ) {
# 		y_pheno <- P.n[p]
# 		x_pheno <- P.n[p+1]
# 		png( paste(PathToPlot,"DRB1_",hap,"_6-BETA.3.",x_pheno,".png",sep=""),height=1200,width=2400,pointsize=36 )
# 		par(mfrow=c(1,2))
# 		XLIM <- extendrange(H.b[,x_pheno],f=.2)
# 		YLIM <- extendrange(H.b[,y_pheno],f=.2)
# 		WTS <- H.f
# 		# WTS <- 1 / rowMeans(H.s[,c(y_pheno,x_pheno)]) 
# 		## Weighted Regression
# 		plot( H.b[,x_pheno],H.b[,y_pheno], pch=21,col=COLS.ph[x_pheno],bg=COLS.ph[y_pheno],lwd=5,cex=WTS/mean(WTS), xlim=XLIM,ylim=YLIM,main=paste("Wtd. Betas of",y_pheno,"vs",x_pheno,"-",hap),xlab=paste("Beta:",x_pheno),ylab=paste("Beta:",y_pheno) )
# 		# print( H.b[,c(y_pheno,x_pheno)])
# 		abline( h=seq(-5,5,.5),v=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
# 		abline( h=0,v=0,lty=1,col="grey50",lwd=1 )
# 		arrows( H.b[,x_pheno]+H.s[,x_pheno],H.b[,y_pheno],H.b[,x_pheno]-H.s[,x_pheno],H.b[,y_pheno], code=3,angle=90,lwd=5,col=COLS.ph[x_pheno] )
# 		arrows( H.b[,x_pheno],H.b[,y_pheno]+H.s[,y_pheno],H.b[,x_pheno],H.b[,y_pheno]-H.s[,y_pheno], code=3,angle=90,lwd=5,col=COLS.ph[y_pheno] )
# 		text( H.b[,x_pheno],H.b[,y_pheno]+.02*diff(YLIM), label=rownames(H.b), col=COLS.ph[y_pheno], pos=4,cex=1.2 )
# 		points( H.b[,x_pheno],H.b[,y_pheno], pch=21,col=COLS.ph[x_pheno],bg=COLS.ph[y_pheno],lwd=5,cex=WTS/mean(WTS) )
# 		MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS )
# 		abline(MOD,lwd=6,lty=2,col="grey25" )
# 		text( quantile(XLIM,.1),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col="grey25",cex=1.2 )
# 		## Haplotype Effect Sizes
# 		XLIM <- c(0,1+nrow(H.b))
# 		YLIM <- extendrange(H.b[,c(y_pheno,x_pheno)],f=.2)
# 		ORDER <- order(H.b[,x_pheno])
# 		plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",main=paste("Effect Sizes for:",hap,"-",y_pheno,"&",x_pheno),ylab="Beta",xlab="Haplotype" )
# 		abline( v=1:nrow(H.b),lty=3,col="grey50",lwd=1 )
# 		points( 1:nrow(H.b),H.b[ORDER,x_pheno], pch=10,cex=log(2*H.f[ORDER]),col=COLS.ph[x_pheno],lwd=5 )
# 		# print( H.b[ORDER,c(y_pheno,x_pheno)])
# 		axis(1,at=1:nrow(H.b),label=rownames(H.b)[ORDER],las=2)
# 		points( 1:nrow(H.b),H.b[ORDER,y_pheno], pch=10,cex=log(2*H.f[ORDER]),col=COLS.ph[y_pheno], lwd=5 )
# 		legend("topleft",pch=10,col=COLS.ph[c(y_pheno,x_pheno)],legend=c(y_pheno,x_pheno),pt.cex=2,pt.lwd=5 )
# 		dev.off()
# 	}
# 	## Compile Outputs
# 	HAP_OUT[[hap]] <- data.frame( HAP=H.n, FREQ=H.f, BETA=H.b, SE=H.s, P=H.p )
# 	write.table( HAP_OUT[[hap]], paste(PathToPlot,"_DRB1_",hap,"_Table.txt",sep=""), col.names=T,row.names=F,quote=F,sep="\t" )
# }


# #############################################################
# ## RE-CREATE VIATTE 2015 PLOTS ##############################
# #############################################################

# ##########################################
# ## FCT: Pull out Beta/P/SE Values for Position & Phenotype
# PULL_BETA <- function( position, pheno ) {
# 	which_pos <- paste( "Pos",position,sep="_" )
# 	which_pos <- gsub("-",".",which_pos,fixed=T)
# 	TEMP.dat <- M.out$Dig_4$AA_DOS$DRB1[[pheno]][[which_pos]]
# 	TEMP <- matrix(unlist( lapply( TEMP.dat, function(x)summary(x)$coefficients[length(coef(x)),] ) ),byrow=T,ncol=4) # ,byrow=F,ncol=length(TEMP.dat) )
# 	rownames(TEMP) <- paste(which_pos,names(TEMP.dat),sep="_")
# 	colnames(TEMP) <- c("BETA","SE","T","P")# colnames(summary(TEMP.dat[[1]])$coefficients) # c( names(coef(TEMP.dat[[1]]))[-length(coef(TEMP.dat[[1]]))], which_pos )
# 	TEMP.maf <- unlist(lapply( TEMP.dat, function(x) mean(x$model[,ncol(x$model)]) )) / 2
# 	OUT <- data.frame( tag=rownames(TEMP), t(sapply(strsplit(rownames(TEMP),"_"),"[",2:3)), MAF=TEMP.maf, TEMP )
# 	colnames(OUT)[1:3] <- c("tag","Pos","AA")
# 	return(OUT)
# }

# Positions <- c(11,71,74)
# TAB_2.list <- list()
# for ( p in 1:8 ) {
# 	pheno <- PHENOS[p]
# 	TEMP.list <- lapply( Positions, function(x) PULL_BETA(x,pheno) )
# 	TAB_2.list[[pheno]] <- Reduce( rbind, TEMP.list )
# 	rownames(TAB_2.list[[pheno]]) <- NULL
# 	# TEMP.arr <- Reduce( rbind, TEMP.list )
# 	# TAB_2.list[[pheno]] <- data.frame( , TEMP.arr )
# 	# colnames(TAB_2.list[[pheno]])[1:2] <- c("Pos","AA")
# }

# ##########################################
# ## FCT: Grouped Haplotype Analysis
#  # Figures (Later)
# HAPLO_GROUP <- function( hap, HAP_GRP, tag ) {
# 	# Set New HAP_GRP variable to be modified in first half of function
# 	HAP_GRP.2 <- HAP_GRP
# 	 # Specify Group Colors
# 	COLS.hap_grp <- c("seagreen2","magenta2","deepskyblue2","chocolate1","brown1","cyan2","blueviolet","yellow1","deeppink2","darkolivegreen1","goldenrod2")
# 	COLS.hap_grp <- COLS.list.2
# 	 # Specify any non-assigned haplotypes
# 	COLS.hap_grp <- c( COLS.hap_grp[1:length(HAP_GRP.2)], "black" )
# 	HAP_GRP.2$Other <- setdiff( rownames(BETAS[[hap]]),Reduce(union,HAP_GRP.2) )
# 	names(COLS.hap_grp) <- names(HAP_GRP.2)
# 	KEY.hap_grp <- data.frame( HAP=unlist(HAP_GRP.2), GRP=rep(names(HAP_GRP.2),lapply(HAP_GRP.2,length)),COL=rep(COLS.hap_grp,lapply(HAP_GRP.2,length)) )
# 	rownames(KEY.hap_grp) <- unlist(HAP_GRP.2)

# 	##########################################
# 	## Figures 2 & 3

# 	## Effect Size for Disease Severity vs Haplotype
# 	 # Size: Haplotype Frequency
# 	 # Color: Group

# 	## Create Plots for Various Phenotypes
# 	for ( p in 1:4 ) {
# 		P.n <- colnames(BETAS[[hap]])
# 		P.num <- length(P.n)
# 		y_pheno <- P.n[p]
# 		x_pheno <- P.n[p+4]
# 		H.n <- rownames(BETAS[[hap]])
# 		H.f <- FREQS[[hap]][H.n]
# 		H.b <- BETAS[[hap]][H.n,]
# 		H.s <- SES[[hap]][H.n,]
# 		H.p <- PS[[hap]][H.n,]
# 		H.c <- as.character(KEY.hap_grp[H.n,"COL"])
# 		## Create Plot/File Location
# 		png( paste(PathToPlot,"HapGroup_",hap,"-",tag,"_1_",x_pheno,".png",sep=""),height=1000,width=2000,pointsize=24 )
# 		par(mfrow=c(1,2))
# 		 # Model Response vs Severity Haplotypes
# 		y_pheno <- PHENOS[p]
# 		x_pheno <- PHENOS[p+4]
# 		XLIM <- extendrange(H.b[,x_pheno],f=.25)
# 		YLIM <- extendrange(H.b[,y_pheno],f=.25)
# 		WTS <- H.f
# 		# WTS <- 1 / rowMeans(H.s[,c(y_pheno,x_pheno)]) 
# 		# WTS.p <- rowMeans(-log10(H.p[,c(y_pheno,x_pheno)]))
# 		plot( H.b[,x_pheno],H.b[,y_pheno], pch=20,col=H.c, lwd=5,cex=WTS/mean(WTS), xlim=XLIM,ylim=YLIM,main=paste("Beta Estimates of",y_pheno,"vs",x_pheno,"-",hap),xlab=paste("Beta:",x_pheno),ylab=paste("Beta:",y_pheno) )
# 		abline( h=seq(-5,5,.5),v=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
# 		abline( h=0,v=0,lty=1,col="grey50",lwd=1 )
# 		arrows( H.b[,x_pheno]+H.s[,x_pheno],H.b[,y_pheno],H.b[,x_pheno]-H.s[,x_pheno],H.b[,y_pheno], code=3,angle=90,lwd=5,col=H.c )
# 		arrows( H.b[,x_pheno],H.b[,y_pheno]+H.s[,y_pheno],H.b[,x_pheno],H.b[,y_pheno]-H.s[,y_pheno], code=3,angle=90,lwd=5,col=H.c )
# 		text( H.b[,x_pheno],H.b[,y_pheno]+.02*diff(YLIM), label=rownames(H.b), col=H.c, pos=4,cex=1.2 )
# 		points( H.b[,x_pheno],H.b[,y_pheno], pch=20,col=H.c, lwd=5,cex=WTS/mean(WTS) )

# 		 # Significance		
# 		MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno] )
# 		abline(MOD,lwd=6,lty=2,col="mediumpurple3" )
# 		text( quantile(XLIM,.1),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col="mediumpurple3",cex=1.2 )
# 		MOD.w <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS )
# 		abline(MOD.w,lwd=6,lty=2,col="steelblue3" )
# 		text( quantile(XLIM,.1),quantile(YLIM,.07), label=paste("(ws)p=",formatC(summary(MOD.w)$coefficients[length(coef(MOD.w)),4],digits=2,format="e"),sep=""), col="steelblue3",cex=1.2 )
# 		# MOD.wp <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS.p )
# 		# abline(MOD.wp,lwd=6,lty=2,col="cadetblue3" )
# 		# text( quantile(XLIM,.1),quantile(YLIM,.12), label=paste("(wp)p=",formatC(summary(MOD.wp)$coefficients[length(coef(MOD.wp)),4],digits=2,format="e"),sep=""), col="cadetblue3",cex=1.2 )

# 		# MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno] )
# 		# abline(MOD,lwd=6,lty=2,col="mediumpurple3" )
# 		# text( quantile(XLIM,.1),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col="mediumpurple3",cex=1.2 )
# 		# MOD.w <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS )
# 		# abline(MOD.w,lwd=6,lty=2,col="cadetblue3" )
# 		# text( quantile(XLIM,.1),quantile(YLIM,.07), label=paste("wp=",formatC(summary(MOD.w)$coefficients[length(coef(MOD.w)),4],digits=2,format="e"),sep=""), col="cadetblue3",cex=1.2 )
# 		 # Haplotype Effect Sizes
# 		XLIM <- c(0,1+nrow(H.b))
# 		YLIM <- extendrange(H.b[,c(y_pheno,x_pheno)])
# 		ORDER <- order(H.b[,x_pheno])
# 		plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",main=paste("Effect Sizes for:",hap),ylab="Beta",xlab="Haplotype" )
# 		abline( h=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
# 		abline( h=0,lty=1,col="grey50",lwd=1 )
# 		abline( v=1:nrow(H.b),lty=3,col="grey50",lwd=1 )
# 		points( 1:nrow(H.b),H.b[ORDER,x_pheno], pch=9,cex=log(2*H.f[ORDER]),col=H.c[ORDER],lwd=5 )
# 		axis(1,at=1:nrow(H.b),label=rownames(H.b)[ORDER],las=2)
# 		points( 1:nrow(H.b),H.b[ORDER,y_pheno], pch=10,cex=log(2*H.f[ORDER]),col=H.c[ORDER], lwd=5 )
# 		legend("topleft",pch=9:10,col="black",legend=c(y_pheno,x_pheno),pt.cex=1.8,pt.lwd=4 )
# 		dev.off()
# 	}

# 	##########################################
# 	## Group Haplotypes ##

# 	## Greate Dosage Array for Patients & Haplotype Groups
# 	HAP_GRP.arr <- matrix(unlist(lapply( HAP_GRP, function(x) rowSums(OUT[[hap]]$HAP[,which(colnames(OUT[[hap]]$HAP)%in%x)]) )),byrow=F,ncol=length(HAP_GRP))
# 	colnames(HAP_GRP.arr) <- names(HAP_GRP)
# 	rownames(HAP_GRP.arr) <- rownames(OUT[[hap]]$HAP)
# 	MG.hap_grp <- merge( HAP_GRP.arr, FT, by.x="row.names",by.y="ID" )
# 	MOD.hap_grp.full <- MOD.hap_grp.ind <- list()
# 	 # Model Phenotypes vs Haplotype Groups
# 	for ( p in grep("ACPA|RF",PHENOS,invert=T) ) {
# 		pheno <- as.character(PH.COV[p,"PHENOS"])
# 		cov <- as.character(PH.COV[p,"COVS"])
# 		formula.full <- as.formula(paste(pheno,"~",cov,"+",paste(names(HAP_GRP),collapse="+")))
# 		MOD.hap_grp.full[[pheno]] <- lm( formula.full, data=MG.hap_grp )
# 		MOD.hap_grp.ind[[pheno]] <- list()
# 		for ( g in 1:length(HAP_GRP) ) {
# 			g_tag <- paste("G",g,sep="")
# 			formula.temp <- as.formula(paste(pheno,"~",cov,"+",g_tag))
# 			MOD.hap_grp.ind[[pheno]][[g_tag]] <- lm( formula.temp, data=MG.hap_grp )
# 		}
# 	}

# 	 # Full Models
# 	# BETAS.hap_grp.full <- lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[,"Estimate"] )
# 	# SES.hap_grp.full <- lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[,"Std. Error"] )
# 	# PS.hap_grp.full <- lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[,"Pr(>|t|)"] )
# 	# OUT.hap_grp.full <- lapply( names(BETAS.hap_grp.full),function(x)data.frame(BETA=BETAS.hap_grp.full[[x]],SE=SES.hap_grp.full[[x]],P=PS.hap_grp.full[[x]]) )
# 	# names(OUT.hap_grp.full) <- names(BETAS.hap_grp.full)
# 	# BETAS.G.hap_grp.full <- matrix(unlist( lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[names(HAP_GRP),"Estimate"] )),byrow=F,nrow=length(HAP_GRP) )
# 	# SES.G.hap_grp.full <- matrix(unlist( lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[names(HAP_GRP),"Std. Error"] )),byrow=F,nrow=length(HAP_GRP) )
# 	# PS.G.hap_grp.full <- matrix(unlist( lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[names(HAP_GRP),"Pr(>|t|)"] )),byrow=F,nrow=length(HAP_GRP) )
# 	# colnames(BETAS.G.hap_grp.full) <- colnames(SES.G.hap_grp.full) <- colnames(PS.G.hap_grp.full) <- names(BETAS.hap_grp.full)
# 	# rownames(BETAS.G.hap_grp.full) <- rownames(SES.G.hap_grp.full) <- rownames(PS.G.hap_grp.full) <- names(HAP_GRP)
#  	# OUT.G.hap_grp.full <- lapply( names(BETAS.G.hap_grp.full),function(x)data.frame(BETA=BETAS.G.hap_grp.full[,x],SE=SES.G.hap_grp.full[,x],P=PS.G.hap_grp.full[,x]) )
# 	# names(OUT.G.hap_grp.full) <- names(BETAS.G.hap_grp.full)
# 	 # Individual Models
# 	BETAS.G.hap_grp.ind <- matrix(unlist( lapply( MOD.hap_grp.ind, function(x) unlist(lapply(names(x),function(y)summary(x[[y]])$coefficients[y,"Estimate"])) )),byrow=F,nrow=length(HAP_GRP) )
# 	SES.G.hap_grp.ind <- matrix(unlist( lapply( MOD.hap_grp.ind, function(x) unlist(lapply(names(x),function(y)summary(x[[y]])$coefficients[y,"Std. Error"])) )),byrow=F,nrow=length(HAP_GRP) )
# 	PS.G.hap_grp.ind <- matrix(unlist( lapply( MOD.hap_grp.ind, function(x) unlist(lapply(names(x),function(y)summary(x[[y]])$coefficients[y,"Pr(>|t|)"])) )),byrow=F,nrow=length(HAP_GRP) )
# 	colnames(BETAS.G.hap_grp.ind) <- colnames(SES.G.hap_grp.ind) <- colnames(PS.G.hap_grp.ind) <- names(MOD.hap_grp.ind)
# 	rownames(BETAS.G.hap_grp.ind) <- rownames(SES.G.hap_grp.ind) <- rownames(PS.G.hap_grp.ind) <- names(HAP_GRP)

# 	## Plot
# 	# H.b <- BETAS.G.hap_grp.full
# 	# H.s <- SES.G.hap_grp.full
# 	# H.p <- PS.G.hap_grp.full
# 	H.b <- BETAS.G.hap_grp.ind
# 	H.s <- SES.G.hap_grp.ind
# 	H.p <- PS.G.hap_grp.ind
# 	P.n <- colnames(H.b)
# 	WTS <- colMeans(HAP_GRP.arr)
# 	# par(mfrow=c(2,2))
# 	for ( p in 1:4 ) {
# 		y_pheno <- P.n[p]
# 		x_pheno <- P.n[p+4]
# 		# WTS <- 1 / rowMeans(H.s[,c(y_pheno,x_pheno)]) 
# 		# WTS.p <- rowMeans(-log10(H.p[,c(y_pheno,x_pheno)]))
# 		XLIM <- extendrange(H.b[,x_pheno],f=.5)
# 		YLIM <- extendrange(H.b[,y_pheno],f=.5)
# 		png( paste(PathToPlot,"HapGroup_",hap,"-",tag,"_2_",x_pheno,".png",sep=""),height=1200,width=1200,pointsize=24 )
# 		plot( H.b[,y_pheno] ~ H.b[,x_pheno], cex=2,pch=20,col=COLS.hap_grp, xlim=XLIM,ylim=YLIM, main=paste("Severity vs Response:",hap,"Haplo Groups"),xlab=x_pheno,ylab=y_pheno )
# 		abline(h=0,v=0,lty=1,col="grey50",lwd=1)
# 		arrows( H.b[,x_pheno], H.b[,y_pheno]+H.s[,y_pheno], H.b[,x_pheno], H.b[,y_pheno]-H.s[,y_pheno], col=COLS.hap_grp, code=3,angle=90,length=.1,lwd=4 )
# 		arrows( H.b[,x_pheno]+H.s[,x_pheno] ,H.b[,y_pheno], H.b[,x_pheno]-H.s[,x_pheno], H.b[,y_pheno], col=COLS.hap_grp, code=3,angle=90,length=.1,lwd=4 )
# 		 # Label Haplotypes in each Group
# 		text( H.b[,x_pheno],H.b[,y_pheno], label=unlist(lapply(HAP_GRP,function(x)paste(x,collapse="\n"))), col=COLS.hap_grp,pos=2,cex=1.2 )
# 		 # Significance		
# 		MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno] )
# 		abline(MOD,lwd=6,lty=2,col="mediumpurple3" )
# 		text( quantile(XLIM,.1),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col="mediumpurple3",cex=1.2 )
# 		MOD.w <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS )
# 		abline(MOD.w,lwd=6,lty=2,col="steelblue3" )
# 		text( quantile(XLIM,.1),quantile(YLIM,.07), label=paste("(ws)p=",formatC(summary(MOD.w)$coefficients[length(coef(MOD.w)),4],digits=2,format="e"),sep=""), col="steelblue3",cex=1.2 )
# 		# MOD.wp <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS.p )
# 		# abline(MOD.wp,lwd=6,lty=2,col="cadetblue3" )
# 		# text( quantile(XLIM,.1),quantile(YLIM,.12), label=paste("(wp)p=",formatC(summary(MOD.wp)$coefficients[length(coef(MOD.wp)),4],digits=2,format="e"),sep=""), col="cadetblue3",cex=1.2 )

# 		# MOD <- lm( H.b[,y_pheno] ~ H.b[,x_pheno] )
# 		# abline(MOD,lwd=6,lty=2,col="mediumpurple3" )
# 		# text( quantile(XLIM,.1),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col="mediumpurple3",cex=1.2 )
# 		# MOD.w <- lm( H.b[,y_pheno] ~ H.b[,x_pheno] )
# 		# abline(MOD.w,lwd=6,lty=2,col="cadetblue3" )
# 		# text( quantile(XLIM,.1),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),"Pr(>|t|)"],digits=2,format="e"),sep=""), col="grey25",cex=1.2 )
# 		dev.off()
# 	}

# 	## Compile Outputs
# 	COMPILE <- list( MOD.f=MOD.hap_grp.full,MOD.i=MOD.hap_grp.ind,BETA=BETAS.G.hap_grp.ind,SE=SES.G.hap_grp.ind,P=PS.G.hap_grp.ind )
# 	return(COMPILE)
# }

# GROUP_BY <- function( hap, col_tag, n_grps ) {
# 	which_cols <- grep( col_tag, colnames(BETAS[[hap]]) )
# 	TEMP <- sort(rowSums(apply( data.frame(BETAS[[hap]][,which_cols]), 2, function(x) (x-mean(x)) / sd(x) )))
# 	QUANTILES <- seq(0,1,length.out=n_grps+1)
# 	TEMP <- data.frame( HAP=names(TEMP), GRP=cut(TEMP, quantile(TEMP,QUANTILES),include.lowest=T,label=paste("G",1:(length(QUANTILES)-1),sep="") ) )
# 	HAP_GRP <- lapply( unique(TEMP$GRP), function(x) as.character(TEMP$HAP[TEMP$GRP==x]) ) 
# 	names(HAP_GRP) <- paste("G",1:length(HAP_GRP),sep="")
# 	return(HAP_GRP)
# }

# ##########################################
# ## RUN GROUPED HAPLOTYPE ANALYSIS
# OUT.haplo_grp <- list()

# ## Viatte 2015 Haplotype Groupings #######
#  # Specify Haplotype Positions (pos 11,71,74)
# hap <- names(OUT)[1]
#  # Specify Groups/Haplotypes
# HAP_GRP <- list()
# HAP_GRP$G4 <- c("SEA","SKR","SRL","LEA")
# HAP_GRP$G3 <- c("SRE","SRA","GRQ","PAA","SKA")
# HAP_GRP$G2 <- c("VEA","DRE","VRE","PRA","LRA")
# HAP_GRP$G1 <- c("VRA","VKA")
#  # Specify Tag for Saving Plots
# tag <- "Viatte"
#  # Run it...
# OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- HAPLO_GROUP( hap, HAP_GRP, tag )

# ## Custom Groupings for pos 11,71,74 #####
# hap <- names(OUT)[1]
# for ( N_Grps in 3:6 ) {
# 	HAP_GRP <- GROUP_BY(hap,"BL",N_Grps)
# 	tag <- paste("BySever",N_Grps,"g",sep="")
# 	OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- HAPLO_GROUP( hap, HAP_GRP, tag )
# }

# ## Custom Groupings for Shared Epitope ###
# hap <- names(OUT)[3]
# for ( N_Grps in 3:6 ) {
# 	HAP_GRP <- GROUP_BY(hap,"BL",N_Grps)
# 	tag <- paste("BySever",N_Grps,"g",sep="")
# 	OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- HAPLO_GROUP( hap, HAP_GRP, tag )
# }

# ## Custom Groupings for pos 11,13,71,74 ##
# hap <- names(OUT)[2]
# for ( N_Grps in 3:6 ) {
# 	HAP_GRP <- GROUP_BY(hap,"BL",N_Grps)
# 	tag <- paste("BySever",N_Grps,"g",sep="")
# 	OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- HAPLO_GROUP( hap, HAP_GRP, tag )
# }

# ## Custom Groupings based only on DAS ##
# for ( hap in names(OUT) ) {
# 	for ( N_Grps in 3:6 ) {
# 		HAP_GRP <- GROUP_BY(hap,"DAS_BL",N_Grps)
# 		tag <- paste("ByDAS",N_Grps,"g",sep="")
# 		OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- try( HAPLO_GROUP( hap, HAP_GRP, tag ), silent=T )
# 	}
# }




# #############################################################
# ## END OF DOC ###############################################
# #############################################################

