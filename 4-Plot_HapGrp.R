## Figure 4 for HLA Manuscript ##
## Collapsed Haplotype Analysis ##
## April 21, 2016 ##
## Kristopher Standish ##

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
COLS.ph <- COLS.list.2[c(6,3,2,5,4,1)][1:length(PHENOS)]
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
PathToAssocP <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160614_HLA_Assoc_0-P_Precise.Rdata" # 20160503_HLA_Assoc_0-P_Precise.Rdata"
PathToAssocB <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160614_HLA_Assoc_0-B_Precise.Rdata" # 20160503_HLA_Assoc_0-B_Precise.Rdata"
PathToAssocHap <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160614_HLA_Assoc_0-Hap_An.Rdata"
PathToHLA <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/"
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_ManuHLA_Fig4/",sep="")
PathToSave <- paste("/Users/kstandis/Data/Janssen/Data/HLA/Association/",DATE,"_HLA_Assoc",sep="")
dir.create( PathToPlot )

## Load Reference Tables
 # Amino Acid Key
AA.key <- read.table("/Users/kstandis/Data/Janssen/Data/HLA/Amino_Acid_Table.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(AA.key) <- AA.key[,1]
 # Viatte HLA-DRB1 (11,13,71,74) Susceptibility Table
VIATTE.2013.l <- readLines("/Users/kstandis/Data/Janssen/Data/HLA/Viatte_2013_DRB1_Types.txt")
VIATTE.2013.1 <- t(sapply(strsplit(VIATTE.2013.l,","),"[",1:5))
VIATTE.2013 <- VIATTE.2013.1[-1,]
colnames(VIATTE.2013) <- c( paste("Pos",VIATTE.2013.1[1,1:4],sep="_"), "OR" )
VIATTE.2013.2 <- apply( VIATTE.2013[,1:4], 2, function(x)AA.key[x,"X1_Let"])
VIATTE.2013 <- data.frame( HAP=apply(VIATTE.2013.2,1,function(x)paste(x,collapse="")), VIATTE.2013, stringsAsFactors=F )
rownames(VIATTE.2013) <- VIATTE.2013$HAP
VIATTE.2013$OR <- as.numeric(VIATTE.2013$OR)
 # Viatte HLA-DRB1 (11,71,74) Joint Erosion Table
VIATTE.2015 <- read.table("Data/Janssen/Data/HLA/Viatte_2015_DRB1_JointErosion_ORs.txt",fill=T)
VIATTE.2015 <- VIATTE.2015[,-c(10,15)]
colnames(VIATTE.2015) <- c("HAP","Pos_11","Pos_71","Pos_74","Fr.HET","Fr.HOM","Fr.HAP","OR.Joint","CI.Joint.L","CI.Joint.U","P.Joint","DEL.Larsen","CI.Larsen.L","CI.Larsen.U","P.Larsen")
VIATTE.2015$HAP <- gsub(",","",gsub( "[a-z]","",VIATTE.2015$HAP ))
VIATTE.2015$CI.Joint.L <- gsub("[Reference]","1",gsub( "(","",VIATTE.2015$CI.Joint.L,fixed=T ),fixed=T)
VIATTE.2015$CI.Joint.U <- gsub("[Reference]","1",gsub( ")","",VIATTE.2015$CI.Joint.L,fixed=T ),fixed=T)
VIATTE.2015$P.Joint <- as.numeric(gsub( "[a-z]","",VIATTE.2015$P.Joint ))
VIATTE.2015$DEL.Larsen <- as.numeric(gsub("p",".",gsub("[[:punct:]]","-", gsub(".","p",VIATTE.2015$DEL.Larsen,fixed=T) )))
VIATTE.2015$CI.Larsen.L <- as.numeric(gsub("p",".",gsub("[[:punct:]]","-", gsub(".","p",gsub("(","",VIATTE.2015$CI.Larsen.L,fixed=T),fixed=T) )))
VIATTE.2015$CI.Larsen.U <- as.numeric(gsub( ")","",VIATTE.2015$CI.Larsen.U,fixed=T ))
VIATTE.2015$P.Larsen <- as.numeric(gsub( "[a-z]","",VIATTE.2015$P.Larsen ))
rownames(VIATTE.2015) <- as.character(VIATTE.2015$HAP)

## Load HLA Association Results
load( PathToAssocP )
summary(P.out)
load( PathToAssocB )
summary(B.out)
load( PathToAssocHap )
load( gsub("Hap_An","Hap_An.Table",PathToAssocHap,fixed=T) )

## Load Janssen HLA Typing Results
 # Types
load( PathToTypes )
TYPES.l <- COMPILE
 # Amino Acids
load( PathToAA )
AA.l <- COMPILE

## Load HLA Types
HLA_AA.l <- read.table(paste(PathToHLA,"20160614_HLA_Assoc_HLA_AA_Table.txt",sep=""),header=T,sep="\t" )
colnames(HLA_AA.l) <- gsub(".","_",colnames(HLA_AA.l),fixed=T)
colnames(HLA_AA.l) <- gsub("-","_",colnames(HLA_AA.l),fixed=T)
HLA_TYP.l <- read.table(paste(PathToHLA,"20160614_HLA_Assoc_HLA_Types_Table.txt",sep=""),header=T,sep="\t" )
colnames(HLA_TYP.l) <- gsub(".","_",colnames(HLA_TYP.l),fixed=T)
colnames(HLA_TYP.l) <- gsub("-","_",colnames(HLA_TYP.l),fixed=T)
HLA_HAP.l <- read.table(paste(PathToHLA,"20160614_HLA_Assoc_HLA_Hap_Table.txt",sep=""),header=T,sep="\t" )
colnames(HLA_HAP.l) <- gsub(".","_",colnames(HLA_HAP.l),fixed=T)
colnames(HLA_HAP.l) <- gsub("-","_",colnames(HLA_HAP.l),fixed=T)

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
FT$RF <- as.numeric(factor(FT$RF=="Positive")) - 1

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
MG.hap <- merge( FT, HLA_HAP.l, by="ID" )

#############################################################
## HAPLOTYPE ANALYSIS (DRB1) ################################
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

##########################################
## Compile Results #######################

## Pull out Relevant Info from Models
PHENOS.temp <- grep("RF|ACPA",PHENOS,invert=T,value=T)
BETAS <- SES <- PS <- FREQS <- list()
BETAS.f <- SES.f <- PS.f <- FREQS.f <- list()
for ( hap in names(OUT) ) {
	# Individual Models
	PS[[hap]] <- matrix( unlist(lapply( PHENOS.temp, function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Pr(>|t|)"] )) )), byrow=F,ncol=length(PHENOS.temp) )
	BETAS[[hap]] <- matrix( unlist(lapply( PHENOS.temp, function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Estimate"] )) )), byrow=F,ncol=length(PHENOS.temp) )
	SES[[hap]] <- matrix( unlist(lapply( PHENOS.temp, function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Std. Error"] )) )), byrow=F,ncol=length(PHENOS.temp) )
	colnames(BETAS[[hap]]) <- colnames(SES[[hap]]) <- colnames(PS[[hap]]) <- PHENOS.temp
	rownames(BETAS[[hap]]) <- rownames(SES[[hap]]) <- rownames(PS[[hap]]) <- names(OUT[[hap]]$MOD.AA[[1]])
	# Full Models
	PS.f[[hap]] <- matrix( unlist(lapply( PHENOS.temp, function(x)summary(OUT[[hap]]$MOD.AA.F[[x]])$coefficients[rownames(PS[[hap]]),"Pr(>|t|)"] )), byrow=F,ncol=length(PHENOS.temp) )
	BETAS.f[[hap]] <- matrix( unlist(lapply( PHENOS.temp, function(x)summary(OUT[[hap]]$MOD.AA.F[[x]])$coefficients[rownames(PS[[hap]]),"Estimate"] )), byrow=F,ncol=length(PHENOS.temp) )
	SES.f[[hap]] <- matrix( unlist(lapply( PHENOS.temp, function(x)summary(OUT[[hap]]$MOD.AA.F[[x]])$coefficients[rownames(PS[[hap]]),"Std. Error"] )), byrow=F,ncol=length(PHENOS.temp) )
	colnames(BETAS.f[[hap]]) <- colnames(SES.f[[hap]]) <- colnames(PS.f[[hap]]) <- PHENOS.temp
	rownames(BETAS.f[[hap]]) <- rownames(SES.f[[hap]]) <- rownames(PS.f[[hap]]) <- rownames(PS[[hap]])
	# Haplotype Frequencies
	FREQS[[hap]] <- colSums( OUT[[hap]]$HAP[ OUT[[hap]]$SAMPS, ] )
}

#############################################################
## RE-CREATE VIATTE 2015 PLOTS ##############################
#############################################################

##########################################
## AMINO ACIDS ###########################

# ## FCT: Pull out Beta/P/SE Values for Position & Phenotype
# PULL_BETA <- function( gene, position, pheno ) {
# 	which_pos <- paste( "Pos",position,sep="_" )
# 	which_pos <- gsub("-",".",which_pos,fixed=T)
# 	TEMP.beta <- B.out$Dig_4$AA_DOS[[gene]][[pheno]][[which_pos]]
# 	TEMP.p <- P.out$Dig_4$AA_DOS[[gene]][[pheno]][[which_pos]]
# 	TEMP.maf.colnames <- grep( paste("Pr4",gene,which_pos,sep="_"), colnames(P.out$Dig_4$AA_DOS.arr[[gene]]),value=T )
# 	TEMP.maf <- apply( P.out$Dig_4$AA_DOS.arr[[gene]][,TEMP.maf.colnames], 2, mean ) / 2
# 	TEMP.freq <- apply( P.out$Dig_4$AA_DOS.arr[[gene]][,TEMP.maf.colnames], 2, sum )
# 	TEMP.pos <- gsub("Pos_","",which_pos)
# 	TEMP.aa <- names(TEMP.beta)
# 	TEMP.tag <- paste( TEMP.pos, TEMP.aa, sep="_" )
# 	OUT <- data.frame( tag=TEMP.tag, Pos=TEMP.pos, AA=TEMP.aa, MAF=TEMP.maf, FREQ=TEMP.freq, Beta=TEMP.beta, P=TEMP.p )
# 	colnames(OUT)[1:3] <- c("tag","Pos","AA")
# 	return(OUT)
# }

# ## FCT: Plot AA Beta Estimates for Severity vs Response
# PLOT_AA_BETA <- function( BETA.aa.list, position ) {
# 	## Plot Parameters
# 	TAB.sev <- BETA.aa.list$DAS_BL_MN[which(BETA.aa.list$DAS_BL_MN$Pos%in%position),]
# 	TAB.resp <- BETA.aa.list$DEL_MNe_MN[which(BETA.aa.list$DEL_MNe_MN$Pos%in%position),]
# 	XLIM <- extendrange( TAB.sev$Beta )
# 	YLIM <- extendrange( TAB.resp$Beta )
# 	COLS.temp <- c(COLS.list.2)[factor(TAB.sev$Pos)]
# 	## Model
# 	MOD <- lm( TAB.resp$Beta ~ TAB.sev$Beta, weights=TAB.sev$MAF )
# 	MOD.p <- summary(MOD)$coefficients[2,"Pr(>|t|)"]
# 	## Plot Results
# 	plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab="Beta: Disease Severity",ylab="Beta: Response",main=paste("DRB1: Pos",position) )
# 	abline( h=seq(-10,10,.5),v=seq(-10,10,.5),lty=3,col="grey50",lwd=1 )
# 	text( TAB.sev$Beta, TAB.resp$Beta, col=COLS.temp, label=TAB.sev$AA, cex=log10(TAB.sev$FREQ) )
# 	abline( MOD, lty=2,lwd=4,col="black" )
# 	text( XLIM[1],YLIM[1], label=paste("p=",formatC(MOD.p,format="e",digits=2),sep=""), pos=4 )
# 	if ( length(position)>1 ) { legend("topleft",legend=unique(TAB.sev$Pos),col=unique(COLS.temp),pch=20,title="Position",ncol=length(unique(TAB.sev$Pos)) ) }
# }


# ## Positions 11,71,74 (DRB1)
#  # Pull Beta Values for DRB1
# gene <- "DRB1"
# Positions <- c(11,71,74)
# BETA.aa <- list()
# for ( p in 1:length(PHENOS) ) {
# 	pheno <- PHENOS[p]
# 	TEMP.list <- lapply( Positions, function(x) PULL_BETA(gene,x,pheno) )
# 	BETA.aa[[pheno]] <- Reduce( rbind, TEMP.list )
# 	rownames(BETA.aa[[pheno]]) <- NULL
# }
#  # Create Plots of AA Betas
# par(mfrow=c(2,2))
# PLOT_AA_BETA( BETA.aa, Positions )
# PLOT_AA_BETA( BETA.aa, Positions[1] )
# PLOT_AA_BETA( BETA.aa, Positions[2] )
# PLOT_AA_BETA( BETA.aa, Positions[3] )

# ## Create Plots of AA Betas

# ## Positions 70:74 (DRB1)
#  # Pull Beta Values for DRB1
# Positions <- 70:74
# Positions <- Positions[which( unlist(apply(PAT_AA$DRB1[,paste("Pos",70:74,sep="_")],2,function(x)length(na.omit(unique(x))) ))>1 )]
# BETA.aa.se <- list()
# for ( p in 1:length(PHENOS) ) {
# 	pheno <- PHENOS[p]
# 	TEMP.list <- lapply( Positions, function(x) PULL_BETA("DRB1",x,pheno) )
# 	BETA.aa.se[[pheno]] <- Reduce( rbind, TEMP.list )
# 	rownames(BETA.aa.se[[pheno]]) <- NULL
# }
# par(mfrow=c(2,3))
# PLOT_AA_BETA( BETA.aa.se, Positions )
# PLOT_AA_BETA( BETA.aa.se, Positions[1] )
# PLOT_AA_BETA( BETA.aa.se, Positions[2] )
# PLOT_AA_BETA( BETA.aa.se, Positions[3] )
# PLOT_AA_BETA( BETA.aa.se, Positions[4] )

##########################################
## PLOT COLLAPSED HAPLOTYPE RESULTS ######

## FCT: Grouped Haplotype Analysis
PLOT_HAPS <- function( hap, tag ) {
	## Merge Clinical & Haplotype Data
	MG.HAP <- merge( OUT[[hap]]$HAP, FT, by.x="row.names",by.y="ID" )
	
	######################################
	## Plot P-Values for Collapsed Haplotypes

	## Pull Association Results
	 # P-Values
	P.comp <- OUT[[hap]]$P
	P.comp.2 <- P.comp[,-grep("ANOVA",colnames(P.comp))]
	 # Betas
	B.comp <- OUT[[hap]]$B
	## Plotting Parameters
	 # Nyholt Correction
	HAPS.uniq <- setdiff( colnames(P.comp), "ANOVA" )
	HAP.com <- OUT[[hap]]$FR.com
	NYH.ph <- NYHOLT( MG.HAP[,PHENOS] )
	NYH.hap <- NYHOLT( MG.HAP[,HAP.com] )
	NYH <- NYH.ph*NYH.hap

	## Plot Haplotype Association P-Values (New Version)
	XLIM <- c(1,ncol(P.comp.2))
	YLIM <- c( 0,-log10(min(P.comp.2)/30) )
	MAIN.1 <- paste("ANOVA:",hap)
	MAIN.2 <- paste("Haplotype Regression:",hap)
	png( paste(PathToPlot,"/DRB1_",hap,"_5A-HapAssoc_Ps.png",sep=""), height=1000,width=2000,pointsize=30 ) # width=2500
	# png( paste(PathToPlot,"_3AB_TYP",gene,".png",sep=""),height=1000,width=2500,pointsize=30 )
	layout( matrix(1:2,byrow=T,ncol=2), width=c(1,4) )
	par(mar=c(7,4,3,2))
	 # ANOVA
	barplot( -log10(P.comp[,"ANOVA"]), beside=T, col=COLS.ph,ylim=YLIM,las=2,border=NA,main="ANOVA",ylab="-log10(p)")
	abline( h=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH.ph)),lty=2,col=COLS.cor,lwd=3 )
	barplot( -log10(P.comp[,"ANOVA"]), beside=T, col=COLS.ph,ylim=YLIM,las=2,border=NA,add=T )
	 # Regression (P-Values)
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,main=MAIN.2,ylab="-log10(p)",xaxt="n",xlab="")
	axis( 1, at=1:XLIM[2], label=colnames(P.comp.2), las=2 )
	abline( h=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH)),lty=2,col=COLS.cor,lwd=4 )
	legend( "topleft", legend=rownames(P.comp.2),title="Phenotype",col=COLS.ph,pch=16,pt.cex=1.5, cex=1.2,ncol=nrow(P.comp.2),bg="white" ) # ncol=ceiling(nrow(TYP)/4),
	points( rep(1:XLIM[2],each=3), -log10(P.comp.2), col=adjustcolor(COLS.ph,.8),pch=16,cex=1.5 )
	dev.off()

	##########################################
	## Plot Response vs Severity Data

	## Pull Phenotype Names/Numbers
	P.n <- grep("RF|ACPA",rownames(B.comp),value=T,invert=T) # colnames(BETAS[[hap]])
	P.num <- length(P.n)

	## Create Plots for Various Phenotypes
	for ( p in grep("_MNe_MN$",P.n) ) {
		
		## Pull Betas & other info
		y_pheno <- P.n[p]
		x_pheno <- P.n[p+P.num/2] # x_pheno <- P.n[p+4]
		H.n <- rownames(BETAS[[hap]])
		H.f <- FREQS[[hap]][H.n]
		H.b <- BETAS[[hap]][H.n,]
		H.s <- SES[[hap]][H.n,]
		H.p <- PS[[hap]][H.n,]
		H.b <- BETAS.f[[hap]][H.n,]
		H.s <- SES.f[[hap]][H.n,]
		H.p <- PS.f[[hap]][H.n,]

		## Create Plot/File Location
		 # Plotting Parameters
		y_lab <- "Beta: Drug Response" # paste("Beta:",x_pheno)
		x_lab <- "Beta: Disease Severity" # paste("Beta:",y_pheno)
		main.1 <- paste( "Effect Sizes of",hap,"Haplotypes:",strsplit(x_pheno,"_")[[1]][1] )
		main.2 <- paste( "Response vs Severity Estimates:",strsplit(x_pheno,"_")[[1]][1] )

		png( paste(PathToPlot,"HapGroup_",hap,"-",tag,"_1_",x_pheno,".png",sep=""),height=1000,width=2000,pointsize=30 )
		par(mfrow=c(1,2))
		## Haplotype Effect Sizes
		XLIM <- c(0,1+nrow(H.b))
		YLIM <- extendrange(H.b[,c(y_pheno,x_pheno)],f=.25)
		ORDER <- order(H.b[,x_pheno])
		plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",main=main.1,ylab="Beta",xlab="" )
		abline( h=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
		abline( h=0,lty=1,col="grey50",lwd=1 )
		abline( v=1:nrow(H.b),lty=3,col="grey50",lwd=1 )
		points( 1:nrow(H.b),H.b[ORDER,x_pheno], pch=16,cex=2*log10(1+H.f[ORDER]),col=adjustcolor(COLS.ph[x_pheno],.5) )
		axis(1,at=1:nrow(H.b),label=rownames(H.b)[ORDER],las=2)
		points( 1:nrow(H.b),H.b[ORDER,y_pheno], pch=16,cex=2*log10(1+H.f[ORDER]),col=adjustcolor(COLS.ph[y_pheno],.5) )
		legend("topright",pch=16,col=adjustcolor(COLS.ph[c(x_pheno,y_pheno)],.8),legend=c(x_pheno,y_pheno),pt.cex=1.8,ncol=2 )
		
		## Model Response vs Severity Haplotypes
		XLIM <- extendrange(H.b[,x_pheno],f=.25)
		YLIM <- extendrange(H.b[,y_pheno],f=.25)
		WTS <- H.f
		# WTS <- 1 / rowMeans(H.s[,c(y_pheno,x_pheno)]) 
		plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,main=main.2,xlab=x_lab,ylab=y_lab )
		abline( h=seq(-5,5,.2),v=seq(-5,5,.2),lty=3,col="grey50",lwd=1 )
		abline( h=0,v=0,lty=1,col="grey50",lwd=1 )
		points( H.b[,x_pheno],H.b[,y_pheno], pch=16,col=adjustcolor(COLS.beta,.5), cex=2*log10(1+WTS) )
		# arrows( H.b[,x_pheno]+H.s[,x_pheno],H.b[,y_pheno],H.b[,x_pheno]-H.s[,x_pheno],H.b[,y_pheno], code=3,angle=90,lwd=5,col=H.c )
		# arrows( H.b[,x_pheno],H.b[,y_pheno]+H.s[,y_pheno],H.b[,x_pheno],H.b[,y_pheno]-H.s[,y_pheno], code=3,angle=90,lwd=5,col=H.c )
		text( H.b[,x_pheno],H.b[,y_pheno]+.02*diff(YLIM), label=rownames(H.b), pos=3,cex=1.2 )
		# points( H.b[,x_pheno],H.b[,y_pheno], pch=20,col=H.c, lwd=5,cex=WTS/mean(WTS) )
		 # Significance
		MOD.w <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS )
		abline(MOD.w,lwd=6,lty=2,col=COLS.beta )
		text( quantile(XLIM,0),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD.w)$coefficients[length(coef(MOD.w)),4],digits=2,format="e"),sep=""),cex=1.2, pos=4 )
		# legend( "topright")
		# MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno] )
		# abline(MOD,lwd=6,lty=2,col=COLS.mods[2] )
		# text( quantile(XLIM,.1),quantile(YLIM,.07), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col=COLS.mods[2],cex=1.2 )
		dev.off()
	}
}

## Custom Groupings (Disease Severity) ###
for ( hap in names(OUT) ) {
	tag <- hap
	PLOT_HAPS( hap, tag )
}

##########################################
## GROUP COLLAPSED HAPLOTYPE #############

## FCT: Group Haplotypes
GROUP_BY <- function( hap, col_tag, n_grps ) {
	which_cols <- grep( col_tag, colnames(BETAS[[hap]]) )
	TEMP <- sort(rowSums(apply( data.frame(BETAS[[hap]][,which_cols]), 2, function(x) (x-mean(x)) / sd(x) )))
	QUANTILES <- seq(0,1,length.out=n_grps+1)
	TEMP <- data.frame( HAP=names(TEMP), GRP=cut(TEMP, quantile(TEMP,QUANTILES),include.lowest=T,label=paste("G",1:(length(QUANTILES)-1),sep="") ) )
	HAP_GRP <- lapply( unique(TEMP$GRP), function(x) as.character(TEMP$HAP[TEMP$GRP==x]) ) 
	names(HAP_GRP) <- paste("G",1:length(HAP_GRP),sep="")
	return(HAP_GRP)
}
# HAP_GRP <- GROUP_BY( hap, "DAS_BL_MN", 4 )

## FCT: Grouped Haplotype Analysis
GROUPED_AN <- function( hap, HAP_GRP, tag ) {
	# Set New HAP_GRP variable to be modified in first half of function
	HAP_GRP.2 <- HAP_GRP
	 # Specify Group Colors
	COLS.hap_grp <- c("seagreen2","magenta2","deepskyblue2","chocolate1","brown1","cyan2","blueviolet","yellow1","deeppink2","darkolivegreen1","goldenrod2")
	COLS.hap_grp <- COLS.list.2
	 # Specify any non-assigned haplotypes
	COLS.hap_grp <- c( COLS.hap_grp[1:length(HAP_GRP.2)], "black" )
	HAP_GRP.2$Other <- setdiff( rownames(BETAS[[hap]]),Reduce(union,HAP_GRP.2) )
	names(COLS.hap_grp) <- names(HAP_GRP.2)
	KEY.hap_grp <- data.frame( HAP=unlist(HAP_GRP.2), GRP=rep(names(HAP_GRP.2),lapply(HAP_GRP.2,length)),COL=rep(COLS.hap_grp,lapply(HAP_GRP.2,length)) )
	rownames(KEY.hap_grp) <- unlist(HAP_GRP.2)

	## Merge Clinical & Haplotype Data
	MG.HAP <- merge( OUT[[hap]]$HAP, FT, by.x="row.names",by.y="ID" )
	
	######################################
	## Plot P-Values for Collapsed Haplotypes

	## Pull Association Results
	 # P-Values
	P.comp <- OUT[[hap]]$P
	 # Betas
	B.comp <- OUT[[hap]]$B
	## Plotting Parameters
	 # Nyholt Correction
	HAPS.uniq <- setdiff( colnames(P.comp), "ANOVA" )
	HAP.com <- OUT[[hap]]$FR.com
	NYH.ph <- NYHOLT( MG.HAP[,PHENOS] )
	NYH.hap <- NYHOLT( MG.HAP[,HAP.com] )
	NYH <- NYH.ph*NYH.hap

	## Plot Haplotype Association P-Values (New Version)
	YLIM <- c( 0,-log10(min(P.comp)/30) )
	png( paste(PathToPlot,"/DRB1_",hap,"_3p-HapAssoc.png",sep=""), height=1000,width=2000,pointsize=30 )
	layout( matrix(1:2,byrow=T,ncol=2), width=c(1,4) )
	par(mar=c(7,4,3,2))
	 # ANOVA
	barplot( -log10(P.comp[,"ANOVA"]), beside=T, col=COLS.ph,ylim=YLIM,las=2,border=NA,main="Haplotype ANOVA",ylab="-log10(p)")
	abline( h=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH.ph)),lty=2,col=COLS.cor,lwd=3 )
	barplot( -log10(P.comp[,"ANOVA"]), beside=T, col=COLS.ph,ylim=YLIM,las=2,border=NA,add=T )
	 # Dosage
	barplot( -log10(P.comp[,HAP.com]), beside=T, col=COLS.ph,ylim=YLIM,las=2,border=NA,main=paste("Haplotype Dosage Regression: HLA-DRB1:",hap),ylab="-log10(p)")
	abline( h=0:20,lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(NYH)),lty=2,col=COLS.cor,lwd=3 )
	legend( "topleft",legend=rownames(P.comp[,HAP.com]),fill=COLS.ph,ncol=nrow(P.comp),cex=.8,border=NA)
	barplot( -log10(P.comp[,HAP.com]), beside=T, col=COLS.ph,ylim=YLIM,las=2,border=NA,add=T )
	dev.off()

	##########################################
	## Plot Response vs Severity Data

	## Pull Phenotype Names/Numbers
	P.n <- grep("RF|ACPA",rownames(B.comp),value=T,invert=T) # colnames(BETAS[[hap]])
	P.num <- length(P.n)

	## Create Plots for Various Phenotypes
	for ( p in grep("_MNe_MN$",P.n) ) {
		
		## Pull Betas & other info
		y_pheno <- P.n[p]
		x_pheno <- P.n[p+P.num/2] # x_pheno <- P.n[p+4]
		H.n <- rownames(BETAS[[hap]])
		H.c <- as.character(KEY.hap_grp[H.n,"COL"])
		H.f <- FREQS[[hap]][H.n]
		H.b <- BETAS[[hap]][H.n,]
		H.s <- SES[[hap]][H.n,]
		H.p <- PS[[hap]][H.n,]
		H.b <- BETAS.f[[hap]][H.n,]
		H.s <- SES.f[[hap]][H.n,]
		H.p <- PS.f[[hap]][H.n,]

		## Create Plot/File Location
		 # Plotting Parameters
		y_lab <- "Beta: Drug Response" # paste("Beta:",x_pheno)
		x_lab <- "Beta: Disease Severity" # paste("Beta:",y_pheno)
		main.1 <- paste( "Effect Sizes of",hap,"Haplotypes:",strsplit(x_pheno,"_")[[1]][1] )
		main.2 <- paste( "Response vs Severity Estimates:",strsplit(x_pheno,"_")[[1]][1] )

		png( paste(PathToPlot,"HapGroup_",hap,"-",tag,"_1_",x_pheno,".png",sep=""),height=1000,width=2000,pointsize=30 )
		par(mfrow=c(1,2))
		## Haplotype Effect Sizes
		XLIM <- c(0,1+nrow(H.b))
		YLIM <- extendrange(H.b[,c(y_pheno,x_pheno)],f=.25)
		ORDER <- order(H.b[,x_pheno])
		plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",main=main.1,ylab="Beta",xlab="" )
		abline( h=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
		abline( h=0,lty=1,col="grey50",lwd=1 )
		abline( v=1:nrow(H.b),lty=3,col="grey50",lwd=1 )
		points( 1:nrow(H.b),H.b[ORDER,x_pheno], pch=9,cex=2*log10(1+H.f[ORDER]),col=H.c[ORDER],lwd=5 )
		axis(1,at=1:nrow(H.b),label=rownames(H.b)[ORDER],las=2)
		points( 1:nrow(H.b),H.b[ORDER,y_pheno], pch=10,cex=2*log10(1+H.f[ORDER]),col=H.c[ORDER], lwd=5 )
		legend("topright",pch=9:10,col="black",legend=c(x_pheno,y_pheno),pt.cex=1.8,pt.lwd=4,ncol=2 )
		
		## Model Response vs Severity Haplotypes
		XLIM <- extendrange(H.b[,x_pheno],f=.25)
		YLIM <- extendrange(H.b[,y_pheno],f=.25)
		WTS <- H.f
		# WTS <- 1 / rowMeans(H.s[,c(y_pheno,x_pheno)]) 
		plot( H.b[,x_pheno],H.b[,y_pheno], pch=20,col=H.c, lwd=5,cex=log10(1+WTS), xlim=XLIM,ylim=YLIM,main=main.2,xlab=x_lab,ylab=y_lab )
		abline( h=seq(-5,5,.5),v=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
		abline( h=0,v=0,lty=1,col="grey50",lwd=1 )
		# arrows( H.b[,x_pheno]+H.s[,x_pheno],H.b[,y_pheno],H.b[,x_pheno]-H.s[,x_pheno],H.b[,y_pheno], code=3,angle=90,lwd=5,col=H.c )
		# arrows( H.b[,x_pheno],H.b[,y_pheno]+H.s[,y_pheno],H.b[,x_pheno],H.b[,y_pheno]-H.s[,y_pheno], code=3,angle=90,lwd=5,col=H.c )
		text( H.b[,x_pheno],H.b[,y_pheno]+.02*diff(YLIM), label=rownames(H.b), col=H.c, pos=4,cex=1.2 )
		points( H.b[,x_pheno],H.b[,y_pheno], pch=20,col=H.c, lwd=5,cex=WTS/mean(WTS) )
		 # Significance
		COLS.mods <- c("black","grey30")
		MOD.w <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS )
		abline(MOD.w,lwd=6,lty=2,col=COLS.mods[1] )
		text( quantile(XLIM,0),quantile(YLIM,.02), label=paste("p(w)=",formatC(summary(MOD.w)$coefficients[length(coef(MOD.w)),4],digits=2,format="e"),sep=""), col=COLS.mods[1],cex=1.2, pos=4 )
		# legend( "topright")
		# MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno] )
		# abline(MOD,lwd=6,lty=2,col=COLS.mods[2] )
		# text( quantile(XLIM,.1),quantile(YLIM,.07), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col=COLS.mods[2],cex=1.2 )
		dev.off()
	}

	##########################################
	## Group Haplotypes ##

	## Greate Dosage Array for Patients & Haplotype Groups
	HAP_GRP.arr <- matrix(unlist(lapply( HAP_GRP, function(x) rowSums(OUT[[hap]]$HAP[,which(colnames(OUT[[hap]]$HAP)%in%x)]) )),byrow=F,ncol=length(HAP_GRP))
	colnames(HAP_GRP.arr) <- names(HAP_GRP)
	rownames(HAP_GRP.arr) <- rownames(OUT[[hap]]$HAP)
	MG.hap_grp <- merge( HAP_GRP.arr, FT, by.x="row.names",by.y="ID" )
	MOD.hap_grp.full <- MOD.hap_grp.ind <- list()
	 # Model Phenotypes vs Haplotype Groups
	for ( p in grep("ACPA|RF",PHENOS,invert=T) ) {
	# for ( p in grep("_MNe_MN$",PHENOS) ) {
		pheno <- as.character(PH.COV[p,"PHENOS"])
		cov <- as.character(PH.COV[p,"COVS"])
		## Model w/ All Hap Groups
		formula.full <- as.formula(paste(pheno,"~",cov,"+",paste(names(HAP_GRP),collapse="+")))
		MOD.hap_grp.full[[pheno]] <- lm( formula.full, data=MG.hap_grp )
		## Model w/ Individual Hap Groups
		MOD.hap_grp.ind[[pheno]] <- list()
		for ( g in 1:length(HAP_GRP) ) {
			g_tag <- paste("G",g,sep="")
			formula.temp <- as.formula(paste(pheno,"~",cov,"+",g_tag))
			MOD.hap_grp.ind[[pheno]][[g_tag]] <- lm( formula.temp, data=MG.hap_grp )
		}
	}

	 # Full Models
	# BETAS.hap_grp.full <- lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[,"Estimate"] )
	# SES.hap_grp.full <- lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[,"Std. Error"] )
	# PS.hap_grp.full <- lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[,"Pr(>|t|)"] )
	# OUT.hap_grp.full <- lapply( names(BETAS.hap_grp.full),function(x)data.frame(BETA=BETAS.hap_grp.full[[x]],SE=SES.hap_grp.full[[x]],P=PS.hap_grp.full[[x]]) )
	# names(OUT.hap_grp.full) <- names(BETAS.hap_grp.full)
	# BETAS.G.hap_grp.full <- matrix(unlist( lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[names(HAP_GRP),"Estimate"] )),byrow=F,nrow=length(HAP_GRP) )
	# SES.G.hap_grp.full <- matrix(unlist( lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[names(HAP_GRP),"Std. Error"] )),byrow=F,nrow=length(HAP_GRP) )
	# PS.G.hap_grp.full <- matrix(unlist( lapply( MOD.hap_grp.full, function(x)summary(x)$coefficients[names(HAP_GRP),"Pr(>|t|)"] )),byrow=F,nrow=length(HAP_GRP) )
	# colnames(BETAS.G.hap_grp.full) <- colnames(SES.G.hap_grp.full) <- colnames(PS.G.hap_grp.full) <- names(BETAS.hap_grp.full)
	# rownames(BETAS.G.hap_grp.full) <- rownames(SES.G.hap_grp.full) <- rownames(PS.G.hap_grp.full) <- names(HAP_GRP)
 	# OUT.G.hap_grp.full <- lapply( names(BETAS.G.hap_grp.full),function(x)data.frame(BETA=BETAS.G.hap_grp.full[,x],SE=SES.G.hap_grp.full[,x],P=PS.G.hap_grp.full[,x]) )
	# names(OUT.G.hap_grp.full) <- names(BETAS.G.hap_grp.full)
	 # Individual Models
	BETAS.G.hap_grp.ind <- matrix(unlist( lapply( MOD.hap_grp.ind, function(x) unlist(lapply(names(x),function(y)summary(x[[y]])$coefficients[y,"Estimate"])) )),byrow=F,nrow=length(HAP_GRP) )
	SES.G.hap_grp.ind <- matrix(unlist( lapply( MOD.hap_grp.ind, function(x) unlist(lapply(names(x),function(y)summary(x[[y]])$coefficients[y,"Std. Error"])) )),byrow=F,nrow=length(HAP_GRP) )
	PS.G.hap_grp.ind <- matrix(unlist( lapply( MOD.hap_grp.ind, function(x) unlist(lapply(names(x),function(y)summary(x[[y]])$coefficients[y,"Pr(>|t|)"])) )),byrow=F,nrow=length(HAP_GRP) )
	colnames(BETAS.G.hap_grp.ind) <- colnames(SES.G.hap_grp.ind) <- colnames(PS.G.hap_grp.ind) <- names(MOD.hap_grp.ind)
	rownames(BETAS.G.hap_grp.ind) <- rownames(SES.G.hap_grp.ind) <- rownames(PS.G.hap_grp.ind) <- names(HAP_GRP)

	## Plot
	# H.b <- BETAS.G.hap_grp.full
	# H.s <- SES.G.hap_grp.full
	# H.p <- PS.G.hap_grp.full
	H.b <- BETAS.G.hap_grp.ind
	H.s <- SES.G.hap_grp.ind
	H.p <- PS.G.hap_grp.ind
	H.f <- colSums(HAP_GRP.arr)
	P.n <- colnames(H.b)
	P.num <- length(P.n)
	WTS <- H.f # colMeans(HAP_GRP.arr)
	# par(mfrow=c(2,2))
	for ( p in grep("_MNe_MN$",PHENOS) ) {
		y_pheno <- P.n[p]
		x_pheno <- P.n[p+P.num/2] # x_pheno <- P.n[p+4]
		# WTS <- 1 / rowMeans(H.s[,c(y_pheno,x_pheno)]) 
		XLIM <- extendrange(H.b[,x_pheno],f=.5)
		YLIM <- extendrange(H.b[,y_pheno],f=.5)
		png( paste(PathToPlot,"HapGroup_",hap,"-",tag,"_2_",x_pheno,".png",sep=""),height=1000,width=1000,pointsize=30 )
		plot( H.b[,y_pheno] ~ H.b[,x_pheno], cex=log10(1+H.f),pch=20,col=COLS.hap_grp, xlim=XLIM,ylim=YLIM, main=paste("Response vs Severity:",hap,"Haplo Groups"),xlab=x_lab,ylab=y_lab ) # ylab=paste("Beta:",y_pheno)
		abline(h=0,v=0,lty=1,col="grey50",lwd=1)
		arrows( H.b[,x_pheno], H.b[,y_pheno]+H.s[,y_pheno], H.b[,x_pheno], H.b[,y_pheno]-H.s[,y_pheno], col=COLS.hap_grp, code=3,angle=90,length=.1,lwd=4 )
		arrows( H.b[,x_pheno]+H.s[,x_pheno] ,H.b[,y_pheno], H.b[,x_pheno]-H.s[,x_pheno], H.b[,y_pheno], col=COLS.hap_grp, code=3,angle=90,length=.1,lwd=4 )
		 # Label Haplotypes in each Group
		text( H.b[,x_pheno],H.b[,y_pheno], label=unlist(lapply(HAP_GRP,function(x)paste(x,collapse="\n"))), col=COLS.hap_grp,pos=2,cex=1.2 )
		 # Significance		
		MOD.w <- lm( H.b[,y_pheno]~H.b[,x_pheno], weights=WTS )
		abline(MOD.w,lwd=6,lty=2,col=COLS.mods[1] )
		text( quantile(XLIM,0),quantile(YLIM,.02), label=paste("p(w)=",formatC(summary(MOD.w)$coefficients[length(coef(MOD.w)),4],digits=2,format="e"),sep=""), col=COLS.mods[1],cex=1.2, pos=4 )
		# MOD <- lm( H.b[,y_pheno]~H.b[,x_pheno] )
		# abline(MOD,lwd=6,lty=2,col=COLS.mods[2] )
		# text( quantile(XLIM,.1),quantile(YLIM,.07), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col=COLS.mods[2],cex=1.2 )
		dev.off()
	}

	## Compile Outputs
	COMPILE <- list( MOD.f=MOD.hap_grp.full,MOD.i=MOD.hap_grp.ind,BETA=BETAS.G.hap_grp.ind,SE=SES.G.hap_grp.ind,P=PS.G.hap_grp.ind, HAP_GRP=HAP_GRP )
	return(COMPILE)
}
GROUPED_AN( hap, HAP_GRP, tag )

## RUN GROUPED HAPLOTYPE ANALYSIS
OUT.haplo_grp <- list()

## Custom Groupings (Disease Severity) ###
for ( hap in names(OUT) ) {
	for ( N_Grps in 4 ) {
		HAP_GRP <- GROUP_BY(hap,"BL",N_Grps)
		tag <- paste("BySever",N_Grps,"g",sep="")
		# OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- try( GROUPED_AN( hap, HAP_GRP, tag ), silent=T )
		PLOT_HAPS( hap, tag )
# 		HAP_GRP <- GROUP_BY(hap,"DAS_BL",N_Grps)
# 		tag <- paste("ByDAS",N_Grps,"g",sep="")
# 		OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- try( GROUPED_AN( hap, HAP_GRP, tag ), silent=T )
	}
}

## Viatte 2015 Haplotype Groupings #######
 # Specify Haplotype Positions (pos 11,71,74)
hap <- names(OUT)[1]
 # Specify Groups/Haplotypes
HAP_GRP <- list()
HAP_GRP$G4 <- c("SEA","SKR","SRL","LEA")
HAP_GRP$G3 <- c("SRE","SRA","GRQ","PAA","SKA")
HAP_GRP$G2 <- c("VEA","DRE","VRE","PRA","LRA")
HAP_GRP$G1 <- c("VRA","VKA")
 # Specify Tag for Saving Plots
tag <- "Viatte"
 # Run it...
OUT.haplo_grp[[paste(hap,tag,sep="_")]] <- GROUPED_AN( hap, HAP_GRP, tag )

VIATTE_vs <- function( xvals, yvals, freqs, names ) {
	plot( yvals ~ xvals, cex=log10(1+freqs[names]), col=adjustcolor(COLS.list.2[1],.7),pch=16 )
	text( yvals ~ xvals, label=names, pos=4 )
	MOD <- lm( yvals ~ xvals, weights=freqs[names] )
	P <- summary(MOD)$coefficients["xvals","Pr(>|t|)"]
	abline(MOD)
	# summary(MOD)
	text( max(xvals,na.rm=T),max(yvals,na.rm=T), paste("p=",formatC(P,3,format="e")),pos=2 )
}

##########################################
## VIATTE RESULTS vs JANSSEN #############

## Viatte 2013 Collapsed Haplotypes ######
 # Positions 11,13,71,74 vs Viatte '13 Susceptibility ORs
VIATTE.2013.haps <- VIATTE.2013$HAP
MY.11137174.haps <- rownames(BETAS$p11137174)
INT.haps <- sort(intersect( VIATTE.2013.haps, MY.11137174.haps ))
par(mfrow=c(1,3))
 # Viatte vs Disease Severity
VIATTE_vs( log2(VIATTE.2013[INT.haps,"OR"]), BETAS$p11137174[INT.haps,"DAS_BL_MN"], FREQS$p11137174, INT.haps )
 # Viatte vs Treatment Response
VIATTE_vs( log2(VIATTE.2013[INT.haps,"OR"]), BETAS$p11137174[INT.haps,"DEL_MNe_MN"], FREQS$p11137174, INT.haps )
 # JJ Severity vs Response
VIATTE_vs( BETAS$p11137174[INT.haps,"DAS_BL_MN"], BETAS$p11137174[INT.haps,"DEL_MNe_MN"], FREQS$p11137174, INT.haps )

## Viatte 2015 Collapsed Haplotypes ######
 # Positions 11,13,71,74 vs Viatte '13 Susceptibility ORs
VIATTE.2015.haps <- VIATTE.2015$HAP
MY.117174.haps <- rownames(BETAS$p117174)
INT.haps <- sort(intersect( VIATTE.2015.haps, MY.117174.haps ))
par(mfrow=c(2,3))
## Viatte Joint Estimates
 # Viatte vs Disease Severity
VIATTE_vs( log2(VIATTE.2015[INT.haps,"OR.Joint"]), BETAS$p117174[INT.haps,"DAS_BL_MN"], FREQS$p117174, INT.haps )
 # Viatte vs Treatment Response
VIATTE_vs( log2(VIATTE.2015[INT.haps,"OR.Joint"]), BETAS$p117174[INT.haps,"DEL_MNe_MN"], FREQS$p117174, INT.haps )
 # JJ Severity vs Response
VIATTE_vs( BETAS$p117174[INT.haps,"DAS_BL_MN"], BETAS$p117174[INT.haps,"DEL_MNe_MN"], FREQS$p117174, INT.haps )
## Viatte Larsen Estimates
 # Viatte vs Disease Severity
VIATTE_vs( VIATTE.2015[INT.haps,"DEL.Larsen"], BETAS$p117174[INT.haps,"DAS_BL_MN"], FREQS$p117174, INT.haps )
 # Viatte vs Treatment Response
VIATTE_vs( VIATTE.2015[INT.haps,"DEL.Larsen"], BETAS$p117174[INT.haps,"DEL_MNe_MN"], FREQS$p117174, INT.haps )
 # JJ Severity vs Response
VIATTE_vs( BETAS$p117174[INT.haps,"DAS_BL_MN"], BETAS$p117174[INT.haps,"DEL_MNe_MN"], FREQS$p117174, INT.haps )



#############################################################
## END OF DOC ###############################################
#############################################################

