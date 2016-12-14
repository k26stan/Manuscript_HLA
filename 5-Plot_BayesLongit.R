## Figure 5 for HLA Manuscript ##
## Bayes Longitudinal Model Results ##
## Compare vs Mean Model ##
## June 15, 2016 ##
## Kristopher Standish ##

## Load Packages
library(nlme)
library(gplots)
library(brms)
library(vioplot)
library(xtable)

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
SHOW_COLORS <- function() { barplot( 1:length(COLS.list.2),col=COLS.list.2,names.arg=1:length(COLS.list.2)) }

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
PathToData <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
PathToAssocP <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160614_HLA_Assoc_0-P_Precise.Rdata" # 20160503_HLA_Assoc_0-P_Precise.Rdata"
PathToAssocB <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160614_HLA_Assoc_0-B_Precise.Rdata" # 20160503_HLA_Assoc_0-B_Precise.Rdata"
PathToAssocM <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160613_HLA_Assoc_0-M_Precise.Rdata" # 20160503_HLA_Assoc_0-B_Precise.Rdata"
PathToAssocHap <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/20160614_HLA_Assoc_0-Hap_An.Rdata"
PathToHLA <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/"
PathToHLAMod.1 <- "/Users/kstandis/Data/Janssen/Data/HLA/Association/"
PathToBayesHLA <- "Data/Janssen/Data/HLA/Association/20160615_BayesLMM_HLA_ModFits.Rdata"
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_ManuHLA_Fig5/",sep="")
PathToSave <- paste("/Users/kstandis/Data/Janssen/Data/HLA/Association/",DATE,"_HLA_Assoc",sep="")
dir.create( PathToPlot )

#####################################
## HLA Reference Tables #############

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

#####################################
## Association Results ##############

## Load HLA Association Results
load( PathToAssocP ) # P.out
summary(P.out)
load( PathToAssocB ) # B.out
summary(B.out)
load( PathToAssocM ) # M.out
summary(M.out)
load( PathToAssocHap ) # OUT
load( gsub("Hap_An","Hap_An.Table",PathToAssocHap,fixed=T) ) # HAP_OUT
load( PathToBayesHLA ) # HLA.mod.fits

## FCT: Load Previously Compiled Models & Data (BRMS)
LOAD_PREV_MODS <- function( Path, Tag ) {
	Temp.files <- grep( "Rdata$",list.files( Path ), value=T )
	Temp.files.2 <- grep( "Model_Summary", Temp.files, value=T,invert=T )
	Temp.files.mods <- grep("Model",Temp.files.2, value=T,invert=F )
	Temp.files.meta <- grep("Model",Temp.files.2, value=T,invert=T )
	# for ( file in Temp.files.meta ) { load(paste(Path,file,sep="")) }
	OUT <- list()
	for ( file in Temp.files.mods ) {
		print(file)
		mod.tag <- gsub( "Rdata.Model.","", gsub(".Rdata","",file, fixed=T),fixed=T)
		print(mod.tag)
		load(paste(Path,file,sep=""))
		OUT[[mod.tag]] <- temp.model
	}
	return(OUT)
}
 # Load HLA Models
load.tag <- "HLA"
if ( !exists("JJ") ) { JJ <- list() }
JJ[[load.tag]] <- LOAD_PREV_MODS( PathToHLAMod.1, load.tag )

#####################################
## Janssen Typing Results ###########

## Load Janssen HLA Typing Results
 # Types
load( PathToTypes )
TYPES.l <- COMPILE
 # Amino Acids
load( PathToAA )
AA.l <- COMPILE

## Load HLA Type Tables
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
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FT.l <- read.table( PathToFT, sep="\t",header=T )

#############################################################
## ORGANIZE DATA ############################################
#############################################################

#########################################
## CLINICAL FULL TABLE ##################

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

#########################################
## CLINICAL TIME-SERIES DATA ############

## Add Column for "TRT" (Treatment: All but Week 0)
TRT <- rep(1,nrow(TAB.l))
TRT[which(TAB.l$WK==0)] <- 0
TAB <- data.frame( TAB.l, TRT )

## Add ID Column w/ Short Hand
TAB.id.short <- sapply( strsplit( TAB$IID, "-" ),"[",1 )
TAB <- data.frame( ID=TAB.id.short, TAB[,-1] )

## Change ACPA to 1/0
TAB$ACPA <- as.numeric(TAB$ACPA=="Positive")
TAB$RF <- as.numeric(TAB$RF=="Positive")
TAB$RF_ACPA <- as.numeric(TAB$RF_ACPA=="Positive")

# ## Add PC's onto Table for later use
# TAB <- merge( TAB, FT[,c("ID_2",paste("PC",1:3,sep=""))], by.x="FID",by.y="ID_2")

## Take out Patients who left before getting DRUG
RM.exit.id <- as.character( FT.l$ID_2[which(FT.l$IN<=4)] )
RM.exit <- which( TAB$FID %in% RM.exit.id )
TAB <- TAB[ -RM.exit, ]
SAMPS <- unique(as.character(TAB$FID))
SAMPS.short <- unique(as.character(TAB$ID))
N.SAMPS <- length(SAMPS)
SAMPS.list <- list()
SAMPS.list$g <- as.character(FT$ID_2[ FT$GRP=="G" & FT$ID_2%in%SAMPS ])
SAMPS.list$p <- as.character(FT$ID_2[ FT$GRP!="G" & FT$ID_2%in%SAMPS ])
SAMPS.list$ne <- as.character(FT$ID_2[ FT$GRP=="P" & FT$ID_2%in%SAMPS ])
SAMPS.list$ee <- as.character(FT$ID_2[ FT$GRP=="PE" & FT$ID_2%in%SAMPS ])
SAMPS.list$sg <- as.character(FT$ID[ FT$GRP=="G" & FT$ID_2%in%SAMPS ])
SAMPS.list$sp <- as.character(FT$ID[ FT$GRP!="G" & FT$ID_2%in%SAMPS ])
SAMPS.list$sne <- as.character(FT$ID[ FT$GRP=="P" & FT$ID_2%in%SAMPS ])
SAMPS.list$see <- as.character(FT$ID[ FT$GRP=="PE" & FT$ID_2%in%SAMPS ])

RESP_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
## FCT: Create Clin Table for each Response Phenotype
MAKE_CLIN_TAB <- function( RESP_PHENO ) {
	## Pull Response Phenotype of Interest
	TEMP <- TAB[ , c(1:grep("PLAC",colnames(TAB)),which(colnames(TAB)==RESP_PHENO)) ]
	## Remove DAS values that are NA
	RM.na <- which(is.na( TEMP[,RESP_PHENO] ))
	if (length(RM.na) > 0 ) { TEMP <- TEMP[-RM.na, ] }
	## Return Table
	return(TEMP)
}

## Make Clinical Table for each Phenotype
RESP_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
CLIN_TABS <- lapply( RESP_PHENOS, function(x)MAKE_CLIN_TAB(x) )
names(CLIN_TABS) <- RESP_PHENOS

#########################################
## SET BRMS MODEL PARAMETERS ############

## Package Parameters
 # From Tutorial: http://www.r-bloggers.com/r-users-will-now-inevitably-become-bayesians/
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

## Model Parameters
Mod.iter <- 3000
Mod.warm <- 600
Mod.chain <- 3
Mod.fam <- "gaussian"

## Take random sample of Patients
N.Samps <- nrow(FT)
Samps <- sample( SAMPS.short, N.Samps )

## Get Models Set Up
JJ.priors <- JJ.form <- JJ.cor <- list()
if ( !exists("JJ") ) { JJ <- JJ.time <- JJ.samps <- list() }

## Set Priors for each Variable
JJ.priors.list <- list()
# JJ.priors.list$Intercept <- set_prior("normal(5,1)", class="b",coef="Intercept")
JJ.priors.list$Intercept <- set_prior("normal(5,1)", class="Intercept")
JJ.priors.list$DRUG <- set_prior("normal(-1,1)", class="b",coef="DRUG")
JJ.priors.list$WK <- set_prior("normal(0,.1)", class="b",coef="WK")
JJ.priors.list$PLAC <- set_prior("normal(0,1)", class="b",coef="PLAC")
JJ.priors.list$TRT <- set_prior("normal(0,1)", class="b",coef="TRT")
JJ.priors.list$ACPA <- set_prior("normal(0,1)", class="b",coef="ACPA")
JJ.priors.list$DRUG_WK <- set_prior("normal(0,.1)", class="b",coef="DRUG:WK")
JJ.priors.list$DRUG_ACPA <- set_prior("normal(0,1)", class="b",coef="DRUG:ACPA")
JJ.priors.list$ID.sd <- set_prior("cauchy(0,2)", class="sd",group="ID")
JJ.priors.list$cor <- set_prior("lkj(1.5)", class="cor")
JJ.priors.list$ar <- set_prior("normal(0,.75)", class="ar")

## HLA Tables
HLA_TYP <- HLA_TYP.l[ which(HLA_TYP.l$ID %in% SAMPS.short), ]
HLA_AA <- HLA_AA.l[ which(HLA_AA.l$ID %in% SAMPS.short), ]
HLA_HAP <- HLA_HAP.l[ which(HLA_HAP.l$ID %in% SAMPS.short), ]

## Which HLA predictors to use?
Mods.hla.dat <- list()
 # Type
temp.types <- grep( "Pr4_DRB",colnames(HLA_TYP),value=T )
temp.types <- grep( "__",temp.types,value=T,invert=T )	
Mods.hla.dat$Type <- HLA_TYP[, c("ID",temp.types) ]
 # AA
temp.which_pos <- c(11,13,70:74)
Mods.hla.dat$AA <- HLA_AA[, c("ID",grep(paste(paste( "Pr4_DRB1_Pos_",temp.which_pos,sep="" ),collapse="|"),colnames(HLA_AA),value=T )) ]
 # Haps
temp.haps.p117174 <- c("DRE","GRQ","LEA","LRA","PAA","PRA","RRA","SEA","SKR","SRA","SRE","SRL","VEA","VKA","VRA","VRE")
Mods.hla.dat$p117174 <- HLA_HAP[, c("ID",temp.haps.p117174) ]
temp.haps.p11137174 <- c("DFRE","GYRQ","LFEA","LFRA","PRAA","PRRA","RSRA","SGRA","SGRE","SGRL","SSEA","SSKR","SSRA","SSRE","VFRA","VHEA","VHKA","VHRA","VHRE")
Mods.hla.dat$p11137174 <- HLA_HAP[, c("ID",temp.haps.p11137174) ]
temp.haps.pSE <- c("DERAA","DRRAA","DRRAL","DRRGQ","QARAA","QKRAA","QKRGR","QRRAA","QRRAE","RRRAA","RRRAE")
Mods.hla.dat$pSE <- HLA_HAP[, c("ID",temp.haps.pSE) ]

## HLA9: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac
MOD.form <- "DAS ~ (1+DRUG+PLAC|ID)+DRUG*(WK+ACPA+HLA)+PLAC"
MOD.priors <- c(JJ.priors.list$Intercept,
	JJ.priors.list$DRUG,
	JJ.priors.list$WK,
	JJ.priors.list$PLAC,
	JJ.priors.list$DRUG_WK,
	JJ.priors.list$ACPA,
	JJ.priors.list$DRUG_ACPA,
	JJ.priors.list$ID.sd,
	JJ.priors.list$cor )

#############################################################
## PLOT BRMS RESULTS ########################################
#############################################################

## FCT: Plot BRMS Model
PLOT_MOD <- function( model, tag, plot_rand ) {
	write(paste(date(),"Model:",tag), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Collect General Model Info
	summ <- summary(model,waic=F)
	m.obs <- summ$nobs
	m.iter <- summ$iter
	m.warm <- summ$warmup
	m.chain <- summ$chains
	# m.waic <- summ$WAIC
	m.pheno <- as.character(model$formula)[2]
	m.covs <- as.character(model$formula)[3]
	m.form <- paste( m.pheno, "~", m.covs )
	RAND <- length(summ$random)>0

	## Collect Model Outputs
	d.prior <- model$prior
	d.post <- posterior_samples(model)
	f.eff <- summ$fixed
	m.covs.f <- rownames(f.eff)
	c.eff <- summ$cor_pars
	s.eff <- summ$spec_pars
	d.post.f.which <- paste("b_",m.covs.f,sep="")
	d.post.c.which <- rownames(c.eff)
	if ( RAND==T ) {
		r.grp <- summ$group
		r.ngrps <- summ$ngrps
		r.eff <- ranef(model)[[1]]
		r.eff.2 <- summ$random
		m.covs.r <- colnames(r.eff)[1:ncol(r.eff)]
		m.samps <- rownames(r.eff)
		n.samps <- length(m.samps)
		all.eff <- list( f.eff, c.eff, Reduce( rbind, r.eff.2 ), s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","R","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		mod.fit$Eff[grep("cor(",rownames(mod.fit),fixed=T)] <- "RC"
		d.post.r.which <- unlist(lapply( m.covs.r, function(x)paste( "sd_",r.grp,"_",x,sep="") ))
		d.post.r2.which <- unlist(lapply( m.covs.r, function(x)paste( "r_",r.grp,"[",m.samps,",",x,"]",sep="") ))
	}else{
		all.eff <- list( f.eff, c.eff, s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
	}
	print("Model Parsed")
	## Compile List & Save
	# out.m <- list( N.Obs=m.obs, N.iter=m.iter, N.warm=m.warm, N.chain=m.chain,
	# 	Pheno=m.pheno, Covs=m.covs, Formula=m.form, Covs.F=m.covs.f )
	# if ( RAND==T ) { out.r <- list( Grp=r.grp, N.grps=r.ngrps, Effects=r.eff, Effects.2=r.eff.2, Samps=m.samps ) }else{ out.r <- "No Random Effects"}
	# out.compile <- list( Model=out.m, Effects=mod.fit, Posterior=d.post, Prior=d.prior, Random=out.r )
	# save( out.compile, file=paste(PathToPlot,"Rdata.Model_Summary.",tag,".Rdata",sep="") )
	# print("Model Summary Saved")

	###################################
	## FIXED EFFECTS PLOTS ############
	 # Set Color Palette
	COLS.list.heat <- c("firebrick3","chocolate2","gold1","springgreen1","steelblue2","slateblue3")
	COLS <- adjustcolor(COLS.list.2,alpha=.6)
	COLS.eff <- COLS[c(1:4,6)]
	names(COLS.eff) <- c("F","C","R","RC","S")

	## Plot Model Summary
	print("Plotting #1 - Model Summary")
	png( paste(PathToPlot,"ModSumm_",tag,".1-PlotFct.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	par(mfrow=c(nrow(mod.fit),2))
	plot( model, N=nrow(mod.fit) )
	dev.off()

	# ## Plot Chains
	#  # (aka) How to Pull all Sampling Data
	# which_vars <- m.covs.f
	# which_vars <- grep( "sd(",rownames(mod.fit),fixed=T,value=T )
	# png( paste(PathToPlot,"ModSumm_",tag,".1-Chains.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	# par( mfrow=c(length(which_vars),1) )
	# for ( v in which_vars ) {
	# 	v.tag <- grep( substr(gsub("(",paste("_",r.grp,"_",sep=""),v,fixed=T),1,10), model$fit@sim$fnames_oi )
	# 	YLIM <- Reduce(range,lapply(model$fit@sim$samples,function(x)range(x[[paste("b",v,sep="_")]])))
	# 	plot( 0,0,type="n",xlim=c(0,m.iter),ylim=YLIM,xlab="Iteration",ylab=v,main="Chains" )
	# 	abline( h=-100:100,lty=3,col="grey50",lwd=1 )
	# 	# SCRAP <- lapply( seq(YLIM[1],YLIM[2],.025),function(x)abline(h=x,col=adjustcolor(COLS.list.2[4],alpha=2*dnorm(x,5,1)),lwd=2 ))
	# 	for ( c in 1:m.chain ) {
	# 		# TEMP <- model$fit@sim$samples[[c]]$b_Intercept
	# 		points( 1:m.iter, model$fit@sim$samples[[c]][[paste("b",v,sep="_")]], type="l",col=adjustcolor(COLS.list.2[c],alpha=.7),lwd=2 )
	# 	}
	# }
	# dev.off()
	
	## Fixed Effect Sizes
	print("Plotting #2 - Effect Sizes")
	YLIM <- extendrange(mod.fit[,"Estimate"], f=.2)
	png( paste(PathToPlot,"ModSumm_",tag,".2-EffSize.png",sep=""),height=1000,width=400+100*nrow(mod.fit),pointsize=26)
	par(mar=c( 7,5,5,3 ))
	# TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n",ylim=YLIM,ylab="Effect Size",main="Effect Size Estimates" )
	TEMP <- 1:nrow(mod.fit)
	plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=range(TEMP),ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
	axis( 1, at=TEMP,label=rownames(mod.fit), las=2 )
	axis( 2, at=seq(-10,10,2), las=2 )
	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
	abline( h=0, lty=1,col="black",lwd=1 )
	 # Plot Prior Distributions
	for ( v in 1:nrow(mod.fit) ) {
		var <- rownames(mod.fit)[v]
		if ( var %in% m.covs.f ) {
			if ( var=="Intercept" ) {
				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
			}else{
				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
			}
			if ( length(temp.priors)==2 ) {
				names(temp.priors) <- c("Mean","SD")
				vioplot( rnorm(1e5,temp.priors["Mean"],temp.priors["SD"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
			}
		}
		# if ( var %in% paste("sd(",m.covs.r,")",sep="") ) {
		# 	temp.priors <- as.numeric(strsplit( gsub(")","",gsub("cauchy(","",d.prior[which(d.prior$class=="sd"&d.prior$group==r.grp)[1],"prior"], fixed=T),fixed=T),"," )[[1]])
		# 	names(temp.priors) <- c("Loc","Scale")
		# 	vioplot( rnorm(1e5,temp.priors["Loc"],temp.priors["Scale"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )
		# }
	}
	 # Plot Posterior Distributions
	for ( v in 1:nrow(mod.fit) ) {
		var <- rownames(mod.fit)[v]
		var.tag <- var
		if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
		if ( grepl("sd(",var,fixed=T) ) { var.tag <- gsub(")","", gsub("sd(",paste("sd_",r.grp,"_",sep=""),var,fixed=T),fixed=T) }
		if ( grepl("cor(",var,fixed=T) ) { var.tag <- gsub(",","_", gsub(")","", gsub("cor(",paste("cor_",r.grp,"_",sep=""),var,fixed=T),fixed=T),fixed=T) }
		if ( var==paste("sigma(",m.pheno,")",sep="") ) { var.tag <- paste("sigma",m.pheno,sep="_") }
		if ( var.tag %in% colnames(d.post) ) {
			vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.eff[mod.fit[v,"Eff"]],add=T,drawRect=F )
		}
	}
	# TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n", add=T)
	arrows( TEMP,mod.fit[,"l.95..CI"],TEMP,mod.fit[,"u.95..CI"],lwd=3,length=0 )
	arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],lwd=3,length=0 )
	# legend("topright",fill=COLS.eff,border=NA,legend=names(COLS.eff),ncol=length(COLS.eff),title="Effect Type",bg="white")
	legend("topright",fill=c(adjustcolor("black",alpha=.2),COLS.eff),border=NA,legend=c("Prior Distribution",names(COLS.eff)),ncol=(1+length(COLS.eff))/2,title="Effect Type",bg="white")
	dev.off()

	## RANDOM EFFECTS PLOTS ###########
	if ( RAND==T & plot_rand!=F ) {
		print("Plotting Random Effects")
		n.covs.r <- length(m.covs.r)

		## Plot Distributions & Shrinkage of Random Effects
		print("Plotting #R1 - Distributions")
		PLOT_SHRINK <- function( var, which_plots ) {
			if ( "TRT"%in%m.covs.f ) { TRT_PBO <- "TRT"
			}else{ TRT_PBO <- "PLAC" }
			## Pull Data
			brm.var <- r.eff[,var] + f.eff[var,"Estimate"]
			if ( var %in% c("Intercept","DRUG") ) {
				if ( var=="Intercept" ) {
					mn.var <- FT[,"DAS_BL_MN"]
					brm.fx <- f.eff[var,"Estimate"]
				}
				if ( var=="DRUG" ) {
					mn.var <- FT[,"DEL_MNe_MN"]
					brm.fx <- f.eff[var,"Estimate"] + f.eff[TRT_PBO,"Estimate"]
					if ( TRT_PBO%in%m.covs.f ) { brm.var <- brm.var + f.eff[TRT_PBO,"Estimate"] }
				}
				mn.mod <- T
				names(mn.var) <- as.character(FT$ID)
			}else{
				mn.var <- brm.var
				mn.mod <- F
				brm.fx <- f.eff[var,"Estimate"]
			}
			mg.var <- merge( data.frame(brm.var),data.frame(mn.var), by="row.names" )
			mg.var <- merge( FT[,c("ID","GRP")], mg.var, by.x="ID",by.y="Row.names" )
			mg.var.ord <- mg.var[ order(mg.var$mn.var), ]
			mn.fx <- mean(mg.var$mn.var)

			## Plot 1: Distributions
			if ( 1 %in% which_plots ) {
				XLIM <- range( c(mg.var$mn.var,mg.var$brm.var) )
				X_BIN <- .25
				BRKS <- seq( floor(XLIM[1]), ceiling(XLIM[2])+X_BIN, X_BIN)
				YLIM <- c( 0, Reduce( max, lapply(list(mg.var$brm.var,mg.var$mn.var),function(x)hist(x,breaks=BRKS,plot=F)$counts) ) )
				hist( mg.var$brm.var, col=COLS[1],border=NA,breaks=BRKS, xlim=XLIM,ylim=YLIM, main=paste(var,"- Mean vs BRMS"),xlab=var,ylab="# Patients" )
				abline( v=brm.fx, col=COLS.list.2[1],lwd=4,lty=1 )
				if ( mn.mod==T ) {
					hist( mg.var$mn.var, col=COLS[2],border=NA,breaks=BRKS, add=T )
					abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=1 )
					legend( "topleft", legend=c("BRMS","Mean"),title="Model",fill=COLS[1:2],border=NA)	
				}			
			}
	
			## Plot 2: Boxplot by Arm
			if ( 2 %in% which_plots ) {
				ARMS <- as.character( unique( mg.var$GRP ) )
				XLIM <- c(0,7)
				YLIM <- range( c(mg.var$mn.var,mg.var$brm.var) )
				plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",xlab="ARM",ylab=var,main=paste(var,"Distribution by Arm") )
				abline( h=-10:10,lty=3,col="grey50",lwd=1 )
				axis( 1,at=seq(1.5,7,2),label=ARMS )
				SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$brm.var[mg.var$GRP==ARMS[x]],at=2*x-1,col=COLS[1],border=NA,add=T) )
				abline( h=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
				if ( mn.mod==T ) {
					SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$mn.var[mg.var$GRP==ARMS[x]],at=2*x,col=COLS[2],border=NA,add=T) )
					abline( h=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
					legend( "topleft", legend=c("BRMS","Mean"),title="Model",fill=COLS[1:2],border=NA)	
					arrows( c(1,3,5),rep(YLIM[1],3),c(2,4,6),rep(YLIM[1],3), length=0,lwd=4,col=COLS.list.2[3:5] )
				}		
			}

			## Plot 3: Shrinkage by Person
			if ( 3 %in% which_plots ) {
				XLIM <- range( c(mg.var.ord$mn.var,mg.var.ord$brm.var) )
				X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5)
				YLIM <- c( 1, n.samps )
				plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=var,ylab="Patient",main=paste(var,"- Estimate Shrinkage") )
				abline( v=X_LNS, lty=3,col="grey50",lwd=1 )
				points( mg.var.ord$brm.var,1:nrow(mg.var.ord), pch=16,cex=1,col=COLS[1] )
				if ( mn.mod==T ) {
					points( mg.var.ord$mn.var,1:nrow(mg.var.ord), pch=16,cex=1,col=COLS[2] )
					arrows( mg.var.ord$mn.var,1:nrow(mg.var.ord),mg.var.ord$brm.var,1:nrow(mg.var.ord), angle=30,length=.2,lwd=2,col=COLS.list.2[3:5][factor(mg.var.ord$GRP)] )
					abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
				}
				abline( v=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
				legend("topleft",col=COLS.list.2[1:5],legend=c("Mod: BRMS","Mod: Mean","Arm: G","Arm: P","Arm: PE"),pch=c(16,16,NA,NA,NA),lty=c(0,0,1,1,1),lwd=2 )
			}

			## Plot 4: MN vs BRMS Model Estimates
			if ( 4 %in% which_plots ) {
				XLIM <- range( mg.var$mn.var )
				YLIM <- range( mg.var$brm.var )
				X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5 )
				Y_LNS <- seq( floor(YLIM[1]-5), ceiling(YLIM[2]+5), .5 )
				plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Mean Estimate",ylab="BRMS Estimate",main=var )
				abline( v=X_LNS,h=Y_LNS, lty=3,col="grey50",lwd=1 )
				abline( 0,1 )
				abline( h=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
				abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
				points( mg.var$brm.var ~ mg.var$mn.var, data=mg.var, pch=16,cex=1,col=COLS[3:5][factor(GRP)] )
				legend( "topleft",pch=16,col=COLS[3:5],legend=levels(factor(mg.var$GRP)))
			}
		}
		 # Plot it
		for ( var in m.covs.r ) {
			if ( var %in% c("Intercept","DRUG") ) {
				png( paste(PathToPlot,"ModSumm_",tag,".R1-Dist.",var,".png",sep=""),height=2000,width=2000,pointsize=32)
				layout( matrix(c(1:4,4,4),byrow=F,ncol=2) )
				PLOT_SHRINK( var, c(1,2,4) )
				PLOT_SHRINK( var, c(3) )
				dev.off()
			}else{
				png( paste(PathToPlot,"ModSumm_",tag,".R1-Dist.",var,".png",sep=""),height=2000,width=2000,pointsize=32)
				layout( matrix(c(1:2,3,3),byrow=F,ncol=2) )
				PLOT_SHRINK( var, c(1,2) )
				PLOT_SHRINK( var, c(3) )
				dev.off()
			}

		}
	
		## Plot Pairs of Random Effects
		panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...) {
			usr <- par("usr"); on.exit(par(usr))
			par(usr = c(0, 1, 0, 1))
			r <- abs(cor(x, y))
			txt <- format(c(r, 0.123456789), digits = digits)[1]
			txt <- paste0(prefix, txt)
			if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
			text(0.5, 0.5, txt, cex = cex.cor * r)
		}
		print("Plotting #R1 - Pairs")
		if ( n.covs.r > 1 ) {
			temp.reff <- t( f.eff[m.covs.r,"Estimate"] + t(r.eff) ) 
			# temp.grp <- FT$GRP[ match(rownames(temp.reff),FT$ID) ] # col=adjustcolor(COLS.list.2[factor(temp.grp)],alpha=.4)
			temp.reff <- merge( FT[,c("ID","GRP")], temp.reff, by.x="ID",by.y="row.names" )
			png( paste(PathToPlot,"ModSumm_",tag,".R1-PairRand.png",sep=""),height=400*n.covs.r,width=400*n.covs.r,pointsize=24+2*n.covs.r)
			# pairs( temp.reff, pch=16,col=adjustcolor(COLS.eff["RC"],alpha=.5),cex=2,upper.panel=panel.cor )
			pairs( temp.reff[,m.covs.r], pch=16,col=COLS[3:5][factor(temp.reff$GRP)],cex=2,upper.panel=panel.cor )
			dev.off()
		}

		## Boxplot of Individual Posterior Distributions
		print("Plotting #R2 - Boxplot")
		# png( paste(PathToPlot,"ModSumm_",tag,".R2-BoxRand.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
		# par(mfrow=c(n.covs.r,1))
		# par(mar=c(6,4,4,1))
		# for ( z in 1:n.covs.r ) {
		# 	z.name <- m.covs.r[z]
		# 	z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
		# 	z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
		# 	z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
		# 	if ( z>0 ) { z.ord <- order( colMeans(z.temp) ) }
		# 	YLIM <- range( z.temp )
		# 	boxplot( z.temp[,z.ord], col=COLS.eff["R"],main=paste("Random",z.name,"by Patient"),ylim=YLIM,ylab=z.name,names=m.samps[z.ord],las=2,pch=16 )
		# 	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
		# 	abline( h=0, lty=1,col="black",lwd=1 )
		# 	boxplot( z.temp[,z.ord], col=COLS.eff["R"],names=m.samps[z.ord],las=2,pch=16,add=T )
		# 	if ( any(z.name %in% c("DRUG","PLAC","TRT")) ) {
		# 		abline(h=c(0,-1),lwd=2,col=COLS.list.2[6:7] )
		# 		text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.5 )
		# 		text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[7],cex=.5 )
		# 	}else{
		# 		legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
		# 	}			
		# }
		# dev.off()
		 # Custom (quicker) "Boxplot" of Individual Posterior Distributions
		png( paste(PathToPlot,"ModSumm_",tag,".R2-BoxRand.2.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
		par(mfrow=c(n.covs.r,1))
		par(mar=c(6,4,4,1))
		for ( z in 1:n.covs.r ) {
			# Pull Data
			z.name <- m.covs.r[z]
			z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
			z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
			z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
			# z.samps.temp <- gsub("r_ID[","", sapply(strsplit(colnames(z.temp),","),"[",1) ,fixed=T)
			z.grp <- FT$GRP[ match(m.samps,as.character(FT$ID)) ]
			which_quants <- c(.025,.25,.5,.75,.975)
			z.quants <- apply( z.temp, 2, function(x)quantile(x,which_quants) )
			if ( z>0 ) { z.ord <- order( colMeans(z.temp) ) }
			# Parameters
			XLIM <- c( 1,ncol(z.temp) )
			YLIM <- range( z.temp )
			# Build Plot
			plot( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp], main=paste("Random",z.name,"by Patient"),xlab="",ylab=z.name,xaxt="n",xlim=XLIM,ylim=YLIM )
			axis( 1, at=1:ncol(z.temp),label=m.samps[z.ord],las=2,cex=.6 )
			abline( h=-10:10,v=1:ncol(z.temp), lty=3,col="grey50",lwd=1 )
			abline( h=0, lty=1,col="black",lwd=1 )
			# Text/Legend
			if ( any(z.name %in% c("DRUG","PLAC","TRT")) ) {
				abline(h=c(0,-1),lwd=2,col=COLS.list.2[6:7] )
				text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.7,srt=90 )
				text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[7],cex=.7,srt=90 )
			}else{
				# legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
				legend( "topleft",pch=c(rep(16,3),16,1),cex=1.5,col=c(COLS.list.2[3:5],rep("black",2)),legend=c("Arm: G","Arm: P","Arm: PE","Median","CI: 95%"),lwd=3 )
			}
			# Points/Arrows
			arrows( 1:ncol(z.temp),z.quants["25%",z.ord],1:ncol(z.temp),z.quants["75%",z.ord],col=COLS[1],lwd=5,code=3,angle=90,length=.1 )
			points( rep(1:ncol(z.temp),2),c(z.quants["2.5%",z.ord],z.quants["97.5%",z.ord]),pch=1,lwd=3,col=COLS.list.2[3:5][z.grp],cex=1.5 )
			points( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp],cex=1.5 )
		}
		dev.off()

		# ## Violin Plot: Prior vs Posterior (Compiled)
		# print("Plotting #R2 - Violin Plot")
		# png( paste(PathToPlot,"ModSumm_",tag,".R2b-ViolRand.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
		# par(mfrow=c(n.covs.r,1))
		# par(mar=c(6,4,4,1))
		# for ( z in 1:n.covs.r ) {
		# 	z.name <- m.covs.r[z]
		# 	z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
		# 	z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
		# 	z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
		# 	if ( z==1 ) { z.ord <- order( colMeans(z.temp) ) }
		# 	XLIM <- c(1,ncol(z.temp))
		# 	YLIM <- extendrange(z.temp,f=.2)
		# 	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="",main="Posterior Population Estimate (by Patient)",xaxt="n" )
		# 	axis( 1, at=1:ncol(z.temp), label=m.samps[z.ord], las=2)
		# 	abline(h=-10:10,lty=3,col="grey50",lwd=1 )
		# 	abline(h=0)
		# 	# Get Prior
		# 	if ( z.name=="Intercept" ) { temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
		# 	}else{ temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==z.name,"prior"], fixed=T),fixed=T),"," )[[1]]) }
		# 	names(temp.priors) <- c("Mean","SD")
		# 	temp.prior.distr <- rnorm(10000,temp.priors["Mean"],temp.priors["SD"])
		# 	# Plot Prior/Posterior
		# 	SCRAP <- lapply( 1:ncol(z.temp), function(x)vioplot(temp.prior.distr,at=x, col=adjustcolor("black",alpha=.2),border=NA,add=T ))
		# 	SCRAP <- lapply( 1:ncol(z.temp), function(x)vioplot(z.temp[,z.ord[x]],at=x, col=COLS.eff["R"],add=T ))
		# 	if ( any(z.name %in% c("DRUG","PLAC")) ) {
		# 		abline(h=c(0,-1),lwd=2,col=COLS.list.2[5:6] )
		# 		text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[5],cex=.5 )
		# 		text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.5 )
		# 	}else{
		# 		legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
		# 	}			
		# }
		# dev.off()

		## Plot Posterior Probabilities for Random Effects
		print("Plotting #R3 - Heatmaps of Posteriors")
		 # Plotting Parameters
		COLS.heat <- colorRampPalette(c(COLS.list.heat[1:3],"white",COLS.list.heat[4:6]))(100)
		BRKS.heat <- seq(0,1,length.out=101)
		n.iter.tot <- m.chain * ( m.iter - m.warm )
		post.prob.brks <- list( Intercept=seq(3.5,7.5,.125),
			DRUG=seq(-3,0,.125),
			PLAC=seq(-2,1,.125),
			TRT=seq(-2,1,.125) )
		 # Calculate/Plot Posteriors at Various Thresholds
		post.prob <- list()
		for ( z in 1:n.covs.r ) {
			z.name <- m.covs.r[z]
			z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
			z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
			z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
			z.grp <- FT$GRP[ match(m.samps,as.character(FT$ID)) ]
			# Calculate Posteriors
			post.prob[[z.name]] <- Reduce( cbind, lapply( post.prob.brks[[z.name]],function(x)apply(z.temp,2,function(y)length(which(y<x))) ) ) / n.iter.tot
			colnames(post.prob[[z.name]]) <- paste("LT",post.prob.brks[[z.name]],sep="_")
			rownames(post.prob[[z.name]]) <- colnames(z.temp)
			if ( z.name=="Intercept" ) { post.prob[[z.name]] <- 1 - post.prob[[z.name]] ; colnames(post.prob[[z.name]]) <- gsub("LT","GT",colnames(post.prob[[z.name]])) }
			# Plot Heatmap
			png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.",z.name,".png",sep=""),height=1200,width=2000+10*n.samps,pointsize=28)
			heatmap.2( t(post.prob[[z.name]]), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5+n.samps/100),main=z.name,ColSideColors=COLS.list.2[3:5][z.grp] )
			dev.off()
		}
		 # Heatmap of All Random Effect Posteriors
		z.grp <- FT$GRP[ match(m.samps,as.character(FT$ID)) ]
		post.prob.temp <- Reduce( cbind, post.prob )
		png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.All.png",sep=""),height=800+400*n.covs.r,width=2000+10*n.samps,pointsize=28)
		heatmap.2( t(post.prob.temp), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5+n.samps/100),rowsep=cumsum(lapply(post.prob,ncol)),ColSideColors=COLS.list.2[3:5][z.grp] )
		dev.off()

		# png( paste(PathToPlot,"ModSumm_",tag,".R3-PostScatter.png",sep=""),height=800+400*n.covs.r,width=2000+10*n.samps,pointsize=28)
		# par(mfrow=c(n.covs.r,1))
		# for ( z in 1:n.covs.r ) {
		# 	z.name <- m.covs.r[z]
		# 	temp.split <- as.numeric(sapply(strsplit(colnames(post.prob[[z.name]]),"_"),"[",2))
		# 	which_cols <- which( temp.split %in% unique(round(2*temp.split)/2) )
		# 	z.ord <- order( post.prob[[z.name]][,which.max(apply(post.prob[[z.name]],2,sd))] )
		# 	COLS.temp <- colorRampPalette(COLS.list.ord)(length(which_cols))
		# 	plot( 0,0,type="n", xlim=c(1,n.samps),ylim=c(0,1) )
		# 	lapply( 1:length(which_cols),function(x)points(1:n.samps,post.prob[[z.name]][z.ord,which_cols[x]],col=adjustcolor(COLS.temp[x],alpha=.7),pch=16,type="o",lwd=2) )
		# }
		# dev.off()
		
		## Plot Distribution (across patients) of Standard Deviations (across iteration) of Random Effects
		# TEST <- lapply(m.covs.r,function(r)unlist(lapply( m.samps, function(s) sd(d.post[,paste("r_ID[",s,",",r,"]",sep="")]) )) )
		# names(TEST) <- m.covs.r

		## Plot DRUG vs Intercept for each Simulation/Patient
		# if ( all(c("DRUG","Intercept")%in%m.covs.r) ) {
		# 	print("Plotting #R4 - Drug v Int per Patient")
		# 	x.name <- "Intercept"
		# 	x.col.1 <- d.post.f.which[ d.post.f.which==paste("b",x.name,sep="_") ] # grep(x.name,d.post.f.which,value=T)
		# 	x.col.2 <- grep( paste(",",x.name,"]",sep=""),d.post.r2.which,value=T)
		# 	x.temp <- d.post[,x.col.1] + d.post[,x.col.2]
		# 	y.name <- "DRUG"
		# 	y.col.1 <- d.post.f.which[ d.post.f.which==paste("b",y.name,sep="_") ] # grep(y.name,d.post.f.which,value=T)
		# 	y.col.2 <- grep( paste(",",y.name,"]",sep=""),d.post.r2.which,value=T)
		# 	y.temp <- d.post[,y.col.1] + d.post[,y.col.2]
		# 	XLIM <- range(x.temp)
		# 	YLIM <- range(y.temp)
		# 	COLS.ind <- colorRampPalette(COLS.list)(ncol(x.temp))
		# 	COLS.ind.2 <- adjustcolor(COLS.ind,alpha=.01)
		# 	png( paste(PathToPlot,"ModSumm_",tag,".R4-Sims.DRvINT.png",sep=""),height=1000,width=1600,pointsize=30)
		# 	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="Drug",main="Individual Simulation Estimates")
		# 	abline(h=-10:10,v=-10:10,lty=3,col="grey50",lwd=1 )
		# 	abline(h=0)
		# 	SCRAP <- lapply( 1:ncol(x.temp),function(i)points(x.temp[,i],y.temp[,i],col=COLS.ind.2[i],pch=16) )
		# 	mod.temp <- lm(unlist(y.temp)~unlist(x.temp))
		# 	abline( mod.temp )
		# 	points( colMeans(x.temp),colMeans(y.temp),pch="+",col=adjustcolor(COLS.ind,alpha=.7),cex=2,lwd=2 )
		# 	dev.off()
		# }

		## Plot a Few Individual Patients' Profiles
		if ( all(c("DRUG","Intercept")%in%m.covs.r) & any(c("PLAC","TRT")%in%m.covs.r) ) {
		# if ( all(c("DRUG","PLAC","Intercept")%in%m.covs.r) ) {
			print("Plotting Individual Patients")
			## PLAC or TRT?
			PBO <- grep("PLAC|TRT",m.covs.r,value=T )
			## Which Patients?
			 # How Many?
			plot.rows <- 4
			plot.cols <- 5
			z.samps.n <- plot.rows * plot.cols
			 # Specify Samples
			# z.samps <- sample(m.samps, z.samps.n )
			z.samps.temp.ord.int <- order( colMeans(d.post[,grep( paste(",Intercept]",sep=""),d.post.r2.which,value=T)]) )
			z.samps.temp.ord.dr <- order( colMeans(d.post[,grep( paste(",DRUG]",sep=""),d.post.r2.which,value=T)]) )
			z.samps.temp <- m.samps[ c( head(z.samps.temp.ord.int,plot.rows),tail(z.samps.temp.ord.int,plot.rows), head(z.samps.temp.ord.dr,plot.rows),tail(z.samps.temp.ord.dr,plot.rows) ) ]
			z.samps <- sample( setdiff(m.samps,z.samps.temp), z.samps.n-length(z.samps.temp) )
			z.samps <- c( head(z.samps.temp,2*plot.rows), z.samps, tail(z.samps.temp,2*plot.rows) )
			# Pull out Sample Data
			z.data <- model$data[model$data[,r.grp]%in%z.samps,]
			z.temp <- predict( model, newdata=z.data )
			z.pred <- cbind( z.data, z.temp )
			z.ranef <- t(f.eff[m.covs.r,"Estimate"] + t(r.eff[z.samps,]))
			if ( PBO=="TRT" ) { z.ranef[,"DRUG"] <- z.ranef[,"DRUG"] + z.ranef[,"TRT"]}

			# print("Plotting #R5 - Predicted Values")
			## FCT: Plot Patient Profile
			sample <- z.samps[1]
			PLOT_IND <- function( sample, which_plots ) {
				## Plot 1: Real vs Predicted Values
				if ( 1 %in% which_plots ) {
					COLS.tru <- COLS.list.2[c(6,2,1)]
					COLS.conf <- COLS.list.2[4]
					sub.tab <- z.pred[ z.pred[,r.grp]==sample, ]
					plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main=sample )
					abline( h=0:10,lty=3,col="grey50",lwd=1)
					points( Estimate ~ WK, data=sub.tab,type="l",lwd=4 )
					polygon( c(sub.tab$WK,rev(sub.tab$WK)), c(sub.tab[,"2.5%ile"],rev(sub.tab[,"97.5%ile"])), col=adjustcolor(COLS.conf,alpha=.2),border=NA )
					points( DAS ~ WK, data=sub.tab,type="p",lwd=3,cex=1.5,pch=c(4,16)[factor(DRUG)],col=COLS.tru[1+rowSums(sub.tab[,c(PBO,"DRUG")])] )
					# pre.post <- cumsum(unlist( FT[FT$ID==sample,c("DAS_BL_MN","DEL_MNe_MN")] ))
					# arrows( c(0,24),pre.post,c(24,100),pre.post, col=COLS.pred[5],lwd=4,lty=2,length=0 )
					pre.post.2 <- z.ranef[sample,c(1,3,2)] + c(0,rep(z.ranef[sample,1],2))
					arrows( c(-5,2,24),pre.post.2,c(100,24,100),pre.post.2, col=COLS.tru,lwd=6,lty=2,length=0 )
					text( rep(100,3),c(9,8.5,8), paste(colnames(z.ranef),round(z.ranef[sample,],2),sep="=")[c(1,3,2)], pos=2,col=COLS.tru )
				}

				## Plot 2: Scatter of Individual Monte Carlo Estimates
				if ( 2 %in% which_plots ) {
					x.name <- "Intercept"
					x.col.1 <- d.post.f.which[ d.post.f.which==paste("b",x.name,sep="_") ] # grep(x.name,d.post.f.which,value=T)
					x.col.2 <- grep( paste(sample,",",x.name,"]",sep=""),d.post.r2.which,value=T)
					x.temp <- d.post[,x.col.1] + d.post[,x.col.2]
					y.name <- "DRUG"
					y.col.1 <- d.post.f.which[ d.post.f.which==paste("b",y.name,sep="_") ] # grep(y.name,d.post.f.which,value=T)
					y.col.2 <- grep( paste(sample,",",y.name,"]",sep=""),d.post.r2.which,value=T)
					y.temp <- d.post[,y.col.1] + d.post[,y.col.2]
					w.name <- PBO
					w.col.1 <- d.post.f.which[ d.post.f.which==paste("b",w.name,sep="_") ] # grep(w.name,d.post.f.which,value=T)
					w.col.2 <- grep( paste(sample,",",w.name,"]",sep=""),d.post.r2.which,value=T)
					w.temp <- d.post[,w.col.1] + d.post[,w.col.2]
					if ( PBO=="TRT" ) { y.temp <- w.temp + y.temp }
					XLIM <- c(1,8) # range(x.temp)
					YLIM <- c(-3,2) # range(y.temp)
					# Plot it
					plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="Drug",main="Individual Simulation Estimates")
					abline(h=-10:10,v=-10:10,lty=3,col="grey50",lwd=1 )
					abline(h=0,v=5)
					SCRAP <- points(x.temp,y.temp,col=adjustcolor(COLS.tru[3],.03),pch=16)
					SCRAP <- points(x.temp,w.temp,col=adjustcolor(COLS.tru[2],.03),pch=16)
					points( mean(x.temp),mean(y.temp),pch="+",col=COLS.tru[3],cex=2,lwd=2 )
					points( mean(x.temp),mean(w.temp),pch="+",col=COLS.tru[2],cex=2,lwd=2 )
					legend( "topleft",col=COLS.tru[2:3],pch=16,legend=c("TRT","DRUG+TRT"),title="Effect Sizes" )
				}					
			}

			for ( samp in z.samps ) {
			# for ( samp in z.samps[seq(1,20,4)] ) {
				png( paste(PathToPlot,"ModSumm_",tag,".Ri.",samp,"-Pred.png",sep=""),height=1000,width=2000,pointsize=30)
				par(mfrow=c(1,2))
				SCRAP <- PLOT_IND( samp, 1:2 )
				dev.off()
			}

		}
	} ## End Random Effects Plots
} ## Close "PLOT_MOD" Function

## FCT: Modify Priors & Formula for various HLA Predictiors
Modify_HLA_Inputs <- function( mod.name ) {
	## Pull Predictor Info
	hla.data <- Mods.hla.dat[[mod.name]]
	hla.predictors <- colnames(hla.data)[-1]
	## Modify Formula
	hla.form.preds <- paste(hla.predictors,collapse="+")
	OUT.form <- gsub( "HLA",hla.form.preds,MOD.form )
	## Modify Priors
	OUT.priors <- MOD.priors
	for ( hla.pred in hla.predictors ) {
		OUT.priors <- rbind( OUT.priors, JJ.priors.list$HLA, JJ.priors.list$DRUG_HLA )
		OUT.priors[which(OUT.priors$coef=="HLA"),"coef"] <- hla.pred
		OUT.priors[which(OUT.priors$coef=="DRUG:HLA"),"coef"] <- paste("DRUG:",hla.pred,sep="")
	}
	## Return Modified Formula/Priors
	OUT <- list( Formula=OUT.form, Priors=OUT.priors, Data=hla.data )
	return(OUT)
}

# ## FCT: Plot Certain Fixed Effects from HLA Models
# PLOT_FIXED <- function( mod.fit.table, eff.tag ) {
# 	mod.fit.table.names <- gsub("_",":",gsub( "Pr4_DRB1_","DRB1*",rownames(mod.fit.table) ))
# 	YLIM <- extendrange(mod.fit.table[,"Estimate"], f=.2)
# 	YLIM <- YLIM - c(.15*diff(YLIM), 0)
# 	TEMP <- 1:nrow(mod.fit.table)
# 	XLIM <- range(TEMP)+c(-.5,.5)
# 	png( paste(PathToPlot,"HLA-",mod,".2-EffSize.",eff.tag,".png",sep=""), height=1000,width=600+45*nrow(mod.fit.table),pointsize=26)
# 	par(mar=c( 8,5,5,3 ))
# 	plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=XLIM,ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
# 	axis( 1, at=TEMP,label=mod.fit.table.names, las=2 )
# 	axis( 2, at=seq(-10,10,2), las=2 )
# 	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
# 	abline( h=0, lty=1,col="black",lwd=1 )
# 	 # Plot Prior Distributions
# 	for ( v in 1:nrow(mod.fit.table) ) {
# 		var <- rownames(mod.fit.table)[v]
# 		if ( var %in% m.covs.f ) {
# 			if ( var=="Intercept" ) {
# 				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
# 			}else{
# 				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
# 			}
# 			if ( length(temp.priors)==2 ) {
# 				names(temp.priors) <- c("Mean","SD")
# 				vioplot( rnorm(1e4,temp.priors["Mean"],temp.priors["SD"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
# 			}
# 		}
# 	}
# 	 # Plot Posterior Distributions
# 	COLS.temp.list <- c( COLS.ph[2:1],COLS.list.2[4] )
# 	COLS.temp <- adjustcolor(COLS.temp.list[1:2][factor(grepl("DRUG",mod.fit.table.names))],.8)
# 	COLS.temp[grep(paste(hla.preds,collapse="|"),rownames(mod.fit.table),invert=T)] <- adjustcolor(COLS.temp.list[3],.8)
# 	for ( v in 1:nrow(mod.fit.table) ) {
# 		var <- rownames(mod.fit.table)[v]
# 		var.tag <- var
# 		if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
# 		if ( var.tag %in% colnames(d.post) ) {
# 			vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.temp[v],add=T,drawRect=F )
# 		}
# 	}
# 	arrows( TEMP,mod.fit.table[,"l.95..CI"],TEMP,mod.fit.table[,"u.95..CI"],lwd=3,length=0 )
# 	arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],lwd=3,length=0 )
# 	legend("topright",fill=c(adjustcolor("black",.2),adjustcolor(COLS.temp.list,.8)),border=NA,legend=c("Prior Distribution","HLA Disease Severity","HLA Drug Response","Clinical Covariate"),ncol=2,title="Effect Type",bg="white")
# 	temp.text <- abs(.5 - round(d.post.prob[rownames(mod.fit.table),"Pr_l0"],2)) + .5
# 	text( TEMP, YLIM[1], temp.text, cex=.9 )
# 	text( XLIM[2], YLIM[1]+.05*diff(YLIM), "Posterior Prob. of Effect", cex=.9, pos=2 )
# 	dev.off()
# }

## FCT: Plot Certain Fixed Effects from HLA Models
PLOT_FIXED_HLA <- function( mod.fit.table.unsrt, eff.tag ) {
	temp.order <- order(grep("DRUG",rownames(mod.fit.table.unsrt),invert=T,value=T))
	temp.order <- temp.order + rep(c(0,max(temp.order)),each=length(temp.order) )
	mod.fit.table <- mod.fit.table.unsrt[ temp.order, ]
	mod.fit.table.names <- gsub("_",":",gsub( "Pr4_DRB1_","DRB1*",rownames(mod.fit.table) ))
	mod.fit.names <- paste(rep(c("","DRUG:"),each=length(mod.fit.table.names)), mod.fit.table.names[temp.order],sep="" )

	YLIM <- extendrange(mod.fit.table[,"Estimate"], f=.2)
	YLIM <- YLIM - c(.1*diff(YLIM), 0)
	YLIM <- c(-2,2)
	TEMP <- rep( 1:(nrow(mod.fit.table)/2), 2 )
	XLIM <- range(TEMP)+c(-.5,.5)
	# png( paste(PathToPlot,"HLA-",m.tag,".2-EffSize.",eff.tag,".png",sep=""), height=1000,width=600+45*nrow(mod.fit.table),pointsize=26)
	if ( any(grepl("DRB",mod.fit.table.names)) ) {
		png( paste(PathToPlot,"HLA-",mod,".2-EffSize.",eff.tag,".png",sep=""), height=1000,width=2000,pointsize=30)
		par(mar=c( 6,4,4,3 ))
	}else{
		png( paste(PathToPlot,"HLA-",mod,".2-EffSize.",eff.tag,".png",sep=""), height=1000,width=2000,pointsize=30)
		par(mar=c( 4,4,4,3 ))
	}
	plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=XLIM,ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
	axis( 1, at=TEMP[1:max(TEMP)],label=mod.fit.table.names[1:max(TEMP)], las=2 )
	axis( 2, at=seq(-2,2,.5), las=2 )
	abline( h=seq(-10,10,.5), lty=3,col="grey50",lwd=1 )
	abline( h=0, lty=1,col="black",lwd=1 )
	 # Plot Prior Distributions
	HLA.prior <- rnorm(1e4,0,.3)
	for ( v in 1:max(TEMP) ) {
		vioplot( HLA.prior, at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
	}
	 # Plot Posterior Distributions
	COLS.temp.list <- c( COLS.ph[2:1] )
	COLS.temp <- rep( adjustcolor(COLS.temp.list[1:2],.4), each=max(TEMP))
	COLS.temp.2 <- rep( COLS.temp.list[1:2], each=max(TEMP))
	for ( v in 1:nrow(mod.fit.table) ) {
		var <- rownames(mod.fit.table)[v]
		var.tag <- var
		if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
		if ( var.tag %in% colnames(d.post) ) {
			vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.temp[v],add=T,drawRect=F )
		}
	}
	arrows( TEMP,mod.fit.table[,"l.95..CI"],TEMP,mod.fit.table[,"u.95..CI"],lwd=3,length=0 )
	arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],lwd=3,length=0 )
	legend("topright",fill=c(adjustcolor("black",.2),adjustcolor(COLS.temp.list,.8)),border=NA,legend=c("Prior Distribution","HLA Disease Severity","HLA Drug Response"),title="Effect Type",bg="white",ncol=3)
	temp.text <- abs(.5 - round(d.post.prob[rownames(mod.fit.table),"Pr_l0"],2)) + .5
	text( TEMP, YLIM[1]+rep(c(.05,.25),each=max(TEMP)), substr(temp.text,2,4), col=COLS.temp.2, cex=.8,srt=0 )
	text( XLIM[2], YLIM[1]+.4, "Posterior Prob. of Effect", cex=.9, pos=2 )
	dev.off()
}

## Compile & Plot Results & Create Tables
tag <- "HLA"
HLA.mod.fits <- HLA.mod.MNvBRMS <- list()
Mod.Names <- sapply(strsplit( names(JJ[[tag]]),"_"),"[",2 )[-2]
which_plots <- 1:3
which_plots <- 1:2
for ( m in 1:length(Mod.Names) ) {
	mod <- Mod.Names[m]
	m.tag <- paste("hla",mod,sep="_")
	model <- JJ[[tag]][[m.tag]]
	summ <- summary(model,waic=F)
	
	## Collect HLA Info
	mod.inputs <- Modify_HLA_Inputs( mod )
	hla.preds <- setdiff( colnames(mod.inputs$Data), "ID" )
	hla.freqs <- colSums( mod.inputs$Data[,-1] )
	if ( mod=="Type" ) {
		gene.tag <- "DRB1"
		hla.tags <- gsub("_",":",gsub("Pr4_DRB1_","",hla.preds ))
	}else{ hla.tags <- hla.preds }
	
	## Collect General Model Info
	m.obs <- summ$nobs
	m.iter <- summ$iter
	m.warm <- summ$warmup
	m.chain <- summ$chains
	# m.waic <- summ$WAIC
	m.pheno <- as.character(model$formula)[2]
	m.covs <- as.character(model$formula)[3]
	m.form <- paste( m.pheno, "~", m.covs )
	RAND <- length(summ$random)>0

	## Collect Model Outputs
	f.eff <- summ$fixed
	c.eff <- summ$cor_pars
	s.eff <- summ$spec_pars
	m.covs.f <- rownames(f.eff)
	all.eff <- list( f.eff, c.eff, s.eff )
	mod.fit <- Reduce( rbind, all.eff )
	mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
	 # Posterior Probabilities
	d.prior <- model$prior
	d.post <- posterior_samples(model)
	d.post.f.which <- paste("b_",m.covs.f,sep="")
	d.post.c.which <- rownames(c.eff)
	d.post.sd <- round(unlist(lapply( d.post.f.which,function(x)sd(d.post[,x]) )),3)
	d.post.prob <- round(unlist(lapply( d.post.f.which,function(x)length(which(d.post[,x]>0)) )) / 7200, 2)
	d.post.prob <- data.frame( F.eff=d.post.f.which, SD=d.post.sd, Pr_l0=d.post.prob, Pr_Abs=.5+abs(.5-d.post.prob) )
	rownames(d.post.prob) <- m.covs.f
	# d.post.prob <- data.frame( d.post.prob, Pr_Abs=.5+abs(.5-round(d.post.prob[,"Pr_l0"],2)) )
	d.post.prob.rank <- order( d.post.prob$Pr_Abs, decreasing=T ) # order( abs(d.post.prob$Pr_l0-.5), decreasing=T )
	out.table <- merge( mod.fit, d.post.prob, by="row.names", all=T )
	rownames(out.table) <- out.table[,"Row.names"]
	out.table <- out.table[ rownames(mod.fit), -1 ]
	HLA.mod.fits[[m.tag]] <- out.table

	## PLOTS ######################

	## Plot Severity vs Response
	if ( 1 %in% which_plots ) {
		COLS.temp <- COLS.list.2[c(1,5)]
		png( paste(PathToPlot,"HLA-",mod,".1-RespVsSev.png",sep=""), height=1000,width=1000,pointsize=30 )
		par(mar=c( 4,4,4,3 ))
		# png( paste(PathToPlot,"HLA-",m.tag,".1-RespVsSev.png",sep=""), height=1200,width=2400,pointsize=32 ) ; par(mfrow=c(1,2))
		## Plot Betas
		XVALS <- f.eff[hla.preds,"Estimate"]
		YVALS <- f.eff[paste("DRUG:",hla.preds,sep=""),"Estimate"]
		XLIM <- extendrange( XVALS )
		YLIM <- extendrange( YVALS )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Disease Severity",ylab="Drug Response",main=paste("HLA Allele Effect Sizes:",m.tag) )
		abline( h=seq(-5,5,.1),v=seq(-5,5,.1),lty=3,col="grey50",lwd=1 )
		abline( h=0,v=0 )
		points( XVALS, YVALS, pch=16,col=adjustcolor(COLS.temp[1],.5),cex=2*log10(1+hla.freqs) )
		MOD <- lm( YVALS ~ XVALS, weights=hla.freqs )
		PVAL <- summary(MOD)$coefficients["XVALS","Pr(>|t|)"]
		abline( MOD, lty=2,lwd=6,col=COLS.temp[1] )
		text( XLIM[1],YLIM[1],paste("p=",formatC(PVAL,3,format="e"),sep=""), pos=4 )
		text( XVALS, YVALS, hla.tags, pos=3,cex=.8 )
		dev.off()
		## Plot Posterior Probs
		png( paste(PathToPlot,"HLA-",mod,".1b-RespVsSev.Post.png",sep=""), height=1000,width=1000,pointsize=30 )
		par(mar=c( 4,4,4,3 ))
		XVALS <- d.post.prob[d.post.prob$F.eff%in%paste("b_",hla.preds,sep=""),"Pr_l0"]
		YVALS <- d.post.prob[d.post.prob$F.eff%in%paste("b_DRUG:",hla.preds,sep=""),"Pr_l0"]
		XLIM <- c(0,1) # extendrange( XVALS )
		YLIM <- c(0,1) # extendrange( YVALS )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Disease Severity",ylab="Drug Response",main=paste("HLA Allele Posterior Probabilities (B>0):",m.tag),xaxt="n",yaxt="n" )
		axis( 1, at=seq(0,1,.1) )
		axis( 2, at=seq(0,1,.1), las=2 )
		abline( h=seq(0,1,.1),v=seq(0,1,.1),lty=3,col="grey50",lwd=1 )
		abline( h=.5,v=.5 )
		points( XVALS, YVALS, pch=16,col=adjustcolor(COLS.temp[2],.5),cex=2*log10(1+hla.freqs) )
		MOD <- lm( YVALS ~ XVALS )
		PVAL <- summary(MOD)$coefficients["XVALS","Pr(>|t|)"]
		abline( MOD, lty=2,lwd=6,col=COLS.temp[2] )
		text( XLIM[1],YLIM[1],paste("p=",formatC(PVAL,3,format="e"),sep=""), pos=4 )
		text( XVALS, YVALS, hla.tags, pos=3,cex=.8 )
		dev.off()
	}
	
	## Fixed Effect Sizes
	if ( 2 %in% which_plots ) {
		#  # Plot top 25 predictors
		# temp.n <- 25
		# mod.fit.f <- mod.fit[m.covs.f,][d.post.prob.rank[1:temp.n],]
		# PLOT_FIXED( mod.fit.f, "top25" )
		#  # Plot Clinical Covariates
		# mod.fit.c <- mod.fit[grep(paste(hla.preds,collapse="|"),rownames(mod.fit),invert=T),]
		# mod.fit.c <- mod.fit.c[intersect(rownames(mod.fit.c),m.covs.f),]
		# PLOT_FIXED( mod.fit.c, "clin" )
		 # Plot HLA Covariates
		mod.fit.h <- mod.fit[grep(paste(hla.preds,collapse="|"),rownames(mod.fit),invert=F),]
		mod.fit.h <- mod.fit.h[intersect(rownames(mod.fit.h),m.covs.f),]
		PLOT_FIXED_HLA( mod.fit.h, "hla" )
	}

	## Plot BRMS vs Mean Models
	if ( 3 %in% which_plots ) {
		COLS.temp <- COLS.ph
		mod.brms <- HLA.mod.fits[[m.tag]]
		if ( mod=="Type" ) {
			mod.mean <- t( B.out$Dig_4$TYP_DOS$DRB1 )
			mod.preds <- grep("^Pr4",rownames(mod.brms),value=T)
			mod.preds.print <- gsub("_",":",gsub("Pr4_DRB1_","",mod.preds ))
			mod.preds.mean <- paste("T",gsub(":","",mod.preds.print),sep="")
			mod.preds.int <- which( mod.preds.mean %in% rownames(mod.mean) )
			mod.preds.mean <- mod.preds.mean[mod.preds.int]
			mod.preds <- mod.preds[mod.preds.int]
			mod.preds.print <- mod.preds.print[mod.preds.int]
			mod.mean.sev <- mod.mean[mod.preds.mean,"DAS_BL_MN"]
			mod.mean.resp <- mod.mean[mod.preds.mean,"DEL_MNe_MN"]
			mod.freq <- hla.freqs[ mod.preds ]
			# For Table
			MN.beta <- c( mod.mean.sev, mod.mean.resp )
			MN.p <- c( t( P.out$Dig_4$TYP_DOS$DRB1 )[mod.preds.mean,c("DAS_BL_MN","DEL_MNe_MN")] )
			MN.se.das <- unlist(lapply( mod.preds.mean, function(x)summary(M.out$Dig_4$TYP_DOS$DRB1$DAS_BL_MN[[x]])$coefficients[x,"Std. Error"] ))
			MN.se.del <- unlist(lapply( mod.preds.mean, function(x)summary(M.out$Dig_4$TYP_DOS$DRB1$DEL_MNe_MN[[x]])$coefficients[x,"Std. Error"] ))
			MN.se <- c( MN.se.das, MN.se.del )
		}else{
			mod.mean <- HAP_OUT[[mod]]
			mod.preds <- rownames( mod.mean )
			mod.preds.print <- mod.preds
			mod.mean.sev <- mod.mean[mod.preds,"BETA.f.DAS_BL_MN"]
			mod.mean.resp <- mod.mean[mod.preds,"BETA.f.DEL_MNe_MN"]
			mod.freq <- mod.mean[mod.preds,"FREQ"]
			# For Table
			MN.beta <- c( mod.mean.sev, mod.mean.resp )
			MN.se <- c( mod.mean[mod.preds,"SE.f.DAS_BL_MN"], mod.mean[mod.preds,"SE.f.DEL_MNe_MN"])
			MN.p <- c( mod.mean[mod.preds,"P.f.DAS_BL_MN"], mod.mean[mod.preds,"P.f.DEL_MNe_MN"])
		}
		mod.brms.sev <- mod.brms[mod.preds,"Estimate"]
		mod.brms.resp <- mod.brms[paste("DRUG:",mod.preds,sep=""),"Estimate"]
		## Plot it
		png( paste(PathToPlot,"HLA-",mod,".3-MeanVsBRMS.png",sep=""), height=1000,width=2000,pointsize=30 )
		par(mfrow=c(1,2))
		 # Severity
		LIM <- extendrange(c(mod.brms.sev,mod.mean.sev))
		plot( 0,0,type="n",xlim=LIM,ylim=LIM,xlab="Mean Model",ylab="BRMS Model",main="Severity Effect Size Estimates" )
		abline(h=seq(-5,5,.5),v=seq(-5,5,.5),lty=3,col="grey50") ; abline(0,1,h=0,v=0)
		points( mod.brms.sev ~ mod.mean.sev, cex=2*log10(1+mod.freq),col=adjustcolor(COLS.temp[1],.6),pch=16 )
		MOD <- lm( mod.brms.sev ~ mod.mean.sev, )
		COR <- cor( mod.brms.sev, mod.mean.sev, method="spearman" )
		abline(MOD,lwd=4,lty=2,col=COLS.temp[1])
		text( LIM[2],LIM[1], paste("R=",round(COR,2),sep=""),pos=2 )
		text( mod.brms.sev ~ mod.mean.sev, label=mod.preds.print,pos=3 )
		 # Response
		LIM <- extendrange(c(mod.brms.resp,mod.mean.resp))
		plot( 0,0,type="n",xlim=LIM,ylim=LIM,xlab="Mean Model",ylab="BRMS Model",main="Response Effect Size Estimates" )
		abline(h=seq(-5,5,.5),v=seq(-5,5,.5),lty=3,col="grey50") ; abline(0,1,h=0,v=0)
		points( mod.brms.resp ~ mod.mean.resp, cex=2*log10(1+mod.freq),col=adjustcolor(COLS.temp[2],.6),pch=16 )
		MOD <- lm( mod.brms.resp ~ mod.mean.resp )
		COR <- cor( mod.brms.resp, mod.mean.resp, method="spearman" )
		abline(MOD,lwd=4,lty=2,col=COLS.temp[2])
		text( LIM[2],LIM[1], paste("R=",round(COR,2),sep=""),pos=2 )
		text( mod.brms.resp ~ mod.mean.resp, label=mod.preds.print,pos=3 )
		dev.off()

		## Combine & Save Table w/ Mean & BRMS Models
		# MN.beta <- c( mod.mean[mod.preds,"BETA.f.DAS_BL_MN"], mod.mean[mod.preds,"BETA.f.DEL_MNe_MN"])
		# MN.se <- c( mod.mean[mod.preds,"SE.f.DAS_BL_MN"], mod.mean[mod.preds,"SE.f.DEL_MNe_MN"])
		# MN.p <- c( mod.mean[mod.preds,"P.f.DAS_BL_MN"], mod.mean[mod.preds,"P.f.DEL_MNe_MN"])
		BRMS.beta <- c( mod.brms[mod.preds,"Estimate"], mod.brms[paste("DRUG:",mod.preds,sep=""),"Estimate"] )
		BRMS.se <- c( mod.brms[mod.preds,"Est.Error"], mod.brms[paste("DRUG:",mod.preds,sep=""),"Est.Error"] )
		BRMS.95l <- c( mod.brms[mod.preds,"l.95..CI"], mod.brms[paste("DRUG:",mod.preds,sep=""),"l.95..CI"] )
		BRMS.95u <- c( mod.brms[mod.preds,"u.95..CI"], mod.brms[paste("DRUG:",mod.preds,sep=""),"u.95..CI"] )
		BRMS.Pr <- c( mod.brms[mod.preds,"Pr_Abs"], mod.brms[paste("DRUG:",mod.preds,sep=""),"Pr_Abs"] )
		MN.FREQ <- data.frame( HAP=mod.preds.print, FREQ=mod.freq ) # mod.mean[mod.preds,c("HAP","FREQ")]
		# out.table <- data.frame( rbind(MN.FREQ,MN.FREQ), MN.beta,MN.se,MN.p, BRMS.beta,BRMS.se,BRMS.95l,BRMS.95u,BRMS.Pr )
		out.table <- data.frame( Effect=rep(c("Disease","Response"),each=nrow(MN.FREQ)), rbind(MN.FREQ,MN.FREQ), MN_SD=paste(round(MN.beta,3)," (",round(MN.se,3),")",sep=""), MN.p, BRMS_SD=paste(round(BRMS.beta,3)," (",round(BRMS.se,3),")",sep=""), BRMS_95=paste(round(BRMS.95l,3),"-",round(BRMS.95u,3),sep=""), BRMS.Pr, stringsAsFactors=F )
		out.table <- out.table[ order(out.table$Effect,out.table$HAP), ]
		out.table.2 <- out.table[ order(out.table$BRMS.Pr,decreasing=T), ]
		out.table.2 <- out.table.2[ order(out.table.2$Effect), ]
		out.table$Effect[which(duplicated(out.table$Effect))] <- ""
		HLA.mod.MNvBRMS[[m.tag]] <- out.table
		write.table( out.table, paste(PathToSave,"_TAB-MNvBRMS.",mod,".txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
		writeLines( print(xtable(out.table), include.rownames=F ), con=paste(PathToSave,"_LaTeX-MNvBRMS.",mod,".txt",sep="") )
		 # Sort by Posterior Prob
		out.table.2$Effect[which(duplicated(out.table.2$Effect))] <- ""
		write.table( out.table.2, paste(PathToSave,"_TAB-MNvBRMS.",mod,".ByPost.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
		writeLines( print(xtable(out.table.2), include.rownames=F ), con=paste(PathToSave,"_LaTeX-MNvBRMS.",mod,".ByPost.txt",sep="") )
	}
}

## Save Lists of Tables
save( HLA.mod.fits, file=paste(PathToSave,"_BayesLMM_HLA_ModFits.Rdata",sep="") )
save( HLA.mod.MNvBRMS, file=paste(PathToSave,"_BayesLMM_HLA_MNvBRMS.Rdata",sep="") )

## Save Tables w/ Fewer Columns
 # hla_Type
temp.write <- HLA.mod.MNvBRMS$hla_Type
temp.write <- temp.write[,c(1:3,6:8)]
temp.write.spl <- sapply( strsplit( temp.write$BRMS_SD, " (", fixed=T ), "[", 1:2 )
temp.write.mn <- as.numeric(temp.write.spl[1,])
temp.write.sd <- as.numeric( gsub(")","",temp.write.spl[2,]) )
temp.write <- data.frame( temp.write[,1:3], Mean=temp.write.mn, St.Dev=temp.write.sd, temp.write[,5:6] )
colnames(temp.write) <- c("Effect","Hap","Freq","Mean","St.Dev","95% CI","Post.Pr.")
write.table( temp.write, paste(PathToSave,"_TAB-MNvBRMS.hla_Type.BRMS.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
writeLines( print(xtable(temp.write), include.rownames=F ), con=paste(PathToSave,"_LaTeX-MNvBRMS.hla_Type.BRMS.txt",sep="") )
 # hla_p117174
temp.write <- HLA.mod.MNvBRMS$hla_p117174
temp.write <- temp.write[,c(1:3,6:8)]
temp.write.spl <- sapply( strsplit( temp.write$BRMS_SD, " (", fixed=T ), "[", 1:2 )
temp.write.mn <- as.numeric(temp.write.spl[1,])
temp.write.sd <- as.numeric( gsub(")","",temp.write.spl[2,]) )
temp.write <- data.frame( temp.write[,1:3], Mean=temp.write.mn, St.Dev=temp.write.sd, temp.write[,5:6] )
colnames(temp.write) <- c("Effect","Hap","Freq","Mean","St.Dev","95% CI","Post.Pr.")
write.table( temp.write, paste(PathToSave,"_TAB-MNvBRMS.hla_p117174.BRMS.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
writeLines( print(xtable(temp.write), include.rownames=F ), con=paste(PathToSave,"_LaTeX-MNvBRMS.hla_p117174.BRMS.txt",sep="") )

## Plot BRMS Model Summary
# PLOT_MOD( JJ[[tag]]$hla_pSE, "pSE", 1 )


#############################################################
## CREATE SUMMARY TABLES OF EFFECTS #########################
#############################################################
























#############################################################
## END OF DOC ###############################################
#############################################################
