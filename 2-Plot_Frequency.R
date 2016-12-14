## Figure 2 for HLA Manuscript ##
## Plot Frequencies of Types, AAs, and Collapsed Haplotypes ##
## April 5, 2016 ##
## Kristopher Standish ##

## New Plan: August 26, 2016 ##
 # 1) 4-digit (to supplemental)
   # Break down by EUR and non-EUR
 # 2) Keep AA plot for DRB1
 # 3) plot B_9, DPB1_9, DRB1_11_71_74
   # Break down by EUR and non-EUR


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
COLS.fr <- COLS.list.2[1]
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
## PLOT HLA TYPE FREQUENCIES ################################
#############################################################

## FCT: Plot Allele Frequency of Each Haplotype
Plot_AF <- function( Samples, tag ) {
	# FREQS <- list()
	# for ( g in 1:N.GENE ) {
	# 	gene <- GENE_LIST[g]
	# 	FREQ.4 <- rowSums( PAT_DOS.4[[gene]][,Samples] ) ; names(FREQ.4)[which(names(FREQ.4)==":")]<-"NA" ; FREQ.4 <- FREQ.4[order(names(FREQ.4))]
	# 	names(FREQ.4) <- gsub(":","",names(FREQ.4),fixed=T)
	# 	FREQ.4 <- FREQ.4[ which(names(FREQ.4)!="NA" & FREQ.4!=0) ]
	# 	FREQ.2 <- rowSums( PAT_DOS.2[[gene]][,Samples] ) ; FREQ.2 <- FREQ.2[order(names(FREQ.2))]
	# 	## Plot it
	# 	# 4 digit
	# 	YLIM <- c( 0, max(FREQ.4) ) * c(1,1.25)
	# 	LINES <- which(!duplicated( substr(names(FREQ.4),1,2) ))
	# 	png( paste(PathToPlot,tag,"-2B-TypeFreq.Pr4.",gene,".png",sep=""), height=800,width=2000,pointsize=30 )
	# 	# par(mfrow=c(2,1))
	# 	par(mar=c(4,4,3,2))
	# 	TEMP.4 <- barplot( FREQ.4, col=COLS.fr,border=NA,ylim=YLIM,
	# 		main=paste("Aggregate Haplotype Count (SOAP-HLA 4-digit): HLA-",gene,sep=""),ylab="# Haplotypes",xlab="Haplotype",las=2 )
	# 	abline( h=seq(0,400,50),lty=3,col="grey50",lwd=1 )
	# 	abline( v=TEMP.4[LINES]-.6,lty=1,col="grey20",lwd=1 )
	# 	TEMP.4 <- barplot( FREQ.4, col=COLS.fr,border=NA,las=2,add=T )
	# 	text( TEMP.4, FREQ.4+YLIM[2]*.1, label=paste("n=",FREQ.4,sep=""),srt=90,cex=.9 )
	# 	dev.off()
	# 	# 2 digit
	# 	YLIM <- c( 0, max(FREQ.2) ) * c(1,1.25)
	# 	png( paste(PathToPlot,tag,"-2B-TypeFreq.Pr2.",gene,".png",sep=""), height=800,width=1600,pointsize=30 )
	# 	par(mar=c(4,4,3,2))
	# 	TEMP.2 <- barplot( FREQ.2, col=COLS.fr,border=NA,ylim=YLIM,
	# 		main=paste("Aggregate Haplotype Count (SOAP-HLA 2-digit): HLA-",gene,sep=""),ylab="# Haplotypes",xlab="Haplotype",las=2 )
	# 	abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
	# 	TEMP.2 <- barplot( FREQ.2, col=COLS.fr,border=NA,las=2,add=T )
	# 	text( TEMP.2, FREQ.2+YLIM[2]*.1, label=paste("n=",FREQ.2,sep=""),srt=90,cex=.9 )
	# 	dev.off()
	# 	FREQS[[gene]] <- list( FREQ.4=FREQ.4, FREQ.2=FREQ.2 )
	# }
	# return(FREQS)
}

## FCT: Plot Allele Frequency broken down by EUR or non-EUR
Plot_AF.EUR <- function( Samples, genes, tag ) {
	FREQS <- list()
	if ( tag=="ART3001" ) {
		Samps.EUR.which <- which( Samples %in% FT$ID[FT$RACE=="WHITE" & FT$ETHN=="NOT HISPANIC OR LATINO"] )
	}else{
		Samps.EUR.which <- which( Samples %in% FT.3002$ID[FT.3002$RACE=="WHITE" & FT.3002$ETHNIC=="NOT HISPANIC OR LATINO"] )
	}
	Samps.EUR <- Samples[Samps.EUR.which]
	N.GENE <- length(genes)
	for ( g in 1:N.GENE ) {
		gene <- genes[g]
		FREQ.4.EUR <- rowSums( PAT_DOS.4[[gene]][,Samples[Samps.EUR.which]] )
		FREQ.4.nEUR <- rowSums( PAT_DOS.4[[gene]][,Samples[-Samps.EUR.which]] )
		FREQ.4 <- data.frame( EUR=FREQ.4.EUR, NON=FREQ.4.nEUR )
		FREQ.4 <- FREQ.4[ rownames(FREQ.4)!=":", ]
		rownames(FREQ.4) <- paste( gene, rownames(FREQ.4), sep="*" )
		rownames(FREQ.4) <- gsub( ":$","", rownames(FREQ.4) )
		FREQ.4 <- FREQ.4[ order(rownames(FREQ.4)), ]
		FREQ.4 <- FREQ.4[ rowSums(FREQ.4)!=0, ]
		FREQ.4.sc <- apply( FREQ.4, 2, function(x)x/sum(x) )

		## Plot it
		# 4 digit
		YLIM <- c( 0, max(FREQ.4.sc) ) * c(1,1.25)
		LINES <- which(!duplicated( substr(rownames(FREQ.4.sc),nchar(gene)+1,nchar(gene)+2) ))
		COLS.fr.2 <- c( BLEND(c("grey30",COLS.fr)), COLS.fr )
		YBRK <- ifelse( YLIM[2]>.1, 5, 2 )
		png( paste(PathToPlot,tag,"-2B-TypeFreq.Pr4.",gene,".png",sep=""), height=800,width=2000,pointsize=30 )
		par(mar=c(7,4,3,2))
		TEMP.4 <- barplot( t(FREQ.4.sc), beside=T,col=COLS.fr.2,border=NA,
			ylim=YLIM,yaxt="n",
			main=paste("Haplotype Frequency (SOAP-HLA 4-digit): HLA-",gene,sep=""),ylab="Percent Haplotypes",xlab="",las=2 )
		axis( 2, at=seq(0,1,YBRK/100), label=paste(seq(0,100,YBRK),"%",sep=""), las=2 )
		abline( h=seq(0,1,YBRK/100),lty=3,col="grey50",lwd=1 )
		TEMP.4 <- barplot( t(FREQ.4.sc), beside=T,col=COLS.fr.2,border=NA,add=T,xaxt="n",yaxt="n" )
		legend( "topright", legend=c("EUR","Non-EUR"), fill=COLS.fr.2,border=NA, bg="white",ncol=2 )		
		dev.off()

		FREQS[[gene]] <- FREQ.4
	}
	return(FREQS)
}

## Run Function on Various Sample Lists
SAMPS.temp <- SAMPS[1:2]
genes <- c("B","DPB1","DRB1")
FREQS <- lapply( names(SAMPS.temp), function(x) Plot_AF.EUR(SAMPS[[x]],genes,x) )
names(FREQS) <- names(SAMPS.temp)

# ## Compare ART3001 to ART3002
# FREQS.DRB1 <- merge( FREQS$ART3001$DRB1, FREQS$ART3002$DRB1, all=T,by="row.names" )
# FREQS.DRB1[is.na(FREQS.DRB1)] <- 0
# rownames(FREQS.DRB1) <- gsub("DRB1","", FREQS.DRB1[,1])
# FREQS.DRB1 <- FREQS.DRB1[,-1]
# colnames(FREQS.DRB1) <- gsub("x","ART3001",colnames(FREQS.DRB1))
# colnames(FREQS.DRB1) <- gsub("y","ART3002",colnames(FREQS.DRB1))
# FREQS.DRB1.sc <- apply( FREQS.DRB1, 2, function(x)x/sum(x) )
# COLS.fr <- COLS.list.2[c(1,7)]
# COLS.temp <- c( BLEND(c("grey30",COLS.fr[1])), COLS.fr[1], BLEND(c("grey30",COLS.fr[2])), COLS.fr[2] )
# # barplot( t(FREQS.DRB1), col=COLS.temp, border=NA, beside=T, names.arg=rownames(FREQS.DRB1),las=2, legend=T )
# barplot( t(FREQS.DRB1.sc), col=COLS.temp, border=NA, beside=T,legend=T,ylab="% of Haplotypes (Within-Group)",main="HLA-DRB1 Frequencies",las=2 )

# MG <- merge( t(PAT_DOS.4$DRB1), FT[,c("ID","RACE","ETHN","COUN")],by.x="row.names",by.y="ID" )
# MG.3002 <- merge( t(PAT_DOS.4$DRB1), FT.3002[,c("ID","RACE","ETHNIC","COUNTRY")],by.x="row.names",by.y="ID" )

# par(mfrow=c(2,2))
# par(mar=c(7,4,3,2))
# barplot(t(prop.table(table(MG$RACE,MG$ETHN))),beside=T,las=2,main="ART3001",ylim=c(0,1),ylab="% Patients",legend=T) ; abline(h=seq(0,1,.1),lty=3)
# barplot(t(prop.table(table(MG.3002$RACE,MG.3002$ETHNIC))),beside=T,las=2,main="ART3002",ylim=c(0,1),ylab="% Patients",legend=T) ; abline(h=seq(0,1,.1),lty=3)

# # par(mfrow=c(1,2))
# # par(mar=c(7,4,3,2))
# barplot(t(prop.table(table(MG$RACE,MG$COUN))),beside=T,las=2,main="ART3001",ylim=c(0,1),ylab="% Patients",legend=T) ; abline(h=seq(0,1,.1),lty=3)
# barplot(t(prop.table(table(MG.3002$RACE,MG.3002$COUNTRY))),beside=T,las=2,main="ART3002",ylim=c(0,1),ylab="% Patients",legend=T) ; abline(h=seq(0,1,.1),lty=3)

# par(mfrow=c(2,1))
# par(mar=c(4,4,1,2))
# barplot( t(FREQS.DRB1.sc[,c(1,3)]), col=COLS.temp[c(1,3)], border=NA, beside=T, las=2,ylab="% of Haplotypes (Within-Group)" )
# barplot( FREQS.DRB1.sc[,1]-FREQS.DRB1.sc[,3], col=COLS.list.2[4],border=NA,las=2 )
# table( MG$COUN,MG$`01:01`,MG$RACE )[,,"WHITE"]
# table( MG.3002$COUNTRY,MG.3002$`01:01`,MG.3002$RACE )[,,"WHITE"]
# table(MG$RACE,MG$COUN,MG$`01:01`)
# table(MG.3002$RACE,MG.3002$COUNTRY,MG.3002$`01:01`)
# table(MG$RACE,MG$`13:02`,MG$ETHN)[,,1:2]
# table(MG.3002$RACE,MG.3002$`13:02`,MG.3002$ETHNIC)[,,1:2]

# MOD.temp <- lm( EUR.ART3001 ~ EUR.ART3002, data=FREQS.DRB1 )
# # par(mfrow=c(1,2))
# layout( matrix(1:2,ncol=2), widths=c(2,3))
# plot( EUR.ART3001 ~ EUR.ART3002, data=FREQS.DRB1 )
# abline(MOD.temp)
# text( EUR.ART3001 ~ EUR.ART3002, data=FREQS.DRB1, labels=rownames(FREQS.DRB1),pos=2 )
# summary(MOD.temp)
# barplot( resid(MOD.temp), las=2 )

# par(mfrow=c(3,1))
# par(mar=c(0,4,4,2))
# barplot( t(FREQS.DRB1.sc), col=COLS.temp, border=NA, beside=T,legend=T,names.arg=rep("",nrow(FREQS.DRB1.sc)),ylab="% of Haplotypes (Within-Group)",main="HLA-DRB1 Frequencies" )
# par(mar=c(0,4,1,2))
# barplot( t(FREQS.DRB1.sc[,c(2,4)]), col=COLS.temp[c(2,4)], border=NA, beside=T,names.arg=rep("",nrow(FREQS.DRB1.sc)),ylab="% of Haplotypes (Within-Group)" )
# par(mar=c(4,4,1,2))
# barplot( t(FREQS.DRB1.sc[,c(1,3)]), col=COLS.temp[c(1,3)], border=NA, beside=T, las=2,ylab="% of Haplotypes (Within-Group)" )

# par(mfrow=c(2,1))
# par(mar=c(4,4,1,2))
# barplot( t(FREQS.DRB1.sc[,c(1,3)]), col=COLS.temp[c(1,3)], border=NA, beside=T, las=2,ylab="% of Haplotypes (Within-Group)" )
# barplot( FREQS.DRB1.sc[,1]-FREQS.DRB1.sc[,3], col=COLS.list.2[4],border=NA,las=2 )
# ## Run Function on Various Sample Lists
# # FREQS <- lapply( names(SAMPS), function(x) Plot_AF(SAMPS[[x]],x) )
# # names(FREQS) <- names(SAMPS)

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
	barplot( t(t(pat_tab.prc[,x])),add=T,xaxt="n",yaxt="n",beside=F,border=NA,col=COLS.AA[pat_tab.names[,x]],space=c(xval,0),width=1)
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
		png( paste(PathToPlot,tag,"-2C-AAfreq_",gene,".png",sep=""), height=800,width=2000,pointsize=30 )
		layout( matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(7,1) ) # layout( matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(3,2), heights=c(1,1) )
		plot( 0,0,type="n",xlim=XLIM,ylim=c(0,1), main=paste("Amino Acid Frequency - HLA-",gene,sep=""),xlab="Amino Acid",ylab="Frequency",xaxt="n")
		axis( 1, at=seq(-1000,1000,20),las=2 )	
		# abline( v=seq(-1000,1000,20),lty=3,col="grey50" )
		abline( h=seq(0,1,.2),lty=3,col="grey50" ) ; abline( h=c(0,1) )
		DWAI <- lapply( which.temp, function(x)allele_bar(pat_tab.prc,pat_tab.names,x) )
		## Amino Acid Key
		Which_AA <- which(names(COLS.AA)%in% pat_aa )
		COLS.aa <- COLS.AA[Which_AA]
		barplot( matrix(rep(1,length(COLS.aa)),ncol=1), beside=F,col=COLS.aa,xaxt="n",yaxt="n",ylab="Amino Acid Residue Key")
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
PLOT_HAP_FREQ <- function( Positions, tag ) {
	# ## Pull together Haplotypes from Specified Positions
	# N.res <- length(Positions)
	# HAP <- apply( PAT_AA$DRB1[,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
	# HAP <- gsub("NA","-",HAP)
	# HAP.uniq <- sort(unique(HAP))
	# HAP.plot <- setdiff( HAP.uniq, paste(rep("-",N.res),collapse=""))
	#  # Get Haplotype Frequencies
	# HAP.freq <- table(HAP)
	# HAP.freq <- HAP.freq[ HAP.plot ]

	# ## Plot Haplotype Frequency
	# YLIM <- c( 0, max(HAP.freq) ) * c(1,1.35)
	# png( paste(PathToPlot,"/DRB1_",tag,"_2D-HapFreq.png",sep=""), height=800,width=1000,pointsize=30 )
	# TEMP <- barplot( HAP.freq,las=2,col=COLS.fr,border=NA,main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype Frequency"),ylab="# Haplotypes",ylim=YLIM)
	# abline(h=seq(0,1000,50),lty=3,col="grey50")
	# barplot( HAP.freq,las=2,col=COLS.fr,border=NA,add=T )
	# text( TEMP, HAP.freq+.13*YLIM[2], label=paste("(n=",HAP.freq,")",sep=""), srt=90 )
	# dev.off()
}

## Function to do Haplotype Analysis of Specified Amino Acid Positions
PLOT_HAP_FREQ.EUR <- function( Positions, tag ) {
	Samples <- SAMPS$ART3001
	Samps.EUR.which <- which( Samples %in% FT$ID[FT$RACE=="WHITE" & FT$ETHN=="NOT HISPANIC OR LATINO"] )
	Samps.EUR <- Samples[Samps.EUR.which]

	## Pull together Haplotypes from Specified Positions
	N.res <- length(Positions)
	HAP <- apply( PAT_AA$DRB1[,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
	HAP <- gsub("NA","-",HAP)
	HAP.uniq <- sort(unique(HAP))
	HAP.plot <- setdiff( HAP.uniq, paste(rep("-",N.res),collapse=""))
	 # Get Haplotype Frequencies
	HAP.samps <- sapply(strsplit( names(HAP), "_" ),"[",1)
	HAP.samps.EUR.which <- which(HAP.samps %in% Samps.EUR )
	
	HAP.freq.EUR <- table( HAP[HAP.samps.EUR.which] )
	HAP.freq.nEUR <- table( HAP[-HAP.samps.EUR.which] )
	HAP.freq <- merge( data.frame(HAP.freq.EUR), data.frame(HAP.freq.nEUR), by="Var1", all=T)
	colnames(HAP.freq) <- c("HAP","EUR","NON")
	rownames(HAP.freq) <- HAP.freq$HAP
	HAP.freq <- HAP.freq[ HAP.plot, ]
	HAP.freq[is.na(HAP.freq)] <- 0
	HAP.freq.sc <- apply( HAP.freq[,-1], 2, function(x)x/sum(x) )

	## Plot it
	# 4 digit
	YLIM <- c( 0, max(HAP.freq.sc) ) * c(1,1.25)
	LINES <- which(!duplicated( substr(rownames(HAP.freq.sc),nchar(gene)+1,nchar(gene)+2) ))
	COLS.fr.2 <- c( BLEND(c("grey30",COLS.fr)), COLS.fr )
	YBRK <- ifelse( YLIM[2]>.1, 5, 2 )
	png( paste(PathToPlot,"/DRB1_",tag,"_2D-HapFreq.png",sep=""), height=800,width=2000,pointsize=30 )
	par(mar=c(4,4,3,2))
	TEMP.4 <- barplot( t(HAP.freq.sc), beside=T,col=COLS.fr.2,border=NA,
		ylim=YLIM,yaxt="n",
		main=paste("Collapsed Haplotype Frequency: Pos",paste(Positions,collapse=",")),ylab="Percent Haplotypes",xlab="",las=2 )
	axis( 2, at=seq(0,1,YBRK/100), label=paste(seq(0,100,YBRK),"%",sep=""), las=2 )
	abline( h=seq(0,1,YBRK/100),lty=3,col="grey50",lwd=1 )
	TEMP.4 <- barplot( t(HAP.freq.sc), beside=T,col=COLS.fr.2,border=NA,add=T,xaxt="n",yaxt="n" )
	legend( "topright", legend=c("EUR","Non-EUR"), fill=COLS.fr.2,border=NA, bg="white",ncol=2 )		
	dev.off()
	## Return Frequencies
	return(HAP.freq)
}


OUT <- list()
##########################################
## POS 11, 71, 74 ##
 # Viatte, et al (2015)
Positions <- c(11,71,74)
tag <- paste(c("p",Positions),collapse="")
OUT[[tag]] <- PLOT_HAP_FREQ.EUR(Positions,tag)
#  # Viatte, et al (2015)
# Positions <- c(11,71,74)
# tag <- paste(c("p",Positions),collapse="")
# OUT[[tag]] <- PLOT_HAP_FREQ(Positions,tag)

##########################################
## POS 11, 13, 71, 74 ##
# Positions <- 70:74
# tag <- "pSE"
# OUT[[tag]] <- PLOT_HAP_FREQ(Positions,tag)

#############################################################
## PLOT AMINO ACID FREQUENCIES (top 5) ######################
#############################################################

## Specify EUR/NON Samples
Samples <- SAMPS$ART3001
Samps.EUR.which <- which( Samples %in% FT$ID[FT$RACE=="WHITE" & FT$ETHN=="NOT HISPANIC OR LATINO"] )
Samps.EUR <- Samples[Samps.EUR.which]
Samps.nEUR <- Samples[-Samps.EUR.which]
Samps.EUR.typ <- paste( rep(Samps.EUR,each=2),rep(1:2,length(Samps.EUR)),sep="_" )
Samps.nEUR.typ <- paste( rep(Samps.nEUR,each=2),rep(1:2,length(Samps.nEUR)),sep="_" )

## Pull out Specific Amino Acids
genes <- c("B","DPB1","DRB1")
AA.cand <- list()
for ( g in 1:3 ) {
	gene <- genes[g]
	if ( gene %in% c("B","DPB1") ) {
		positions <- 9
	}else{
		positions <- c(11,71,74)
	}
	position.labs <- paste("Pos",positions,sep="_")
	# Pull out Patient Data for Gene
	pat_aa.EUR <- data.frame( EUR=PAT_AA[[gene]][Samps.EUR.typ,position.labs] )
	pat_aa.nEUR <- data.frame( NON=PAT_AA[[gene]][Samps.nEUR.typ,position.labs] )
	# Return Patient AA Freqs
	AA.cand[[gene]] <- list( EUR=pat_aa.EUR, NON=pat_aa.nEUR )

}

## Get Counts for each AA and each Population
AA.cand.tab <- lapply( AA.cand, function(x)lapply(x,function(y)apply(y,2,table)) )
 # Compile into Tables
AA.cand.fr <- list()
AA.cand.fr$B_9 <- merge( AA.cand.tab$B$EUR,AA.cand.tab$B$NON, by="row.names", all=T )
AA.cand.fr$DPB1_9 <- merge( AA.cand.tab$DPB1$EUR,AA.cand.tab$DPB1$NON, by="row.names", all=T )
AA.cand.fr$DRB1_11 <- merge( data.frame(AA.cand.tab$DRB1$EUR$EUR.Pos_11),data.frame(AA.cand.tab$DRB1$NON$NON.Pos_11), by="Var1", all=T )
AA.cand.fr$DRB1_71 <- merge( data.frame(AA.cand.tab$DRB1$EUR$EUR.Pos_71),data.frame(AA.cand.tab$DRB1$NON$NON.Pos_71), by="Var1", all=T )
AA.cand.fr$DRB1_74 <- merge( data.frame(AA.cand.tab$DRB1$EUR$EUR.Pos_74),data.frame(AA.cand.tab$DRB1$NON$NON.Pos_74), by="Var1", all=T )

## Convert to % of population
AA.cand.sc <- list()
for ( tag in names(AA.cand.fr) ) {
	AA.cand.fr[[tag]][is.na(AA.cand.fr[[tag]])] <- 0
	AA.cand.sc[[tag]] <- apply(AA.cand.fr[[tag]][,-1],2,function(x)x/sum(x) )
	rownames(AA.cand.sc[[tag]]) <- AA.cand.fr[[tag]][,1]
}

## Plot it!!
COLS.fr.2 <- c( BLEND(c("grey30",COLS.fr)), COLS.fr )
png( paste(PathToPlot,tag,"-2E-Cand.AAfreq.png",sep=""), height=800,width=2000,pointsize=36 )
layout( matrix(1:5,ncol=5), widths=sapply(AA.cand.sc,nrow)+1 )
for ( x in names(AA.cand.sc) ) {
	ylab <- ifelse( x==names(AA.cand.sc)[1], "Percent Haplotypes", "" )
	temp.split <- strsplit(x,"_")[[1]]
	barplot( t(AA.cand.sc[[x]]),beside=T, col=COLS.fr.2,border=NA, main=paste("HLA",temp.split[1],sep="-"),xlab=paste("Pos",temp.split[2]),ylab=ylab,yaxt="n" )
	axis( 2, at=seq(0,1,.1), label=paste(seq(0,100,10),"%",sep=""), las=2 )
	abline( h=seq(0,1,.1), lty=3,col="grey50",lwd=1 )
	barplot( t(AA.cand.sc[[x]]),beside=T, col=COLS.fr.2,border=NA, main=paste("HLA",temp.split[1],sep="-"),xlab=paste("Pos",temp.split[2]),ylab=ylab,yaxt="n", add=T )
}
legend( "topright", legend=c("EUR","Non-EUR"), fill=COLS.fr.2,border=NA, bg="white",ncol=1 )
dev.off()






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
