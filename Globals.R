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
