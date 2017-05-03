setwd('./')
#Install required libraries
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
if("scales" %in% rownames(installed.packages()) == FALSE) {install.packages("scales")}
if("grid" %in% rownames(installed.packages()) == FALSE) {install.packages("grid")}
if("gtable" %in% rownames(installed.packages()) == FALSE) {install.packages("gtable")}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}
if("shape" %in% rownames(installed.packages()) == FALSE) {install.packages("shape")}
if("VennDiagram" %in% rownames(installed.packages()) == FALSE) {install.packages("VennDiagram")}
if("gridBase" %in% rownames(installed.packages()) == FALSE) {install.packages("gridBase")}
if("lattice" %in% rownames(installed.packages()) == FALSE) {install.packages("lattice")}
if("gplots" %in% rownames(installed.packages()) == FALSE) {install.packages("gplot")}
if("lattice" %in% rownames(installed.packages()) == FALSE) {install.packages("lattice")}
if("pracma" %in% rownames(installed.packages()) == FALSE) {install.packages("pracma")}

#Load libraries
require(RColorBrewer)
require(ggplot2)
require(gplots)
require(scales)
require(grid)
require(gtable)
require(gridExtra)
require(shape)
require(reshape2)
require(VennDiagram)
require(gridBase)
require(lattice)

setwd(dir="data")

#functions
source("../R/functions.R")
####################################
#         Summary figure          #
###################################
source("../R/summary_figure_drawing.R")

####################################
#         Fold-enrichment         #
###################################
source("../R/fold_enrichment_analysis_per_dataset.R")
source("../R/fold_enrichment_score_calibration_analysis.R")

####################################
#           Correlation analysis       #
###################################
source("../R/correlation_analysis.R")



