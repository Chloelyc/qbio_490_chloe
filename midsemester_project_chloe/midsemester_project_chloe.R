# create the directory and set working directory
dir.create("/Users/chloelyc/Desktop/QBIO490/qbio_490_chloe/midsemester_project_chloe/outputs")
setwd("/Users/chloelyc/Desktop/QBIO490/qbio_490_chloe/midsemester_project_chloe/outputs")

# install and load all the necessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
BiocManager::install("maftools")
library(maftools)
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
install.packages("ggplot2")
library(ggplots2)
install.packages("survival")
install.packages("survminer")
library(survival)
library(survminer)
BiocManager::install("maftools")
library(maftools)


# Download, query and prepare all the datasets
# perpare the clinical data
clinical_query <- GDCquery(project="TCGA-BRCA", data.category="Clinical", file.type = "xml")
GDCdownload(clinical_query)
clinical <- GDCprepare_clinic(clinical_query, clinical.info="patient")
colnames(clinical)[1] <- "Tumor_Sample_Barcode"

# create age category in clinical data frame
clinical$age_category <- ifelse(clinical$age_at_initial_pathologic_diagnosis<50, "young", "old")
age_mask <- ifelse(clinical$age_category == "young", T, F)
young_clinical <- clinical[age_mask,] # subset out the young clinical data frame
old_clinical <- clinical[!age_mask,] # subset out the old clinical data frame

# prepare the MAF data
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf <- GDCprepare(maf_query)

# create the MAF object
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

# read the patient drug and radiation data
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")




# creating the coOncoplot
# create the young and old maf onject
young_mask <- ifelse(maf_object@clinical.data$age_category=="young", T, F)
young_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[young_mask]
young_maf <- subsetMaf(maf = maf_object, tsb = young_patient_barcodes)
old_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[!young_mask]
old_maf <- subsetMaf(maf = maf_object, tsb = old_patient_barcodes)

# plotting the oncoplot
jpeg("./coOncoplot.jpg")
coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "young", 
           m2Name = "old")
dev.off()


# lollipopPlot2
jpeg("./lollipopPlot2.jpg")
lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              gene = "CDH1",
              m1_name = "young", 
              m2_name = "old")
dev.off()

# MAF survival plot
# creating the survival time and death event for the two maf_objects
old_maf@clinical.data$survival_time <- ifelse(is.na(old_maf@clinical.data$days_to_death)==FALSE, old_maf@clinical.data$days_to_death, old_maf@clinical.data$days_to_last_followup)
young_maf@clinical.data$survival_time <- ifelse(is.na(young_maf@clinical.data$days_to_death)==FALSE, young_maf@clinical.data$days_to_death, young_maf@clinical.data$days_to_last_followup)
old_maf@clinical.data$death_event <- ifelse(old_maf@clinical.data$vital_status=="Dead", T, F)
young_maf@clinical.data$death_event <- ifelse(young_maf@clinical.data$vital_status=="Dead", T, F)

# plot the survival plots
jpeg("./survival plot1.jpg")
mafSurvival(maf = old_maf,
            genes = "CDH1",
            time = "survival_time", 
            Status = "death_event", 
            isTCGA = TRUE)
dev.off()

jpeg("./survival plot2.jpg")
mafSurvival(maf = young_maf,
            genes = "CDH1", 
            time = "survival_time", 
            Status = "death_event", 
            isTCGA = TRUE)
dev.off()

# boxplot
clinical$age_category <- ifelse(clinical$age_at_initial_pathologic_diagnosis<50, "young", "old")
age_mask <- ifelse(clinical$age_category=="young", TRUE, FALSE)
young_clinical <- clinical[age_mask, ]
old_clinical <- clinical[!age_mask, ]
# calculate the days of radiation
clinical.rad$days_of_radiation <- clinical.rad$days_to_radiation_therapy_end
- clinical.rad$days_to_radiation_therapy_start
# filter out the NA
clinical.rad <- clinical.rad[!is.na(clinical.rad$days_of_radiation), ]
# filter out the duplicate
clinical.rad <- clinical.rad[!duplicated(clinical.rad$bcr_patient_barcode),]
# filter out the outliers
too_large <- ifelse(clinical.rad$days_of_radiation>500|clinical.rad$days_of_radiation<0,F,T)
clinical.rad <- clinical.rad[too_large,]

# correspond the patients in clinical with clinical.rad
young_rad_mask <- ifelse(unique(young_clinical)$Tumor_Sample_Barcode %in% clinical.rad$bcr_patient_barcode==T, T, F)
old_rad_mask <- ifelse(unique(old_clinical)$Tumor_Sample_Barcode %in% clinical.rad$bcr_patient_barcode==T, T, F)

# creating boxplot for young patients days of radiation
jpeg("./boxplot_young.jpg")
boxplot(clinical.rad$days_of_radiation[young_rad_mask])
dev.off()

# creating boxplot for old patients days of radiation
jpeg("./boxplot_old.jpg")
boxplot(clinical.rad$days_of_radiation[old_rad_mask])
dev.off()

# printing the statistics of days of radiation
summary(clinical.rad$days_of_radiation[young_rad_mask])
summary(clinical.rad$days_of_radiation[old_rad_mask])

