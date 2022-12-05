
knitr::opts_knit$set(root.dir = normalizePath("/Users/chloelyc/Desktop/QBIO490/qbio_490_chloe/final_project"))

map_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "MAP3K1"]
gata_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "GATA3"]
cdh_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "CDH1"]
tp_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "TP53"]

map_mask <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode%in%levels(factor(map_data)),T,F)
map_clinical <- maf_object@clinical.data[map_mask,]
maptp <- intersect(map_data, tp_data)
map_clinical$mutation <- ifelse(map_clinical$Tumor_Sample_Barcode%in%maptp, "MAP3K1 and TP53", "MAP3K1")

surv_object_map <- Surv(time = as.numeric(map_clinical$survival_time), event = as.logical(map_clinical$death_event))

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
map_fit <- surv_fit( surv_object_map ~ map_clinical$mutation,
                     data = map_clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_map = ggsurvplot(map_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_map = survplot_map$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))


KM_plot_map


gata_mask <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode%in%levels(factor(gata_data)),T,F)
gata_clinical <- maf_object@clinical.data[gata_mask,]
gatatp <- intersect(gata_data, tp_data)
gata_clinical$mutation <- ifelse(gata_clinical$Tumor_Sample_Barcode%in%gatatp, "GATA3 and TP53", "GATA3")

surv_object_gata <- Surv(time = as.numeric(gata_clinical$survival_time), event = as.logical(gata_clinical$death_event))

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
gata_fit <- surv_fit( surv_object_gata ~ gata_clinical$mutation,
                     data = gata_clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_gata = ggsurvplot(gata_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_gata = survplot_gata$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))


KM_plot_gata


cdh_mask <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode%in%levels(factor(cdh_data)),T,F)
cdh_clinical <- maf_object@clinical.data[cdh_mask,]
cdhtp <- intersect(cdh_data, tp_data)
cdh_clinical$mutation <- ifelse(cdh_clinical$Tumor_Sample_Barcode%in%cdhtp, "CDH1 and TP53", "CDH1")

surv_object_cdh <- Surv(time = as.numeric(cdh_clinical$survival_time), event = as.logical(cdh_clinical$death_event))

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
cdh_fit <- surv_fit( surv_object_cdh ~ cdh_clinical$mutation,
                      data = cdh_clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_cdh = ggsurvplot(cdh_fit, 
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                           legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_cdh = survplot_cdh$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))


KM_plot_cdh


