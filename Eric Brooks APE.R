


# 1. Load Library----------------------------------------------------------------------------------------------------------------------------------------

library("rmarkdown")
library("jsonlite")
library("readr")
library("dplyr")
library("sesame")
# sesameDataCache() # only on new install of sesame
library("BiocParallel")
library("tidyr")
library("ggplot2")
library("pals")
library("tibble")
library("SummarizedExperiment")
library("DelayedArray")
library("stats")
library("randomForest")
library("parallel")
library("tidyverse")
library("DNAcopy")
library("knitr")
library("ggrepel")
library("data.table")
library("impute")
library("factoextra")
library("NbClust")
library("HDF5Array")

# # Edit MulticoreParam to SnowParam to parallelize on windows, leave MulticoreParam if on mac
trace(DML, edit=TRUE)

# 2. Read in Data Using Sesame Pipeline--------------------------------------------------------------------------------------------------------------------

# # read in methylation data
# # cpg beta values from raw IDAT files
betas <- openSesame("./APE/TARGET-AML/data types/dna methylation by platform/HM450",BPPARAM = BiocParallel::SnowParam(12))
closeAllConnections()

# # Can read in from CSV after writing betas to csv to save time


# # SigDef Files from IDAT
sdfs <- openSesame("./APE/TARGET-AML/data types/dna methylation by platform/HM450", func = NULL, BPPARAM = BiocParallel:: SnowParam(12))
closeAllConnections()

# # preprocessing, qualityMask, inferInfiniumChannel, dyeBiasNL, pOOBAH, noob
sdfs <- openSesame(sdfs, prep = "QCDPB", func = NULL)

# !!OPTION - Read in betas_filtered and betas_no_xy_no_fibro-----------------------------------------------
# betas <- fread("./APE/TARGET-AML/data files/betas.csv")
# betas <- as.data.frame(betas)
# rownames(betas) <- betas[,1]
# betas <- betas[,-1]


# 3. Sample Record Masked Intensities -------------------------------------------------------------------------------------------------------------------

joined_record_MI <- fromJSON("./APE/TARGET-AML/data files/joined_record_masked_intensities.json", flatten = T)




# 4. Reorganize Data-------------------------------------------------------------------------------------------------------------------------------------

# # rename IDs by matching file_id to sample_id in record file
ind_sdfs_full <- match(names(sdfs), joined_record_MI$file_name)

names(sdfs) <- joined_record_MI$sample_id[ind_sdfs_full]



# # Only need to do if reading betas in with Sesame
# ind_betas_full <- match(colnames(betas), joined_record_MI$file_name)
# colnames(betas) <- joined_record_MI$sample_id[ind_betas_full]


joined_record_MI_clean_full <- as.data.frame(joined_record_MI[,-c(1,3,4,5,7,8)])



# # Write to CSV
# fwrite(as.data.frame(betas), "./APE/TARGET-AML/data files/betas.csv", row.names = T, col.names = T)


betas <- as.data.frame(betas)


joined_record_MI_clean_matched <- match(colnames(betas), joined_record_MI_clean_full$sample_id)

joined_record_MI_clean_matched_df <- joined_record_MI_clean_full[joined_record_MI_clean_matched,]


remove(joined_record_MI_clean_full)
remove(joined_record_MI)
remove(joined_record_MI_clean_matched)
remove(ind_betas_full)
remove(ind_sdfs_full)


# 5. Quality Control --------------------------------------------------------------------------------------------------------------------


qcs <- openSesame(sdfs, prep ="", func = sesameQC_calcStats)
qcs[[1]]

qcs_det <- openSesame(sdfs, prep ="", func = sesameQC_calcStats, funs = "detection")
qcs_det[[1]]

sesameQC_getStats(qcs_det[[1]], "frac_dt")
qcs_det_df <- head(do.call(rbind, lapply(qcs_det, as.data.frame)))


qcs_int <- openSesame(sdfs, prep = "", func = sesameQC_calcStats, funs = c("intensity"))

qcs_int[[1]]

sesameQC_rankStats(qcs_int[[1]], platform =  "EPIC")

gc()

# QC Plots

sesameQC_plotBar(lapply(sdfs, sesameQC_calcStats, "detection"))

sesameQC_plotRedGrnQQ(sdfs[[1]])

sesameQC_plotIntensVsBetas(sdfs[[1]])

head(getAFs(sdfs[[1]]))

head(getAFTypeIbySumAlleles(sdfs[[1]]))

sesameQC_plotHeatSNPs(sdfs)

gc()

rm(qcs)
rm(qcs_det)
rm(qcs_det_df)
rm(qcs_int)


# 6. Inferences ---------------------------------------------------------------------------------------------------------------------------------------------

inferSex(sdfs[[1]], platform = "HM450")
inferSexKaryotypes(sdfs[[1]])

inferEthnicity(sdfs[[1]])

sdf = sdfs[[1]]
segs <- cnSegmentation(sdf, sdfs.normal = sdfs)
visualizeSegments(segs)

sdf212 = sdfs[[212]]
segs212 <- cnSegmentation(sdf212, sdfs.normal = sdfs)
visualizeSegments(segs212)
gc()

rm(sdfs)
rm(sdf)
rm(sdf212)

# 7. Data filtering for Qulcore -------------------------------------------------------------------------------------------------------------------------

JR_full <- joined_record_MI_clean_matched_df %>%
  separate(sample_id,c("project", "tumor_code", "case_identifier", "tissue_code" ), sep = "-", remove = FALSE)


JR_full <- JR_full %>% mutate(tumor_type =
                                case_when(tumor_code == 20 ~ "AML",
                                          tumor_code == 21 ~ "AML-IF")
)


JR_full <- JR_full[,-c(2,3,5)]

JR_full <- JR_full[, c(1,2,8,3,4,5,6,7)]

rownames(JR_full) <- JR_full[,1]

JR_full <- JR_full[,-1]


# write to CSV
# fwrite(JR_full, "./APE/TARGET-AML/data files/JR_full.csv", row.names = T, col.names = T)


betas_df_t <- as.data.frame(t(betas))


betas_merge <- as.data.frame(t(merge(JR_full, betas_df_t,
                     by = 'row.names', all = TRUE)))

colnames(betas_merge) <- betas_merge[1,]

betas_merge <- betas_merge[-1,]

# remove rows with all NA values, every column is NA
betas_no_na_col <- betas_merge[rowSums(is.na(betas_merge)) != ncol(betas_merge), ]


# Remove columns and rows with more than 10% NA
betas_filtered <- betas_no_na_col[which(rowMeans(!is.na(betas_no_na_col)) > 0.9), which(colMeans(!is.na(betas_no_na_col)) > 0.9)]

betas_filtered <- as.data.frame(betas_filtered)

# write to CSV
# fwrite(betas_filtered, "./APE/TARGET-AML/data files/betas_filtered.csv", row.names = T, col.names = T)

betas_filtered_no_metadata <- betas_filtered[-c(1:7),]

betas_filtered_metadata <- betas_filtered[c(1:7),]



#load HM450k CpG map
hm450 <- fread("./APE/TARGET-AML/data files/hm450_map.csv", data.table = F)

hm450 <- as.data.frame(hm450)

colnames(hm450) <- hm450[8,]

rownames(hm450) <- hm450[,1]

hm450 <- hm450[-c(1:8),-1]


# #Merge betas and hm450 map

map_merge <- merge(betas_filtered_no_metadata, hm450,
                          by = 'row.names', all = TRUE)

#filter map_merge
map_merge <- map_merge[which(rowMeans(!is.na(map_merge)) > 0.9),]

rownames(map_merge) <- map_merge[,1]

map_merge <- map_merge[,-1]

map_merge_with_metadata <- betas_filtered_metadata %>%
   bind_rows(map_merge)

# remove CpGs on X and Y chrom
map_no_XY<-map_merge[!(map_merge$CHR=="X" | map_merge$CHR=="Y"),]

map_no_XY <- betas_filtered_metadata %>%
   bind_rows(map_no_XY)

# remove extra clinical data
betas_no_XY <- map_no_XY[,-c(309:340)]

# Remove Fibroblasts
betas_no_XY_no_fibro <-betas_no_XY[,!betas_no_XY[3,]== "Fibroblasts from Bone Marrow Normal"]

betas_no_XY_no_fibro <- betas_no_XY_no_fibro[which(rowMeans(!is.na(betas_no_XY_no_fibro)) > 0.9),]


fwrite(map_merge_with_metadata, "./APE/TARGET-AML/data files/mapped_cpg_betas.csv", col.names =  T, row.names = T)
fwrite(map_no_XY, "./APE/TARGET-AML/data files/mapped_cpg_betas_no_xy.csv", col.names =  T, row.names = T)
fwrite(betas_no_XY_no_fibro, "./APE/TARGET-AML/data files/mapped_cpg_betas_no_xy_no_fibro.csv", col.names =  T, row.names = T)


betas_no_XY_no_fibro_metadata <- betas_no_XY_no_fibro[c(1:7),]
gc()

remove(joined_record_MI_clean_matched_df)
remove(JR_full)
rm(map_merge)
rm(betas_no_XY)
rm(betas_filtered_no_metadata)
remove(betas_merge)
remove(betas_no_na_col)
rm(map_merge_with_metadata)
rm(map_no_XY)
remove(betas_df_t)
remove(betas)


# !! OPTION: Read in Files created in sections 1-7----------------------------------------------------------------------------------------------------------------------
# # can read in csv in place of completing previous steps


# betas_filtered <- fread("./APE/TARGET-AML/data files/betas_filtered.csv")
# betas_filtered <- as.data.frame(betas_filtered)
# rownames(betas_filtered) <- betas_filtered[,1]
# betas_filtered <- betas_filtered[,-1]

# betas_filtered_metadata <- betas_filtered[c(1:7),]

# # load HM450k CpG map
# hm450 <- fread("./APE/TARGET-AML/data files/hm450_map.csv", data.table = F)
# 
# hm450 <- as.data.frame(hm450)
# 
# colnames(hm450) <- hm450[8,]
# 
# rownames(hm450) <- hm450[,1]
# 
# hm450 <- hm450[-c(1:8),-1]


# betas_no_XY_no_fibro <- fread("./APE/TARGET-AML/data files/mapped_cpg_betas_no_xy_no_fibro.csv")
# betas_no_XY_no_fibro <- as.data.frame(betas_no_XY_no_fibro)
# rownames(betas_no_XY_no_fibro) <- betas_no_XY_no_fibro[,1]


# betas_no_XY_no_fibro <- betas_no_XY_no_fibro[,-1]
# betas_no_XY_no_fibro_metadata <- betas_no_XY_no_fibro[c(1:7),]



# 8. Find Number of Clusters ----------------------------------------------------------------------------------------------------------------------
Nb <- betas_no_XY_no_fibro[c(8:1008),]
Nb <- as.data.frame(sapply(Nb, as.numeric)) #<- sapply is here
Nb <- scale(Nb)
Nb <- na.omit(Nb)
set.seed(1234)
pdf(file="./APE/TARGET-AML/data files/cluster plots.pdf")
# Elbow method
plot(fviz_nbclust(Nb, FUNcluster = kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method"))


# Silhouette method
plot(fviz_nbclust(Nb, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method"))

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
plot(fviz_nbclust(Nb, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+ labs(subtitle = "Gap statistic method"))


dev.off()

# number of unique sample types
nsample_type <- apply(betas_filtered[3,], 1, function(x)length(unique(x)))

rm(Nb)

# ^^use this csv to do analysis in qlucore



# 9. Summarized Experiment Object Creation ----------------------------------------------------------------------------------

JRunlist <- as.data.frame(t(betas_no_XY_no_fibro[1:7,]))

data <- as.matrix(betas_no_XY_no_fibro[-c(1:7),], rownames.force = T)

se_data_full <- SummarizedExperiment(assays = data, colData = JRunlist)

names(se_data_full@assays@data@listData)[[1]] <- "betas"

assay(se_data_full)

colData(se_data_full)

cd <- as.data.frame(colData(se_data_full)); rownames(cd) = NULL
cd

fwrite(cd, "./APE/TARGET-AML/data files/clinical_patient_data.csv", row.names = T, col.names = T)
gc()


se_ok_full = (checkLevels(assay(se_data_full), colData(se_data_full)$case_identifier) 
              & checkLevels(assay(se_data_full), colData(se_data_full)$tumor_type) 
              & checkLevels(assay(se_data_full), colData(se_data_full)$sample_type) 
              & checkLevels(assay(se_data_full), colData(se_data_full)$gender) 
              & checkLevels(assay(se_data_full), colData(se_data_full)$race) 
              & checkLevels(assay(se_data_full), colData(se_data_full)$ethnicity) 
              & checkLevels(assay(se_data_full), colData(se_data_full)$vital_status))

sum(se_ok_full)

se_data_full <- se_data_full[se_ok_full,]


colData(se_data_full)$tumor_type <- relevel(factor(colData(se_data_full)$tumor_type), "AML")
colData(se_data_full)$sample_type <- relevel(factor(colData(se_data_full)$sample_type), "Bone Marrow Normal")
colData(se_data_full)$gender <- relevel(factor(colData(se_data_full)$gender), "female")
colData(se_data_full)$race <- relevel(factor(colData(se_data_full)$race), "white")
colData(se_data_full)$ethnicity <- relevel(factor(colData(se_data_full)$ethnicity), "not hispanic or latino")
colData(se_data_full)$vital_status <- relevel(factor(colData(se_data_full)$vital_status), "Alive")


colData(se_data_full)






rm(se_ok_full)
rm(JRunlist)
rm(cd)
rm(data)


gc()
closeAllConnections()


# 10. DIfferential Methylation by Loci ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

smry_full <- DML(se_data_full, ~tumor_type + race + vital_status + sample_type + ethnicity, mc.cores = 12)
closeAllConnections()

smry_full[[1]] 

dmContrasts <- dmContrasts(smry_full)

#Test interpretation
test_result_full = summaryExtractTest(smry_full)
colnames(test_result_full)
head(test_result_full)
dim(test_result_full)
gc()

fwrite(test_result_full, "./APE/TARGET-AML/data files/test_result_full.csv", row.names = T, col.names = T)


# 11. Filtering DML Test Results--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_result_filter <- test_result_full %>% filter(Pval_tumor_typeAML.IF
< 0.0000000942 | Pval_raceasian < 0.0000000942 | Pval_raceblack.or.african.american < 0.0000000942 | Pval_racenative.hawaiian.or.other.pacific.islander < 0.0000000942 |
  Pval_raceother < 0.0000000942 | Pval_raceUnknown < 0.0000000942 | Pval_vital_statusDead < 0.0000000942 | Pval_sample_typeBlood.Derived.Normal < 0.0000000942 |
  Pval_sample_typePrimary.Blood.Derived.Cancer...Bone.Marrow < 0.0000000942 | Pval_sample_typePrimary.Blood.Derived.Cancer...Peripheral.Blood < 0.0000000942 | 
  Pval_sample_typeRecurrent.Blood.Derived.Cancer...Bone.Marrow < 0.0000000942 | Pval_sample_typeRecurrent.Blood.Derived.Cancer...Peripheral.Blood < 0.0000000942 | 
  Pval_ethnicityhispanic.or.latino < 0.0000000942 | Pval_ethnicityUnknown < 0.0000000942 | FPval_tumor_type < 0.0000000942 | FPval_race < 0.0000000942 | 
  FPval_vital_status < 0.0000000942 | FPval_sample_type < 0.0000000942 | FPval_ethnicity < 0.0000000942)

fwrite(test_result_filter, "./APE/TARGET-AML/data files/test_result_filter.csv", row.names = T, col.names = T)

test_result_fp_filter <- test_result_full %>% filter(FPval_tumor_type < 0.0000000942 | FPval_race < 0.0000000942 | FPval_vital_status < 0.0000000942 | 
                                                       FPval_sample_type < 0.0000000942 | FPval_ethnicity < 0.0000000942)
fwrite(test_result_fp_filter, "./APE/TARGET-AML/data files/test_result_fp_filter.csv", row.names = T, col.names = T)


# !!Option: Read in test results and filtered test results-----------------------------------------------------------------------------------------

# test_result_full <- fread("./APE/TARGET-AML/data files/test_result_full.csv")
# 
# test_result_full <- as.data.frame(test_result_full)
# rownames(test_result_full) <- test_result_full[,1]
# 
# test_result_filter <- fread("./APE/TARGET-AML/data/files/test_result_filter.csv")
# 
# test_result_filter <- as.data.frame(test_result_filter)
# rownames(test_result_filter) <- test_result_filter[,1]
# 
# test_result_fp_filter <- fread("./APE/TARGET-AML/data files/test_result_fp_filter.csv")
# 
# test_result_fp_filter <- as.data.frame(test_result_fp_filter)
# rownames(test_result_fp_filter) <- test_result_fp_filter[,1]


# 12. Goodness of Fit-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 

race_test = test_result_full %>% dplyr::filter(FPval_race < 0.0000000942, Eff_race > 0.1) %>%
  select(FPval_race, Eff_race)

sample_type_test = test_result_full %>% dplyr::filter(FPval_sample_type < 0.0000000942, Eff_sample_type > 0.1) %>%
  select(FPval_sample_type, Eff_sample_type)

ethnicity_test = test_result_full %>% dplyr::filter(FPval_ethnicity < 0.0000000942, Eff_ethnicity > 0.1) %>%
  select(FPval_ethnicity, Eff_ethnicity)

vital_status_test = test_result_full %>% dplyr::filter(FPval_vital_status < 0.0000000942, Eff_vital_status > 0.1) %>%
  select(FPval_vital_status, Eff_vital_status)

tumor_type_test = test_result_full %>% dplyr::filter(FPval_tumor_type < 0.0000000942, Eff_tumor_type >0.1) %>%
  select(FPval_tumor_type, Eff_tumor_type)

fwrite(race_test, "./APE/TARGET-AML/data files/race_specific_cpgs.csv", row.names = T, col.names = T)
fwrite(sample_type_test, "./APE/TARGET-AML/data files/sample_type_specific_cpgs.csv", row.names = T, col.names = T)
fwrite(ethnicity_test, "./APE/TARGET-AML/data files/ethnicity_specific_cpgs.csv", row.names = T, col.names = T)
fwrite(vital_status_test, "./APE/TARGET-AML/data files/vital_status_specific_cpgs.csv", row.names = T, col.names = T)
fwrite(tumor_type_test, "./APE/TARGET-AML/data files/tumor_type_specific_cpgs.csv", row.names = T, col.names = T)

rm(race_test)
rm(sample_type_test)
rm(ethnicity_test)
rm(vital_status_test)
rm(tumor_type_test)


test_result_full %>%
  mutate(race_specific =
           ifelse(FPval_race < 0.0000000942 & Eff_race > 0.1, TRUE, FALSE)) %>%
  mutate(sample_type_specific =
           ifelse(FPval_sample_type < 0.0000000942 & Eff_sample_type > 0.1, TRUE, FALSE)) %>%
  mutate(ethnicity_specific = 
           ifelse(FPval_ethnicity < 0.0000000942 & Eff_ethnicity > 0.1, TRUE, FALSE)) %>%
  mutate(vital_status_specific =
           ifelse(FPval_vital_status < 0.0000000942 & Eff_vital_status > 0.1, TRUE, FALSE))%>%
  mutate(tumor_type_specific =
           ifelse(FPval_tumor_type < 0.0000000942 & Eff_tumor_type > 0.1, TRUE, FALSE))%>%
  select(race_specific, sample_type_specific, ethnicity_specific, vital_status_specific, tumor_type_specific) %>% table



# 13. DML:Pairwise comparison-----------------------------------------------------------------------------------------------------------------------------------------
# start pdf device
pdf(file="./APE/TARGET-AML/data files/pairwise comparison plots.pdf")

plot(ggplot(test_result_filter) + geom_point(aes(Est_raceblack.or.african.american, -log10(Pval_raceblack.or.african.american))))


plot(ggplot(test_result_filter) + geom_point(aes(Est_raceasian, -log10(Pval_raceasian))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_racenative.hawaiian.or.other.pacific.islander, -log10(Pval_racenative.hawaiian.or.other.pacific.islander))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_raceUnknown, -log10(Pval_raceUnknown))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_raceother, -log10(Pval_raceother))))



plot(ggplot(test_result_filter) + geom_point(aes(Est_sample_typeBlood.Derived.Normal, -log10(Pval_sample_typeBlood.Derived.Normal))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_sample_typePrimary.Blood.Derived.Cancer...Bone.Marrow, -log10(Pval_sample_typePrimary.Blood.Derived.Cancer...Bone.Marrow))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_sample_typeRecurrent.Blood.Derived.Cancer...Bone.Marrow, -log10(Pval_sample_typeRecurrent.Blood.Derived.Cancer...Bone.Marrow))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_sample_typePrimary.Blood.Derived.Cancer...Peripheral.Blood, -log10(Pval_sample_typePrimary.Blood.Derived.Cancer...Peripheral.Blood))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_sample_typeRecurrent.Blood.Derived.Cancer...Peripheral.Blood, -log10(Pval_sample_typeRecurrent.Blood.Derived.Cancer...Peripheral.Blood))))




plot(ggplot(test_result_filter) + geom_point(aes(Est_tumor_typeAML.IF, -log10(Pval_tumor_typeAML.IF))))




plot(ggplot(test_result_filter) + geom_point(aes(Est_ethnicityUnknown, -log10(Pval_ethnicityUnknown))))

plot(ggplot(test_result_filter) + geom_point(aes(Est_ethnicityhispanic.or.latino, -log10(Pval_ethnicityhispanic.or.latino))))




plot(ggplot(test_result_filter) + geom_point(aes(Est_vital_statusDead, -log10(Pval_vital_statusDead))))


# close graphics device
dev.off()

gc()





# 14. Visualization----------------------------------------------------------------------------------------------------------------------

visualizeRegion('chr19', 10260000, 10380000, betas_no_XY_no_fibro[-c(1:7),] , platform = "HM450")

visualizeGene('DNMT1', betas_no_XY_no_fibro[-c(1:7),] , platform = "HM450")

visualizeProbes(c("cg02382400", "cg03738669"), betas_no_XY_no_fibro[-c(1:7),] , platform = "HM450")



# 15. Final Data Analysis --------------------------------------------------------------------------------------------------------------


classifier <- as.data.frame(fread("./APE/TARGET-AML/data files/classifier_noxy_cpg_no_fibro.csv", header = F))
rownames(classifier) <- classifier[,1]
colnames(classifier) <- c("CpGs of Interest")


CpGs <- merge(classifier, hm450, by = 'row.names', all = F)

genes <- cbind(CpGs$`CpGs of Interest`, CpGs$CHR, CpGs$UCSC_RefGene_Name)

colnames(genes) <- c('CpGs of Interest','CHR','UCSC_RefGene_Name')

genes <- as.data.frame(genes)

rownames(genes) <- genes[,1]
genes <- genes[,-1]

betas_of_interest <- merge(genes, betas_no_XY_no_fibro, by = 'row.names', all = T)

rownames(betas_of_interest) <- betas_of_interest[,1]
betas_of_interest <- betas_of_interest[,-1]

clin_data3 <- cbind(UCSC_RefGene_Name = NA, betas_no_XY_no_fibro_metadata)
clin_data3 <- cbind(CHR = NA, clin_data3)

CpGs_of_Interest <- clin_data3 %>% bind_rows(betas_of_interest)






ncase_ID <- apply(betas_filtered_metadata[1,], 1, function(x)length(unique(x)))
ncase_ID <- as.data.frame(ncase_ID)
ncase_ID <- setDT(ncase_ID, keep.rownames = TRUE)[]
colnames(ncase_ID) <- c("value", "count")

clin_data4 <- betas_filtered_metadata
colnames(clin_data4) <- clin_data4[1,]
clin_data4 <- clin_data4[-1,]
clin_data4 <- clin_data4[-2,]

clin_data5 <- clin_data4[!duplicated(colnames(clin_data4))]


clin_t <- as.data.frame(t(clin_data5))
counts <- rbind(ncase_ID, 
  aggregate(data.frame(count = clin_t$tumor_type), list(value = clin_t$tumor_type), length),
  aggregate(data.frame(count = clin_t$race), list(value = clin_t$race), length), 
  aggregate(data.frame(count = clin_t$gender), list(value = clin_t$gender), length),
  aggregate(data.frame(count = clin_t$ethnicity), list(value = clin_t$ethnicity), length),
  aggregate(data.frame(count = clin_t$vital_status), list(value = clin_t$vital_status), length))

fwrite(CpGs_of_Interest, "./APE/TARGET-AML/data files/cpgs_of_interest.csv", row.names = T, col.names = T)
fwrite(counts, "./APE/TARGET-AML/data files/sample_counts.csv", row.names = T, col.names = T)

rm(genes)
rm(classifier)
rm(CpGs)
rm(betas_of_interest)
rm(clin_data3)
rm(clin_data4)
rm(clin_data5)
rm(ncase_ID)
rm(clin_t)

gc()

