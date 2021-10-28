library(here)

maindir= here()
resdir <- file.path(maindir,"output","bulk_ChIPseq")
if(!dir.exists(resdir)) dir.create(resdir)
source(file.path(maindir,"Scripts","global_var.R"))
library(dplyr)

outputDir = file.path(resdir,"ChIPreChIP_in_vitro")

list_peaks_reChIP = list()
files = c(
    "MM468/ChIPreChIP/fisher_pvalue_K27_K4_peak_0.001_0.15.csv",
    "BT20/ChIPreChIP/fisher_pvalue_K4_K27_peak_0.001_0.15.csv",
    "BT20/ChIPreChIP/fisher_pvalue_K27_K4_peak_0.001_0.15.csv",
    "HCC38/ChIPreChIP/fisher_pvalue_K4_K27_peak_0.001_0.15.csv",
    "HCC38/ChIPreChIP/fisher_pvalue_K27_K4_peak_0.001_0.15.csv",
    "MM231/ChIPreChIP/fisher_pvalue_K27_K4_3m_peak_0.001_0.15.csv",
    "MM231/ChIPreChIP/fisher_pvalue_K27_K4_10m_peak_0.001_0.15.csv")
names(files) =  c("MM468_K27_K4","BT20_K4_K27", "BT20_K27_K4",
                  "HCC38_K4_K27", "HCC38_K27_K4",
                  "MM231_K27_K4_3m", "MM231_K27_K4_10m")

for(i in seq_along(files)){
    file = files[i]
    name_file = names(files)[i]
    list_peaks_reChIP[[name_file]] = read.csv(file.path(resdir,file), row.names = 1)
    list_peaks_reChIP[[name_file]]$Sample = name_file
    list_peaks_reChIP[[name_file]]$Type = "K27/K4"
}

reChIP_df = do.call("rbind", list_peaks_reChIP)

list_peaks_reChIP_IgG = list()
files = c(
    "MM468/ChIPreChIP/fisher_pvalue_K27_IgG_peak_0.001_0.15.csv",
    "BT20/ChIPreChIP/fisher_pvalue_K4_IgG_peak_0.001_0.15.csv",
    "BT20/ChIPreChIP/fisher_pvalue_K27_IgG_peak_0.001_0.15.csv",
    "HCC38/ChIPreChIP/fisher_pvalue_K4_IgG_peak_0.001_0.15.csv",
    "HCC38/ChIPreChIP/fisher_pvalue_K27_IgG_peak_0.001_0.15.csv",
    "MM231/ChIPreChIP/fisher_pvalue_K27_IgG_3m_peak_0.001_0.15.csv",
    "MM231/ChIPreChIP/fisher_pvalue_K27_IgG_10m_peak_0.001_0.15.csv")

names(files) =  c("MM468_K27_IgG","BT20_K4_IgG", "BT20_K27_IgG",
                  "HCC38_K4_IgG", "HCC38_K27_IgG",
                  "MM231_K27_IgG_3m", "MM231_K27_IgG_10m")

for(i in seq_along(files)){
    file = files[i]
    name_file = names(files)[i]
    list_peaks_reChIP_IgG[[name_file]] = read.csv(file.path(resdir,file), row.names = 1)
    list_peaks_reChIP_IgG[[name_file]]$Sample = name_file
    list_peaks_reChIP_IgG[[name_file]]$Type = "IgG"
}

reChIP_IgG_df = do.call("rbind", list_peaks_reChIP_IgG)

reChIP = rbind(reChIP_df, reChIP_IgG_df)
reChIP$Sample.Grouped =  ""
reChIP$Sample.Grouped[reChIP$Sample %in% c("MM468_K27_K4","MM468_K27_IgG")] =  "MM468_K27"
reChIP$Sample.Grouped[reChIP$Sample %in% c("BT20_K4_K27","BT20_K4_IgG")] =  "BT20_K4"
reChIP$Sample.Grouped[reChIP$Sample %in% c("BT20_K27_K4","BT20_K27_IgG")] =  "BT20_K27s"
reChIP$Sample.Grouped[reChIP$Sample %in% c("HCC38_K4_K27","HCC38_K4_IgG")] =  "HCC38_K4"
reChIP$Sample.Grouped[reChIP$Sample %in% c("HCC38_K27_K4","HCC38_K27_IgG")] =  "HCC38_K27"
reChIP$Sample.Grouped[reChIP$Sample %in% c("MM231_K27_K4_3m","MM231_K27_IgG_3m")] =  "MM231_K27_3m"
reChIP$Sample.Grouped[reChIP$Sample %in% c("MM231_K27_K4_10m","MM231_K27_IgG_10m")] =  "MM231_K27_10m"

table(reChIP$Sample.Grouped)

# Barplot of the number of Bivalent peaks in K27/K4 reChIP compared to IgG reChIP
reChIP %>% filter(Significant == TRUE) %>% group_by(Sample, Sample.Grouped, Type) %>% summarise(number_bivalent_peaks = n())

png(file.path(outputDir, "Barplot_number_true_false_positive_bivalents.png"), width = 2000, height = 1600, res =200)
reChIP %>% filter(Significant == TRUE) %>% group_by(Sample, Sample.Grouped, Type) %>% summarise(number_bivalent_peaks = n()) %>%
    ggplot(aes(x = Sample.Grouped, y = number_bivalent_peaks, fill = Type)) +
    geom_bar(stat = "identity") + theme_classic() +
    theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90))
dev.off()
    

# Number of genes in common - Venn Diagrams
list_bivalent_peaks_reChIP = sapply(list_peaks_reChIP, function(x){
    x = x %>% filter(Significant == TRUE) 
    unique(unlist(str_split(x$gene, pattern = ",")))
})

# All together, selecting samples with most peaks per cell line
png(file.path(outputDir, "Venn_diagram_5.png"), width = 2000, height = 1600, res =200)
names(list_bivalent_peaks_reChIP)
plot(venn(list_bivalent_peaks_reChIP[c("MM468_K27_K4", "BT20_K4_K27" ,"HCC38_K4_K27", "HCC38_K27_K4",
                                          "MM231_K27_K4_10m")]))
dev.off()

png(file.path(outputDir, "Venn_diagram_4.png"), width = 2000, height = 1600, res =200)
names(list_bivalent_peaks_reChIP)
plot(venn(list_bivalent_peaks_reChIP[c("MM468_K27_K4", "BT20_K4_K27" ,"HCC38_K4_K27",
                                       "MM231_K27_K4_10m")]))
dev.off()

# Intersection between one way (K27->K4) & the other (K4 -> K27) in each cell line
png(file.path(outputDir, "Venn_diagram_HCC38.png"), width = 2000, height = 1600, res =200)
names(list_bivalent_peaks_reChIP)
plot(venn(list_bivalent_peaks_reChIP[c("HCC38_K4_K27","HCC38_K27_K4")]))
dev.off()

png(file.path(outputDir, "Venn_diagram_BT20.png"), width = 2000, height = 1600, res =200)
names(list_bivalent_peaks_reChIP)
plot(venn(list_bivalent_peaks_reChIP[c("BT20_K4_K27","BT20_K27_K4")]))
dev.off()
    
png(file.path(outputDir, "Venn_diagram_MM231.png"), width = 2000, height = 1600, res =200)
names(list_bivalent_peaks_reChIP)
plot(venn(list_bivalent_peaks_reChIP[c("MM231_K27_K4_10m","MM231_K27_K4_3m")]))
dev.off()


# Correlation of the pathways
pathways_list = list()
files = c(
    "MM468/ChIPreChIP/Enrichment_test_Bivalent_pathways_K27_K4.xlsx",
    "BT20/ChIPreChIP/Enrichment_test_Bivalent_pathways_K4_K27.xlsx",
    "BT20/ChIPreChIP/Enrichment_test_Bivalent_pathways_K27_K4.xlsx",
    "HCC38/ChIPreChIP/Enrichment_test_Bivalent_pathways_K4_K27.xlsx",
    "HCC38/ChIPreChIP/Enrichment_test_Bivalent_pathways_K27_K4.xlsx",
    "MM231/ChIPreChIP/Enrichment_test_Bivalent_pathways_K27_K4_3m.xlsx",
    "MM231/ChIPreChIP/Enrichment_test_Bivalent_pathways_K27_K4_10m.xlsx")
names(files) =  c("MM468_K27_K4","BT20_K4_K27", "BT20_K27_K4",
                  "HCC38_K4_K27", "HCC38_K27_K4",
                  "MM231_K27_K4_3m", "MM231_K27_K4_10m")

for(i in seq_along(files)){
    file = files[i]
    name_file = names(files)[i]
    pathways_list[[name_file]] = readxl::read_xlsx(file.path(resdir,file))
    pathways_list[[name_file]] = pathways_list[[name_file]] %>% filter(Class %in% c("c2_curated", "c5_GO", "hallmark"))
    pathways_list[[name_file]]$Sample = name_file
    pathways_list[[name_file]]$Type = "K27/K4"
}

pathways = do.call("rbind", pathways_list)
pathways$`q-value` = as.numeric(pathways$`q-value`)

# Choice : remove MM231_K27_K4_3m & BT20_K27_K4 which have < 500 targets
# pathways = pathways %>% filter(Sample != "MM231_K27_K4_3m" & Sample != "BT20_K27_K4")

# All three categories
top100 = pathways %>% filter(Nb_of_genes < 500) %>% group_by(Gene.Set) %>% summarise(`q-value` = mean(`q-value`)) %>%
    dplyr::arrange(`q-value`) %>% head(100)

pathways_top100 = pathways[which(pathways$Gene.Set %in% top100$Gene.Set),]
pathways_top100$mean_qvalue = top100$`q-value`[match(pathways_top100$Gene.Set, top100$Gene.Set)]

pdf(file.path(outputDir, "GSA_top10_all.pdf"))
pathways_top100 %>% mutate("-log10(q.value)" = -log10(`q-value`)) %>%
    arrange(mean_qvalue) %>% head(7*10) %>%
    ggplot(aes(x = forcats::fct_reorder(Gene.Set, mean_qvalue), y = `-log10(q.value)`)) +
    geom_point() + theme_classic() +
    geom_point(aes(y = -log10(mean_qvalue)), size = 2, shape = 24, col = "red") +
    theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90)) +
    xlab("")
dev.off()

pdf(file.path(outputDir, "GSA_top20_all.pdf"))
pathways_top100 %>% mutate("-log10(q.value)" = -log10(`q-value`)) %>%
    arrange(mean_qvalue) %>% head(7*20) %>%
    ggplot(aes(x = forcats::fct_reorder(Gene.Set, mean_qvalue, .desc=TRUE), y = `-log10(q.value)`)) +
    geom_point() + theme_classic() +
    geom_point(aes(y = -log10(mean_qvalue)), size = 2, shape = 24, col = "red") +
    theme(text = element_text(size = 12), axis.text.x = element_text(angle = 90)) +
    xlab("") + rotate()
dev.off()

pdf(file.path(outputDir, "GSA_top50_all.pdf"))
pathways_top100 %>% mutate("-log10(q.value)" = -log10(`q-value`)) %>%
    arrange(mean_qvalue) %>% head(7*50) %>%
    ggplot(aes(x = forcats::fct_reorder(Gene.Set, mean_qvalue, .desc=TRUE), y = `-log10(q.value)`)) +
    geom_point() + theme_classic() +
    geom_point(aes(y = -log10(mean_qvalue)), size = 2, shape = 24, col = "red") +
    theme(text = element_text(size = 12), axis.text.x = element_text(angle = 90)) +
    xlab("") + rotate()
dev.off()


pdf(file.path(outputDir, "GSA_top100_all.pdf"))
pathways_top100 %>% mutate("-log10(q.value)" = -log10(`q-value`)) %>%
    arrange(mean_qvalue) %>% head(7*100) %>%
    ggplot(aes(x = forcats::fct_reorder(Gene.Set, mean_qvalue, .desc=TRUE), y = `-log10(q.value)`)) +
    geom_point() + theme_classic() +
    geom_point(aes(y = -log10(mean_qvalue)), size = 2, shape = 24, col = "red") +
    theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
    xlab("") + rotate()
dev.off()

# Hallmark
pathways_h = pathways %>% filter(Class %in% "hallmark")


# All three categories
top100 = pathways_h %>% filter(Nb_of_genes < 500) %>% group_by(Gene.Set) %>% summarise(`q-value` = mean(`q-value`)) %>%
    dplyr::arrange(`q-value`) %>% head(100)

pathways_top100 = pathways_h[which(pathways_h$Gene.Set %in% top100$Gene.Set),]
pathways_top100$mean_qvalue = top100$`q-value`[match(pathways_top100$Gene.Set, top100$Gene.Set)]

pdf(file.path(outputDir, "GSA_top10_all.pdf"))
pathways_top100 %>% mutate("-log10(q.value)" = -log10(`q-value`)) %>%
    arrange(mean_qvalue) %>% head(7*10) %>%
    ggplot(aes(x = forcats::fct_reorder(Gene.Set, mean_qvalue), y = `-log10(q.value)`)) +
    geom_point() + theme_classic() +
    geom_point(aes(y = -log10(mean_qvalue)), size = 2, shape = 24, col = "red") +
    theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90)) +
    xlab("")
dev.off()

pdf(file.path(outputDir, "GSA_top20_all.pdf"))
pathways_top100 %>% mutate("-log10(q.value)" = -log10(`q-value`)) %>%
    arrange(mean_qvalue) %>% head(7*20) %>%
    ggplot(aes(x = forcats::fct_reorder(Gene.Set, mean_qvalue, .desc=TRUE), y = `-log10(q.value)`)) +
    geom_point() + theme_classic() +
    geom_point(aes(y = -log10(mean_qvalue)), size = 2, shape = 24, col = "red") +
    theme(text = element_text(size = 12), axis.text.x = element_text(angle = 90)) +
    xlab("") + rotate()
dev.off()


