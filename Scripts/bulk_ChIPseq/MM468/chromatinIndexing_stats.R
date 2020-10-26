library(here)

maindir= here()
source(file.path(maindir, "Scripts", "global_var.R"))
source(file.path(maindir, "Scripts", "functions.R"))
inputdir = file.path(maindir,"input","bulk_ChIPseq","MM468")
outputdir = file.path(maindir,"output","bulk_ChIPseq","MM468","ChromatinIndexing")
if(!dir.exists(outputdir)) dir.create(outputdir)

raw_qc=read.csv(file.path(inputdir,"chromatin_indexing_qc.csv"))

raw_qc = raw_qc %>% mutate(Total.reads = Total.reads/2)

stats = read.csv(file.path(inputdir,"chromatin_indexing_demultiplex_stats_bowtie2.csv"),
                 header = F,col.names = c("Sample","index","Pool","index_reads"))

stats = unique(stats)
stats = stats %>% group_by(Pool) %>% mutate(Pool_reads = sum(index_reads)) %>%
  ungroup %>% mutate(index_percent = 100*index_reads / Pool_reads)

stats$Sample = gsub("tumor_","",stats$Sample)
stats$Sample = gsub("MM468BC_","",stats$Sample)
stats = left_join(stats,raw_qc %>% dplyr::select(Sample.ID,Total.reads), by=c("Pool"="Sample.ID"))

colnames(stats)[which(colnames(stats)=="Total.reads")] = "Pool_total_reads"
stats = stats %>% dplyr::mutate(Pool_percent_demult = 100*Pool_reads / Pool_total_reads)
stats = stats %>% dplyr::filter(Pool %in% c("D276C01", "D276C02", "D276C03", "D276C04"))

for(i in paste0("D276C0",1:4)){
  tmp = stats %>% filter(Pool == i)
  tmp = tmp[order(sapply(tmp$index, function(x) which(x == c(
    "A02","C02","D02","E02","F02","G02")))),]
  tmp$index  = as_factor(tmp$index)
  pdf(file.path(outputdir,paste0(i,".pie.pdf")))
  print(pie(tmp$index_percent,radius = 0.5, col = rainbow_hcl(6), cex = 1,
            labels = paste0(tmp$index," - ",tmp$Sample),main = i))
  dev.off()
}

for(i in c(1,3)){
  input = paste0("D276C0",c(i))
  IP = paste0("D276C0",c(i+1))
  tmp_input = stats %>% filter(Pool %in% input)
  tmp_IP = stats %>% filter(Pool %in% IP)
  tmp_input = tmp_input[order(sapply(tmp_input$index, function(x) which(x == c(
    "A02","C02","D02","E02","F02","G02")))),]
  tmp_IP = tmp_IP[order(sapply(tmp_IP$index, function(x) which(x == c(
    "A02","C02","D02","E02","F02","G02")))),]
  
  tmp_input$index  = as_factor(tmp_input$index)
  tmp_IP$index  = as_factor(tmp_IP$index)
  
  vec = tmp_IP$index_percent/tmp_input$index_percent
  png(file.path(outputdir,paste0("ratio_",input,"_",IP,".pie.png")),width = 1000, height = 1000,res = 150)
  par(mar=c(4, 4, 2, 2) + 0.1)
  print(barplot(vec, col = rainbow_hcl(6),las=2))
  dev.off()
}

# Plot enrichment 
for(j in c(1,3)){
  stats. = stats %>% filter(.data$Pool %in% c(paste0("D276C0",j),paste0("D276C0",j+1)))
  if(j==3){
    stats. <- stats. %>% dplyr::arrange(index)
  }
  enr = c()
  for(i in seq(1,nrow(stats.),2) ){
    
    enr = c(enr, stats.$index_reads[i+1] / stats.$index_reads[i])
  }
  enr = as_tibble(enr)
  
  if(j == 1){
    enr$Sample = stats.$Sample[grep("K27",stats.$Sample)]
    enr$Sample = factor(enr$Sample,levels = c("DMSO_D33_K27","GSKJ4_D33_K27", "UNC_D33_K27", 
                                              "UNC_5FU_D33_K27", "5FU_D33_K27"))
  } else{
    enr$Sample = paste0(stats.$Sample[grep("K4",stats.$Sample)],"_", stats.$index[grep("K4",stats.$Sample)])
    enr$Sample = factor(enr$Sample,levels = c("DMSO_D33_K4_A02","DMSO_D33_K4_C02", "5FU_D33_K4_D02", 
                                              "5FU_D33_K4_E02", "5FU_D33_K4_G02"))
  }
  pool = ifelse(j==1,"K27","K4")
  png(file.path(outputdir,paste0("enrichment_pool_",pool,".png")), width = 800,height = 800,res = 250)
  print(ggplot(enr) + geom_point(aes(x=Sample, y = log2(value))) + theme_classic() +
    theme(text = element_text(angle=90)))
  dev.off()
  
  png(file.path(outputdir,paste0("enrichment_pool_",pool,"_wolegend.png")), width = 800,height = 800,res = 250)
  print(ggplot(enr) + geom_point(aes(x=Sample, y = log2(value))) + theme_classic() +
    theme(text = element_blank()))
  dev.off()
  
  write.csv(enr, file = file.path(outputdir,paste0("ratio_ip_input_",pool,".csv")),
            row.names = F)
}

p1 <- ggplot(stats, aes(x=Sample, y=index_percent,fill=index)) + geom_bar(stat = "identity",width = 0.5) + 
  theme_bw() +
  theme(text= element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =20),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  labs(title = "Demultiplexing efficiency per sample (in % of Pool)") + 
  scale_y_continuous(name = "Percentage indexed reads (%)", limits=c(0, 100)) +
  geom_hline(yintercept = 75, linetype="dashed") 

ggsave(file.path(outputdir,"Efficiency_demultiplex_statistics_percent.png"), plot=p1, width=16, height=6)

p2 <-ggplot(stats, aes(x=Sample, y=index_reads,fill=index)) + geom_bar(stat = "identity",width = 0.5) + 
  theme_bw() +
  theme(text= element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =20),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  labs(title = "Demultiplexing efficiency (Reads)") + 
  scale_y_continuous(name = "Reads") +
  geom_hline(yintercept = 75, linetype="dashed") 

ggsave(file.path(outputdir,"Efficiency_demultiplex_statistics_reads.png"), plot=p2, width=16, height=6)


p3 <-ggplot(stats %>% group_by(index) %>% summarise(TotalReads = sum(index_reads)),
           aes(x=index, y=TotalReads,fill=index)) + geom_bar(stat = "identity",width = 0.5) + 
  theme_bw() +
  theme(text= element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =20),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  labs(title = "Demultiplexing efficiency by index") + 
  scale_y_continuous(name = "Total Reads") +
  geom_hline(yintercept = 75, linetype="dashed") 

ggsave(file.path(outputdir, "Efficiency_demultiplex_statistics_per_index.png"), plot=p3, width=12, height=8)


p4 <-ggplot(stats, aes(x=Pool, y=Pool_percent_demult,fill=Pool)) + geom_bar(stat = "identity",width = 0.5) + 
  theme_bw() +
  theme(text= element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =20),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  labs(title = "Demultiplexing efficiency") + 
  scale_y_continuous(name = "Percentage indexed reads (%)", limits=c(0, 100)) +
  geom_hline(yintercept = 75, linetype="dashed") 

ggsave(file.path(outputdir,"Efficiency_demultiplex_statistics_per_pool.png"), plot=p4, width=14, height=8)

raw_qc$Sample.fraction.... = as.numeric(gsub("%","",as.character(raw_qc$Sample.fraction....)))
p5 <- ggplot(raw_qc, aes(x=Sample.ID, y=Sample.fraction....,fill=Sample.ID)) + geom_bar(stat = "identity",width = 0.5) + 
  theme_bw() +
  theme(text= element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =20),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  labs(title = "Sample representation in sequencing") + 
  scale_y_continuous(name = "Sample representation (%)", limits=c(0, 100)) +
  geom_hline(yintercept = 75, linetype="dashed") 

ggsave(file.path(outputdir,"Sample_fraction.png"), plot=p5, width=14, height=8)

