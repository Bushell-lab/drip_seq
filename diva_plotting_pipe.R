# libraries and functions ####
library(ggplot2)
library(dplyr)
library(gplots)
library(zoo)



cov_filter <- function(jeremy, window, asi_assign=FALSE, bed_file=""){
  bed_table <- read.table(bed_file, header=FALSE, stringsAsFactors=FALSE)
  listchr <- bed$V1
  listco1 <-bed$V2
  listco2 <- bed$V3
  listsite <- bed$V4
  list_activity <- bed$V5
  filtered = NULL
  for (i in 1:99){
    outputintermediate <- jeremy %>%
      filter(chromosome == listchr[i]) %>%
      filter(coordinate > (listco1[i]+4)-window) %>%
      filter(coordinate < (listco2[i]-4)+window) %>%
      mutate(coverage = readcount) %>%
      mutate(site = listsite[i]) %>%
      mutate(around_break = (coordinate - (listco1[i]+4))) %>%
      mutate(activity = list_activity[i]) %>%
      select("around_break", "coverage", "site")
    filtered <- rbind(outputintermediate, filtered)
  }
  if (asi_assign == TRUE){
    hr_sites <- bed_file[(bed_file$V5=="HR"), bed$V4]
    nhej_sites <- bed_file[(bed_file$V5=="NHEJ"), bed$V4]
    low_sites <- bed_file[(bed_file$V6=="Low"), bed$V4]
    med_sites <- bed_file[(bed_file$V6=="Med"), bed$V]
    high_sites <- bed_file[(bed_file$V6=="High"), bed$V4]
    gen_sites <- bed_file[(bed_file$V7=="Genic"), bed$V4]
    ngen_sites <- bed_file[(bed_file$V7=="Non_Genic"), bed$V4]
    filtered[(filtered$site %in% low_sites),"Transcriptional_Activity"]<-"Low"
    filtered[(filtered$site %in% med_sites),"Transcriptional_Activity"]<-"Med"
    filtered[(filtered$site %in% high_sites),"Transcriptional_Activity"]<-"High"
    filtered[(filtered$site %in% gen_sites),"Location"]<-"Genic"
    filtered[(filtered$site %in% ngen_sites),"Location"]<-"Non-genic"
    filtered[(filtered$site %in% hr_sites),"Repair_pathway"]<-"HR"
    filtered[(filtered$site %in% nhej_sites),"Repair_pathway"]<-"NHEJ"
  }
  return(filtered)
}

feature_plot <- function(filtered_cov, feature="Transcriptional_Activity", roll=0){
  featured <- filtered_cov %>%
    group_by_("around_break", feature) %>%
    summarise(coverage = mean(coverage)) %>%
    select_("around_break", feature, "coverage")
  colnames(featured)[2] <- c("feature")
  if (roll>0){
    featured <- featured %>%
      group_by(feature) %>%
      mutate(coverage = rollapply(coverage, roll, FUN=mean, partial=TRUE)) %>%
      select(around_break, feature, coverage)
  }
  return(featured)
}

Pile_Peaker <- function(dfdf, range, asi_assign=FALSE){
  final_out = NULL
  for (i in 1:99){
    outputer <- dfdf %>%
      filter(site == i) %>%
      filter(between(around_break, -(range), range)) %>%
      summarise(peak = mean(coverage)) %>%
      mutate(site = i)
    final_out <- rbind(outputer, final_out)
  }
  if(asi_assign==TRUE){
    hr_sites <- bed_file[(bed_file$V5=="HR"), bed$V4]
    nhej_sites <- bed_file[(bed_file$V5=="NHEJ"), bed$V4]
    low_sites <- bed_file[(bed_file$V6=="Low"), bed$V4]
    med_sites <- bed_file[(bed_file$V6=="Med"), bed$V]
    high_sites <- bed_file[(bed_file$V6=="High"), bed$V4]
    gen_sites <- bed_file[(bed_file$V7=="Genic"), bed$V4]
    ngen_sites <- bed_file[(bed_file$V7=="Non_Genic"), bed$V4]
    final_out[(final_out$site %in% low_sites),"Transcriptional_Activity"]<-"Low"
    final_out[(final_out$site %in% med_sites),"Transcriptional_Activity"]<-"Med"
    final_out[(final_out$site %in% high_sites),"Transcriptional_Activity"]<-"High"
    final_out[(final_out$site %in% gen_sites),"Location"]<-"Genic"
    final_out[(final_out$site %in% ngen_sites),"Location"]<-"Non-genic"
    final_out[(final_out$site %in% hr_sites),"Repair_pathway"]<-"HR"
    final_out[(final_out$site %in% nhej_sites),"Repair_pathway"]<-"NHEJ"
  }
  return(final_out)
}

chipseq_bin <- function(jo, bins=10){
  secfun <- function(sumlist){
    reter <- mean(as.numeric(sumlist))
    return(reter)
  }
  sumfunc <- function(sumframe){
    teret <- lapply(sumframe[1:5], secfun)
    return(teret)
  }
  split_list=NULL
  for (i in 1:99){
    inter <- split(filter(jo, site==i), cut(filter(jo, site==i)$around_break, bins))
    split_list <- rbind(split_list, lapply(inter, FUN=sumfunc))
  }
  output_df <- as.data.frame(t(sapply(split_list, unlist)))
  return(output_df)
}

bin_stat <- function(nicola, bin_width=2000){
  df_size <- nrow(nicola)
  bin_number <- round((df_size/bin_width), digits = 0)
  bin_end = min(nicola$around_break)
  decile = 10
  stats_results <- NULL
  bin_range <- NULL
  for (i in 1:bin_number){
    if(i == round(bin_number*(decile/100))){
      cat(decile, "%  complete\n")
      decile = decile + 10
    } 
    bin_start = bin_end
    bin_end = bin_end + bin_width
    binned_df <- nicola %>%
      filter(between(around_break, bin_start, bin_end)) %>%
      select(coverage, feature)
    bin_range <- c(bin_range, paste(bin_start, bin_end, sep = "-"))
    stats_results <- c(stats_results, (wilcox.test(binned_df$coverage ~ binned_df$feature, paired = T, alternative = c("greater")))$p.value)
  }
  final <- data.frame(bin_range, stats_results)
  return(final)
}

matrix_make <- function(boris){
  listactivity <- boris$activity
  new_df <- unique(boris$around_break)
  for (i in 1:99){
    inter <- filter(arrange(boris, activity), activity==listactivity[i])[2]
    colnames(inter) <- listactivity[i]
    new_df <- cbind(new_df, inter)
  }
  rownames(new_df) <- new_df[,"new_df"]
  new_df[,"new_df"]=NULL
  final <-data.matrix(t(new_df))
  return(final)
}

peak_stats <- function(peak_frame){
  final_result = NULL
  stats_temp <- wilcox.test(filter(peak_frame, Transcriptional_Activity %in% c("High", "Low"))$peak ~ filter(peak_frame, Transcriptional_Activity %in% c("High", "Low"))$Transcriptional_Activity, paired = F, alternative = c("greater"))
  final_result <- rbind(final_result, cbind("High vs Low Transcription", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(filter(peak_frame, Repair_pathway %in% c("NHEJ", "HR"))$peak ~ filter(peak_frame, Repair_pathway %in% c("NHEJ", "HR"))$Repair_pathway, paired = F, alternative = c("greater"))
  final_result <- rbind(final_result, cbind("NHEJ vs HR Repair", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(peak_frame$peak ~ peak_frame$Location, paired = F, alternative = c("greater"))
  final_result <- rbind(final_result, cbind("Genic vs Non-genic", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(filter(peak_frame, Transcriptional_Activity %in% c("High", "Low"))$peak ~ filter(peak_frame, Transcriptional_Activity %in% c("High", "Low"))$Transcriptional_Activity, paired = F, alternative = c("less"))
  final_result <- rbind(final_result, cbind("High vs Low Transcription", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(filter(peak_frame, Repair_pathway %in% c("NHEJ", "HR"))$peak ~ filter(peak_frame, Repair_pathway %in% c("NHEJ", "HR"))$Repair_pathway, paired = F, alternative = c("less"))
  final_result <- rbind(final_result, cbind("NHEJ vs HR Repair", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(peak_frame$peak ~ peak_frame$Location, paired = F, alternative = c("less"))
  final_result <- rbind(final_result, cbind("Genic vs Non-genic", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(filter(peak_frame, Transcriptional_Activity %in% c("High", "Low"))$peak ~ filter(peak_frame, Transcriptional_Activity %in% c("High", "Low"))$Transcriptional_Activity, paired = F, alternative = c("two.sided"))
  final_result <- rbind(final_result, cbind("High vs Low Transcription", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(filter(peak_frame, Repair_pathway %in% c("NHEJ", "HR"))$peak ~ filter(peak_frame, Repair_pathway %in% c("NHEJ", "HR"))$Repair_pathway, paired = F, alternative = c("two.sided"))
  final_result <- rbind(final_result, cbind("NHEJ vs HR Repair", stats_temp$alternative, stats_temp$p.value))
  stats_temp <- wilcox.test(peak_frame$peak ~ peak_frame$Location, paired = F, alternative = c("two.sided"))
  final_result <- rbind(final_result, cbind("Genic vs Non-genic", stats_temp$alternative, stats_temp$p.value))
  colnames(final_result) <- c("comparison", "alternative", "p_value")
  return(final_result)
}



# filter gcov file for breaksites and assign features to each site####

args <- commandArgs()
setwd(args[6])

cat(args[7], "started\n")
master <- cov_filter(read.table(paste0(args[7], "_log.cov"), quote="", col.names=c("chromosome", "coordinate", "readcount"),  fill=TRUE, stringsAsFactors = FALSE), window=5000, asi_assign=TRUE, bed_file=args[8])
cat(args[7], "imported\n")

# plotting ####

plot_temp <- ggplot() +
  geom_line(data=(feature_plot(master, feature="Transcriptional_Activity", roll=200)), aes(x=around_break, y=coverage, col=feature), size = 5) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 1.5),
    axis.line = element_line(colour = 'black', size = 1.5)
  ) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), alpha=0.3, fill="grey") +
  geom_vline(xintercept = 0, size=1, alpha = 1, linetype="dashed") +
  scale_x_continuous(name = "No. of nucleotides from breaksite") +
  scale_y_continuous(name = "Log2 fold change (cut/uncut)") +
  scale_colour_manual(values = c("blue", "grey", "cyan2")) +
  ggsave(filename = paste0(args[7], "_metagene_tran.png"), width = 10, height = 10, dpi=1000, bg="transparent")
cat(args[7], "activity metagene plotted\n")

plot_temp <- ggplot() +
  geom_line(data=(feature_plot(master, feature="Repair_pathway", roll=200)), aes(x=around_break, y=coverage, col=feature), size = 5) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 1.5),
    axis.line = element_line(colour = 'black', size = 1.5)
  ) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), alpha=0.3, fill="grey") +
  geom_vline(xintercept = 0, size=1, alpha = 1, linetype="dashed") +
  scale_x_continuous(name = "No. of nucleotides from breaksite") +
  scale_y_continuous(name = "Log2 fold change (cut/uncut)") +
  scale_colour_manual(values = c("#ffcd3a", "#5ed4ff", "grey")) +
  ggsave(filename = paste0(args[7], "_metagene_rep.png"), width = 10, height = 10, dpi=1000, bg="transparent")
cat(args[7], "repair metagene plotted\n")

plot_temp <- ggplot() +
  geom_line(data=(feature_plot(master, feature="Location", roll=200)), aes(x=around_break, y=coverage, col=feature), size = 5) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 1.5),
    axis.line = element_line(colour = 'black', size = 1.5)
  ) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), alpha=0.3, fill="grey") +
  geom_vline(xintercept = 0, size=1, alpha = 1, linetype="dashed") +
  scale_x_continuous(name = "No. of nucleotides from breaksite") +
  scale_y_continuous(name = "Log2 fold change (cut/uncut)") +
  scale_colour_manual(values = c("green1", "forestgreen")) +
  ggsave(filename = paste0(args[7], "_metagene_loc.png"), width = 10, height = 10, dpi=1000, bg="transparent")
cat(args[7], "location metagene plotted\n")

save.image(file="master_image.RData")

# box-plotting ####

positions = c("Low", "Med", "High")
positions2 = c("HR", "NHEJ", "Unknown")

cat(args[7], "peaking\n")
master_peaks <- Pile_Peaker(master, range=500, asi_assign=TRUE)
write.csv(master_peaks, file=paste0(args[7], "_peaks.csv"), row.names=FALSE)
cat(args[7], "peaked\n")

plot_temp <- ggplot(master_peaks, aes(x=factor(Transcriptional_Activity, levels = positions), y=peak)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_y_continuous() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(colour = 'black', size = 2.5), axis.text = element_text(size =35), axis.title = element_blank())+
  stat_boxplot(geom ='errorbar', lwd = 3.7) +
  geom_boxplot(lwd=3.5) +
  ggsave(filename = paste0(args[7], "_boxplot_tran.png"), width = 14.5, height = 10, dpi=1000, bg="transparent")
cat(args[7], "activity box plotted\n")

plot_temp <- ggplot(master_peaks, aes(x=Location, y=peak)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_y_continuous() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(colour = 'black', size = 2.5), axis.text = element_text(size =35), axis.title = element_blank())+
  stat_boxplot(geom ='errorbar', lwd = 3.7) +
  geom_boxplot(lwd=3.5) +
  ggsave(filename = paste0(args[7], "_boxplot_loc.png"), width = 10, height = 10, dpi=1000, bg="transparent")
cat(args[7], "location box plotted\n")

plot_temp <- ggplot(master_peaks, aes(x=Repair_pathway, y=peak)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_y_continuous() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(colour = 'black', size = 2.5), axis.text = element_text(size =35), axis.title = element_blank())+
  stat_boxplot(geom ='errorbar', lwd = 3.7) +
  geom_boxplot(lwd=3.5) +
  ggsave(filename = paste0(args[7], "_boxplot_rep.png"), width = 10, height = 10, dpi=1000, bg="transparent")
cat(args[7], "repair box plotted\n")

# stat testing ####

cat(args[7], "testing peaks\n")
write.csv(peak_stats(master_peaks), file=paste0(args[7], "_stats.csv"), row.names = FALSE)
cat(args[7], "peaks statted\n")

# heatmapping ####

cat(args[7], "heatmapping\n")
hmcols2<- c((colorRampPalette(c("blue","#00b2ff"))(35)),(colorRampPalette(c("#00b2ff","white"))(35)),(colorRampPalette(c("white","#ff6100"))(35)),(colorRampPalette(c("#ff6100","red"))(35)))
act_temp <- matrix_make(chipseq_bin(master, bins=500))

png(file = paste0(args[7], "_heatmap_tran.png") ,width=1000,height=3000,units="px")
heatmap.2(act_temp[order(as.integer(rownames(act_temp)), decreasing=TRUE),], trace = "none", Colv = FALSE, Rowv = FALSE, scale = "row", density.info = "none", labCol = "", col=hmcols2, dendrogram = "none")
dev.off()

# correlation ####

cat(args[7], "correlating\n")
master_peaks$Relative_Activity <- log2(master_peaks$activity/max(master_peaks$activity))
master_peaks$Relative_Peak <- master_peaks$peak/max(master_peaks$peak)

plot_temp <- ggplot(master_peaks, aes(x=Relative_Peak, y=(Relative_Activity))) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 1.5),
    axis.line = element_line(colour = 'black', size = 1.5)
  ) +
  geom_point(size = 5, aes(col=Transcriptional_Activity))+
  #geom_text(aes(label=site), size = 7) +
  geom_smooth(method='lm',formula=y~x, size = 3, colour = "black") +
  scale_colour_manual(values = c("blue", "cyan2", "slategray2")) +
  ggsave(filename = paste0(args[7], "_correlation_tran.png"), width = 10, height = 10, dpi=1000, bg="transparent")
cat(args[7], "activity correlation plotted\n")

cat(args[7], "testing correlation\n")
cor_data =NULL
cor_temp <- cor.test(x=master_peaks[-5,]$Relative_Activity, y=master_peaks[-5,]$Relative_Peak, method = c("pearson"))
cor_data <- rbind(cor_data, cbind(cor_temp$data.name, cor_temp$estimate, cor_temp$p.value))
colnames(cor_data) <- c("comparison", "correlation", "p_value")
write.csv(cor_data, file=paste0(args[7], "_correlation.csv"), row.names=FALSE)
cat(args[7], "correlation tested\n")

cat(args[7], "finished\n")

