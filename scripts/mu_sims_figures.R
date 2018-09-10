library(tidyverse)  # tidy data functions
library(ggthemr)  # more ggplot themes
library(seqinr)  # load fasta files
library(biomaRt)  # dn ds values
library(patchwork)  # combine plots
# packages not used for paper
library(scatterplot3d)  # 3d scatterplot
library(splines)  # bs, b-spline regression
library(ggrepel)  # geom_text_repel
library(SDMTools)  # wt.sd, wt.mean

# # # set default ggthemr
# theme_set(theme_classic(base_size=14))
fresh_edit <- ggthemr("fresh")$palette
fresh_edit$text <- list(inner = '#000000', outer = '#000000')
fresh_edit$line <- list(inner = '#000000', outer = '#000000')
fresh_edit$gridline <- '#FFFFFF'
ggthemr(fresh_edit, text_size=16)

# # # helper functions
process_file <- function(file_name_list) {
  final_file <- NULL
  for (file_name in file_name_list) {
    print(file_name)
    temp_file <- read_tsv(file_name)
    if (is.null(file_name)) {
      final_file <- temp_file
    } else {
      final_file <- rbind(final_file, temp_file)
    }
  }
  # fix floating point issue
  final_file <- final_file %>% 
    mutate(mut_scaled=round(mut_scaled*10)/10)
  final_file <- final_file %>% 
    mutate(constr_cond=factor(
      constr_cond, levels=c("neutral", "grantham", "identity"),
      labels=c("Non-coding", "Intermediate", "Protein-coding")))
  # remove intermediate sims
  final_file <- final_file %>%
    filter(constr_cond %in% c("Non-coding", "Protein-coding")) %>%
    mutate(constr_cond=droplevels(constr_cond))
  return(final_file)
}

# https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
# following method of Cochran (1977) as recommended by Gatz and Smith (1994) for estimate of se of wmean
se_wmean_cochran <- function(x, w) {
  n = length(w)
  w_bar = mean(w)
  x_wbar = weighted.mean(x, w)
  x_wmean = n/((n-1)*sum(w)^2)*(sum((w*x-w_bar*x_wbar)^2)-
    2*x_wbar*sum((w-w_bar)*(w*x-w_bar*x_wbar))+x_wbar^2*sum((w-w_bar)^2))
  return(x_wmean)
}

get_mut_mean <- function(temp_mut) {
  temp_mut_mean <- temp_mut %>% 
    group_by(constr_cond, mut_scaled) %>% 
    summarise(mut_sd=sd(mut_mean),
              esr_sd=sd(esr_mean),
              mut_sdw=wt.sd(mut_mean, seq_length), 
              esr_sdw=wt.sd(esr_mean, seq_length), 
              mut_sew=se_wmean_cochran(mut_mean, seq_length), 
              esr_sew=se_wmean_cochran(esr_mean, seq_length),
              mut_mean=wt.mean(mut_mean, seq_length),  # mean(mut_mean*seq_length)/mean(seq_length),
              esr_mean=wt.mean(esr_mean, seq_length))  # mean(esr_mean*seq_length)/mean(seq_length))
  return(temp_mut_mean)
}

extract_pfvalue <- function(model) {
  temp <- summary(model)
  pvalue <- 1-pf(temp$fstatistic[1], temp$fstatistic[2], temp$fstatistic[3])
  if (pvalue < .Machine$double.eps) {
    pvalue <- "< 2.2e-16"
  } else {
    pvalue <- paste(c("= ", round(pvalue, digits=2)), collapse="")
  }
  return(pvalue)
}

process_kmers <- function(exon_kmers) {
  exon_kmers_temp <- exon_kmers %>% gather(motif, count, -c(1:8)) %>% 
    group_by(constr_cond, mut_scaled, motif) %>% summarise(count=sum(count)) %>% ungroup() %>% 
    group_by(constr_cond, mut_scaled) %>% mutate(prop=count/sum(count)) %>% 
    mutate(CpG=grepl("CG", motif)) %>% ungroup()
  return(exon_kmers_temp)
}

process_delta <- function(exon_kmers_temp) {
  final <- unique(exon_kmers_temp[c("motif", "CpG")])
  final$real_v_identity <- filter(exon_kmers_temp, mut_scaled==0, constr_cond=="Protein-coding")$prop/
    filter(exon_kmers_temp, mut_scaled==max(mut_scaled), constr_cond=="Protein-coding")$prop
  final$real_v_none <- filter(exon_kmers_temp, mut_scaled==0, constr_cond=="Non-coding")$prop/
    filter(exon_kmers_temp, mut_scaled==max(mut_scaled), constr_cond=="Non-coding")$prop
  final$identity_v_none <- filter(exon_kmers_temp, mut_scaled==max(mut_scaled), constr_cond=="Protein-coding")$prop/
    filter(exon_kmers_temp, mut_scaled==max(mut_scaled), constr_cond=="Non-coding")$prop
  final$real_rank <- rank(filter(exon_kmers_temp, mut_scaled==0, constr_cond=="Non-coding")$prop)
  final$none_rank <- rank(filter(exon_kmers_temp, mut_scaled==max(mut_scaled), constr_cond=="Non-coding")$prop)
  final$identity_rank <- rank(filter(exon_kmers_temp, mut_scaled==max(mut_scaled), constr_cond=="Protein-coding")$prop)
  return(as_tibble(final))
}

plot_scale <- 0.8
plot_scale_alt <- 0.9
plot_scale_alt_2 <- 1.0
plot_text <- theme(
  axis.text.x=element_text(size=14), 
  axis.text.y=element_text(size=14),
  legend.text=element_text(size=14))
plot_text_alt <- theme(
  axis.text.x=element_text(size=14), 
  axis.text.y=element_text(size=14),
  legend.text=element_text(size=14))
legend_constraint <- "Constraint"  # "Simulation type"

# # # input datasets
# 4-mers, use ~1/5 of dataset
# for computational convenience
exon_track_mut <- paste(
  "../results/simulations/exon_", c(200*c(0:106)), "_", c(200*c(1:106), 21328), "_track_mut.txt", sep="")
exon_track_kmers1 <- paste(
  "../results/simulations/exon_", c(200*c(0:106)), "_", c(200*c(1:106), 21328), "_track_kmers1.txt", sep="")
exon_track_kmers2 <- paste(
  "../results/simulations/exon_", c(200*c(0:106)), "_", c(200*c(1:106), 21328), "_track_kmers2.txt", sep="")
exon_track_kmers3 <- paste(
  "../results/simulations/exon_", c(200*c(0:106)), "_", c(200*c(1:106), 21328), "_track_kmers3.txt", sep="")
exon_track_kmers4 <- paste(
  "../results/simulations/exon_", c(200*c(0:106)), "_", c(200*c(1:106), 21328), "_track_kmers4.txt", sep="")[1:20]
exon_track_kmers6 <- paste(
  "../results/simulations/exon_", c(200*c(0:106)), "_", c(200*c(1:106), 21328), "_track_kmers6.txt", sep="")

intron_track_mut <- paste(
  "../results/simulations/intron_", c(200*c(0:93)), "_", c(200*c(1:93), 18730), "_track_mut.txt", sep="")
intron_track_kmers2 <- paste(
  "../results/simulations/intron_", c(200*c(0:93)), "_", c(200*c(1:93), 18730), "_track_kmers2.txt", sep="")
intron_track_kmers6 <- paste(
  "../results/simulations/intron_", c(200*c(0:93)), "_", c(200*c(1:93), 18730), "_track_kmers6.txt", sep="")

random_track_mut <- paste(
  "../results/simulations/random_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_mut.txt", sep="")
random_track_kmers1 <- paste(
  "../results/simulations/random_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers1.txt", sep="")
random_track_kmers2 <- paste(
  "../results/simulations/random_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers2.txt", sep="")
random_track_kmers3 <- paste(
  "../results/simulations/random_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers3.txt", sep="")
random_track_kmers4 <- paste(
  "../results/simulations/random_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers4.txt", sep="")[1:20]
random_track_kmers6 <- paste(
  "../results/simulations/random_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers6.txt", sep="")

# mutation bias datasets
random_scaled_0_track_mut <- paste(
  "../results/simulations/random_scaled_0_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_mut.txt", sep="")
random_scaled_0_track_kmers1 <- paste(
  "../results/simulations/random_scaled_0_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers1.txt", sep="")
random_scaled_0_track_kmers2 <- paste(
  "../results/simulations/random_scaled_0_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers2.txt", sep="")
random_scaled_0_track_kmers3 <- paste(
  "../results/simulations/random_scaled_0_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers3.txt", sep="")
random_scaled_0_track_kmers4 <- paste(
  "../results/simulations/random_scaled_0_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers4.txt", sep="")[1:20]

random_scaled_50_track_mut <- paste(
  "../results/simulations/random_scaled_50_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_mut.txt", sep="")
random_scaled_50_track_kmers1 <- paste(
  "../results/simulations/random_scaled_50_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers1.txt", sep="")
random_scaled_50_track_kmers2 <- paste(
  "../results/simulations/random_scaled_50_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers2.txt", sep="")
random_scaled_50_track_kmers3 <- paste(
  "../results/simulations/random_scaled_50_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers3.txt", sep="")
random_scaled_50_track_kmers4 <- paste(
  "../results/simulations/random_scaled_50_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers4.txt", sep="")[1:20]

random_scaled_100_track_mut <- paste(
  "../results/simulations/random_scaled_100_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_mut.txt", sep="")
random_scaled_100_track_kmers1 <- paste(
  "../results/simulations/random_scaled_100_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers1.txt", sep="")
random_scaled_100_track_kmers2 <- paste(
  "../results/simulations/random_scaled_100_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers2.txt", sep="")
random_scaled_100_track_kmers3 <- paste(
  "../results/simulations/random_scaled_100_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers3.txt", sep="")
random_scaled_100_track_kmers4 <- paste(
  "../results/simulations/random_scaled_100_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers4.txt", sep="")[1:20]

random_scaled_200_track_mut <- paste(
  "../results/simulations/random_scaled_200_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_mut.txt", sep="")
random_scaled_200_track_kmers1 <- paste(
  "../results/simulations/random_scaled_200_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers1.txt", sep="")
random_scaled_200_track_kmers2 <- paste(
  "../results/simulations/random_scaled_200_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers2.txt", sep="")
random_scaled_200_track_kmers3 <- paste(
  "../results/simulations/random_scaled_200_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers3.txt", sep="")
random_scaled_200_track_kmers4 <- paste(
  "../results/simulations/random_scaled_200_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers4.txt", sep="")[1:20]

# # # exon simulations
exon_mut <- process_file(exon_track_mut)
write_tsv(exon_mut, "../results/figures/exon_sims-mut_esr_all.txt")
exon_mut_mean <- get_mut_mean(exon_mut)
write_tsv(exon_mut_mean, "../results/figures/exon_sims-mut_esr_mean.txt")
intron_mut <- process_file(intron_track_mut)
write_tsv(intron_mut, "../results/figures/intron_sims-mut_esr_all.txt")

ggplot(filter(exon_mut, mut_scaled==max(mut_scaled)), 
       aes(x=constr_cond, y=mut_mean, color=constr_cond)) + plot_text + 
  geom_violin() + xlab(legend_constraint) + ylab("Final ERM score") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  geom_violin(data=filter(intron_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=3, y=mut_mean, linetype="Introns"), color="#65ADC2") +  # 4
  geom_violin(data=filter(exon_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=4, y=mut_mean, linetype="Exons"), color="#000000") +  # 5
  scale_linetype_manual(name="Genomic type", values=c(2, 2), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(aspect.ratio=1)
ggsave("../results/figures/exon_sims-mut_violin.pdf", scale=plot_scale)

ggplot(filter(exon_mut, mut_scaled==max(mut_scaled)), 
       aes(x=constr_cond, y=mut_mean, color=constr_cond)) + plot_text + 
  geom_violin() + xlab(legend_constraint) + ylab("Final ERM score") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  geom_violin(data=filter(intron_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=3, y=mut_mean, linetype="Introns"), color="#65ADC2") +  # 4
  geom_violin(data=filter(exon_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=4, y=mut_mean, linetype="Exons"), color="#000000") +  # 5
  scale_linetype_manual(name="Genomic type", values=c(2, 2), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(aspect.ratio=1)
ggsave("../results/figures/exon_sims-mut_violin.pdf", scale=plot_scale)

ggplot(filter(exon_mut, mut_scaled==max(mut_scaled)), 
       aes(x=constr_cond, y=esr_mean, color=constr_cond)) + plot_text + 
  geom_violin() + xlab(legend_constraint) + ylab("Final EI score") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  geom_violin(data=filter(intron_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=3, y=esr_mean, linetype="Introns"), color="#65ADC2") +  # 4
  geom_violin(data=filter(exon_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=4, y=esr_mean, linetype="Exons"), color="#000000") +  # 5
  scale_linetype_manual(name="Genomic type", values=c(2, 2), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(aspect.ratio=1)
ggsave("../results/figures/exon_sims-esr_violin.pdf", scale=plot_scale)

# variance in difference in ei score
exon_mut_var <- as_tibble(cbind(
  filter(exon_mut, mut_scaled==0 & constr_cond=="Protein-coding")$run_number, 
  filter(exon_mut, mut_scaled==0 & constr_cond=="Protein-coding")$esr_mean - 
    filter(exon_mut, mut_scaled==max(mut_scaled) & constr_cond=="Non-coding")$esr_mean, 
  filter(exon_mut, mut_scaled==max(mut_scaled) & constr_cond=="Protein-coding")$esr_mean - 
    filter(exon_mut, mut_scaled==max(mut_scaled) & constr_cond=="Non-coding")$esr_mean))
names(exon_mut_var) <- c("run_number", "real_neutral", "identity_neutral")
exon_mut_var <- exon_mut_var %>% 
  mutate(ei_protein_prop=identity_neutral/real_neutral)
ggplot(exon_mut_var, aes(ei_protein_prop)) + plot_text + 
  geom_density() + xlim(c(0, 1)) + theme(aspect.ratio=1) + 
  xlab("Protein-coding - Non-coding)/(Real - Non-coding)")
ggsave("../results/figures/exon_sims-esr_protein_prop.pdf", scale=plot_scale)
mean(exon_mut_var$ei_protein_prop)
median(exon_mut_var$ei_protein_prop)

# mean traj in mut rate and ei score
exon_mean <- exon_mut %>% 
  filter(mut_scaled==0, constr_cond=="Non-coding") %>% 
  summarise(mut_sdw=wt.sd(mut_mean, seq_length), 
            esr_sdw=wt.sd(esr_mean, seq_length), 
            mut_mean=wt.mean(mut_mean, seq_length),  # mean(mut_mean*seq_length)/mean(seq_length), 
            esr_mean=wt.mean(esr_mean, seq_length))  # mean(esr_mean*seq_length)/mean(seq_length))
write_tsv(exon_mean, "../results/figures/exon_sims-exon_mean.txt")

intron_mut <- process_file(intron_track_mut)
intron_mean <- intron_mut %>% 
  filter(mut_scaled==0, constr_cond=="Non-coding") %>% 
  summarise(mut_sdw=wt.sd(mut_mean, seq_length), 
            esr_sdw=wt.sd(esr_mean, seq_length), 
            mut_mean=wt.mean(mut_mean, seq_length),  # mean(mut_mean*seq_length)/mean(seq_length), 
            esr_mean=wt.mean(esr_mean, seq_length))  # mean(esr_mean*seq_length)/mean(seq_length))
write_tsv(intron_mean, "../results/figures/intron_sims-intron_mean.txt")

# https://stackoverflow.com/questions/39119917/how-to-add-a-legend-to-hline
plot_temp_a <- ggplot(exon_mut_mean, aes(x=mut_scaled, y=mut_mean, color=constr_cond)) + plot_text + 
  geom_line() + xlab("Mutations per base") + ylab("Mean ERM rate") + 
  # geom_ribbon(aes(x=mut_scaled, ymin=mut_mean-1*mut_sdw, ymax=mut_mean+1*mut_sdw, fill=constr_cond), 
  #   color=NA, linetype="dashed", alpha=0.2) + 
  geom_hline(aes(yintercept=intron_mean$mut_mean, linetype="Introns"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_mean$mut_mean, linetype="Exons"), color="#000000") + 
  # geom_ribbon(aes(x=mut_scaled, ymin=intron_mean$mut_mean-1*intron_mean$mut_sdw, 
  #   ymax=intron_mean$mut_mean+1*intron_mean$mut_sdw, linetype="Introns"), 
  #   fill="#65ADC2", color=NA, linetype="dashed", alpha=0.2) + 
  # geom_ribbon(aes(x=mut_scaled, ymin=exon_mean$mut_mean-1*exon_mean$mut_sdw, 
  #   ymax=exon_mean$mut_mean+1*exon_mean$mut_sdw, linetype="Exons"), 
  #   fill="#000000", color=NA, linetype="dashed", alpha=0.2) + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1), fill=FALSE) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + theme(aspect.ratio=1)
plot_temp_a
ggsave("../results/figures/exon_sims-mut_mean.pdf", scale=plot_scale)

plot_temp_b <- ggplot(exon_mut_mean, aes(x=mut_scaled, y=esr_mean, color=constr_cond)) + plot_text + 
  geom_line() + xlab("Mutations per base") + ylab("Mean EI score") + 
  # geom_ribbon(aes(x=mut_scaled, ymin=esr_mean-1*esr_sdw, ymax=esr_mean+1*esr_sdw, fill=constr_cond), 
  #   color=NA, linetype="dashed", alpha=0.2) + 
  geom_hline(aes(yintercept=intron_mean$esr_mean, linetype="Introns"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_mean$esr_mean, linetype="Exons"), color="#000000") + 
  # geom_ribbon(aes(x=mut_scaled, ymin=intron_mean$esr_mean-1*intron_mean$esr_sdw, 
  #   ymax=intron_mean$esr_mean+1*intron_mean$esr_sdw, linetype="Introns"), 
  #   fill="#65ADC2", color=NA, linetype="dashed", alpha=0.2) + 
  # geom_ribbon(aes(x=mut_scaled, ymin=exon_mean$esr_mean-1*exon_mean$esr_sdw, 
  #   ymax=exon_mean$esr_mean+1*exon_mean$esr_sdw, linetype="Exons"), 
  #   fill="#000000", color=NA, linetype="dashed", alpha=0.2) + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1), fill=FALSE) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + theme(aspect.ratio=1)
plot_temp_b
ggsave("../results/figures/exon_sims-esr_mean.pdf", scale=plot_scale)

plot_temp <- plot_temp_a + plot_text_alt + guides(color=FALSE, linetype=FALSE) + theme(aspect.ratio=2) + 
  plot_temp_b + plot_text_alt + theme(aspect.ratio=2)
plot_temp
ggsave("../results/figures/exon_sims-mut_mean-esr_mean.pdf", plot=plot_temp, scale=plot_scale_alt_2)

# proportion of CpG
exon_kmers <- process_file(exon_track_kmers2)
exon_kmers_mean <- process_kmers(exon_kmers)
write_tsv(exon_kmers_mean, "../results/figures/exon_sims-CpG_mean.txt")

intron_kmers <- process_file(intron_track_kmers2)
intron_kmers_mean <- process_kmers(intron_kmers)
write_tsv(intron_kmers_mean, "../results/figures/intron_sims-CpG_mean.txt")

exon_CpG <- filter(exon_kmers_mean, CpG, mut_scaled==0, constr_cond=="Non-coding")["prop"]
intron_CpG <- filter(intron_kmers_mean, CpG, mut_scaled==0, constr_cond=="Non-coding")["prop"]
write_tsv(exon_CpG, "../results/figures/exon_sims-exon_CpG.txt")
write_tsv(intron_CpG, "../results/figures/exon_sims-intron_CpG.txt")

ggplot(filter(exon_kmers_mean, CpG), 
       aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=constr_cond)) + plot_text + 
  geom_line() + xlab("Mutations per base") + ylab("Mean CpG content") + 
  geom_hline(aes(yintercept=intron_CpG$prop, linetype="Exons"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_CpG$prop, linetype="Introns"), color="#000000") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + theme(aspect.ratio=1)
ggsave("../results/figures/exon_sims-CpG_mean.pdf", scale=plot_scale)

# mut and esr motif correlations
exon_motif <- exon_kmers %>% 
  gather(motif, count, -c(1:8)) %>% 
  group_by(constr_cond, mut_scaled) %>% 
  mutate(prop=count/sum(count)) %>% 
  mutate(CpG=grepl("CG", motif)) %>% ungroup()

exon_motif <- left_join(exon_motif, exon_mut) %>% drop_na()
exon_motif_corr <- exon_motif %>% 
  group_by(motif, constr_cond) %>% 
  summarise(cor_mut=cor(prop, mut_mean, method="pearson"), 
            cor_esr=cor(prop, esr_mean, method="pearson"))
write_tsv(exon_motif_corr, "../results/figures/exon_sims-mut_esr_motif_corr.txt")

ggplot(exon_motif_corr, aes(x=cor_mut, y=cor_esr, color=constr_cond, label=motif)) + plot_text + 
  geom_point() + geom_text(hjust=-0, vjust=-0.3, show.legend=FALSE) + 
  ylim(c(-0.15, 0.44)) + xlim(c(-0.15, 0.64)) + 
  labs(x="Correlation with ERM rate", y="Correlation with EI score", color=legend_constraint) + 
  theme(aspect.ratio=1)
ggsave("../results/figures/exon_sims-mut_esr_motif_corr.pdf", scale=plot_scale)

# exon_motif <- NULL
# exon_motif_corr <- NULL
# exon_mut <- NULL
# intron_mut <- NULL

# # # # random simulations
# random_mut <- process_file(random_track_mut)
# random_mut_mean <- get_mut_mean(random_mut)
# write_tsv(random_mut_mean, "../results/figures/random_sims-mut_esr_mean.txt")

# ggplot(random_mut_mean, aes(x=mut_scaled, y=mut_mean, color=constr_cond)) + plot_text + 
#   geom_line() + xlab("Mutations per base") + ylab("Mean ERM rate") + 
#   geom_hline(aes(yintercept=intron_mean$mut_mean, linetype="Introns"), color="#65ADC2") + 
#   geom_hline(aes(yintercept=exon_mean$mut_mean, linetype="Exons"), color="#000000") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_manual(name="Genomic mean", values=c(2, 2), 
#                         guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-mut_mean.pdf", scale=plot_scale)

# ggplot(random_mut_mean, aes(x=mut_scaled, y=esr_mean, color=constr_cond)) + plot_text + 
#   geom_line() + xlab("Mutations per base") + ylab("Mean EI score") + 
#   geom_hline(aes(yintercept=intron_mean$esr_mean, linetype="Exons"), color="#65ADC2") + 
#   geom_hline(aes(yintercept=exon_mean$esr_mean, linetype="Introns"), color="#000000") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_manual(name="Genomic mean", values=c(2, 2), 
#                         guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-esr_mean.pdf", scale=plot_scale)

# # proportion of CpG
# random_kmers <- process_file(random_track_kmers2)
# random_kmers_mean <- process_kmers(random_kmers)
# write_tsv(random_kmers_mean, "../results/figures/random_sims-CpG_mean.txt")

# ggplot(filter(random_kmers_mean, CpG), 
#        aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=constr_cond)) + plot_text + 
#   geom_line() + xlab("Mutations per base") + ylab("Mean CpG content") + 
#   geom_hline(aes(yintercept=intron_CpG$prop, linetype="Exons"), color="#65ADC2") + 
#   geom_hline(aes(yintercept=exon_CpG$prop, linetype="Introns"), color="#000000") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_manual(name="Genomic mean", values=c(2, 2), 
#                         guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-CpG_mean.pdf", scale=plot_scale)
# random_mut <- NULL

# # # # mutation bias mut esr
# random_scaled_temp <- process_file(random_scaled_0_track_mut)
# random_scaled_0_mut <- get_mut_mean(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_50_track_mut)
# random_scaled_50_mut <- get_mut_mean(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_100_track_mut)
# random_scaled_100_mut <- get_mut_mean(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_200_track_mut)
# random_scaled_200_mut <- get_mut_mean(random_scaled_temp)
# write_tsv(random_scaled_0_mut, "../results/figures/random_sims-scaled_0-mut_esr_mean.txt")
# write_tsv(random_scaled_50_mut, "../results/figures/random_sims-scaled_50-mut_esr_mean.txt")
# write_tsv(random_scaled_100_mut, "../results/figures/random_sims-scaled_100-mut_esr_mean.txt")
# write_tsv(random_scaled_200_mut, "../results/figures/random_sims-scaled_200-mut_esr_mean.txt")
# random_scaled_temp <- NULL

# random_scaled_all <- rbind(
#   random_scaled_0_mut %>% mutate(scaling_factor="0.0 (None)"), 
#   random_scaled_50_mut %>% mutate(scaling_factor="0.5"), 
#   random_scaled_100_mut %>% mutate(scaling_factor="1.0 (Original)"), 
#   random_scaled_200_mut %>% mutate(scaling_factor="2.0")
# )

# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-mut_esr_mean.txt")

# ggplot(random_scaled_all, 
#        aes(x=mut_scaled, y=mut_mean, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
#   geom_line() + xlab("Mutations per base") + ylab("Mean ERM rate") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_discrete("Mutational bias") + 
#   facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/random_sims-mut_bias-simsall_mut.pdf", scale=plot_scale_alt)

# ggplot(random_scaled_all, 
#        aes(x=mut_scaled, y=esr_mean, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
#   geom_line() + xlab("Mutations per base") + ylab("Mean EI score") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_discrete("Mutational bias") + 
#   facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/random_sims-mut_bias-simsall_esr.pdf", scale=plot_scale_alt)

# # mutation bias motif variance
# # kmers1
# random_scaled_temp <- process_file(random_scaled_0_track_kmers1)
# random_scaled_0_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_50_track_kmers1)
# random_scaled_50_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_100_track_kmers1)
# random_scaled_100_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_200_track_kmers1)
# random_scaled_200_kmers <- process_kmers(random_scaled_temp)
# write_tsv(random_scaled_0_kmers, "../results/figures/random_sims-scaled_0-kmers1_mean.txt")
# write_tsv(random_scaled_50_kmers, "../results/figures/random_sims-scaled_50-kmers1_mean.txt")
# write_tsv(random_scaled_100_kmers, "../results/figures/random_sims-scaled_100-kmers1_mean.txt")
# write_tsv(random_scaled_200_kmers, "../results/figures/random_sims-scaled_200-kmers1_mean.txt")
# random_scaled_temp <- NULL

# random_scaled_all <- rbind(
#   random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
#   random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
#   random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
#   random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
# )
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers1_mean.txt")

# ggplot(random_scaled_all, 
#        aes(x=mut_scaled, y=prop, group=interaction(scaling_factor, constr_cond, motif))) + plot_text + 
#   geom_line() + facet_wrap(~interaction(scaling_factor, constr_cond)) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-mut_bias-kmers1_simsall_motif.pdf", scale=plot_scale)

# random_scaled_var <- random_scaled_all %>% 
#   group_by(mut_scaled, scaling_factor, constr_cond) %>% 
#   summarise(prop_var=var(prop)) %>% ungroup()
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers1_motif_var.txt")

# ggplot(random_scaled_var, 
#        aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
#   geom_line() + xlab("Mutations per base") + ylab("Variance in nucleotide frequencies") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_discrete("Mutational bias") + 
#   facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/random_sims-mut_bias-kmers1_simsall_prop_var.pdf", scale=plot_scale_alt)
 
# # kmers2
# random_scaled_temp <- process_file(random_scaled_0_track_kmers2)
# random_scaled_0_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_50_track_kmers2)
# random_scaled_50_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_100_track_kmers2)
# random_scaled_100_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_200_track_kmers2)
# random_scaled_200_kmers <- process_kmers(random_scaled_temp)
# write_tsv(random_scaled_0_kmers, "../results/figures/random_sims-scaled_0-kmers2_mean.txt")
# write_tsv(random_scaled_50_kmers, "../results/figures/random_sims-scaled_50-kmers2_mean.txt")
# write_tsv(random_scaled_100_kmers, "../results/figures/random_sims-scaled_100-kmers2_mean.txt")
# write_tsv(random_scaled_200_kmers, "../results/figures/random_sims-scaled_200-kmers2_mean.txt")
# random_scaled_temp <- NULL

# random_scaled_all <- rbind(
#   random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
#   random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
#   random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
#   random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
# )
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers2_mean.txt")

# ggplot(random_scaled_all, 
#        aes(x=mut_scaled, y=prop, group=interaction(scaling_factor, constr_cond, motif))) + plot_text + 
#   geom_line() + facet_wrap(~interaction(scaling_factor, constr_cond)) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-mut_bias-kmers2_simsall_motif.pdf", scale=plot_scale)

# random_scaled_var <- random_scaled_all %>% 
#   group_by(mut_scaled, scaling_factor, constr_cond) %>% 
#   summarise(prop_var=var(prop)) %>% ungroup()
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers2_motif_var.txt")

# ggplot(random_scaled_var, 
#        aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
#   geom_line() + xlab("Mutations per base") + ylab("Variance in 2-mer frequencies") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_discrete("Mutational bias") + 
#   facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/random_sims-mut_bias-kmers2_simsall_prop_var.pdf", scale=plot_scale_alt)

# # proportion of CpG
# ggplot(filter(random_scaled_all, CpG), 
#        aes(x=mut_scaled, y=prop, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
#   geom_line() + xlab("Mutations per base") + ylab("Mean CpG content") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_discrete("Mutational bias") + 
#   facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/random_sims-CpG_prop-kmers2_simsall_mean.pdf", scale=plot_scale_alt)
 
# # kmers3
# random_scaled_temp <- process_file(random_scaled_0_track_kmers3)
# random_scaled_0_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_50_track_kmers3)
# random_scaled_50_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_100_track_kmers3)
# random_scaled_100_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_200_track_kmers3)
# random_scaled_200_kmers <- process_kmers(random_scaled_temp)
# write_tsv(random_scaled_0_kmers, "../results/figures/random_sims-scaled_0-kmers3_mean.txt")
# write_tsv(random_scaled_50_kmers, "../results/figures/random_sims-scaled_50-kmers3_mean.txt")
# write_tsv(random_scaled_100_kmers, "../results/figures/random_sims-scaled_100-kmers3_mean.txt")
# write_tsv(random_scaled_200_kmers, "../results/figures/random_sims-scaled_200-kmers3_mean.txt")
# random_scaled_temp <- NULL

# random_scaled_all <- rbind(
#   random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
#   random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
#   random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
#   random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
# )
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers3_mean.txt")

# ggplot(random_scaled_all, 
#        aes(x=mut_scaled, y=prop, group=interaction(scaling_factor, constr_cond, motif))) + plot_text + 
#   geom_line() + facet_wrap(~interaction(scaling_factor, constr_cond)) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-mut_bias-kmers3_simsall_motif.pdf", scale=plot_scale)

# random_scaled_var <- random_scaled_all %>% 
#   group_by(mut_scaled, scaling_factor, constr_cond) %>% 
#   summarise(prop_var=var(prop)) %>% ungroup()
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers3_motif_var.txt")

# ggplot(random_scaled_var, 
#        aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
#   geom_line() + xlab("Mutations per base") + ylab("Variance in 3-mer frequencies") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_discrete("Mutational bias") + 
#   facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/random_sims-mut_bias-kmers3_simsall_prop_var.pdf", scale=plot_scale_alt)
 
# # kmers4
# random_scaled_temp <- process_file(random_scaled_0_track_kmers4)
# random_scaled_0_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_50_track_kmers4)
# random_scaled_50_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_100_track_kmers4)
# random_scaled_100_kmers <- process_kmers(random_scaled_temp)
# random_scaled_temp <- process_file(random_scaled_200_track_kmers4)
# random_scaled_200_kmers <- process_kmers(random_scaled_temp)
# write_tsv(random_scaled_0_kmers, "../results/figures/random_sims-scaled_0-kmers4_mean.txt")
# write_tsv(random_scaled_50_kmers, "../results/figures/random_sims-scaled_50-kmers4_mean.txt")
# write_tsv(random_scaled_100_kmers, "../results/figures/random_sims-scaled_100-kmers4_mean.txt")
# write_tsv(random_scaled_200_kmers, "../results/figures/random_sims-scaled_200-kmers4_mean.txt")
# random_scaled_temp <- NULL

# random_scaled_all <- rbind(
#   random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
#   random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
#   random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
#   random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
# )
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers4_mean.txt")

# ggplot(random_scaled_all, 
#        aes(x=mut_scaled, y=prop, group=interaction(scaling_factor, constr_cond, motif))) + plot_text + 
#   geom_line() + facet_wrap(~interaction(scaling_factor, constr_cond)) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-mut_bias-kmers4_simsall_motif.pdf", scale=plot_scale)

# random_scaled_var <- random_scaled_all %>% 
#   group_by(mut_scaled, scaling_factor, constr_cond) %>% 
#   summarise(prop_var=var(prop)) %>% ungroup()
# write_tsv(random_scaled_all, "../results/figures/random_sims-mut_bias-kmers4_motif_var.txt")

# ggplot(random_scaled_var, 
#        aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
#   geom_line() + xlab("Mutations per base") + ylab("Variance in 4-mer frequencies") + 
#   scale_color_discrete(legend_constraint) + 
#   guides(color=guide_legend(order=-1)) + 
#   scale_linetype_discrete("Mutational bias") + 
#   facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/random_sims-mut_bias-kmers4_simsall_prop_var.pdf", scale=plot_scale_alt)
# # rm(list=ls()[grepl("random_scaled", ls())])

# # # # mutation rate vs equilibrium frequency
# random_track_kmers6 <- paste(
#   "../results/simulations/random_", c(200*c(0:24)), "_", c(200*c(1:25)), "_track_kmers6.txt", sep="")
# hexamers_mut <- read_tsv("../data/hexamers_ERVaveraged.txt") %>% 
#   dplyr::rename(motif=Id, mut_rate=Avg.ht)
# random_scaled_temp <- process_file(random_track_kmers6)
# random_scaled_kmers <- process_kmers(random_scaled_temp) %>% 
#   filter(mut_scaled==max(mut_scaled))
# random_scaled_final <- full_join(
#   random_scaled_kmers, hexamers_mut, by="motif")
# write_tsv(random_scaled_final, "../results/figures/random_sims-mut_rate-prop-kmers6.txt")

# ggplot(filter(random_scaled_final, constr_cond=="Protein-coding")) + 
#   geom_point(aes(x=mut_rate, y=prop, color=CpG)) + plot_text_alt + 
#   xlab("Hexamer ERM rate") + ylab("Equil prop (random sim)") + 
#   labs(color="CpG in hexamer") + theme(aspect.ratio=1)
# ggsave("../results/figures/random_sims-mut_rate-prop-kmers6.pdf", scale=plot_scale_alt)

# delta kmer equilibrium tracking
random_kmers <- process_file(random_track_kmers3)
random_kmers_mean <- process_kmers(random_kmers)

random_kmers_temp <- random_kmers_mean %>% 
  group_by(constr_cond, motif) %>% 
  mutate(prop.diff = prop-lag(prop)) %>% 
  mutate(prop.diff.diff = prop.diff-lag(prop.diff))

random_kmers_temp_mean <- random_kmers_temp %>% 
  group_by(constr_cond, mut_scaled) %>% 
  summarise(prop.diff.mean=mean(abs(prop.diff)), prop.diff.sd=sd(abs(prop.diff)))

ggplot(filter(random_kmers_temp, constr_cond=="Non-coding")) + plot_text + 
  geom_line(aes(x=mut_scaled, y=abs(prop.diff), group=motif)) + scale_y_sqrt() + 
  geom_text(data=filter(random_kmers_temp, constr_cond=="Non-coding" & mut_scaled==1.0), 
    mapping=aes(x=mut_scaled, y=abs(prop.diff), label=motif)) + 
  stat_summary(aes(x=mut_scaled, y=abs(prop.diff), group=1), fun.y=mean, colour="red", geom="line", group=1) # + 
  # geom_ribbon(data=filter(random_kmers_temp_mean, constr_cond=="Non-coding"), 
  #   aes(x=mut_scaled, ymin=prop.diff.mean-1.96*prop.diff.sd, ymax=prop.diff.mean+1.96*prop.diff.sd), 
  #   fill="#d73027", color=NA, size=2/3, linetype="dashed", alpha=0.2)
ggsave("../results/figures/random_sims-abs_delta_motif_non-coding.pdf", scale=plot_scale)

ggplot(filter(random_kmers_temp, constr_cond=="Protein-coding")) + plot_text + 
  geom_line(aes(x=mut_scaled, y=abs(prop.diff), group=motif)) + scale_y_sqrt() + 
  geom_text(data=filter(random_kmers_temp, constr_cond=="Protein-coding" & mut_scaled==1.0), 
    mapping=aes(x=mut_scaled, y=abs(prop.diff), label=motif)) + 
  stat_summary(aes(x=mut_scaled, y=abs(prop.diff), group=1), fun.y=mean, colour="red", geom="line", group=1) # + 
  # geom_ribbon(data=filter(random_kmers_temp_mean, constr_cond=="Non-coding"), 
  #   aes(x=mut_scaled, ymin=prop.diff.mean-1.96*prop.diff.sd, ymax=prop.diff.mean+1.96*prop.diff.sd), 
  #   fill="#d73027", color=NA, size=2/3, linetype="dashed", alpha=0.2)
ggsave("../results/figures/random_sims-abs_delta_motif_protein-coding.pdf", scale=plot_scale)

# # # # save values of lines
# esr_a <- filter(exon_mut_mean, constr_cond=="Protein-coding", mut_scaled==max(mut_scaled))$esr_mean
# esr_c <- filter(exon_mut_mean, constr_cond=="Non-coding", mut_scaled==max(mut_scaled))$esr_mean
# esr_d <- exon_mean$esr_mean
# esr_e <- intron_mean$esr_mean

# mut_a <- filter(exon_mut_mean, constr_cond=="Protein-coding", mut_scaled==max(mut_scaled))$mut_mean
# mut_c <- filter(exon_mut_mean, constr_cond=="Non-coding", mut_scaled==max(mut_scaled))$mut_mean
# mut_d <- exon_mean$mut_mean
# mut_e <- intron_mean$mut_mean

# equil_vals <- as_tibble(
#   as.list(c("esr_real_exon"=esr_d, "esr_real_intron"=esr_e, 
#             "esr_sim_identity"=esr_a, "esr_sim_none"=esr_c, 
#             "mut_real_exon"=mut_d, "mut_real_intron"=mut_e, 
#             "mut_sim_identity"=mut_a, "mut_sim_none"=mut_c,
#             "esr_real-esr_identity"=esr_d-esr_a,
#             "esr_identity-esr_none"=esr_a-esr_c,
#             "esr_real-esr_none"=esr_d-esr_c,
#             "esr_noncoding_prop"=(esr_d-esr_a)/(esr_d-esr_c),
#             "esr_coding_prop"=(esr_a-esr_c)/(esr_d-esr_c)))) %>% 
#   gather()
# write_tsv(equil_vals, "../results/figures/exon_sims-equil_vals-exon_mut_esr_mean.txt")
# exon_mut_mean <- NULL

# # # # kmer plots
# # kmers1
# exon_kmers <- process_file(exon_track_kmers1)
# exon_kmers_mean <- process_kmers(exon_kmers)
# write_tsv(exon_kmers_mean, "../results/figures/kmer_motifs-kmers1_mean.txt")
# ggplot(exon_kmers_mean, 
#        aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=CpG)) + plot_text_alt + 
#   geom_line() + facet_grid(vars(constr_cond), rows=NULL) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/kmer_motifs-kmers1_mean.pdf", scale=plot_scale_alt)

# # kmers2
# exon_kmers <- process_file(exon_track_kmers2)
# exon_kmers_mean <- process_kmers(exon_kmers)
# write_tsv(exon_kmers_mean, "../results/figures/kmer_motifs-kmers2_mean.txt")
# ggplot(exon_kmers_mean, 
#        aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=CpG)) + plot_text_alt + 
#   geom_line() + facet_grid(vars(constr_cond), rows=NULL) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/kmer_motifs-kmers2_mean.pdf", scale=plot_scale_alt)

# # kmers3
# exon_kmers <- process_file(exon_track_kmers3)
# exon_kmers_mean <- process_kmers(exon_kmers)
# write_tsv(exon_kmers_mean, "../results/figures/kmer_motifs-kmers3_mean.txt")
# ggplot(exon_kmers_mean, 
#        aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=CpG)) + plot_text_alt + 
#   geom_line() + facet_grid(vars(constr_cond), rows=NULL) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/kmer_motifs-kmers3_mean.pdf", scale=plot_scale_alt)

# # kmers4
# exon_kmers <- process_file(exon_track_kmers4)  # ~1/5 subset
# exon_kmers_mean <- process_kmers(exon_kmers)
# write_tsv(exon_kmers_mean, "../results/figures/kmer_motifs-kmers4_mean.txt")
# ggplot(exon_kmers_mean, 
#        aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=CpG)) + plot_text_alt + 
#   geom_line() + facet_grid(vars(constr_cond), rows=NULL) + 
#   xlab("Mutations per base") + ylab("Proportion of motifs") + 
#   theme(aspect.ratio=2) + theme(panel.spacing=unit(0.5, "cm"))
# ggsave("../results/figures/kmer_motifs-kmers4_mean.pdf", scale=plot_scale_alt)
# exon_kmers <- NULL
# exon_kmers_mean <- NULL

# # # # dN/dS gene-level analysis
# # get uniprot ids
# # https://www.biostars.org/p/151249/ Accessed 2018-07-10
# gene_mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", 
#                      host="grch37.ensembl.org", version="Ensembl Genes 93")
# hgnc_uniprot <- as_tibble(getBM(attributes=c(
#   "hgnc_symbol", "ensembl_gene_id", "ensembl_peptide_id", "uniprotswissprot"), mart=gene_mart))

# uniprot_exon <- 
#   as_tibble(names(read.fasta(
#     "../data/hg19-unipAliSwissprot-cds_genes-names.txt"))) %>% 
#   mutate(value=gsub("hg19_unipAliSwissprot_", "", value)) %>% 
#   mutate(value=gsub("-.*", "", value)) %>% 
#   dplyr::rename(uniprotswissprot=value) %>% 
#   mutate(run_number=row_number()) %>% 
#   dplyr::select(run_number, uniprotswissprot)

# exon_mut <- read_tsv("../results/figures/exon_sims-mut_esr_all.txt")
# exon_mut <- exon_mut %>% 
#   filter(mut_scaled==0, constr_cond=="Non-coding") %>% 
#   dplyr::select(c("run_number", "seq_length", "mut_mean", "esr_mean"))
# dnds_hgnc <- left_join(uniprot_exon, hgnc_uniprot, by="uniprotswissprot")
# dnds_hgnc <- left_join(dnds_hgnc, exon_mut, by="run_number")
# exon_mut <- NULL

# # get dn, ds, and dn/ds, mmusculus (mouse), mmulatta (rhesus macaques), ptroglodytes (chimpanzees)
# # https://www.biostars.org/p/95252/ Accessed 2017-07-10
# my_attributes <- c(
#   "mmusculus_homolog_canonical_transcript_protein", "mmusculus_homolog_dn", "mmusculus_homolog_ds", 
#   "mmulatta_homolog_canonical_transcript_protein", "mmulatta_homolog_dn", "mmulatta_homolog_ds",
#   "ptroglodytes_homolog_canonical_transcript_protein", "ptroglodytes_homolog_dn", "ptroglodytes_homolog_ds")
# my_rename <- c(
#   "mmusculus_homolog_protein", "mmusculus_homolog_dn", "mmusculus_homolog_ds", 
#   "mmulatta_homolog_protein", "mmulatta_homolog_dn", "mmulatta_homolog_ds", 
#   "ptroglodytes_homolog_protein", "ptroglodytes_homolog_dn", "ptroglodytes_homolog_ds")

# dnds_exon <- getBM(
#   attributes=my_attributes, filters="hgnc_symbol", 
#   values=unique(dnds_hgnc$hgnc_symbol), mart=gene_mart)
# names(dnds_exon) <- my_rename
# dnds_exon <- as_tibble(dnds_exon) %>% 
#   mutate(mmusculus_homolog_dnds=mmusculus_homolog_dn/mmusculus_homolog_ds, 
#          mmulatta_homolog_dnds=mmulatta_homolog_dn/mmulatta_homolog_ds, 
#          ptroglodytes_homolog_dnds=ptroglodytes_homolog_dn/ptroglodytes_homolog_ds)
# write_tsv(dnds_exon, "../results/figures/dnds_exon-species_all.txt")

# dnds_mmusculus <- dnds_exon %>% 
#   dplyr::rename(ensembl_peptide_id=mmusculus_homolog_protein) %>% 
#   right_join(dnds_hgnc, by="ensembl_peptide_id") %>% 
#   dplyr::select(hgnc_symbol, seq_length, mut_mean, esr_mean, 
#                 mmusculus_homolog_dn, mmusculus_homolog_ds, mmusculus_homolog_dnds) %>% 
#   unique() %>% drop_na()
# write_tsv(dnds_mmusculus, "../results/figures/dnds_exon-mmusculus.txt")

# dnds_mmulatta <- dnds_exon %>% 
#   dplyr::rename(ensembl_peptide_id=mmulatta_homolog_protein) %>% 
#   right_join(dnds_hgnc, by="ensembl_peptide_id") %>% 
#   dplyr::select(hgnc_symbol, seq_length, mut_mean, esr_mean, 
#                 mmulatta_homolog_dn, mmulatta_homolog_ds, mmulatta_homolog_dnds) %>% 
#   unique() %>% drop_na()
# write_tsv(dnds_mmulatta, "../results/figures/dnds_exon-mmulatta.txt")

# dnds_ptroglodytes <- dnds_exon %>% 
#   dplyr::rename(ensembl_peptide_id=ptroglodytes_homolog_protein) %>% 
#   right_join(dnds_hgnc, by="ensembl_peptide_id") %>% 
#   dplyr::select(hgnc_symbol, seq_length, mut_mean, esr_mean, 
#                 ptroglodytes_homolog_dn, ptroglodytes_homolog_ds, ptroglodytes_homolog_dnds) %>% 
#   unique() %>% drop_na()
# write_tsv(dnds_ptroglodytes, "../results/figures/dnds_exon-ptroglodytes.txt")

# # heuristics for extremely large or small dn/ds: 
# # https://www.biostars.org/p/84428/, (Villanueva-Canas et al. 2013)
# # dnds mmusculus
# temp <- filter(dnds_mmusculus, 
#                is.finite(mmusculus_homolog_dnds) & mmusculus_homolog_dnds != 0)
# temp <- filter(temp, mmusculus_homolog_dnds<=10 & 
#                  mmusculus_homolog_dn <= 2 & mmusculus_homolog_ds <= 2)

# model <- lm(mut_mean~log10(mmusculus_homolog_dnds), data=temp)
# sink()  # make sure everything is clear
# sink("../results/figures/dnds_exon-mmusculus-lm_mut.txt")
# print(model)
# print(summary(model))
# sink()

# ggplot(temp, aes(x=mmusculus_homolog_dnds, y=mut_mean)) + plot_text + 
#   geom_point(alpha=.5) + geom_density2d(color="black", linetype="dashed") + scale_x_log10() + 
#   xlab(bquote(italic("Homo sapiens")~"-"~italic("Mus musculus")~" dN/dS")) + ylab("Mean ERM rate") + 
#   geom_line(aes(y=fitted(model)), color="black", size=1) + 
#   ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
#   theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmusculus-plot_mut.pdf", scale=plot_scale)

# ggplot(temp, aes(x=cut_number(mmusculus_homolog_dnds, 5), y=mut_mean)) +plot_text + 
#   geom_boxplot(alpha=.5) + xlab(bquote(italic("Homo sapiens")~"-"~italic("Mus musculus")~" dN/dS")) + 
#   ylab("Mean ERM rate") + theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmusculus-boxplot_mut.pdf", scale=plot_scale)

# model <- lm(esr_mean~log10(mmusculus_homolog_dnds), data=temp)
# sink("../results/figures/dnds_exon-mmusculus-lm_esr.txt")
# print(model)
# print(summary(model))
# sink()
# ggplot(temp, aes(x=mmusculus_homolog_dnds, y=esr_mean)) + plot_text + 
#   geom_point(alpha=.5) + geom_density2d(color="black", linetype="dashed") + scale_x_log10() + 
#   xlab(bquote(italic("Homo sapiens")~"-"~italic("Mus musculus")~" dN/dS")) + ylab("Mean EI score") + 
#   geom_line(aes(y=fitted(model)), color="black", size=1) + 
#   ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
#   theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmusculus-plot_esr.pdf", scale=plot_scale)

# ggplot(temp, aes(x=cut_number(mmusculus_homolog_dnds, 5), y=esr_mean)) + plot_text + 
#   geom_boxplot(alpha=.5) + xlab(bquote(italic("Homo sapiens")~"-"~italic("Mus musculus")~" dN/dS")) + 
#   ylab("Mean EI score") + theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmusculus-boxplot_esr.pdf", scale=plot_scale)

# # dnds mmulatta
# temp <- filter(dnds_mmulatta, 
#                is.finite(mmulatta_homolog_dnds) & mmulatta_homolog_dnds != 0)
# temp <- filter(temp, mmulatta_homolog_dnds<=10 & 
#                  mmulatta_homolog_dn <= 2 & mmulatta_homolog_ds <= 2)

# model <- lm(mut_mean~log10(mmulatta_homolog_dnds), data=temp)
# sink("../results/figures/dnds_exon-mmulatta-lm_mut.txt")
# print(model)
# print(summary(model))
# sink()
# ggplot(temp, aes(x=mmulatta_homolog_dnds, y=mut_mean)) + plot_text + 
#   geom_point(alpha=.5) + geom_density2d(color="black", linetype="dashed") + scale_x_log10() + 
#   xlab(bquote(italic("Homo sapiens")~"-"~italic("Macaca mulatta")~" dN/dS")) + ylab("Mean ERM rate") + 
#   geom_line(aes(y=fitted(model)), color="black", size=1) + 
#   ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
#   theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmulatta-plot_mut.pdf", scale=plot_scale)

# ggplot(temp, aes(x=cut_number(mmulatta_homolog_dnds, 5), y=mut_mean)) + plot_text + 
#   geom_boxplot(alpha=.5) + xlab(bquote(italic("Homo sapiens")~"-"~italic("Macaca mulatta")~" dN/dS")) + 
#   ylab("Mean ERM rate") + theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmulatta-boxplot_mut.pdf", scale=plot_scale)

# model <- lm(esr_mean~log10(mmulatta_homolog_dnds), data=temp)
# sink("../results/figures/dnds_exon-mmulatta-lm_esr.txt")
# print(model)
# print(summary(model))
# sink()
# ggplot(temp, aes(x=mmulatta_homolog_dnds, y=esr_mean)) + plot_text + 
#   geom_point(alpha=.5) + geom_density2d(color="black", linetype="dashed") + scale_x_log10() + 
#   xlab(bquote(italic("Homo sapiens")~"-"~italic("Macaca mulatta")~" dN/dS")) + ylab("Mean EI score") + 
#   geom_line(aes(y=fitted(model)), color="black", size=1) + 
#   ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
#   theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmulatta-plot_esr.pdf", scale=plot_scale)

# ggplot(temp, aes(x=cut_number(mmulatta_homolog_dnds, 5), y=esr_mean)) + plot_text + 
#   geom_boxplot(alpha=.5) + xlab(bquote(italic("Homo sapiens")~"-"~italic("Macaca mulatta")~" dN/dS")) + 
#   ylab("Mean EI score") + theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-mmulatta-boxplot_esr.pdf", scale=plot_scale)

# # dnds ptroglodytes
# temp <- filter(dnds_ptroglodytes, 
#                is.finite(ptroglodytes_homolog_dnds) & ptroglodytes_homolog_dnds != 0)
# temp <- filter(temp, ptroglodytes_homolog_dnds<=10 & 
#                  ptroglodytes_homolog_dn <= 2 & ptroglodytes_homolog_ds <= 2)

# model <- lm(mut_mean~log10(ptroglodytes_homolog_dnds), data=temp)
# sink("../results/figures/dnds_exon-ptroglodytes-lm_mut.txt")
# print(model)
# print(summary(model))
# sink()
# ggplot(temp, aes(x=ptroglodytes_homolog_dnds, y=mut_mean)) + plot_text + plot_text + 
#   geom_point(alpha=.5) + geom_density2d(color="black", linetype="dashed") + scale_x_log10() + 
#   xlab(bquote(italic("Homo sapiens")~"-"~italic("Pan troglodytes")~" dN/dS")) + ylab("Mean ERM rate") + 
#   geom_line(aes(y=fitted(model)), color="black", size=1) + 
#   ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
#   theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-ptroglodytes-plot_mut.pdf", scale=plot_scale)

# ggplot(temp, aes(x=cut_number(ptroglodytes_homolog_dnds, 5), y=mut_mean)) + plot_text + 
#   geom_boxplot(alpha=.5) + xlab(bquote(italic("Homo sapiens")~"-"~italic("Pan troglodytes")~" dN/dS")) + 
#   ylab("Mean ERM rate") + theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-ptroglodytes-boxplot_mut.pdf", scale=plot_scale)

# model <- lm(esr_mean~log10(ptroglodytes_homolog_dnds), data=temp)
# sink("../results/figures/dnds_exon-ptroglodytes-lm_esr.txt")
# print(model)
# print(summary(model))
# sink()
# ggplot(temp, aes(x=ptroglodytes_homolog_dnds, y=esr_mean)) + plot_text + 
#   geom_point(alpha=.5) + geom_density2d(color="black", linetype="dashed") + scale_x_log10() + 
#   xlab(bquote(italic("Homo sapiens")~"-"~italic("Pan troglodytes")~" dN/dS")) + ylab("Mean EI score") + 
#   geom_line(aes(y=fitted(model)), color="black", size=1) + 
#   ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
#   theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-ptroglodytes-plot_esr.pdf", scale=plot_scale)

# ggplot(temp, aes(x=cut_number(ptroglodytes_homolog_dnds, 5), y=esr_mean)) + plot_text + 
#   geom_boxplot(alpha=.5) + xlab(bquote(italic("Homo sapiens")~"-"~italic("Pan troglodytes")~" dN/dS")) + 
#   ylab("Mean EI score") + theme(aspect.ratio=1)
# ggsave("../results/figures/dnds_exon-ptroglodytes-boxplot_esr.pdf", scale=plot_scale)
# # seq length binning

# # # # clear objects
# gene_mart <- NULL; uniprot_exon <- NULL; hgnc_uniprot <- NULL; temp <- NULL; model <- NULL
# rm(list=ls()[grepl("dnds", ls())])

# # # exon vs intron enrichment
exon_kmers_temp <- process_file(exon_track_kmers6) %>% 
  filter(mut_scaled==0, constr_cond=="Non-coding") %>% 
  process_kmers()
intron_kmers_temp <- process_file(intron_track_kmers6) %>% 
  filter(mut_scaled==0, constr_cond=="Non-coding") %>% 
  process_kmers()
exon_kmers_temp$count_intron <- intron_kmers_temp$count
exon_kmers_temp$prop_intron <- intron_kmers_temp$prop
exon_kmers_temp <- exon_kmers_temp %>% 
  mutate(exon_v_intron=prop/prop_intron) %>% 
  dplyr::rename(count_exon=count, prop_exon=prop)
 
rescue_ese_hexamers <- read_tsv(
  "../data/ESE.txt", col_names="motif")
hexamer_erv_scores <- read_tsv(
  "../data/chasin_ESS_ESE_numbers.txt", col_names=c("motif", "esr_score", "esr_type"))
hexamer_erv_scores <- hexamer_erv_scores %>% 
  mutate(rescue_ese=ifelse(motif %in% rescue_ese_hexamers$motif, "Y", "N"))
exon_kmers_esr <- full_join(exon_kmers_temp, hexamer_erv_scores, by="motif")
write_tsv(exon_kmers_esr, "../results/figures/motif_enrich-kmers6_exon_intron.txt")

ggplot(exon_kmers_esr, aes(x=prop_intron, y=exon_v_intron, color=CpG)) + plot_text + scale_y_log10() + 
  geom_point(alpha=0.5) + xlab("Proportion in introns") + ylab("Ratio in exons vs introns") + theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-prop_intron.pdf", scale=plot_scale)

ggplot(exon_kmers_esr, aes(x=prop_exon, y=exon_v_intron, color=CpG)) + plot_text + scale_y_log10() + 
  geom_point(alpha=0.5) + xlab("Proportion in exons") + ylab("Ratio in exons vs introns") + theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-prop_exon.pdf", scale=plot_scale)

ggplot(exon_kmers_esr, aes(x=prop_intron, y=prop_exon, color=CpG)) + plot_text + scale_y_log10() + 
  geom_point(alpha=0.5) + xlab("Proportion in introns") + ylab("Proportion in exons") + theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-prop_intron-prop_exon.pdf", scale=plot_scale)

ggplot(exon_kmers_esr, aes(y=prop_intron, x=esr_score, color=CpG)) + plot_text + geom_point(alpha=0.5) + scale_y_log10() + 
  xlab("Hexamer EI score") + ylab("Proportion in introns") + theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-esr_score-prop_intron.pdf", scale=plot_scale)

ggplot(exon_kmers_esr, aes(y=prop_exon, x=esr_score, color=CpG)) + plot_text + geom_point(alpha=0.5) + scale_y_log10() + 
  xlab("Hexamer EI score") + ylab("Proportion in exons") + theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-esr_score-prop_exon.pdf", scale=plot_scale)

# consider esr_score only
model <- lm(log10(exon_v_intron)~esr_score, data=exon_kmers_esr)
sink("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_short.txt")
print(model)
print(summary(model))
sink()
ggplot(exon_kmers_esr, aes(y=exon_v_intron, x=esr_score, color=CpG)) + plot_text + 
  geom_point(alpha=0.5) + scale_y_log10() + 
  xlab("Hexamer EI score") + ylab("Ratio in exon vs intron") + 
  geom_line(aes(y=10**fitted(model)), color="black", size=1) + 
  ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
  theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_short.pdf", scale=plot_scale)

# # 
# consider interaction between esr_score and CpG
model <- lm(log10(exon_v_intron)~esr_score + factor(CpG) + esr_score*factor(CpG), data=exon_kmers_esr)
sink("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_inter.txt")
print(model)
print(summary(model))
sink()
ggplot(exon_kmers_esr, aes(y=exon_v_intron, x=esr_score, color=CpG, group=CpG)) + plot_text + 
  geom_point(alpha=0.5) + scale_y_log10() + 
  xlab("Hexamer EI score") + ylab("Ratio in exon vs intron") + 
  # geom_line(aes(y=10**fitted(model)), size=2.2, alpha=0.5, color="white") + 
  geom_line(aes(y=10**fitted(model)), size=1) + 
  labs(color="CpG in hexamer") + # guides(color=F) + 
  ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
  theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_inter.pdf", scale=plot_scale)

# consider esr_score and CpG with no interaction
model <- lm(log10(exon_v_intron)~esr_score + factor(CpG), data=exon_kmers_esr)
sink("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_long.txt")
print(model)
print(summary(model))
sink()
ggplot(exon_kmers_esr, aes(y=exon_v_intron, x=esr_score, color=CpG, group=CpG)) + plot_text + 
  geom_point(alpha=0.5) + scale_y_log10() + labs(color="CpG in hexamer") + # guides(color=F) + 
  xlab("Hexamer EI score") + ylab("Ratio in exon vs intron") + 
  # geom_line(aes(y=10**fitted(model)), size=2, color="white") + 
  geom_line(aes(y=10**fitted(model)), size=1) + 
  ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
  theme(aspect.ratio=1)
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_long.pdf", scale=plot_scale)

model <- lm(log10(exon_v_intron)~esr_score, data=filter(exon_kmers_esr, CpG))
sink("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_CpGyes.txt")
print(model)
print(summary(model))
sink()

model <- lm(log10(exon_v_intron)~esr_score, data=filter(exon_kmers_esr, !CpG))
sink("../results/figures//motif_enrich-kmers6-exon_intron_ratio-esr_score-lm_CpGno.txt")
print(model)
print(summary(model))
sink()

# 3d plots, 6mers
hexamers_mut <- read_tsv("../data/hexamers_ERVaveraged.txt") %>% 
  dplyr::rename(motif=Id, mut_rate=Avg.ht)
exon_kmers_mut <- full_join(exon_kmers_esr, hexamers_mut, by="motif") %>% 
  mutate(log_exon_v_intron=log10(exon_v_intron))
write_tsv(exon_kmers_mut, "../results/figures/motif_enrich-kmers6_exon_intron_mut.txt")

colors <- c("#65ADC280", "#00000080")  # alpha=80
colors <- colors[as.numeric(exon_kmers_mut$CpG)+1]

pdf("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-mut_rate-3d.pdf")
p <- scatterplot3d(exon_kmers_mut[c("esr_score", "mut_rate", "exon_v_intron")], color=colors, type="h", 
                   xlab="Hexamer EI score", ylab="Mean ERM rate", zlab="Ratio in exon vs intron")
legend(p$xyz.convert(-1.5, 0.03, 35), legend=c("TRUE", "FALSE"), col=c("#000000", "#65ADC2"), pch=1, title="CpG in hexamer", cex=1)
dev.off()

pdf("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-mut_rate-3d_log10.pdf")
p <- scatterplot3d(exon_kmers_mut[c("esr_score", "mut_rate", "log_exon_v_intron")], color=colors, # type="h", 
                   xlab="Hexamer EI score", ylab="Mean ERM rate", zlab="Log10 ratio in exon v intron")
legend(p$xyz.convert(-1.5, 0.03, 1.5), legend=c("TRUE", "FALSE"), col=c("#000000", "#65ADC2"), pch=1, title="CpG in hexamer", cex=1)
dev.off()

# 2d plots, 6mers
# equation from above, with mut_rate bins
# # 
model <- lm(log10(exon_v_intron) ~ esr_score, data=exon_kmers_mut)
ggplot(exon_kmers_mut, aes(y=exon_v_intron, x=esr_score, color=cut_number(round(mut_rate, 3), 5))) + plot_text + 
  geom_point(alpha=0.5) + scale_y_log10() + 
  scale_color_manual(values=c('#d7191c','#fdae61','#fee090','#abd9e9','#2c7bb6')) + 
  xlab("Hexamer EI score") + ylab("Ratio in exon vs intron") + labs(color="ERM rate quintile") + 
  geom_line(aes(y=10**fitted(model)), color="black", size=1) + 
  ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
  theme(aspect.ratio=1) + guides(color = guide_legend(override.aes = list(alpha=1)))
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-mut_rate-lm.pdf", scale=plot_scale)

ggplot(exon_kmers_mut, aes(y=exon_v_intron, x=esr_score, color=cut_number(round(mut_rate, 3), 5), shape=CpG)) + plot_text + 
  geom_point(alpha=0.5) + scale_y_log10() + 
  scale_shape_manual(values=c(1, 4), name="CpG in hexamer") + 
  scale_color_manual(values=c('#d7191c','#fdae61','#fee090','#abd9e9','#2c7bb6')) + 
  xlab("Hexamer EI score") + ylab("Ratio in exon vs intron") + labs(color="ERM rate quintile") + 
  geom_line(aes(y=10**fitted(model)), color="black", size=1) + 
  ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
  theme(aspect.ratio=1) + guides(color = guide_legend(override.aes = list(alpha=1)))
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-mut_rate-lm_old.pdf", scale=plot_scale)

# predict on esr_score, mut_rate, CpG, and interactions
model <- lm(log10(exon_v_intron) ~ esr_score + mut_rate + factor(CpG) +
              esr_score*factor(CpG) + mut_rate*factor(CpG) + esr_score*mut_rate, data=exon_kmers_mut)
sink("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-mut_rate-lm.txt")
print(model)
print(summary(model))
sink()

ggplot(exon_kmers_mut, aes(y=exon_v_intron, x=esr_score, color=cut_number(round(mut_rate, 3), 5), shape=CpG)) + plot_text + 
  geom_point(aes(fitted(model)), alpha=0.5) + scale_y_log10() + 
  scale_shape_manual(values=c(1, 4), name="CpG in hexamer") + 
  scale_color_manual(values=c('#d7191c','#fdae61','#fee090','#abd9e9','#2c7bb6')) +
  xlab("Hexamer EI score") + ylab("Ratio in exon vs intron") + labs(color="ERM rate quintile") + 
  theme(aspect.ratio=1) + guides(color = guide_legend(override.aes = list(alpha=1)))
ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-mut_rate-lm_predict.pdf", scale=plot_scale)

# compare against protein-coding/non-coding from random

random_kmers <- process_file(random_track_kmers6)
random_kmers_mean <- process_kmers(random_kmers)
random_kmers_enrich <- cbind(
  filter(random_kmers_mean, constr_cond=="Protein-coding" & mut_scaled==10)[c("motif", "prop")], 
  filter(random_kmers_mean, constr_cond=="Non-coding" & mut_scaled==10)[c("prop")])
names(random_kmers_enrich) <- c("motif", "prop_protein", "prop_noncoding")
random_kmers_enrich <- as_tibble(random_kmers_enrich) %>% 
  mutate(protein_v_noncoding=prop_protein/prop_noncoding)
kmers_enrich <- full_join(exon_kmers_mut, random_kmers_enrich)
write_tsv(kmers_enrich, "../results/figures/motif_enrich-kmers6_random_protein_noncoding.txt")

model <- lm(log10(exon_v_intron) ~ log10(protein_v_noncoding), data=kmers_enrich)
sink("../results/figures/motif_enrich-kmers6_random_protein_noncoding-lm.txt")
print(model)
print(summary(model))
sink()

ggplot(kmers_enrich) + plot_text + 
  geom_point(aes(x=protein_v_noncoding, y=exon_v_intron, color=cut_width(esr_score, 0.3))) + scale_x_log10() + scale_y_log10() + 
  xlab("Ratio in protein-coding vs non-coding") + ylab("Ratio in exon vs intron") + 
  # geom_line(aes(y=10**fitted(model)), color="black", size=1) + 

  # 
  # 
  # 

  ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
  theme(aspect.ratio=1) 
ggsave("../results/figures/motif_enrich-kmers6_random_protein_noncoding-lm.pdf")

model <- lm(log10(protein_v_noncoding) ~ esr_score, data=kmers_enrich)
ggplot(kmers_enrich, aes(y=protein_v_noncoding, x=esr_score, color=cut_number(round(mut_rate, 3), 5))) + plot_text + 
  geom_point(alpha=0.5) + scale_y_log10() + 
  scale_color_manual(values=c('#d7191c','#fdae61','#fee090','#abd9e9','#2c7bb6')) + 
  xlab("Hexamer EI score") + ylab("Ratio in exon vs intron") + labs(color="ERM rate quintile") + 
  geom_line(aes(y=10**fitted(model)), color="black", size=1) + 
  ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) + 
  theme(aspect.ratio=1) + guides(color = guide_legend(override.aes = list(alpha=1)))
# ggsave("../results/figures/motif_enrich-kmers6-exon_intron_ratio-esr_score-mut_rate-lm.pdf", scale=plot_scale)

model <- lm(esr_score ~ log10(protein_v_noncoding) + log10(exon_v_intron), data=kmers_enrich)
anova(model)
model <- lm(log10(exon_v_intron) ~ log10(protein_v_noncoding) + esr_score, data=kmers_enrich)
anova(model)
model <- lm(log10(exon_v_intron) ~ log10(protein_v_noncoding) + esr_score + mut_rate, data=kmers_enrich)
anova(model)

model <- lm(mut_rate  ~ log10(protein_v_noncoding), data=kmers_enrich)
anova(model)
model <- lm(esr_score ~ log10(protein_v_noncoding), data=kmers_enrich)
anova(model)

model <- lm(log10(exon_v_intron) ~ log10(protein_v_noncoding) + esr_score, data=kmers_enrich)
anova(model)

# 
#
# 

# # # # replot of AA dimers
# aa_dimer_temp <- read_tsv("../data/pair_all_data.txt")
# aa_dimer_temp <- aa_dimer_temp %>% 
#   mutate(emphasis=(Human_EI_Avg>0.3 & Human_enrichment>1)) %>% 
#   mutate(bin_human_ei=cut(Human_EI_Avg, breaks=c(-10:10)/10)) %>% 
#   mutate(bin_human_ei_alt=cut(Human_EI_Avg, breaks=c(-5:5)/5))
# aa_dimer <- aa_dimer_temp %>% 
#   gather(key=enrichment_key, value=enrichment_value, Bacteria_enrichment, Human_enrichment)
# aa_dimer <- aa_dimer %>% 
#   mutate(emphasis=(Human_EI_Avg>0.3 & enrichment_key=="Human_enrichment" & enrichment_value>1)) %>% 
#   mutate(bin_human_ei=cut(Human_EI_Avg, breaks=c(-10:10)/10)) %>% 
#   mutate(bin_human_ei_alt=cut(Human_EI_Avg, breaks=c(-5:5)/5))

# sink("../results/figures/dimer_aa-human_enrichment_long.txt")
# model <- lm(enrichment_value~Human_EI_Avg, 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~ordered(bin_human_ei), 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~bin_human_ei, 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~ordered(bin_human_ei), 
#             data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~ordered(bin_human_ei_alt), 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~bin_human_ei_alt, 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~ordered(bin_human_ei_alt), 
#             data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~poly(Human_EI_Avg, 5), 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# print(model)
# print(summary(model))

# model <- lm(enrichment_value~poly(Human_EI_Avg, 5), 
#             data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment"))
# print(model)
# print(summary(model))

# model <- lm((Human_enrichment-Bacteria_enrichment)~ordered(bin_human_ei), 
#             data=aa_dimer_temp)
# print(model)
# print(summary(model))

# model <- lm((Human_enrichment-Bacteria_enrichment)~bin_human_ei, 
#             data=aa_dimer_temp)
# print(model)
# print(summary(model))

# model <- lm((Human_enrichment-Bacteria_enrichment)~ordered(bin_human_ei_alt), 
#             data=aa_dimer_temp)
# print(model)
# print(summary(model))

# model <- lm((Human_enrichment-Bacteria_enrichment)~bin_human_ei_alt, 
#             data=aa_dimer_temp)
# print(model)
# print(summary(model))

# model <- lm((Human_enrichment-Bacteria_enrichment)~poly(Human_EI_Avg, 5), 
#             data=aa_dimer_temp)
# print(model)
# print(summary(model))
# sink()

# # three main models
# sink("../results/figures/dimer_aa-human_enrichment_short.txt")
# model <- function(d) {
#   if (d > 0) {
#     return(lm(enrichment_value~poly(Human_EI_Avg, d), 
#               data=filter(aa_dimer, enrichment_key=="Human_enrichment")))
#   } else {
#     return(lm(enrichment_value~1, 
#               data=filter(aa_dimer, enrichment_key=="Human_enrichment")))
#   }
# }
# print("Human enrichment")
# print("ANOVA")
# print(anova(model(0), model(1), model(2), model(3), model(4), model(5)))  # quadratic 2
# print("AIC")
# print(AIC(model(0), model(1), model(2), model(3), model(4), model(5)))  # quintic
# print("BIC")
# print(BIC(model(0), model(1), model(2), model(3), model(4), model(5)))  # quadratic
# print(model(3))
# print(summary(model(3)))
# print(model(5))
# print(summary(model(5)))


# model <- function(d) {
#   if (d > 0) {
#     return(lm(enrichment_value~poly(Human_EI_Avg, d), 
#               data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment")))
#   } else {
#     return(lm(enrichment_value~1, 
#               data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment")))
#   }
# }
# print("Bacteria enrichment")
# print("ANOVA")
# print(anova(model(0), model(1), model(2), model(3), model(4), model(5)))  # constant 0
# print("AIC")
# print(AIC(model(0), model(1), model(2), model(3), model(4), model(5)))  # linear
# print("BIC")
# print(BIC(model(0), model(1), model(2), model(3), model(4), model(5)))  # constant
# print(model(3))
# print(summary(model(3)))
# print(model(5))
# print(summary(model(5)))

# model <- function(d) {
#   if (d > 0) {
#     return(lm((Human_enrichment-Bacteria_enrichment)~poly(Human_EI_Avg, d), 
#               data=aa_dimer_temp))
#   } else {
#     return(lm((Human_enrichment-Bacteria_enrichment)~1, 
#               data=aa_dimer_temp))
#   }
# }
# print("Delta enrichment")
# print("ANOVA")
# print(anova(model(0), model(1), model(2), model(3), model(4), model(5)))  # cubic 3
# print("AIC")
# print(AIC(model(0), model(1), model(2), model(3), model(4), model(5)))  # cubic
# print("BIC")
# print(BIC(model(0), model(1), model(2), model(3), model(4), model(5)))  # quadratic
# print(model(3))
# print(summary(model(3)))
# print(model(5))
# print(summary(model(5)))
# sink()

# # dimer plots with prediction
# example_line <- filter(aa_dimer, AA_seq=="ED")
# ggplot(aa_dimer) + plot_text + 
#   geom_point(aes(x=Human_EI_Avg, y=enrichment_value, color=emphasis, shape=enrichment_key)) + 
#   geom_text(data=filter(aa_dimer, emphasis), 
#                   mapping=aes(x=Human_EI_Avg, y=enrichment_value, label=AA_seq), hjust=0.5, vjust=-0.5) + guides(color=FALSE) + 
#   geom_hline(yintercept=1, color="black") + geom_vline(xintercept=0.3, color="black") + 
#   xlab("Average EI score in humans") + ylab("Enrichment (vs. single AA frequencies)") + 
#   scale_shape_manual(name="Organism", values=c(1, 16), labels=c("Bacteria?", "Human")) + theme(aspect.ratio=1)
# ggsave("../results/figures/dimer_aa-human_enrichment.pdf", scale=plot_scale)

# model_1 <- lm(enrichment_value~ordered(bin_human_ei), 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# model_2 <- lm(enrichment_value~ordered(bin_human_ei), 
#             data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment"))
# ggplot(aa_dimer) + plot_text + 
#   geom_point(aes(x=Human_EI_Avg, y=enrichment_value, color=emphasis, shape=enrichment_key)) + 
#   geom_text(data=filter(aa_dimer, emphasis), 
#                   mapping=aes(x=Human_EI_Avg, y=enrichment_value, label=AA_seq), hjust=0.5, vjust=-0.5) + guides(color=FALSE) + 
#   geom_hline(yintercept=1, color="black") + geom_vline(xintercept=0.3, color="black") + 
#   xlab("Average EI score in humans") + ylab("Enrichment (vs. single AA frequencies)") + 
#   scale_shape_manual(name="Organism", values=c(1, 16), labels=c("Bacteria?", "Human")) + 
#   geom_point(data=filter(aa_dimer, enrichment_key=="Human_enrichment"), 
#              mapping=aes(x=Human_EI_Avg, y=fitted(model_1)), color="red") + 
#   geom_point(data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment"), 
#              mapping=aes(x=Human_EI_Avg, y=fitted(model_2)), color="purple") + theme(aspect.ratio=1)
# ggsave("../results/figures/dimer_aa-human_enrichment-predict.pdf", scale=plot_scale)

# # # 
# model_1 <- lm(enrichment_value~poly(Human_EI_Avg, 3),  # poly(Human_EI_Avg, 2), 
#             data=filter(aa_dimer, enrichment_key=="Human_enrichment"))
# model_2 <- lm(enrichment_value~poly(Human_EI_Avg, 3),  # 1,
#             data=filter(aa_dimer, enrichment_key=="Bacteria_enrichment"))

# xpred_1 <- as_tibble(seq(min(aa_dimer_temp$Human_EI_Avg), max(aa_dimer_temp$Human_EI_Avg), 0.001))
# names(xpred_1) <- c("Human_EI_Avg")
# xpred_1$pred_mean <- predict(model_1, xpred_1, se.fit=T)$fit
# xpred_1$pred_se <- predict(model_1, xpred_1, se.fit=T)$se.fit

# xpred_2 <- as_tibble(seq(min(aa_dimer_temp$Human_EI_Avg), max(aa_dimer_temp$Human_EI_Avg), 0.001))
# names(xpred_2) <- c("Human_EI_Avg")
# xpred_2$pred_mean <- predict(model_2, xpred_2, se.fit=T)$fit
# xpred_2$pred_se <- predict(model_2, xpred_2, se.fit=T)$se.fit

# ggplot(aa_dimer) + plot_text + 
#   geom_segment(data=filter(aa_dimer_temp, Human_EI_Avg>0.3), mapping=aes(
#     x=Human_EI_Avg, xend=Human_EI_Avg, y=Human_enrichment, yend=Bacteria_enrichment), linetype="dotted") + 
#   geom_point(aes(x=Human_EI_Avg, y=enrichment_value, size=emphasis, color=enrichment_key), alpha=0.8) + 
#   scale_color_manual(name="Organism", labels=c("Bacteria", "Humans"), values=c("#74add1", "#f46d43")) + 
#   scale_size_manual(values=c(1.5, 3)) + guides(size=FALSE) + 
#   geom_hline(yintercept=1, color="black", linetype="dashed") + geom_vline(xintercept=0.3, color="black", linetype="dashed") + 
#   geom_line(data=xpred_1, aes(x=Human_EI_Avg, y=pred_mean), color="#d73027", size=1) + 
#   geom_ribbon(data=xpred_1, aes(x=Human_EI_Avg, ymin=pred_mean-1.96*pred_se, ymax=pred_mean+1.96*pred_se), 
#     fill="#d73027", color=NA, size=2/3, linetype="dashed", alpha=0.2) + 
#   geom_line(data=xpred_2, aes(x=Human_EI_Avg, y=pred_mean), color="#4575b4", size=1) + 
#     geom_ribbon(data=xpred_2, aes(x=Human_EI_Avg, ymin=pred_mean-1.96*pred_se, ymax=pred_mean+1.96*pred_se), 
#     fill="#4575b4", color=NA, size=2/3, linetype="dashed", alpha=0.2) + 
#   geom_text(data=filter(aa_dimer, emphasis), 
#     mapping=aes(x=Human_EI_Avg, y=enrichment_value, label=AA_seq), hjust=-0, vjust=-0.4) + 
#   xlab("Mean EI score in humans") + ylab("AA pair Obs/Exp ratio") + theme(aspect.ratio=1) + 
#   scale_x_continuous(breaks=c(-0.3, 0, 0.3, 0.6), lim=c(-0.35, 0.65)) + 
#   scale_y_continuous(breaks=c(0.8, 1.0, 1.2, 1.4, 1.6))
# ggsave("../results/figures/dimer_aa-human_enrichment-predict_poly.pdf", scale=plot_scale)

# ggplot(aa_dimer) + plot_text + 
#   geom_boxplot(aes(x=bin_human_ei, y=enrichment_value, fill=enrichment_key)) + theme(aspect.ratio=1) + 
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# ggsave("../results/figures/dimer_aa-human_enrich_boxplot.pdf", scale=plot_scale)

# ggplot(aa_dimer) + plot_text + 
#   geom_boxplot(aes(x=bin_human_ei_alt, y=enrichment_value, fill=enrichment_key)) + theme(aspect.ratio=1) + 
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# ggsave("../results/figures/dimer_aa-human_enrich_boxplot_alt.pdf", scale=plot_scale)

# ggplot(aa_dimer_temp) + plot_text + 
#   geom_boxplot(aes(x=bin_human_ei, y=Human_enrichment-Bacteria_enrichment)) + theme(aspect.ratio=1) + 
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# ggsave("../results/figures/dimer_aa-human_enrich_boxplot_diff.pdf", scale=plot_scale)

# ggplot(aa_dimer_temp) + plot_text + 
#   geom_boxplot(aes(x=bin_human_ei_alt, y=Human_enrichment-Bacteria_enrichment)) + theme(aspect.ratio=1) + 
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# ggsave("../results/figures/dimer_aa-human_enrich_boxplot_diff_alt.pdf", scale=plot_scale)

# ggplot(aa_dimer_temp) + plot_text + 
#   geom_point(aes(x=Human_EI_Avg, y=Human_enrichment-Bacteria_enrichment, color=emphasis)) + theme(aspect.ratio=1) + 
#   geom_text(data=filter(aa_dimer_temp, emphasis), mapping=aes(x=Human_EI_Avg, y=Human_enrichment-Bacteria_enrichment, label=AA_seq), 
#     hjust=0.5, vjust=-0.5) + guides(color=FALSE)
# ggsave("../results/figures/dimer_aa-human_enrichment_diff.pdf", scale=plot_scale)

# model <- lm((Human_enrichment-Bacteria_enrichment)~ordered(bin_human_ei), data=aa_dimer_temp)
# ggplot(aa_dimer_temp) + plot_text + 
#   geom_point(aes(x=Human_EI_Avg, y=Human_enrichment-Bacteria_enrichment, color=emphasis)) + 
#   geom_point(aes(x=Human_EI_Avg, y=fitted(model)), color="red") + theme(aspect.ratio=1) + 
#   geom_hline(yintercept=0, color="black") + geom_vline(xintercept=0.3, color="black") + 
#   geom_text(data=filter(aa_dimer_temp, emphasis), mapping=aes(x=Human_EI_Avg, y=Human_enrichment-Bacteria_enrichment, label=AA_seq), 
#     hjust=0.5, vjust=-0.5) + guides(color=FALSE)
# ggsave("../results/figures/dimer_aa-human_enrichment_diff-predict.pdf", scale=plot_scale)

# # # 
# model <- lm((Human_enrichment-Bacteria_enrichment)~poly(Human_EI_Avg, 3), data=aa_dimer_temp)
# xpred <- as_tibble(seq(min(aa_dimer_temp$Human_EI_Avg), max(aa_dimer_temp$Human_EI_Avg), 0.001))
# names(xpred) <- c("Human_EI_Avg")
# xpred$pred_mean <- predict(model, xpred, se.fit=T)$fit
# xpred$pred_se <- predict(model, xpred, se.fit=T)$se.fit
# ggplot(aa_dimer_temp) + plot_text + 
#   geom_point(aes(x=Human_EI_Avg, y=Human_enrichment-Bacteria_enrichment, size=emphasis), alpha=0.8) + 
#   scale_size_manual(values=c(1.5, 3)) + guides(size=FALSE) + 
#   geom_hline(yintercept=0, color="black", linetype="dashed") + geom_vline(xintercept=0.3, color="black", linetype="dashed") + 
#   geom_line(data=xpred, aes(x=Human_EI_Avg, y=pred_mean), color="red", size=1) + 
#   geom_ribbon(data=xpred, aes(x=Human_EI_Avg, ymin=pred_mean-1.96*pred_se, ymax=pred_mean+1.96*pred_se), 
#     fill="red", color=NA, size=2/3, linetype="dashed", alpha=0.2) + 
#   geom_text(data=filter(aa_dimer_temp, emphasis), 
#     mapping=aes(x=Human_EI_Avg, y=Human_enrichment-Bacteria_enrichment, label=AA_seq), hjust=-0, vjust=-0.4) + 
#   guides(color=FALSE) + theme(aspect.ratio=1) + theme(aspect.ratio=1) + 
#   xlab("Mean EI score in humans") + ylab(bquote("Human Obs/Exp ratio - Bacteria Obs/Exp ratio")) + 
#   scale_x_continuous(breaks=c(-0.3, 0.0, 0.3, 0.6)) + 
#   scale_y_continuous(breaks=c(-0.4, -0.2, 0.0, 0.2, 0.4, 0.6))
# ggsave("../results/figures/dimer_aa-human_enrichment_diff-predict_poly.pdf", scale=plot_scale)

# # # # replot of delta esr
# erv_hep_chasin <- read_tsv("../data/ERVheps_chasin.txt")
# model <- lm(changewt_mtESRav ~ ERVrate, data=erv_hep_chasin)
# ggplot(erv_hep_chasin, aes(x=ERVrate, y=changewt_mtESRav)) + plot_text + 
#   geom_point(alpha=0.5) + 
#   geom_line(aes(y=fitted(model)), color="black", size=1) + 
#   xlab("ERM rate") + ylab(expression(Delta*"EI (wt-mt)")) + 
#   ggtitle(label="", subtitle=bquote(italic(r)~"² = "~.(round(summary(model)$r.squared, digits=4))~", "~italic(p)~"-value "~.(extract_pfvalue(model)))) +   
#   theme(aspect.ratio=1)
# ggsave("../results/figures/erv_hep_chasin-scatterplot_lm.pdf", scale=plot_scale)

# # # # replot of AA singlets
# aa_singlets <- read_tsv("../data/single_aa_scores_sem.txt") %>% 
#   arrange(Avg_EI_score)
# ggplot(aa_singlets) + plot_text + 
#   geom_bar(aes(x=reorder(AA_seq, sort(Avg_EI_score)), y=Avg_EI_score), stat="identity", fill="grey50", color=NA) + 
#   theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Amino Acid") + ylab("Mean EI score") + 
#   geom_errorbar(aes(x=reorder(AA_seq, sort(Avg_EI_score)), ymin=Avg_EI_score-EI_score_SEM, 
#     ymax=Avg_EI_score+EI_score_SEM), width=.2, position=position_dodge(.9), color="black")
# ggsave("../results/figures/singlet_aa-human_ranked_barplot.pdf", scale=0.7)  # plot_scale)

# # # # replot of SS score plot
# ss_scored <- as_tibble(read.table("../data/ss_scored_for_r.txt", header=FALSE))
# ss_scored <- ss_scored %>% mutate(V18_V19=V18-V19)
# model <- lm(V7 ~ poly(V18_V19, 2), data=ss_scored)
# xpred <- as_tibble(seq(min(ss_scored$V18_V19), max(ss_scored$V18_V19), 0.001))
# names(xpred) <- c("V18_V19")
# xpred$pred_mean <- predict(model, xpred, se.fit=T)$fit
# xpred$pred_se <- predict(model, xpred, se.fit=T)$se.fit
# ggplot(ss_scored) + # plot_text + 
#   geom_line(data=xpred, aes(x=V18_V19, y=pred_mean), color="grey50", size=1) + 
#   geom_ribbon(data=xpred, aes(x=V18_V19, ymin=pred_mean-1.96*pred_se, ymax=pred_mean+1.96*pred_se), 
#     fill="grey50", color=NA, size=2/3, linetype="dashed", alpha=0.2) + 
#   # geom_hline(yintercept=0, color="black", linetype="dashed") + 
#   # geom_vline(xintercept=0, color="black", linetype="dashed") + 
#   geom_point(aes(x=V18_V19, y=V7, color=V12), size=2) + 
#   xlab(expression(Delta*"SS score (wt-mt)")) + ylab("M/W splice ratio") + labs(color="Splice site") + 
#   theme(aspect.ratio=0.9, legend.position=c(0.2, 0.2), axis.text.x=element_text(angle=45, hjust=1), 
#     axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_rect(color = "black", size=1, fill=NA))
# ggsave("../results/figures/ssscore-scatterplot_lm.pdf", scale=0.6)
