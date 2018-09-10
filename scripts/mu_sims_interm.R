library(tidyverse)  # tidy data functions
library(ggthemr)  # more ggplot themes
library(seqinr)  # load fasta files
library(biomaRt)  # dn ds values
# packages not used for paper
library(scatterplot3d)  # 3d scatterplot
library(splines)  # bs, b-spline regression
library(ggrepel)  # geom_text_repel

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
    # filter(constr_cond %in% c("Non-coding", "Protein-coding")) %>%
    mutate(constr_cond=droplevels(constr_cond))
  return(final_file)
}

get_mut_mean <- function(temp_mut) {
  temp_mut_mean <- temp_mut %>% 
    group_by(constr_cond, mut_scaled) %>% 
    summarise(mut_mean=mean(mut_mean*seq_length)/mean(seq_length), 
              esr_mean=mean(esr_mean*seq_length)/mean(seq_length))
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
plot_scale_alt <- 1.0
plot_text <- theme(
  axis.text.x=element_text(size=14), 
  axis.text.y=element_text(size=14),
  legend.text=element_text(size=14))
plot_text_alt <- theme(
  axis.text.x=element_text(size=12), 
  axis.text.y=element_text(size=12),
  legend.text=element_text(size=12))
legend_constraint <- "Constraint"  # "Simulation type"
manual_color <- scale_color_manual(
  name=legend_constraint, 
  values=c("#000000", "#E84646", "#65ADC2"))

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
write_tsv(exon_mut, "../results/figures/supp_interm_figures/exon_sims-mut_esr_all.txt")
exon_mut_mean <- get_mut_mean(exon_mut)
write_tsv(exon_mut_mean, "../results/figures/supp_interm_figures/exon_sims-mut_esr_mean.txt")
intron_mut <- process_file(intron_track_mut)
write_tsv(intron_mut, "../results/figures/supp_interm_figures/intron_sims-mut_esr_all.txt")

ggplot(filter(exon_mut, mut_scaled==max(mut_scaled)), 
       aes(x=constr_cond, y=mut_mean, color=constr_cond)) + plot_text +  
  geom_violin() + xlab(legend_constraint) + ylab("Final ERM score") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  geom_violin(data=filter(intron_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=4, y=mut_mean, linetype="Introns"), color="#65ADC2") +  # 4
  geom_violin(data=filter(exon_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=5, y=mut_mean, linetype="Exons"), color="#000000") +  # 5
  scale_linetype_manual(name="Genomic type", values=c(2, 2), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/exon_sims-mut_violin.pdf", scale=plot_scale)

ggplot(filter(exon_mut, mut_scaled==max(mut_scaled)), 
       aes(x=constr_cond, y=mut_mean, color=constr_cond)) + plot_text +  
  geom_violin() + xlab(legend_constraint) + ylab("Final ERM score") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  geom_violin(data=filter(intron_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=4, y=mut_mean, linetype="Introns"), color="#65ADC2") +  # 4
  geom_violin(data=filter(exon_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=5, y=mut_mean, linetype="Exons"), color="#000000") +  # 5
  scale_linetype_manual(name="Genomic type", values=c(2, 2), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/exon_sims-mut_violin.pdf", scale=plot_scale)

ggplot(filter(exon_mut, mut_scaled==max(mut_scaled)), 
       aes(x=constr_cond, y=esr_mean, color=constr_cond)) + plot_text +  
  geom_violin() + xlab(legend_constraint) + ylab("Final EI score") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  geom_violin(data=filter(intron_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=4, y=esr_mean, linetype="Introns"), color="#65ADC2") +  # 4
  geom_violin(data=filter(exon_mut, mut_scaled==0 & constr_cond=="Non-coding"), 
              mapping=aes(x=5, y=esr_mean, linetype="Exons"), color="#000000") +  # 5
  scale_linetype_manual(name="Genomic type", values=c(2, 2), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/exon_sims-esr_violin.pdf", scale=plot_scale)

# mean traj in mut rate and ei score
exon_mean <- exon_mut %>% 
  filter(mut_scaled==0, constr_cond=="Non-coding") %>% 
  summarise(mut_mean=mean(mut_mean*seq_length/mean(seq_length)), 
            esr_mean=mean(esr_mean*seq_length/mean(seq_length)))
write_tsv(exon_mean, "../results/figures/supp_interm_figures/exon_sims-exon_mean.txt")

intron_mut <- process_file(intron_track_mut)
intron_mean <- intron_mut %>% 
  filter(mut_scaled==0, constr_cond=="Non-coding") %>% 
  summarise(mut_mean=mean(mut_mean*seq_length)/mean(seq_length), 
            esr_mean=mean(esr_mean*seq_length)/mean(seq_length))
write_tsv(intron_mean, "../results/figures/supp_interm_figures/intron_sims-intron_mean.txt")

# https://stackoverflow.com/questions/39119917/how-to-add-a-legend-to-hline
ggplot(exon_mut_mean, aes(x=mut_scaled, y=mut_mean, color=constr_cond)) + plot_text +  
  geom_line() + xlab("Mutations per base") + ylab("Mean ERM rate") + 
  geom_hline(aes(yintercept=intron_mean$mut_mean, linetype="Introns"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_mean$mut_mean, linetype="Exons"), color="#000000") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/exon_sims-mut_mean.pdf", scale=plot_scale)

ggplot(exon_mut_mean, aes(x=mut_scaled, y=esr_mean, color=constr_cond)) + plot_text +  
  geom_line() + xlab("Mutations per base") + ylab("Mean EI score") + 
  geom_hline(aes(yintercept=intron_mean$esr_mean, linetype="Introns"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_mean$esr_mean, linetype="Exons"), color="#000000") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/exon_sims-esr_mean.pdf", scale=plot_scale)

# proportion of CpG
exon_kmers <- process_file(exon_track_kmers2)
exon_kmers_mean <- process_kmers(exon_kmers)
write_tsv(exon_kmers_mean, "../results/figures/supp_interm_figures/exon_sims-CpG_mean.txt")

intron_kmers <- process_file(intron_track_kmers2)
intron_kmers_mean <- process_kmers(intron_kmers)
write_tsv(intron_kmers_mean, "../results/figures/supp_interm_figures/intron_sims-CpG_mean.txt")

exon_CpG <- filter(exon_kmers_mean, CpG, mut_scaled==0, constr_cond=="Non-coding")["prop"]
intron_CpG <- filter(intron_kmers_mean, CpG, mut_scaled==0, constr_cond=="Non-coding")["prop"]
write_tsv(exon_CpG, "../results/figures/supp_interm_figures/exon_sims-exon_CpG.txt")
write_tsv(intron_CpG, "../results/figures/supp_interm_figures/exon_sims-intron_CpG.txt")

ggplot(filter(exon_kmers_mean, CpG), 
       aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=constr_cond)) + plot_text +  
  geom_line() + xlab("Mutations per base") + ylab("Mean CpG content") + 
  geom_hline(aes(yintercept=intron_CpG$prop, linetype="Exons"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_CpG$prop, linetype="Introns"), color="#000000") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/exon_sims-CpG_mean.pdf", scale=plot_scale)

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
write_tsv(exon_motif_corr, "../results/figures/supp_interm_figures/exon_sims-mut_esr_motif_corr.txt")

ggplot(exon_motif_corr, aes(x=cor_mut, y=cor_esr, color=constr_cond, label=motif)) + plot_text +  
  geom_point() + geom_text(hjust=-0, vjust=-0.3, show.legend=FALSE) + 
  ylim(c(-0.15, 0.44)) + xlim(c(-0.15, 0.64)) + 
  labs(x="Correlation with ERM rate", y="Correlation with EI score", color=legend_constraint) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/exon_sims-mut_esr_motif_corr.pdf", scale=plot_scale)

exon_motif <- NULL
exon_motif_corr <- NULL
exon_mut <- NULL
intron_mut <- NULL

# # # random simulations
random_mut <- process_file(random_track_mut)
random_mut_mean <- get_mut_mean(random_mut)
write_tsv(random_mut_mean, "../results/figures/supp_interm_figures/random_sims-mut_esr_mean.txt")

ggplot(random_mut_mean, aes(x=mut_scaled, y=mut_mean, color=constr_cond)) + plot_text +  
  geom_line() + xlab("Mutations per base") + ylab("Mean ERM rate") + 
  geom_hline(aes(yintercept=intron_mean$mut_mean, linetype="Introns"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_mean$mut_mean, linetype="Exons"), color="#000000") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/random_sims-mut_mean.pdf", scale=plot_scale)

ggplot(random_mut_mean, aes(x=mut_scaled, y=esr_mean, color=constr_cond)) + plot_text +  
  geom_line() + xlab("Mutations per base") + ylab("Mean EI score") + 
  geom_hline(aes(yintercept=intron_mean$esr_mean, linetype="Exons"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_mean$esr_mean, linetype="Introns"), color="#000000") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
  theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/random_sims-esr_mean.pdf", scale=plot_scale)

# proportion of CpG
random_kmers <- process_file(random_track_kmers2)
random_kmers_mean <- process_kmers(random_kmers)
write_tsv(random_kmers_mean, "../results/figures/supp_interm_figures/random_sims-CpG_mean.txt")

ggplot(filter(random_kmers_mean, CpG), 
       aes(x=mut_scaled, y=prop, group=interaction(constr_cond, motif), color=constr_cond)) + plot_text +  
  geom_line() + xlab("Mutations per base") + ylab("Mean CpG content") + 
  geom_hline(aes(yintercept=intron_CpG$prop, linetype="Exons"), color="#65ADC2") + 
  geom_hline(aes(yintercept=exon_CpG$prop, linetype="Introns"), color="#000000") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_manual(name="Genomic mean", values=c(2, 2), labels=c("Introns", "Exons"), 
                        guide=guide_legend(override.aes=list(color=c("#65ADC2", "#000000")))) + 
theme(aspect.ratio=1) + manual_color
ggsave("../results/figures/supp_interm_figures/random_sims-CpG_mean.pdf", scale=plot_scale)
random_mut <- NULL

# # # mutation bias mut esr
random_scaled_temp <- process_file(random_scaled_0_track_mut)
random_scaled_0_mut <- get_mut_mean(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_50_track_mut)
random_scaled_50_mut <- get_mut_mean(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_100_track_mut)
random_scaled_100_mut <- get_mut_mean(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_200_track_mut)
random_scaled_200_mut <- get_mut_mean(random_scaled_temp)
write_tsv(random_scaled_0_mut, "../results/figures/supp_interm_figures/random_sims-scaled_0-mut_esr_mean.txt")
write_tsv(random_scaled_50_mut, "../results/figures/supp_interm_figures/random_sims-scaled_50-mut_esr_mean.txt")
write_tsv(random_scaled_100_mut, "../results/figures/supp_interm_figures/random_sims-scaled_100-mut_esr_mean.txt")
write_tsv(random_scaled_200_mut, "../results/figures/supp_interm_figures/random_sims-scaled_200-mut_esr_mean.txt")
random_scaled_temp <- NULL

random_scaled_all <- rbind(
  random_scaled_0_mut %>% mutate(scaling_factor="0.0 (None)"), 
  random_scaled_50_mut %>% mutate(scaling_factor="0.5"), 
  random_scaled_100_mut %>% mutate(scaling_factor="1.0 (Original)"), 
  random_scaled_200_mut %>% mutate(scaling_factor="2.0")
)

write_tsv(random_scaled_all, "../results/figures/supp_interm_figures/random_sims-mut_bias-mut_esr_mean.txt")

ggplot(random_scaled_all, 
       aes(x=mut_scaled, y=mut_mean, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
  geom_line() + xlab("Mutations per base") + ylab("Mean ERM rate") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_discrete("Mutational bias") + 
  facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
  theme(aspect.ratio=2) + manual_color + theme(panel.spacing=unit(0.5, "cm"))
ggsave("../results/figures/supp_interm_figures/random_sims-mut_bias-simsall_mut.pdf", scale=plot_scale_alt)

ggplot(random_scaled_all, 
       aes(x=mut_scaled, y=esr_mean, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
  geom_line() + xlab("Mutations per base") + ylab("Mean EI score") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_discrete("Mutational bias") + 
  facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
  theme(aspect.ratio=2) + manual_color + theme(panel.spacing=unit(0.5, "cm"))
ggsave("../results/figures/supp_interm_figures/random_sims-mut_bias-simsall_esr.pdf", scale=plot_scale_alt)

# mutation bias motif variance
# kmers1
random_scaled_temp <- process_file(random_scaled_0_track_kmers1)
random_scaled_0_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_50_track_kmers1)
random_scaled_50_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_100_track_kmers1)
random_scaled_100_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_200_track_kmers1)
random_scaled_200_kmers <- process_kmers(random_scaled_temp)
write_tsv(random_scaled_0_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_0-kmers1_mean.txt")
write_tsv(random_scaled_50_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_50-kmers1_mean.txt")
write_tsv(random_scaled_100_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_100-kmers1_mean.txt")
write_tsv(random_scaled_200_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_200-kmers1_mean.txt")
random_scaled_temp <- NULL

random_scaled_all <- rbind(
  random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
  random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
  random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
  random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
)

random_scaled_var <- random_scaled_all %>% 
  group_by(mut_scaled, scaling_factor, constr_cond) %>% 
  summarise(prop_var=var(prop)) %>% ungroup()
write_tsv(random_scaled_all, "../results/figures/supp_interm_figures/random_sims-mut_bias-kmers1_motif_var.txt")

ggplot(random_scaled_var, 
       aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
  geom_line() + xlab("Mutations per base") + ylab("Variance in nucleotide frequencies") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_discrete("Mutational bias") + 
  facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
  theme(aspect.ratio=2) + manual_color + theme(panel.spacing=unit(0.5, "cm"))
ggsave("../results/figures/supp_interm_figures/random_sims-mut_bias-kmers1_simsall_prop_var.pdf", scale=plot_scale_alt)
 
# kmers2
random_scaled_temp <- process_file(random_scaled_0_track_kmers2)
random_scaled_0_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_50_track_kmers2)
random_scaled_50_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_100_track_kmers2)
random_scaled_100_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_200_track_kmers2)
random_scaled_200_kmers <- process_kmers(random_scaled_temp)
write_tsv(random_scaled_0_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_0-kmers2_mean.txt")
write_tsv(random_scaled_50_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_50-kmers2_mean.txt")
write_tsv(random_scaled_100_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_100-kmers2_mean.txt")
write_tsv(random_scaled_200_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_200-kmers2_mean.txt")
random_scaled_temp <- NULL

random_scaled_all <- rbind(
  random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
  random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
  random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
  random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
)

random_scaled_var <- random_scaled_all %>% 
  group_by(mut_scaled, scaling_factor, constr_cond) %>% 
  summarise(prop_var=var(prop)) %>% ungroup()
write_tsv(random_scaled_all, "../results/figures/supp_interm_figures/random_sims-mut_bias-kmers2_motif_var.txt")

ggplot(random_scaled_var, 
       aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
  geom_line() + xlab("Mutations per base") + ylab("Variance in 2-mer frequencies") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_discrete("Mutational bias") + 
  facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
  theme(aspect.ratio=2) + manual_color + theme(panel.spacing=unit(0.5, "cm"))
ggsave("../results/figures/supp_interm_figures/random_sims-mut_bias-kmers2_simsall_prop_var.pdf", scale=plot_scale_alt)

# proportion of CpG
ggplot(filter(random_scaled_all, CpG), 
       aes(x=mut_scaled, y=prop, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
  geom_line() + xlab("Mutations per base") + ylab("Mean CpG content") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_discrete("Mutational bias") + 
  facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
  theme(aspect.ratio=2) + manual_color + theme(panel.spacing=unit(0.5, "cm"))
ggsave("../results/figures/supp_interm_figures/random_sims-CpG_prop-kmers2_simsall_mean.pdf", scale=plot_scale_alt)
 
# kmers3
random_scaled_temp <- process_file(random_scaled_0_track_kmers3)
random_scaled_0_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_50_track_kmers3)
random_scaled_50_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_100_track_kmers3)
random_scaled_100_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_200_track_kmers3)
random_scaled_200_kmers <- process_kmers(random_scaled_temp)
write_tsv(random_scaled_0_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_0-kmers3_mean.txt")
write_tsv(random_scaled_50_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_50-kmers3_mean.txt")
write_tsv(random_scaled_100_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_100-kmers3_mean.txt")
write_tsv(random_scaled_200_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_200-kmers3_mean.txt")
random_scaled_temp <- NULL

random_scaled_all <- rbind(
  random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
  random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
  random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
  random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
)

random_scaled_var <- random_scaled_all %>% 
  group_by(mut_scaled, scaling_factor, constr_cond) %>% 
  summarise(prop_var=var(prop)) %>% ungroup()
write_tsv(random_scaled_all, "../results/figures/supp_interm_figures/random_sims-mut_bias-kmers3_motif_var.txt")

ggplot(random_scaled_var, 
       aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
  geom_line() + xlab("Mutations per base") + ylab("Variance in 3-mer frequencies") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_discrete("Mutational bias") + 
  facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
  theme(aspect.ratio=2) + manual_color + theme(panel.spacing=unit(0.5, "cm"))
ggsave("../results/figures/supp_interm_figures/random_sims-mut_bias-kmers3_simsall_prop_var.pdf", scale=plot_scale_alt)
 
# kmers4
random_scaled_temp <- process_file(random_scaled_0_track_kmers4)
random_scaled_0_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_50_track_kmers4)
random_scaled_50_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_100_track_kmers4)
random_scaled_100_kmers <- process_kmers(random_scaled_temp)
random_scaled_temp <- process_file(random_scaled_200_track_kmers4)
random_scaled_200_kmers <- process_kmers(random_scaled_temp)
write_tsv(random_scaled_0_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_0-kmers4_mean.txt")
write_tsv(random_scaled_50_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_50-kmers4_mean.txt")
write_tsv(random_scaled_100_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_100-kmers4_mean.txt")
write_tsv(random_scaled_200_kmers, "../results/figures/supp_interm_figures/random_sims-scaled_200-kmers4_mean.txt")
random_scaled_temp <- NULL

random_scaled_all <- rbind(
  random_scaled_0_kmers %>% mutate(scaling_factor="0.0 (None)"), 
  random_scaled_50_kmers %>% mutate(scaling_factor="0.5"), 
  random_scaled_100_kmers %>% mutate(scaling_factor="1.0 (Original)"), 
  random_scaled_200_kmers %>% mutate(scaling_factor="2.0")
)

random_scaled_var <- random_scaled_all %>% 
  group_by(mut_scaled, scaling_factor, constr_cond) %>% 
  summarise(prop_var=var(prop)) %>% ungroup()
write_tsv(random_scaled_all, "../results/figures/supp_interm_figures/random_sims-mut_bias-kmers4_motif_var.txt")

ggplot(random_scaled_var, 
       aes(x=mut_scaled, y=prop_var, color=constr_cond, linetype=as.factor(scaling_factor))) + plot_text_alt + 
  geom_line() + xlab("Mutations per base") + ylab("Variance in 4-mer frequencies") + 
  scale_color_discrete(legend_constraint) + 
  guides(color=guide_legend(order=-1)) + 
  scale_linetype_discrete("Mutational bias") + 
  facet_wrap(~constr_cond) + theme(strip.text.x = element_blank()) + 
  theme(aspect.ratio=2) + manual_color + theme(panel.spacing=unit(0.5, "cm"))
ggsave("../results/figures/supp_interm_figures/random_sims-mut_bias-kmers4_simsall_prop_var.pdf", scale=plot_scale_alt)
