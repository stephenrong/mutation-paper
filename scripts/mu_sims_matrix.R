library(tidyverse)
library(ggthemr)
ggthemr("fresh")

# Note 2018-05-21: 
# A 7mer mutation and its reverse complement are missing from the dled table, 
# CGTATCG > CGTCTCG and CGATACG > CGAGACG. We thus used the corresponding
# 5mer mutation rate in place of the 7mer mutation rate for these 7mers,
# GTATC > GTCTC (0.000782534) and GATAC > GAGAC (0.000782534).

# # # 5mers
pentamer_mu <- read_tsv("../data/108290-5mers.txt")

pentamer_mu_edit <- pentamer_mu %>% 
  separate(Type, into=c("wt", "mt"), sep="_") %>% 
  separate(Motif, into=c("motif_forward", "motif_revcomp"), sep="[(]") %>% 
  mutate(motif_revcomp=gsub("[)]", "", motif_revcomp))

pentamer_mu_forward <- dplyr::select(pentamer_mu_edit, -c(motif_revcomp)) %>% 
  dplyr::rename(wtMotif=motif_forward) %>% 
  mutate(mtMotif=ifelse(substr(wtMotif, 3, 3)==substr(wt, 1, 1), 
                        paste(substr(wtMotif, 1, 2), substr(mt, 1, 1), substr(wtMotif, 4, 5), sep=""),
                        paste(substr(wtMotif, 1, 2), substr(mt, 2, 2), substr(wtMotif, 4, 5), sep=""))) %>% 
  dplyr::select(c(wt, mt, wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate))

pentamer_mu_revcomp <- dplyr::select(pentamer_mu_edit, -c(motif_forward)) %>% 
  dplyr::rename(wtMotif=motif_revcomp) %>%
  mutate(mtMotif=ifelse(substr(wtMotif, 3, 3)==substr(wt, 1, 1),
                        paste(substr(wtMotif, 1, 2), substr(mt, 1, 1), substr(wtMotif, 4, 5), sep=""),
                        paste(substr(wtMotif, 1, 2), substr(mt, 2, 2), substr(wtMotif, 4, 5), sep=""))) %>%
  dplyr::select(c(wt, mt, wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate))

pentamer_mu_original <- rbind(
  pentamer_mu_forward %>% mutate(k = row_number()), 
  pentamer_mu_revcomp %>% mutate(k = row_number())) %>% 
  arrange(by=k) %>% dplyr::select(-c(k))

pentamer_mu_final <- pentamer_mu_original %>% mutate(mtMotif=gsub("~", "", mtMotif))
write_tsv(pentamer_mu_final, "../data/108290-5mers-edit.txt")

# save matrix
pentamer_mu_matrix <- pentamer_mu_final %>% 
  mutate(ERV_rel_rate=(ERV_rel_rate)) %>% 
  dplyr::select(wtMotif, mtMotif, ERV_rel_rate)
write_tsv(pentamer_mu_matrix, "../data/mu-matrix-5mers.txt")

# # # 7mers
heptamer_mu <- read_tsv("../data/108290-7mers.txt")

heptamer_mu_edit <- heptamer_mu %>% 
  separate(Type, into=c("wt", "mt"), sep="_") %>% 
  separate(Motif, into=c("motif_forward", "motif_revcomp"), sep="[(]") %>% 
  mutate(motif_revcomp=gsub("[)]", "", motif_revcomp))
  
heptamer_mu_forward <- dplyr::select(heptamer_mu_edit, -c(motif_revcomp)) %>% 
  dplyr::rename(wtMotif=motif_forward) %>% 
  mutate(mtMotif=ifelse(substr(wtMotif, 4, 4)==substr(wt, 1, 1), 
                        paste(substr(wtMotif, 1, 3), substr(mt, 1, 1), substr(wtMotif, 5, 7), sep=""),
                        paste(substr(wtMotif, 1, 3), substr(mt, 2, 2), substr(wtMotif, 5, 7), sep=""))) %>% 
  dplyr::select(c(wt, mt, wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate, 
                  nERVs_mask, nMotifs_mask, ERV_rel_rate_mask, nMAC10, MAC10_rel_rate))
heptamer_mu_forward <- rbind(
  heptamer_mu_forward, c("AT", "CG", "CGTATCG", "CGTCTCG", NA, NA, "0.000782534", NA, NA, NA, NA, NA))

heptamer_mu_revcomp <- dplyr::select(heptamer_mu_edit, -c(motif_forward)) %>% 
  dplyr::rename(wtMotif=motif_revcomp) %>% 
  mutate(mtMotif=ifelse(substr(wtMotif, 4, 4)==substr(wt, 1, 1), 
                        paste(substr(wtMotif, 1, 3), substr(mt, 1, 1), substr(wtMotif, 5, 7), sep=""),
                        paste(substr(wtMotif, 1, 3), substr(mt, 2, 2), substr(wtMotif, 5, 7), sep=""))) %>% 
  dplyr::select(c(wt, mt, wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate, 
                  nERVs_mask, nMotifs_mask, ERV_rel_rate_mask, nMAC10, MAC10_rel_rate))
heptamer_mu_revcomp <- rbind(
  heptamer_mu_revcomp, c("AT", "CG", "CGATACG", "CGAGACG", NA, NA, "0.000782534", NA, NA, NA, NA, NA))

heptamer_mu_original <- rbind(
  heptamer_mu_forward %>% mutate(k = row_number()), 
  heptamer_mu_revcomp %>% mutate(k = row_number())) %>% 
  arrange(by=k) %>% dplyr::select(-c(k))
heptamer_mu_final <- heptamer_mu_original %>% 
  mutate(mtMotif=gsub("~", "", mtMotif))

# join 5mer mutation rate
pentamer_mu_final_temp <- pentamer_mu_final %>% 
  dplyr::select(c(wtMotif, mtMotif, ERV_rel_rate)) %>%
  dplyr::rename(wtMotif_5mer=wtMotif, 
                mtMotif_5mer=mtMotif, 
                ERV_rel_rate_5mer=ERV_rel_rate)

heptamer_mu_final <- heptamer_mu_final %>% 
  mutate(wtMotif_5mer=substr(wtMotif, 2, 6), 
         mtMotif_5mer=substr(mtMotif, 2, 6)) %>% 
  mutate(ERV_rel_rate=as.numeric(ERV_rel_rate)) %>%
  full_join(pentamer_mu_final_temp)
write_tsv(heptamer_mu_final, "../data/108290-7mers-edit.txt")

# save matrix
heptamer_mu_matrix <- heptamer_mu_final %>% 
  dplyr::select(wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate)
write_tsv(heptamer_mu_matrix, "../data/mu-matrix-7mers.txt")

# plot matrix
heptamer_mu_final$row_number <- as.numeric(row.names(heptamer_mu_final))
heptamer_mu_final <- heptamer_mu_final %>% 
  mutate(mutate_CpG_CpA=(grepl("CG", substr(wtMotif, 3, 4)) & (grepl("CA", substr(mtMotif, 3, 4))))) %>% 
  mutate(mutate_CpG_TpG=(grepl("CG", substr(wtMotif, 4, 5)) & (grepl("TG", substr(mtMotif, 4, 5))))) %>% 
  mutate(mutate_CpG = mutate_CpG_CpA | mutate_CpG_TpG)

# ggplot(heptamer_mu_final) + 
#   geom_point(aes(x=as.numeric(ERV_rel_rate), y=as.numeric(ERV_rel_rate_5mer))) + 
#   xlab("ERV relative 7mer mutation rate") + ylab("ERV relative 5mer mutation rate")
# ggsave("../results/figures/mu_matrix-7mer_v_5mer.pdf")

# ggplot(heptamer_mu_final) + 
#   geom_point(aes(x=as.numeric(row_number), y=as.numeric(ERV_rel_rate), color=mutate_CpG)) + 
#   xlab("7mer row index") + ylab("ERV relative 7mer mutation rate") + labs(color="CpG transition")
# ggsave("../results/figures/mu_matrix-CpG_index.pdf")

# mutation bias
mu_mean <- mean(heptamer_mu_final$ERV_rel_rate)
mu_mean

heptamer_mu_scaled_0 <- heptamer_mu_final %>% 
  mutate(ERV_rel_rate = 10**(0*log10(ERV_rel_rate))) %>% 
  mutate(ERV_rel_rate = round(mu_mean*ERV_rel_rate/mean(ERV_rel_rate), 9))
ggplot(heptamer_mu_scaled_0) + 
  geom_point(aes(x=as.numeric(row_number), y=as.numeric(ERV_rel_rate), color=mutate_CpG)) + 
  xlab("7mer row index") + ylab("ERV relative 7mer mutation rate") + labs(color="CpG transition")
ggsave("../results/figures/mu_matrix-scaled_0.pdf")
mean(heptamer_mu_scaled_0$ERV_rel_rate)

heptamer_mu_scaled_50 <- heptamer_mu_final %>% 
  mutate(ERV_rel_rate = 10**(0.5*log10(ERV_rel_rate))) %>% 
  mutate(ERV_rel_rate = round(mu_mean*ERV_rel_rate/mean(ERV_rel_rate), 9))
ggplot(heptamer_mu_scaled_50) + 
  geom_point(aes(x=as.numeric(row_number), y=as.numeric(ERV_rel_rate), color=mutate_CpG)) + 
  xlab("7mer row index") + ylab("ERV relative 7mer mutation rate") + labs(color="CpG transition")
ggsave("../results/figures/mu_matrix-scaled_50.pdf")
mean(heptamer_mu_scaled_50$ERV_rel_rate)

heptamer_mu_scaled_100 <- heptamer_mu_final %>% 
  mutate(ERV_rel_rate = 10**(1*log10(ERV_rel_rate))) %>% 
  mutate(ERV_rel_rate = round(mu_mean*ERV_rel_rate/mean(ERV_rel_rate), 9))
ggplot(heptamer_mu_scaled_100) + 
  geom_point(aes(x=as.numeric(row_number), y=as.numeric(ERV_rel_rate), color=mutate_CpG)) + 
  xlab("7mer row index") + ylab("ERV relative 7mer mutation rate") + labs(color="CpG transition")
ggsave("../results/figures/mu_matrix-scaled_100.pdf")
mean(heptamer_mu_scaled_100$ERV_rel_rate)

heptamer_mu_scaled_200 <- heptamer_mu_final %>% 
  mutate(ERV_rel_rate = 10**(2*log10(ERV_rel_rate))) %>% 
  mutate(ERV_rel_rate = round(mu_mean*ERV_rel_rate/mean(ERV_rel_rate), 9))
ggplot(heptamer_mu_scaled_200) + 
  geom_point(aes(x=as.numeric(row_number), y=as.numeric(ERV_rel_rate), color=mutate_CpG)) + 
  xlab("7mer row index") + ylab("ERV relative 7mer mutation rate") + labs(color="CpG transition")
ggsave("../results/figures/mu_matrix-scaled_200.pdf")
mean(heptamer_mu_scaled_200$ERV_rel_rate)

# save matrix
heptamer_mu_matrix_0 <- heptamer_mu_scaled_0 %>% 
  dplyr::select(wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate)
write_tsv(heptamer_mu_matrix_0, "../data/mu-matrix-7mers-scaled_0.txt")

heptamer_mu_matrix_50 <- heptamer_mu_scaled_50 %>% 
  dplyr::select(wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate)
write_tsv(heptamer_mu_matrix_50, "../data/mu-matrix-7mers-scaled_50.txt")

heptamer_mu_matrix_100 <- heptamer_mu_scaled_100 %>% 
  dplyr::select(wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate)
write_tsv(heptamer_mu_matrix_100, "../data/mu-matrix-7mers-scaled_100.txt")

heptamer_mu_matrix_200 <- heptamer_mu_scaled_200 %>% 
  dplyr::select(wtMotif, mtMotif, nERVs, nMotifs, ERV_rel_rate)
write_tsv(heptamer_mu_matrix_200, "../data/mu-matrix-7mers-scaled_200.txt")
