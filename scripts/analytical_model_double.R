library(tidyverse)
# library(MCMCpack)
# library(ramify)

sample_M_double <- function(N, M_mu, M_sigma) {
  # log-normally distributed entries
  M_prime <- matrix(nrow=N, ncol=1)
  for (i in 1:N) {
    M_prime[i, 1] <- exp(rnorm(1, M_mu, M_sigma))
  }
  M <- matrix(nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in 1:N) { 
      if (i != j) {
        # M[i, j] <- M_prime[j, 1]*pmax(
        #   0, exp(rnorm(1, M_mu, M_sigma)))
        M[i, j] <- M_prime[i, 1]*pmax(
          0, exp(rnorm(1, M_mu, M_sigma)))
      }
    }
  }
  # M <- t(M)
  for (i in 1:N) {
    # M[i, i] <- -sum(M[, i], na.rm=T)
    M[i, i] <- -sum(M[i, ], na.rm=T)
  }
  return(M)
}

# test sampling
M_prime <- sample_M_double(10, 0, 1)
M_prime[M_prime<0] <- 0
rowSums(M_prime)
range(M_prime[M_prime>0])
max(range(M_prime[M_prime>0]))/
  min(range(M_prime[M_prime>0]))
plot(sort(rowSums(M_prime)/min(rowSums(M_prime))))

hexamer_erv_motif <- read_tsv("../data/hexamers_ERVaveraged.txt")
names(hexamer_erv_motif) <- c("wtMotif", "ERV_rel_rate")
plot(sort(hexamer_erv_motif$ERV_rel_rate/
            min(hexamer_erv_motif$ERV_rel_rate)))

# remaining parameters
sample_Q <- function(N, Q_alpha) {
  # symmetric dirichlet distributed
  if (Q_alpha == Inf) { 
    Q <- rep(1/N, N)
  } else {
    Q <- t(rdirichlet(1, rep(Q_alpha, N)))
  }
  return(Q)
}

get_equilibrium_frequencies <- function(N, M, S, Q) {
  E_N <- as.matrix(c(rep(0, N-1), 1))
  # M_N <- M; M_N[N, ] <- rep(1, N)
  M_N <- M; M_N[, N] <- rep(1, N)
  Q_S <- S*Q
  # P_0 <- solve(M_N) %*% E_N
  P_0 <- t(t(E_N) %*% solve(M_N))
  # analytic derivation
  P <- (1-S)*P_0 + Q_S
  return(P)
}

# constants
N = 20  # 50
M_mu = 0

# initialize
output_list_v2 <- tibble(
  "S"=as.numeric(), "M_sigma"=as.numeric(), "Q_alpha"=as.numeric(), "iter_index"=as.numeric(), 
  "motif_index"=as.numeric(), "motif_exon"=as.numeric(), "motif_intron"=as.numeric(), 
  "motif_ratio"=as.numeric(), "motif_mutrate"=as.numeric(), "motif_ratio_max"=as.numeric())

# iterate
for (S in c(0, 0.05, 0.1, 0.2, 0.4, 0.8)) {
  for (M_sigma in c(0:100)/100) {
  # for (M_sigma in 2**(c(0:40)/5)) {
    for (Q_alpha in c(Inf)) {  # , 10, 1)) { 
      for (i in 1:1) {
        print(c(S, M_sigma, Q_alpha, i))

        # sample parameters
        M <- sample_M_double(N, 0, M_sigma)
        Q <- sample_Q(N, Q_alpha)
        
        # equilibrium frequencies
        P_exon <- get_equilibrium_frequencies(N, M, S, Q)
        P_intron <- get_equilibrium_frequencies(N, M, 0, Q)
        
        # append to output
        for (j in 1:length(P_exon)) {
          output_list_v2 <- bind_rows(output_list_v2, c(
            "S"=S, "M_sigma"=M_sigma, "Q_alpha"=Q_alpha, "iter_index"=i, "motif_index"=j, 
            "motif_exon"=P_exon[j], "motif_intron"=P_intron[j], "motif_ratio"=P_exon[j]/P_intron[j],
            "motif_mutrate"=M[j, j], "mutrate_ratio_max"=max(M[M > 0])/min(M[M > 0])))
        }
      }
    }
  }
}

write_tsv(output_list_v2, "../results/figures/output_list_v2_double.txt")
output_list_v2 <- read_tsv("../results/figures/output_list_v2_double.txt")

# # # original
palette_blue_red <- c("#297bb0", "#a9d4d9", "#f5e090", "#faae61", "#d71d1c")

output_list_v2_temp <- filter(output_list_v2, S!=1) %>% 
  group_by(M_sigma, S, Q_alpha) %>% summarise(
    motif_ratio_0 = quantile(motif_ratio, probs=0/100, na.rm=T), 
    motif_ratio_25 = quantile(motif_ratio, probs=25/100, na.rm=T),
    motif_ratio_50 = quantile(motif_ratio, probs=50/100, na.rm=T),
    motif_ratio_75 = quantile(motif_ratio, probs=75/100, na.rm=T),
    motif_ratio_100 = quantile(motif_ratio, probs=100/100, na.rm=T)) %>% 
  gather(motif_ratio_summary, motif_ratio, motif_ratio_0, 
    motif_ratio_25, motif_ratio_50, motif_ratio_75, 
    motif_ratio_100) %>% 
  mutate(motif_ratio_summary = as.factor(as.numeric(gsub(
    "motif_ratio_", "", motif_ratio_summary))))

output_list_v2_mean <- output_list_v2_temp %>% 
  dplyr::select(M_sigma, S, Q_alpha, motif_ratio_summary, motif_ratio)
write_tsv(output_list_v2_mean, "../results/figures/output_list_v2_mean_double.txt")
output_list_v2_mean <- read_tsv("../results/figures/output_list_v2_mean_double.txt")

# # # revision
output_list_v2_exon_temp <- filter(output_list_v2, S!=1) %>% 
  group_by(M_sigma, S, Q_alpha) %>% summarise(
    motif_exon_0 = quantile(motif_exon, probs=0/100, na.rm=T), 
    motif_exon_25 = quantile(motif_exon, probs=25/100, na.rm=T),
    motif_exon_50 = quantile(motif_exon, probs=50/100, na.rm=T),
    motif_exon_75 = quantile(motif_exon, probs=75/100, na.rm=T),
    motif_exon_100 = quantile(motif_exon, probs=100/100, na.rm=T)) %>% 
  gather(motif_exon_summary, motif_exon, motif_exon_0, 
    motif_exon_25, motif_exon_50, motif_exon_75, 
    motif_exon_100) %>% 
  mutate(motif_exon_summary = as.factor(as.numeric(gsub(
    "motif_exon_", "", motif_exon_summary))))

output_list_v2_exon_mean <- output_list_v2_exon_temp %>% 
  dplyr::select(M_sigma, S, Q_alpha, motif_exon_summary, motif_exon)
write_tsv(output_list_v2_exon_mean, "../results/figures/output_list_v2_exon_mean_double.txt")
output_list_v2_exon_mean <- read_tsv("../results/figures/output_list_v2_exon_mean_double.txt")

# # # revision
output_list_v2_intron_temp <- filter(output_list_v2, S!=1) %>% 
  group_by(M_sigma, S, Q_alpha) %>% summarise(
    motif_intron_0 = quantile(motif_intron, probs=0/100, na.rm=T), 
    motif_intron_25 = quantile(motif_intron, probs=25/100, na.rm=T),
    motif_intron_50 = quantile(motif_intron, probs=50/100, na.rm=T),
    motif_intron_75 = quantile(motif_intron, probs=75/100, na.rm=T),
    motif_intron_100 = quantile(motif_intron, probs=100/100, na.rm=T)) %>% 
  gather(motif_intron_summary, motif_intron, motif_intron_0, 
    motif_intron_25, motif_intron_50, motif_intron_75, 
    motif_intron_100) %>% 
  mutate(motif_intron_summary = as.factor(as.numeric(gsub(
    "motif_intron_", "", motif_intron_summary))))

output_list_v2_intron_mean <- output_list_v2_intron_temp %>% 
  dplyr::select(M_sigma, S, Q_alpha, motif_intron_summary, motif_intron)
write_tsv(output_list_v2_intron_mean, "../results/figures/output_list_v2_intron_mean_double.txt")
output_list_v2_intron_mean <- read_tsv("../results/figures/output_list_v2_intron_mean_double.txt")

# plots
ggplot(output_list_v2_mean) + 
  geom_line(aes(x=M_sigma, y=log2(motif_ratio), color=as.factor(motif_ratio_summary)), 
               position="identity") + theme(aspect.ratio=3) + 
  facet_grid(cols=vars(S)) + 
  xlab("Increasing mutational bias (mb)") + ylab("Exon vs intron motif ratio (log2)") + 
  scale_color_manual(name="Quantile", values=palette_blue_red) + 
  scale_x_continuous(breaks=c(0, 0.5, 1.0), limits=c(-0.1, 1.1))
ggsave("../results/figures/analytical_model_uniform_mean_double_transpose.pdf", width=5, height=5)

ggplot(output_list_v2_exon_mean) + 
  geom_line(aes(x=M_sigma, y=log2(motif_exon), color=as.factor(motif_exon_summary)), 
               position="identity") + theme(aspect.ratio=3) + 
  facet_grid(cols=vars(S)) + 
  xlab("Increasing mutational bias (mb)") + ylab("Exon motif frequency (log2)") + 
  scale_color_manual(name="Quantile", values=palette_blue_red) + 
  scale_x_continuous(breaks=c(0, 0.5, 1.0), limits=c(-0.1, 1.1)) + 
  ylim(log2(min(c(output_list_v2_exon_mean$motif_exon, output_list_v2_intron_mean$motif_intron))), 
    log2(max(c(output_list_v2_exon_mean$motif_exon, output_list_v2_intron_mean$motif_intron))))
ggsave("../results/figures/analytical_model_uniform_mean_double_transpose_exon.pdf", width=5, height=5)

ggplot(output_list_v2_intron_mean) + 
  geom_line(aes(x=M_sigma, y=log2(motif_intron), color=as.factor(motif_intron_summary)), 
               position="identity") + theme(aspect.ratio=3) + 
  facet_grid(cols=vars(S)) + 
  xlab("Increasing mutational bias (mb)") + ylab("Intron motif frequency (log2)") + 
  scale_color_manual(name="Quantile", values=palette_blue_red) + 
  scale_x_continuous(breaks=c(0, 0.5, 1.0), limits=c(-0.1, 1.1)) + 
  ylim(log2(min(c(output_list_v2_exon_mean$motif_exon, output_list_v2_intron_mean$motif_intron))), 
    log2(max(c(output_list_v2_exon_mean$motif_exon, output_list_v2_intron_mean$motif_intron))))
ggsave("../results/figures/analytical_model_uniform_mean_double_transpose_intron.pdf", width=5, height=5)
