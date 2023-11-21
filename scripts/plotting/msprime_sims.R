source("~/sherlock/oak/stabilizing_selection/scripts/plotting/shared.R")
setwd("~/Documents/selection_Genetics_2024/figures/")

sims = read_delim("~/sherlock/scratch/stabilizing_selection/data_analysis/data/msprime_sims/prob_derived_linked.txt",
                  col_names = FALSE, delim = ' ')
sims = sims %>% mutate(ld = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999))
sims = sims %>% pivot_longer(-ld, values_to = "prob") %>% 
  mutate(allele_counts = as.numeric(substr(name, 2, 4))) %>%
  mutate(allele_frequency = allele_counts / 200)
ggplot(sims %>% filter(ld >= 0.7) %>% mutate(ld = as.factor(ld)), 
       aes(allele_frequency, prob, color=ld)) + 
  geom_point(show.legend=FALSE, size=point_size) + universal_theme + 
  labs(x = "derived allele frequency at tag variant", 
       y = "probability of positive LD with derived allele\nat causal variant",
       title = "B) Probability that two derived alleles are in positive LD") +
  scale_color_manual(values=c("#FBB13C", "#F0774B", "#E43D59", "#D90368"))
ggsave("plots/ld_B.pdf", height=85, width=85, units="mm")

freqs_x = read_csv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/msprime_sims/linked_freqs_1.txt",
                   col_names = c("freq"))
freqs_y = read_csv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/msprime_sims/linked_freqs_2.txt",
                   col_names = c("freq"))
tib = tibble(x = freqs_x$freq, y = freqs_y$freq)
ggplot(tib, aes(x, y)) + 
  geom_point(alpha = 0.03, color="#104F55", size=point_size) + 
  universal_theme + 
  labs(x = "minor allele frequency at tag variant", 
       y = "minor allele frequency at causal variant",
       title = "A) Minor allele frequencies at variants in LD")
ggsave("plots/ld_A.pdf", height=85, width=85, units="mm")
