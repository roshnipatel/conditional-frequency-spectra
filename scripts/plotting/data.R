source("~/sherlock/oak/stabilizing_selection/scripts/plotting/shared.R")
setwd("~/Documents/selection_Genetics_2024/figures/")

#################################### Fig 0A ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/pooled_all_gwas_DAF.txt")
control = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/pooled_all_control_DAF.txt")
tib = bind_rows(gwas %>% select(UKB_WBfreq) %>% mutate(type="gwas"),
                control %>% select(UKB_WBfreq) %>% mutate(type="matched"))
tib = tib %>% mutate(bin = cut(UKB_WBfreq, 100)) %>% 
  group_by(bin, type) %>% 
  summarize(n = n(), freq = min(UKB_WBfreq))
tib = tib %>% group_by(type) %>% mutate(prob = n / sum(n)) %>%
  filter(prob > 0.0005)
ggplot(tib, aes(freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme +
  labs(title = "A) Empirical frequency distributions",
       x = "frequency in UKB White British",
       y = "probability") + pop_color + scale_y_log10()
ggsave("plots/0A.pdf", height=56, width=56, units="mm")

ggplot(tib, aes(freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme +
  labs(title = "Derived allele frequency distribution",
       x = "frequency in UKB White British",
       y = "probability") + pop_color + scale_y_log10()
ggsave("plots/tweet_2.pdf", height=56, width=56, units="mm")

################################### Fig 4AB ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/binary0_joint_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/binary0_joint_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "A) Mean YRI frequency at variants associated with quantitative\ntraits",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/4A.pdf", height=85, width=85, units="mm")

ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "Frequency spectrum in YRI, conditional on frequencies in the\nGWAS cohort (UKB White British)",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/tweet_6.pdf", height=85, width=85, units="mm")

ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, CHBfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=CHBfreq_lower, ymax=CHBfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, CHBfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=CHBfreq_lower, ymax=CHBfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "B) Mean CHB frequency at variants associated with quantitative\ntraits",
       x = "UKB White British frequency",
       y = "mean frequency in CHB")
ggsave("plots/4B.pdf", height=85, width=85, units="mm")

#################################### Fig 4C ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait90_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait90_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "C) Height, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/4C.pdf", height=56, width=56, units="mm")

#################################### Fig 4D ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait96_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait96_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "D) Trunk mass, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/4D.pdf", height=56, width=56, units="mm")

#################################### Fig 4E ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/binary1_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/binary1_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "E) Complex diseases, risk-increasing",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/4E.pdf", height=56, width=56, units="mm")

#################################### Fig 4F ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait90_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait90_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "F) Height, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/4F.pdf", height=56, width=56, units="mm")

#################################### Fig 4G ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait96_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait96_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "G) Trunk mass, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/4G.pdf", height=56, width=56, units="mm")

#################################### Fig 4H ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/binary1_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/binary1_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "H) Complex diseases, risk-decreasing",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/4H.pdf", height=56, width=56, units="mm")

################################### Fig 5CD ####################################
dummy = tibble(freq = seq(0, 1, 0.01)) %>% mutate(het = 2 * freq * (1 - freq))
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/all_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/all_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_line(data=dummy, mapping=aes(freq, het), linetype="dotted") +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, CHBhet_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=CHBhet_lower, ymax=CHBhet_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, CHBhet_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=CHBhet_lower, ymax=CHBhet_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "C) Heterozygosity in CHB",
       x = "UKB White British frequency",
       y = "mean heterozygosity in CHB")
ggsave("plots/5C.pdf", height=85, width=85, units="mm")
ggplot() +
  geom_line(data=dummy, mapping=aes(freq, het), linetype="dotted") +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIhet_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIhet_lower, ymax=YRIhet_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIhet_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIhet_lower, ymax=YRIhet_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "D) Heterozygosity in YRI",
       x = "UKB White British frequency",
       y = "mean heterozygosity in YRI")
ggsave("plots/5D.pdf", height=85, width=85, units="mm")

################################## SUPPLEMENT ##################################
##################################### CpG ######################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/noCpG_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/noCpG_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, CHBfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=CHBfreq_lower, ymax=CHBfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, CHBfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=CHBfreq_lower, ymax=CHBfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "A) Mean CHB frequency at all trait-associated variants,\nexcluding transitions",
       x = "UKB White British frequency",
       y = "mean frequency in CHB")
ggsave("plots/CpG_A.pdf", height=85, width=85, units="mm")
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "B) Mean YRI frequency at all trait-associated variants,\nexcluding transitions",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/CpG_B.pdf", height=85, width=85, units="mm")

################################### CHB PDF ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/all_gwas_CHB_cfs.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/all_matched_CHB_cfs.txt")
gwas = gwas %>% mutate(CHB_freq = CHB_bin / 206) %>% 
  group_by(UKB_sparse_bin) %>% 
  mutate(gwas = UKB_WBfreq / sum(UKB_WBfreq)) %>%
  select(-UKB_WBfreq, -CHB_bin)
matched = matched %>% mutate(CHB_freq = CHB_bin / 206) %>% 
  group_by(UKB_sparse_bin) %>% 
  mutate(matched = UKB_WBfreq / sum(UKB_WBfreq)) %>%
  select(-UKB_WBfreq, -CHB_bin)
tib = full_join(gwas, matched, by=c("UKB_sparse_bin", "CHB_freq")) %>% 
  pivot_longer(-c(UKB_sparse_bin, CHB_freq), names_to = "type", values_to = "prob") %>%
  mutate(prob = replace_na(prob, 0))

a = ggplot(tib %>% filter(UKB_sparse_bin == 1), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "A) UKB WB frequency: [0.01, 0.056)")
b = ggplot(tib %>% filter(UKB_sparse_bin == 2), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "B) UKB WB frequency: [0.056, 0.12)")
c = ggplot(tib %>% filter(UKB_sparse_bin == 3), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "C) UKB WB frequency: [0.12, 0.19)")
d = ggplot(tib %>% filter(UKB_sparse_bin == 4), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "D) UKB WB frequency: [0.19, 0.26)")
e = ggplot(tib %>% filter(UKB_sparse_bin == 5), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "E) UKB WB frequency: [0.26, 0.32)")
f = ggplot(tib %>% filter(UKB_sparse_bin == 6), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "F) UKB WB frequency: [0.32, 0.40)")
g = ggplot(tib %>% filter(UKB_sparse_bin == 7), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "G) UKB WB frequency: [0.40, 0.48)")
h = ggplot(tib %>% filter(UKB_sparse_bin == 8), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "H) UKB WB frequency: [0.48, 0.58)")
i = ggplot(tib %>% filter(UKB_sparse_bin == 9), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "I) UKB WB frequency: [0.58, 0.70)")
j = ggplot(tib %>% filter(UKB_sparse_bin == 10), aes(CHB_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "CHB frequency", y = "probability", title = "J) UKB WB frequency: [0.70, 1)")
ggarrange(a, b, c, d, e, f, g, h, i, j, nrow=4, ncol=3)
ggsave("plots/CHB_pdf.pdf", height=168, width=170, units="mm")

################################### YRI PDF ####################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/all_gwas_YRI_cfs.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/all_matched_YRI_cfs.txt")
gwas = gwas %>% mutate(YRI_freq = YRI_bin / 206) %>% 
  group_by(UKB_sparse_bin) %>% 
  mutate(gwas = UKB_WBfreq / sum(UKB_WBfreq)) %>%
  select(-UKB_WBfreq, -YRI_bin)
matched = matched %>% mutate(YRI_freq = YRI_bin / 206) %>% 
  group_by(UKB_sparse_bin) %>% 
  mutate(matched = UKB_WBfreq / sum(UKB_WBfreq)) %>%
  select(-UKB_WBfreq, -YRI_bin)
tib = full_join(gwas, matched, by=c("UKB_sparse_bin", "YRI_freq")) %>% 
  pivot_longer(-c(UKB_sparse_bin, YRI_freq), names_to = "type", values_to = "prob") %>%
  mutate(prob = replace_na(prob, 0))

a = ggplot(tib %>% filter(UKB_sparse_bin == 1), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "A) UKB WB frequency: [0.01, 0.056)")
b = ggplot(tib %>% filter(UKB_sparse_bin == 2), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "B) UKB WB frequency: [0.056, 0.12)")
c = ggplot(tib %>% filter(UKB_sparse_bin == 3), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "C) UKB WB frequency: [0.12, 0.19)")
d = ggplot(tib %>% filter(UKB_sparse_bin == 4), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "D) UKB WB frequency: [0.19, 0.26)")
e = ggplot(tib %>% filter(UKB_sparse_bin == 5), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "E) UKB WB frequency: [0.26, 0.32)")
f = ggplot(tib %>% filter(UKB_sparse_bin == 6), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "F) UKB WB frequency: [0.32, 0.40)")
g = ggplot(tib %>% filter(UKB_sparse_bin == 7), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "G) UKB WB frequency: [0.40, 0.48)")
h = ggplot(tib %>% filter(UKB_sparse_bin == 8), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "H) UKB WB frequency: [0.48, 0.58)")
i = ggplot(tib %>% filter(UKB_sparse_bin == 9), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "I) UKB WB frequency: [0.58, 0.70)")
j = ggplot(tib %>% filter(UKB_sparse_bin == 10), aes(YRI_freq, prob, color=type)) + 
  geom_line(show.legend=FALSE) + universal_theme + pop_color +
  labs(x = "YRI frequency", y = "probability", title = "J) UKB WB frequency: [0.70, 1)")
ggarrange(a, b, c, d, e, f, g, h, i, j, nrow=4, ncol=3)
ggsave("plots/YRI_pdf.pdf", height=168, width=170, units="mm")

################################ More traits! ##################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait63_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait63_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "A) Mean corpuscular hemoglobin,\npositive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_A.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait63_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait63_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "D) Mean corpuscular hemoglobin,\nnegative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_D.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait77_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait77_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "B) Platelet crit, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_B.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait77_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait77_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "E) Platelet crit, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_E.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait7_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait7_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "C) Alkaline phosphatase, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_C.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait7_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait7_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "F) Alkaline phosphatase, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_F.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait16_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait16_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "G) Basal metabolic rate, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_G.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait16_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait16_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "J) Basal metabolic rate, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_J.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait51_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait51_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "H) IGF-1, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_H.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait51_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait51_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "K) IGF-1, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_K.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait35_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait35_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "I) Gamma glutyltransferase,\npositive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_I.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait35_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait35_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "L) Gamma glutyltransferase,\nnegative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits_L.pdf", height=56, width=56, units="mm")

################################ More traits 2 #################################
gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait27_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait27_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "A) Creatinine, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_A.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait27_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait27_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "D) Creatinine, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_D.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait86_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait86_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "B) SHBG, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_B.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait86_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait86_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "E) SHBG, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_E.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait93_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait93_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "C) Triglycerides, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_C.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait93_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait93_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "F) Triglycerides, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_F.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait30_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait30_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "G) FEV1, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_G.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait30_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait30_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "J) FEV1, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_J.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait23_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait23_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "H) Calcium, positive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_H.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait23_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait23_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "K) Calcium, negative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_K.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait42_positive_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait42_positive_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "I) Heel bone mineral density,\npositive effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_I.pdf", height=56, width=56, units="mm")

gwas = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait42_negative_gwas_summary.txt")
matched = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/empirical/trait42_negative_matched_summary.txt")
gwas = gwas %>% mutate(bin = seq(1, 10), type="gwas")
matched = matched %>% mutate(bin = seq(0, 10), type="matched") %>% 
  select(-UKB_WBfreq_median) %>% filter(bin > 0) %>% 
  inner_join(gwas %>% select(bin, UKB_WBfreq_median)) %>%
  mutate(UKB_WBfreq_median = UKB_WBfreq_median + 0.01)
ggplot() +
  geom_point(data=gwas, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=gwas, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  geom_point(data=matched, mapping=aes(UKB_WBfreq_median, YRIfreq_mean, color=type), show.legend=FALSE) + 
  geom_errorbar(data=matched, mapping=aes(UKB_WBfreq_median, ymin=YRIfreq_lower, ymax=YRIfreq_upper, color=type), show.legend=FALSE, width = 0) +
  universal_theme + pop_color +
  labs(title = "L) Heel bone mineral density,\nnegative effect",
       x = "UKB White British frequency",
       y = "mean frequency in YRI")
ggsave("plots/more_traits2_L.pdf", height=56, width=56, units="mm")

######################### Old functions for loading dfs ########################
load_yri_mean = function(filepath) {
  yri = read_tsv(filepath)
  yri = yri %>% 
    mutate(bin = ntile(UKB_WBfreq, 10))
  bin_bounds = yri %>% group_by(bin) %>% 
    summarize(min = min(UKB_WBfreq), 
              max = max(UKB_WBfreq), 
              median = median(UKB_WBfreq))
  yri_control = yri %>% group_by(bin) %>% 
    summarize_all(mean) %>% 
    select(contains('.'), bin) %>% 
    pivot_longer(-bin, names_to="YRIfreq", values_to="prob") %>% 
    group_by(bin) %>%
    summarize(matched = sum(as.numeric(YRIfreq) * prob))
  yri = yri %>% 
    group_by(bin) %>% 
    summarize(gwas = mean(YRIfreq)) %>% 
    inner_join(bin_bounds) %>%
    inner_join(yri_control) %>%
    pivot_longer(-c(bin, min, max, median), 
                 names_to = "type", values_to = "YRIfreq") %>%
    arrange(desc(type))
  return(yri)
}

load_chb_mean = function(filepath) {
  chb = read_tsv(filepath)
  chb = chb %>% 
    mutate(bin = ntile(UKB_WBfreq, 10))
  bin_bounds = chb %>% group_by(bin) %>% 
    summarize(min = min(UKB_WBfreq), 
              max = max(UKB_WBfreq), 
              median = median(UKB_WBfreq))
  chb_control = chb %>% group_by(bin) %>% 
    summarize_all(mean) %>% 
    select(contains('.'), bin) %>% 
    pivot_longer(-bin, names_to="CHBfreq", values_to="prob") %>% 
    group_by(bin) %>%
    summarize(matched = sum(as.numeric(CHBfreq) * prob))
  chb = chb %>% 
    group_by(bin) %>% 
    summarize(gwas = mean(CHBfreq)) %>% 
    inner_join(bin_bounds) %>%
    inner_join(chb_control) %>%
    pivot_longer(-c(bin, min, max, median), 
                 names_to = "type", values_to = "CHBfreq") %>%
    arrange(desc(type))
  return(chb)
}

load_yri_het = function(filepath) {
  yri = read_tsv(filepath)
  yri = yri %>% 
    mutate(bin = ntile(UKB_WBfreq, 10))
  bin_bounds = yri %>% group_by(bin) %>% 
    summarize(min = min(UKB_WBfreq), 
              max = max(UKB_WBfreq), 
              median = median(UKB_WBfreq))
  yri_control = yri %>% group_by(bin) %>% 
    summarize_all(mean) %>% 
    select(contains('.'), bin) %>% 
    pivot_longer(-bin, names_to="YRIfreq", values_to="prob") %>% 
    group_by(bin) %>%
    summarize(matched = sum(2 * (1 - as.numeric(YRIfreq)) * as.numeric(YRIfreq) * prob))
  yri = yri %>% 
    group_by(bin) %>% 
    summarize(gwas = mean(2 * YRIfreq * (1 - YRIfreq))) %>% 
    inner_join(bin_bounds) %>%
    inner_join(yri_control) %>%
    pivot_longer(-c(bin, min, max, median), 
                 names_to = "type", values_to = "YRIhet") %>%
    arrange(desc(type))
  return(yri)
}

load_chb_het = function(filepath) {
  chb = read_tsv(filepath)
  chb = chb %>% 
    mutate(bin = ntile(UKB_WBfreq, 10))
  bin_bounds = chb %>% group_by(bin) %>% 
    summarize(min = min(UKB_WBfreq), 
              max = max(UKB_WBfreq), 
              median = median(UKB_WBfreq))
  chb_control = chb %>% group_by(bin) %>% 
    summarize_all(mean) %>% 
    select(contains('.'), bin) %>% 
    pivot_longer(-bin, names_to="CHBfreq", values_to="prob") %>% 
    group_by(bin) %>%
    summarize(matched = sum(2 * (1 - as.numeric(CHBfreq)) * as.numeric(CHBfreq) * prob))
  chb = chb %>% 
    group_by(bin) %>% 
    summarize(gwas = mean(2 * CHBfreq * (1 - CHBfreq))) %>% 
    inner_join(bin_bounds) %>%
    inner_join(chb_control) %>%
    pivot_longer(-c(bin, min, max, median), 
                 names_to = "type", values_to = "CHBhet") %>%
    arrange(desc(type))
  return(chb)
}

