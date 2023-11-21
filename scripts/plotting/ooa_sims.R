source("~/sherlock/oak/stabilizing_selection/scripts/plotting/shared.R")
setwd("~/Documents/selection_Genetics_2024/figures/")

ancestral = bind_rows(read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h5e6_s-1.0e-10_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "stabilizing, hs=-5e-4", count = seq(1, 47442)),
                      read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s-1.0e-3_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "negative, hs=-5e-4", count = seq(1, 47442)),
                      read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s+1.0e-3_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "positive, hs=+5e-4", count = seq(1, 47442)),
                      read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s0.0_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "0", count = seq(1, 47442)))
ancestral = ancestral %>% mutate(freq = count / 47442)

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-jouganous_wo_migration/")
dtwf_ceu_anc = load_ooa_dist(dtwf_filepath_prefix, "ceu_conditional_ancestor.txt")
dtwf_ooa_ceu = load_ooa_dist(dtwf_filepath_prefix, "ooa_conditional_ceu.txt")
dtwf_ceu_ooa = load_ooa_dist(dtwf_filepath_prefix, "ceu_conditional_ooa.txt")
dtwf_anc_ceu = load_ooa_dist(dtwf_filepath_prefix, "ancestor_conditional_ceu.txt")
dtwf_anc_yri = load_ooa_dist(dtwf_filepath_prefix, "ancestor_conditional_yri.txt")
dtwf_yri_anc = load_ooa_dist(dtwf_filepath_prefix, "yri_conditional_ancestor.txt")
dtwf_ceu_yri = load_ooa_dist(dtwf_filepath_prefix, "ceu_conditional_yri.txt")
dtwf_yri_ceu = load_ooa_dist(dtwf_filepath_prefix, "yri_conditional_ceu.txt")
dtwf_chb_ceu = load_ooa_dist(dtwf_filepath_prefix, "chb_conditional_ceu.txt")
dtwf_ceu = load_ooa_marginal(dtwf_filepath_prefix, "ceu_marginal.txt")

sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/jouganous_wo_migration/")
sims_ceu_anc = load_ooa_dist(sims_filepath_prefix, "ceu_conditional_ancestor.txt")
sims_ooa_ceu = load_ooa_dist(sims_filepath_prefix, "ooa_conditional_ceu.txt")
sims_ceu_ooa = load_ooa_dist(sims_filepath_prefix, "ceu_conditional_ooa.txt")
sims_anc_ceu = load_ooa_dist(sims_filepath_prefix, "ancestor_conditional_ceu.txt")
sims_anc_yri = load_ooa_dist(sims_filepath_prefix, "ancestor_conditional_yri.txt")
sims_yri_anc = load_ooa_dist(sims_filepath_prefix, "yri_conditional_ancestor.txt")
sims_ceu_yri = load_ooa_dist(sims_filepath_prefix, "ceu_conditional_yri.txt")
sims_yri_ceu = load_ooa_dist(sims_filepath_prefix, "yri_conditional_ceu.txt")
sims_chb_ceu = load_ooa_dist(sims_filepath_prefix, "chb_conditional_ceu.txt")
sims_ceu = load_ooa_marginal(sims_filepath_prefix, "ceu_marginal.txt")
sims_ooa = load_ooa_marginal(sims_filepath_prefix, "ooa_marginal.txt")
sims_yri = load_ooa_marginal(sims_filepath_prefix, "yri_marginal.txt")

################################ CREATE FILTERS ################################
ooa_filter = sims_ooa %>% filter(prob >= 1/3e9) %>% transmute(ooAfreq = freq, s)
ceu_filter = sims_ceu %>% 
  filter(freq > 0, freq < 1) %>%
  group_by(s) %>% 
  mutate(prob = prob / sum(prob)) %>% 
  filter(prob >= 1/3e9) %>% 
  transmute(CEUfreq = freq, s)

################################## Fig 3B ######################################
dtwf_ceu_cumulative = dtwf_ceu %>% 
  filter(freq > 0, freq < 1) %>%
  group_by(s) %>% 
  mutate(prob = prob / sum(prob)) %>% 
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
sims_ceu_cumulative = sims_ceu %>% 
  filter(freq > 0, freq < 1) %>%
  group_by(s) %>% 
  mutate(prob = prob / sum(prob)) %>% 
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_ceu_cumulative, 
            mapping=aes(freq, cum_prob, color=s),
            size=line_size, show.legend = FALSE) + 
  # geom_point(data=sims_ceu_cumulative,
  #            mapping=aes(freq, cum_prob, color=s),
  #            size=point_size, shape=1,
  #            show.legend = FALSE) +
  universal_theme + selection_color + 
  labs(x = "CEU frequency", 
       y = "cumulative probability",
       title = "B) Marginal distribution in CEU")
ggsave("plots/3B.pdf", height=56, width=56, units="mm")

################################## Fig 3C ######################################
dtwf_chb_ceu_mean = dtwf_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  summarize(mean_chb = sum(as.numeric(CHBfreq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s), CEUfreq < 1)
sims_chb_ceu_mean = sims_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  summarize(mean_chb = sum(as.numeric(CHBfreq) * prob)) %>%
  inner_join(ceu_filter) %>% 
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s), CEUfreq < 1)
ggplot() + 
  geom_line(data=dtwf_chb_ceu_mean,
            mapping=aes(CEUfreq, mean_chb, color=s),
            show.legend = FALSE,
            size=line_size) +
  # geom_point(data=sims_chb_ceu_mean,
  #            mapping=aes(CEUfreq, mean_chb, color=s),
  #            show.legend = FALSE,
  #            shape=1, size=point_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "CEU frequency", 
       y = "expected frequency in CHB",
       title = "C)")
ggsave("plots/3C.pdf", height=56, width=56, units="mm")

################################## Fig 3D ######################################
dtwf_yri_ceu_mean = dtwf_yri_ceu %>% group_by(CEUfreq, s) %>%
  summarize(mean_yri = sum(as.numeric(YRIfreq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s), CEUfreq < 1)
sims_yri_ceu_mean = sims_yri_ceu %>% group_by(CEUfreq, s) %>%
  summarize(mean_yri = sum(as.numeric(YRIfreq) * prob)) %>%
  inner_join(ceu_filter) %>% 
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s), CEUfreq < 1)
ggplot() + 
  geom_line(data=dtwf_yri_ceu_mean, 
            mapping=aes(CEUfreq, mean_yri, color=s),
            show.legend = FALSE,
            size=line_size) +
  # geom_point(data=sims_yri_ceu_mean,
  #            mapping=aes(CEUfreq, mean_yri, color=s),
  #            show.legend = FALSE,
  #            shape=1, size=point_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "CEU frequency", 
       y = "expected frequency in YRI",
       title = "D)")
ggsave("plots/3D.pdf", height=56, width=56, units="mm")

################################### Fig 5A #####################################
dtwf_chb_ceu_port = dtwf_chb_ceu %>% 
  group_by(s, CEUfreq) %>% 
  mutate(CHBfreq = as.numeric(CHBfreq)) %>%
  summarize(expected_het = sum(prob * 2 * CHBfreq * (1 - CHBfreq))) %>% 
  mutate(CEUfreq = as.numeric(CEUfreq)) %>%
  filter(CEUfreq > 0, CEUfreq < 1) %>%
  mutate(CEU_het = 2 * CEUfreq * (1 - CEUfreq))

ggplot(dtwf_chb_ceu_port) + 
  geom_line(mapping=aes(CEUfreq, CEU_het), 
            linetype="dotted", show.legend = FALSE) +
  geom_line(mapping=aes(CEUfreq, expected_het, color=s), 
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "CEU frequency", y = "expected CHB heterozygosity", 
       title = "A) Expected CHB heterozygosity")
ggsave("plots/5A.pdf", height=85, width=85, units="mm")

################################### Fig 5B #####################################
dtwf_yri_ceu_port = dtwf_yri_ceu %>% 
  group_by(s, CEUfreq) %>% 
  mutate(YRIfreq = as.numeric(YRIfreq)) %>%
  summarize(expected_het = sum(prob * 2 * YRIfreq * (1 - YRIfreq))) %>% 
  mutate(CEUfreq = as.numeric(CEUfreq)) %>%
  filter(CEUfreq > 0, CEUfreq < 1) %>%
  mutate(CEU_het = 2 * CEUfreq * (1 - CEUfreq))

ggplot(dtwf_yri_ceu_port) + 
  geom_line(mapping=aes(CEUfreq, CEU_het), 
            linetype="dotted", show.legend = FALSE) +
  geom_line(mapping=aes(CEUfreq, expected_het, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "CEU frequency", y = "expected YRI heterozygosity", 
       title = "B) Expected YRI heterozygosity")
ggsave("plots/5B.pdf", height=85, width=85, units="mm")

################################## SUPPLEMENT ##################################
########################### ooA forward + backward #############################
### CEU | ancestor (D)
dtwf_ceu_anc_mean = dtwf_ceu_anc %>% 
  group_by(ancestral_freq, s) %>% 
  summarize(mean_ceu = sum(as.numeric(CEUfreq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq)) %>% 
  filter(!grepl("e-5", s))
sims_ceu_anc_mean = sims_ceu_anc %>% 
  group_by(ancestral_freq, s) %>% 
  summarize(mean_ceu = sum(as.numeric(CEUfreq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_ceu_anc_mean,
            mapping=aes(ancestral_freq, mean_ceu, color=s),
            show.legend = FALSE,
            size=line_size) +
  # geom_point(data=sims_ceu_anc_mean,
  #            mapping=aes(ancestral_freq, mean_ceu, color=s),
  #            show.legend = FALSE,
  #            shape=1, size=point_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-YRI ancestor", 
       y = "expected frequency in CEU",
       title = "D) Forward transition: CEU-YRI ancestor")
ggsave("plots/ooa_forward_backward_D.pdf", height=85, width=85, units="mm")

### CEU | ooA (B)
dtwf_ceu_ooa_mean = dtwf_ceu_ooa %>% 
  group_by(ooAfreq, s) %>% 
  summarize(mean_ceu = sum(as.numeric(CEUfreq) * prob)) %>%
  mutate(ooAfreq = as.numeric(ooAfreq)) %>% 
  filter(!grepl("e-5", s))
sims_ceu_ooa_mean = sims_ceu_ooa %>% 
  group_by(ooAfreq, s) %>% 
  summarize(mean_ceu = sum(as.numeric(CEUfreq) * prob)) %>%
  mutate(ooAfreq = as.numeric(ooAfreq)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_ceu_ooa_mean,
            mapping=aes(ooAfreq, mean_ceu, color=s),
            show.legend = FALSE,
            size=line_size) +
  # geom_point(data=sims_ceu_ooa_mean,
  #            mapping=aes(ooAfreq, mean_ceu, color=s),
  #            show.legend = FALSE,
  #            shape=1, size=point_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-CHB ancestor", 
       y = "expected frequency in CEU",
       title = "B) Forward transition: CEU-CHB ancestor")
ggsave("plots/ooa_forward_backward_B.pdf", height=85, width=85, units="mm")

### ancestor | CEU (C)
dtwf_anc_ceu_mean = dtwf_anc_ceu %>% 
  group_by(CEUfreq, s) %>% filter(CEUfreq < 1) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s))
sims_anc_ceu_mean = sims_anc_ceu %>% 
  group_by(CEUfreq, s) %>% filter(CEUfreq < 1) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_anc_ceu_mean,
            mapping=aes(CEUfreq, mean_anc, color=s),
            show.legend = FALSE,
            size=line_size) +
  # geom_point(data=sims_anc_ceu_mean,
  #            mapping=aes(CEUfreq, mean_anc, color=s),
  #            show.legend = FALSE,
  #            shape=1, size=point_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU", 
       y = "expected frequency in CEU-YRI ancestor",
       title = "C) Backward transition: CEU-YRI ancestor")
ggsave("plots/ooa_forward_backward_C.pdf", height=85, width=85, units="mm")

### ooA | CEU (A)
dtwf_ooa_ceu_mean = dtwf_ooa_ceu %>% 
  group_by(CEUfreq, s) %>% filter(CEUfreq < 1) %>%
  summarize(mean_ooa = sum(as.numeric(ooAfreq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s))
sims_ooa_ceu_mean = sims_ooa_ceu %>% 
  group_by(CEUfreq, s) %>% filter(CEUfreq < 1) %>%
  summarize(mean_ooa = sum(as.numeric(ooAfreq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_ooa_ceu_mean,
            mapping=aes(CEUfreq, mean_ooa, color=s),
            show.legend = FALSE,
            size=line_size) +
  # geom_point(data=sims_ooa_ceu_mean,
  #            mapping=aes(CEUfreq, mean_ooa, color=s),
  #            show.legend = FALSE,
  #            shape=1, size=point_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU", 
       y = "expected frequency in CEU-CHB ancestor",
       title = "A) Backward transition: CEU-CHB ancestor")
ggsave("plots/ooa_forward_backward_A.pdf", height=85, width=85, units="mm")

################################ ooA fb PDFs ###################################
dtwf_anc_ceu_low = dtwf_anc_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq)) %>% 
  filter(as.numeric(CEUfreq) > 0.05,
         as.numeric(CEUfreq) < 0.06) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_anc_ceu_low, 
            mapping=aes(ancestral_freq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-YRI ancestor", 
       y = "cumulative probability",
       title = "A) CEU frequency = 0.05")
ggsave("plots/ooa_fb_A.pdf", height=56, width=56, units="mm")

dtwf_anc_ceu_mid = dtwf_anc_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq)) %>% 
  filter(as.numeric(CEUfreq) > 0.3,
         as.numeric(CEUfreq) < 0.31) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_anc_ceu_mid, 
            mapping=aes(ancestral_freq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-YRI ancestor", 
       y = "cumulative probability",
       title = "B) CEU frequency = 0.3")
ggsave("plots/ooa_fb_B.pdf", height=56, width=56, units="mm")

dtwf_anc_ceu_high = dtwf_anc_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq)) %>% 
  filter(as.numeric(CEUfreq) > 0.6,
         as.numeric(CEUfreq) < 0.61) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_anc_ceu_high, 
            mapping=aes(ancestral_freq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-YRI ancestor", 
       y = "cumulative probability",
       title = "C) CEU frequency = 0.6")
ggsave("plots/ooa_fb_C.pdf", height=56, width=56, units="mm")

dtwf_ooa_ceu_low = dtwf_ooa_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(ooAfreq = as.numeric(ooAfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.05,
         as.numeric(CEUfreq) < 0.06) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_ooa_ceu_low, 
            mapping=aes(ooAfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-CHB ancestor", 
       y = "cumulative probability",
       title = "D) CEU frequency = 0.05")
ggsave("plots/ooa_fb_D.pdf", height=56, width=56, units="mm")

dtwf_ooa_ceu_mid = dtwf_ooa_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(ooAfreq = as.numeric(ooAfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.3,
         as.numeric(CEUfreq) < 0.31) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_ooa_ceu_mid, 
            mapping=aes(ooAfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-CHB ancestor", 
       y = "cumulative probability",
       title = "E) CEU frequency = 0.3")
ggsave("plots/ooa_fb_E.pdf", height=56, width=56, units="mm")

dtwf_ooa_ceu_high = dtwf_ooa_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(ooAfreq = as.numeric(ooAfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.6,
         as.numeric(CEUfreq) < 0.61) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_ooa_ceu_high, 
            mapping=aes(ooAfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CEU-CHB ancestor", 
       y = "cumulative probability",
       title = "F) CEU frequency = 0.6")
ggsave("plots/ooa_fb_F.pdf", height=56, width=56, units="mm")

########################### ooA conditional PDFs ###############################
dtwf_yri_ceu_low = dtwf_yri_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(YRIfreq = as.numeric(YRIfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.05,
         as.numeric(CEUfreq) < 0.06) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_yri_ceu_low, 
            mapping=aes(YRIfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in YRI", 
       y = "cumulative probability",
       title = "A) CEU frequency = 0.05")
ggsave("plots/ooa_cond_pdf_A.pdf", height=56, width=56, units="mm")

dtwf_yri_ceu_mid = dtwf_yri_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(YRIfreq = as.numeric(YRIfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.3,
         as.numeric(CEUfreq) < 0.31) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_yri_ceu_mid, 
            mapping=aes(YRIfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in YRI", 
       y = "cumulative probability",
       title = "B) CEU frequency = 0.3")
ggsave("plots/ooa_cond_pdf_B.pdf", height=56, width=56, units="mm")

dtwf_yri_ceu_high = dtwf_yri_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(YRIfreq = as.numeric(YRIfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.6,
         as.numeric(CEUfreq) < 0.61) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_yri_ceu_high, 
            mapping=aes(YRIfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in YRI", 
       y = "cumulative probability",
       title = "C) CEU frequency = 0.6")
ggsave("plots/ooa_cond_pdf_C.pdf", height=56, width=56, units="mm")

dtwf_chb_ceu_low = dtwf_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(CHBfreq = as.numeric(CHBfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.05,
         as.numeric(CEUfreq) < 0.06) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_chb_ceu_low, 
            mapping=aes(CHBfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CHB", 
       y = "cumulative probability",
       title = "D) CEU frequency = 0.05")
ggsave("plots/ooa_cond_pdf_D.pdf", height=56, width=56, units="mm")

dtwf_chb_ceu_mid = dtwf_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(CHBfreq = as.numeric(CHBfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.3,
         as.numeric(CEUfreq) < 0.31) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_chb_ceu_mid, 
            mapping=aes(CHBfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CHB", 
       y = "cumulative probability",
       title = "E) CEU frequency = 0.3")
ggsave("plots/ooa_cond_pdf_E.pdf", height=56, width=56, units="mm")

dtwf_chb_ceu_high = dtwf_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(CHBfreq = as.numeric(CHBfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.6,
         as.numeric(CEUfreq) < 0.61) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-5", s))
ggplot() + 
  geom_line(data=dtwf_chb_ceu_high, 
            mapping=aes(CHBfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CHB", 
       y = "cumulative probability",
       title = "F) CEU frequency = 0.6")
ggsave("plots/ooa_cond_pdf_F.pdf", height=56, width=56, units="mm")

#################### ooA conditional PDFs - weak selection #####################
dtwf_yri_ceu_low = dtwf_yri_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(YRIfreq = as.numeric(YRIfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.05,
         as.numeric(CEUfreq) < 0.06) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-4", s))
ggplot() + 
  geom_line(data=dtwf_yri_ceu_low, 
            mapping=aes(YRIfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in YRI", 
       y = "cumulative probability",
       title = "A) CEU frequency = 0.05")
ggsave("plots/ooa_cond_pdf_weak_A.pdf", height=56, width=56, units="mm")

dtwf_yri_ceu_mid = dtwf_yri_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(YRIfreq = as.numeric(YRIfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.3,
         as.numeric(CEUfreq) < 0.31) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-4", s))
ggplot() + 
  geom_line(data=dtwf_yri_ceu_mid, 
            mapping=aes(YRIfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in YRI", 
       y = "cumulative probability",
       title = "B) CEU frequency = 0.3")
ggsave("plots/ooa_cond_pdf_weak_B.pdf", height=56, width=56, units="mm")

dtwf_yri_ceu_high = dtwf_yri_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(YRIfreq = as.numeric(YRIfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.6,
         as.numeric(CEUfreq) < 0.61) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-4", s))
ggplot() + 
  geom_line(data=dtwf_yri_ceu_high, 
            mapping=aes(YRIfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in YRI", 
       y = "cumulative probability",
       title = "C) CEU frequency = 0.6")
ggsave("plots/ooa_cond_pdf_weak_C.pdf", height=56, width=56, units="mm")

dtwf_chb_ceu_low = dtwf_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(CHBfreq = as.numeric(CHBfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.05,
         as.numeric(CEUfreq) < 0.06) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-4", s))
ggplot() + 
  geom_line(data=dtwf_chb_ceu_low, 
            mapping=aes(CHBfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CHB", 
       y = "cumulative probability",
       title = "D) CEU frequency = 0.05")
ggsave("plots/ooa_cond_pdf_weak_D.pdf", height=56, width=56, units="mm")

dtwf_chb_ceu_mid = dtwf_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(CHBfreq = as.numeric(CHBfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.3,
         as.numeric(CEUfreq) < 0.31) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-4", s))
ggplot() + 
  geom_line(data=dtwf_chb_ceu_mid, 
            mapping=aes(CHBfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CHB", 
       y = "cumulative probability",
       title = "E) CEU frequency = 0.3")
ggsave("plots/ooa_cond_pdf_weak_E.pdf", height=56, width=56, units="mm")

dtwf_chb_ceu_high = dtwf_chb_ceu %>% 
  group_by(CEUfreq, s) %>%
  mutate(CHBfreq = as.numeric(CHBfreq)) %>% 
  filter(as.numeric(CEUfreq) > 0.6,
         as.numeric(CEUfreq) < 0.61) %>%
  mutate(cum_prob = cumsum(prob)) %>% 
  filter(!grepl("e-4", s))
ggplot() + 
  geom_line(data=dtwf_chb_ceu_high, 
            mapping=aes(CHBfreq, cum_prob, color=s),
            show.legend = FALSE,
            size=line_size) +
  universal_theme + selection_color + ylim(0, 1.01) +
  labs(x = "frequency in CHB", 
       y = "cumulative probability",
       title = "F) CEU frequency = 0.6")
ggsave("plots/ooa_cond_pdf_weak_F.pdf", height=56, width=56, units="mm")

#################################### KL ########################################
yri_selection_KL = dtwf_yri_ceu %>%
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate_at(vars(-YRIfreq, -CEUfreq, -`0`), ~ . * (log(.) - log(`0`))) %>% 
  summarise_at(vars(-YRIfreq, -`0`), sum) %>%
  pivot_longer(-CEUfreq, names_to="s", values_to="KL") %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(CEUfreq >= 0.01, CEUfreq <= 0.99)

ggplot() + 
  geom_line(data=yri_selection_KL, 
            mapping=aes(CEUfreq, KL, color=s),
            show.legend = FALSE,
            size=line_size) + 
  universal_theme + selection_color + scale_y_log10() +
  labs(title = "A)", 
       x = "frequency in CEU", 
       y = "Kullback-Leibler divergence from neutrality")
ggsave("plots/KL_A.pdf", width=85, height=85, units="mm")

chb_selection_KL = dtwf_chb_ceu %>%
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate_at(vars(-CHBfreq, -CEUfreq, -`0`), ~ . * (log(.) - log(`0`))) %>% 
  summarise_at(vars(-CHBfreq, -`0`), sum) %>%
  pivot_longer(-CEUfreq, names_to="s", values_to="KL") %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  filter(CEUfreq >= 0.01, CEUfreq <= 0.99)

ggplot() + 
  geom_line(data=chb_selection_KL, 
            mapping=aes(CEUfreq, KL, color=s),
            show.legend = FALSE,
            size=line_size) + 
  universal_theme + selection_color + scale_y_log10() +
  labs(title = "B)", 
       x = "frequency in CEU", 
       y = "Kullback-Leibler divergence from neutrality")
ggsave("plots/KL_B.pdf", width=85, height=85, units="mm")

yri_strong_stabilizing_KL = yri_ceu %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
  pivot_wider(names_from = s, values_from = prob) %>% 
  pivot_longer(cols=-c(YRIfreq, CEUfreq, `stabilizing, hs=-5e-4`), 
               names_to = "s", values_to = "prob") %>%
  group_by(CEUfreq, s) %>% 
  summarize(KL=sum(`stabilizing, hs=-5e-4` * 
                     (log(`stabilizing, hs=-5e-4`) - log(prob)))) %>%
  mutate(n_snp = 1.96 / KL) %>% filter(grepl("neg", s)) %>%
  group_by(CEUfreq) %>% summarize(max_snp = max(n_snp)) %>%
  mutate(s = "stabilizing, hs=-5e-4")

yri_moderate_stabilizing_KL = yri_ceu %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
  pivot_wider(names_from = s, values_from = prob) %>% 
  pivot_longer(cols=-c(YRIfreq, CEUfreq, `stabilizing, hs=-2.5e-5`), 
               names_to = "s", values_to = "prob") %>%
  group_by(CEUfreq, s) %>% 
  summarize(KL=sum(`stabilizing, hs=-2.5e-5` * 
                     (log(`stabilizing, hs=-2.5e-5`) - log(prob)))) %>%
  mutate(n_snp = 1.96 / KL) %>% filter(grepl("neg", s)) %>%
  group_by(CEUfreq) %>% summarize(max_snp = max(n_snp)) %>%
  mutate(s = "stabilizing, hs=-2.5e-5")

yri_weak_stabilizing_KL = yri_ceu %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
  pivot_wider(names_from = s, values_from = prob) %>% 
  pivot_longer(cols=-c(YRIfreq, CEUfreq, `stabilizing, hs=-5e-5`), 
               names_to = "s", values_to = "prob") %>%
  group_by(CEUfreq, s) %>% 
  summarize(KL=sum(`stabilizing, hs=-5e-5` * 
                     (log(`stabilizing, hs=-5e-5`) - log(prob)))) %>%
  mutate(n_snp = 1.96 / KL) %>% filter(grepl("neg", s)) %>%
  group_by(CEUfreq) %>% summarize(max_snp = max(n_snp)) %>%
  mutate(s = "stabilizing, hs=-5e-5")

yri_stabilizing_KL = bind_rows(yri_strong_stabilizing_KL, 
                               yri_moderate_stabilizing_KL,
                               yri_weak_stabilizing_KL) %>%
  mutate(CEUfreq = as.numeric(CEUfreq))

ggplot() + 
  geom_point(data=yri_stabilizing_KL, 
             mapping=aes(CEUfreq, max_snp, color=s),
             show.legend = FALSE) + 
  universal_theme + selection_color + scale_y_log10() +
  labs(title = "B)", 
       x = "CEU frequency", 
       y = "number of SNPs to reject negative selection")
ggsave("plots/5B.pdf", width=85, height=85, units="mm")

##################### Unused KL stuff ##########################################
chb_selection_KL = chb_ceu %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate_at(vars(-CHBfreq, -CEUfreq, -`0`), ~ . * (log10(.) - log10(`0`))) %>% 
  summarise_at(vars(-CHBfreq, -`0`), sum) %>%
  pivot_longer(-CEUfreq, names_to="selection", values_to="KL") %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% 
  inner_join(ceu_filter) %>% mutate(n_snp = 2/KL) %>% 
  filter(CEUfreq >= 0.01, CEUfreq <= 0.99)

ggplot() + 
  geom_point(data=chb_selection_KL, 
             mapping=aes(CEUfreq, n_snp, color=selection),
             show.legend = FALSE) + 
  universal_theme + selection_color + scale_y_log10() +
  labs(title = "B) CHB", 
       x = "CEU frequency", 
       y = "number of SNPs to reject neutrality")
ggsave("plots/KLB.pdf", width=85, height=85, units="mm")

################################ Unused stuff ##################################
ceu_ooa_forward = ceu_ooa %>% group_by(ooAfreq, s) %>%
  summarize(mean_ceu = sum(as.numeric(CEUfreq) * prob)) %>%
  filter(!grepl("e-5", s)) %>%
  mutate(ooAfreq = as.numeric(ooAfreq)) %>% 
  inner_join(ooa_filter)
ggplot(ceu_ooa_forward, aes(ooAfreq, mean_ceu, color=s)) + 
  geom_point(show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in CEU-CHB ancestor", 
       y = "expected frequency in CEU",
       title = "B)")
ggsave("plots/3B.pdf", height=85, width=85, units="mm")

ceu_anc_forward = ceu_anc %>% group_by(ancestral_freq, s) %>%
  summarize(mean_ceu = sum(as.numeric(CEUfreq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq)) %>% 
  filter(!grepl("e-5", s))
ggplot(ceu_anc_forward, aes(ancestral_freq, mean_ceu, color=s)) + 
  geom_point(show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestor", 
       y = "expected frequency in CEU",
       title = "D)")
ggsave("plots/3D.pdf", height=85, width=85, units="mm")

################################ Unused PDFs ###################################
ggplot(yri, aes(freq, prob, color=s)) + geom_line(show.legend=FALSE) +
  universal_theme + selection_color + scale_y_log10()
ggplot(ceu, aes(freq, prob, color=s)) + geom_line(show.legend=FALSE) +
  universal_theme + selection_color + scale_y_log10()
ggplot(ooa, aes(freq, prob, color=s)) + geom_line(show.legend=FALSE) +
  universal_theme + selection_color + scale_y_log10()

ceu_forward = ceu_anc %>% group_by(ancestral_freq, s) %>%
  summarize(mean_ceu = sum(as.numeric(CEUfreq) * prob)) %>%
  ungroup %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))
ggplot(ceu_forward, aes(ancestral_freq, mean_ceu, color=s)) + 
  geom_point(size=obj_size, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestral population", 
       y = "expected frequency in CEU")

yri_ceu_summary = yri_ceu %>% group_by(CEUfreq, s) %>%
  summarize(mean_yri = sum(as.numeric(YRIfreq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq))
ggplot(yri_ceu_summary, aes(CEUfreq, mean_yri, color=s)) + 
  geom_point(size=obj_size, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "CEU frequency", 
       y = "expected YRI frequency")

chb_ceu_summary = dtwf_chb_ceu %>% group_by(CEUfreq, s) %>%
  summarize(mean_chb = sum(as.numeric(CHBfreq) * prob)) %>%
  mutate(CEUfreq = as.numeric(CEUfreq)) %>% filter(!grepl("e-5", s))
ggplot(chb_ceu_summary, aes(CEUfreq, mean_chb, color=s)) + 
  geom_point(size=point_size, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "CEU frequency", 
       y = "expected CHB frequency")

################################ CHB/YRI PDF ###################################
chb_cond = chb_ceu %>% group_by(CEUfreq, s) %>% 
  mutate(CHBfreq = as.numeric(CHBfreq)) %>%
  mutate(cum_prob = cumsum(prob))

shared_gg = c(list(geom_line(mapping=aes(CHBfreq, cum_prob, color=s), 
                             show.legend = FALSE),
                   labs(x = "CHB frequency", y = "probability")), 
              selection_color, universal_theme)

ggplot(chb_cond %>% filter(CEUfreq == 0.1)) + 
  labs(title="D) CEU frequency = 0.1") + shared_gg
ggsave("plots/3D.pdf", height=40, width=40, units="mm")

ggplot(chb_cond %>% filter(CEUfreq == 0.3)) + 
  labs(title="E) CEU frequency = 0.3") + shared_gg
ggsave("plots/3E.pdf", height=40, width=40, units="mm")

ggplot(chb_cond %>% filter(CEUfreq == 0.5)) + 
  labs(title="F) CEU frequency = 0.5") + shared_gg
ggsave("plots/3F.pdf", height=40, width=40, units="mm")

ggplot(chb_cond %>% filter(CEUfreq == 0.8)) + 
  labs(title="G) CEU frequency = 0.8") + shared_gg
ggsave("plots/3G.pdf", height=40, width=40, units="mm")

yri_cond = yri_ceu %>% group_by(CEUfreq, s) %>% 
  mutate(cum_prob = cumsum(prob)) %>% 
  mutate(YRIfreq = as.numeric(YRIfreq))

shared_gg = c(list(geom_line(mapping=aes(YRIfreq, cum_prob, color=s), show.legend = FALSE),
                   labs(x = "YRI frequency", y = "probability"), ylim(0, 1)), 
              selection_color, universal_theme)

ggplot(yri_cond %>% filter(CEUfreq == 0.1)) + 
  labs(title="H) CEU frequency = 0.1") + shared_gg
ggsave("plots/3H.pdf", height=40, width=40, units="mm")

ggplot(yri_cond %>% filter(CEUfreq == 0.3)) + 
  labs(title="I) CEU frequency = 0.3") + shared_gg
ggsave("plots/3I.pdf", height=40, width=40, units="mm")

ggplot(yri_cond %>% filter(CEUfreq == 0.5)) + 
  labs(title="J) CEU frequency = 0.5") + shared_gg
ggsave("plots/3J.pdf", height=40, width=40, units="mm")

ggplot(yri_cond %>% filter(CEUfreq == 0.8)) + 
  labs(title="K) CEU frequency = 0.8") + shared_gg
ggsave("plots/3K.pdf", height=40, width=40, units="mm")

############################## expected diff ###################################
yri_exp_diff = yri_ceu %>% group_by(CEUfreq, s) %>% 
  summarize(exp_sq_diff = sum(prob * (as.numeric(YRIfreq) - as.numeric(CEUfreq)) ^ 2)) %>%
  ungroup %>% mutate(CEUfreq = as.numeric(CEUfreq))

ggplot(yri_exp_diff, aes(CEUfreq, exp_sq_diff, color=s)) + 
  geom_point(show.legend = FALSE) +
  universal_theme + selection_color + 
  labs(y = "expected square difference in frequency",
       x = "CEU frequency",
       title = "A) YRI | CEU")
ggsave("plots/expdiffA.pdf", width=85, height=85, units="mm")

chb_exp_diff = chb_ceu %>% group_by(CEUfreq, s) %>% 
  summarize(exp_sq_diff = sum(prob * (as.numeric(CHBfreq) - as.numeric(CEUfreq)) ^ 2)) %>%
  ungroup %>% mutate(CEUfreq = as.numeric(CEUfreq))

ggplot(chb_exp_diff, aes(CEUfreq, exp_sq_diff, color=s)) + 
  geom_point(show.legend = FALSE) +
  universal_theme + selection_color + 
  labs(y = "expected square difference in frequency",
       x = "CEU frequency",
       title = "B) CHB | CEU")
ggsave("plots/expdiffB.pdf", width=85, height=85, units="mm")