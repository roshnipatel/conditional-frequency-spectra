source("~/sherlock/oak/stabilizing_selection/scripts/plotting/shared.R")
setwd("~/Documents/selection_Genetics_2024/figures/")

################################### Fig 1D #####################################
sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
sims_p1_anc = load_two_pop_dist(sims_filepath_prefix, "p1_conditional_ancestor.txt")
sims_p1_anc_mean = sims_p1_anc %>% group_by(ancestral_freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
dtwf_p1_anc = load_two_pop_dist(dtwf_filepath_prefix, "p1_conditional_ancestor.txt")
dtwf_p1_anc_mean = dtwf_p1_anc %>% 
  filter(ancestral_freq < 1) %>%
  group_by(ancestral_freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

ggplot() + 
  geom_line(data=dtwf_p1_anc_mean, mapping=aes(ancestral_freq, mean_p1, color=s),
            size=line_size, show.legend = FALSE) +
  # geom_point(data=sims_p1_anc_mean, mapping=aes(ancestral_freq, mean_p1, color=s),
  #            size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestor", y = "expected descendant\nfrequency",
       title = "D)")
ggsave("plots/1D.pdf", height=42, width=42, units="mm")

################################### Fig 2B #####################################
sims_filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/"
sims_anc_p1 = load_two_pop_dist(sims_filepath_prefix, "ancestor_conditional_p1.txt")
sims_anc_p1_mean = sims_anc_p1 %>% 
  filter(p1freq < 1) %>%
  group_by(p1freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p1freq = as.numeric(p1freq))

dtwf_filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0/"
dtwf_anc_p1 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p1.txt")
dtwf_anc_p1_mean = dtwf_anc_p1 %>% 
  filter(p1freq < 1) %>%
  group_by(p1freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p1freq = as.numeric(p1freq))

ggplot() + 
  geom_line(data=dtwf_anc_p1_mean, 
            mapping=aes(p1freq, mean_anc, color=s, linetype=s),
            size=line_size, show.legend = FALSE) +
  # geom_point(data=sims_anc_p1_mean, mapping=aes(p1freq, mean_anc, color=s),
  #            size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color +
  scale_linetype_manual(values=c("0" = "solid",
                                 "positive, hs=+5e-4" = "35",
                                 "stabilizing, hs=-5e-4" = "solid",
                                 "negative, hs=-5e-4" = "1331")) +
  labs(x = "descendant frequency", y = "expected frequency in ancestor",
       title = "B) Equilibrium demography")
ggsave("plots/2B.pdf", height=55, width=55, units="mm")

################################### Fig 2C #####################################
sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
sims_anc_p2 = load_two_pop_dist(sims_filepath_prefix, "ancestor_conditional_p2.txt")
sims_anc_p2_mean = sims_anc_p2 %>% 
  filter(p2freq < 1) %>%
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq))

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
dtwf_anc_p2_mean = dtwf_anc_p2 %>% 
  filter(p2freq < 1) %>%
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq))

ggplot() + 
  geom_line(data=dtwf_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
            size=line_size, show.legend = FALSE) +
  # geom_point(data=sims_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
  #            size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "descendant frequency", y = "expected frequency in ancestor",
       title = "C) Bottleneck")
ggsave("plots/2C.pdf", height=55, width=55, units="mm")

################################### Fig 2D #####################################
sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0.001/")
sims_anc_p2 = load_two_pop_dist(sims_filepath_prefix, "ancestor_conditional_p2.txt")
sims_anc_p2_mean = sims_anc_p2 %>% 
  filter(p2freq < 1) %>%
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq))

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0.001/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
dtwf_anc_p2_mean = dtwf_anc_p2 %>% 
  filter(p2freq < 1) %>%
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq))

ggplot() + 
  geom_line(data=dtwf_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
            size=line_size, show.legend = FALSE) +
  # geom_point(data=sims_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
  #            size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "descendant frequency", y = "expected frequency in ancestor",
       title = "D) Bottleneck + exponential growth")
ggsave("plots/2D.pdf", height=55, width=55, units="mm")


################################### Fig 2E #####################################
sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.001/")
sims_anc_p2 = load_two_pop_dist(sims_filepath_prefix, "ancestor_conditional_p2.txt")
sims_anc_p2_mean = sims_anc_p2 %>% 
  filter(p2freq < 1) %>%
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq))

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.001/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
dtwf_anc_p2_mean = dtwf_anc_p2 %>% 
  filter(p2freq < 1) %>%
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq))

ggplot() + 
  geom_line(data=dtwf_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
            size=line_size, show.legend = FALSE) +
  # geom_point(data=sims_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
  #            size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "descendant frequency", y = "expected frequency in ancestor",
       title = "E) Exponential growth")
ggsave("plots/2E.pdf", height=55, width=55, units="mm")

# Needs to be the last panel in this file bc it uses files loaded by other panels
################################### Fig 2A #####################################
backward_dtwf_ex = dtwf_anc_p1 %>% 
  mutate(ancestral_freq = as.numeric(ancestral_freq),
         p1freq = as.numeric(p1freq)) %>%
  filter(grepl("tive", s)) %>% 
  filter(p1freq > 0.49, p1freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 50)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=backward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s, linetype=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  theme(axis.text = element_blank()) +
  scale_linetype_manual(values=c("0" = "solid",
                                 "positive, hs=+5e-4" = "35",
                                 "stabilizing, hs=-5e-4" = "solid",
                                 "negative, hs=-5e-4" = "1331")) +
  labs(x = "frequency in\nancestor", 
       y = "probability")
ggsave("plots/2A-1.pdf", height=30, width=30, units="mm")

forward_dtwf_ex = dtwf_p1_anc %>% 
  mutate(p1freq = as.numeric(p1freq),
         ancestral_freq = as.numeric(ancestral_freq)) %>%
  filter(grepl("tive", s)) %>% 
  filter(p1freq > 0.49, p1freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 50)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=forward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  theme(axis.text = element_blank()) +
  labs(x = "frequency in\nancestor", 
       y = "likelihood")
ggsave("plots/2A-2.pdf", height=30, width=30, units="mm")

ancestral = bind_rows(read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s-1.0e-3_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "negative, hs=-5e-4", count = seq(1, 47442)),
                      read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s+1.0e-3_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "positive, hs=+5e-4", count = seq(1, 47442)))
ancestral = ancestral %>% mutate(freq = count / 47442)
ggplot(ancestral, aes(freq, prob, color=s)) + 
  geom_line(show.legend=FALSE) + 
  universal_theme + selection_color + 
  theme(axis.text = element_blank()) +
  scale_y_log10(limits=c(1e-27, 1)) +
  labs(x = "frequency in\nancestor", y = "log-scaled\nprobability")
ggsave("plots/2A-3.pdf", height=30, width=30, units="mm")

################################# SUPPLEMENT ###################################
################################ Two-pop forward ###############################
### Bottleneck (A)
sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
sims_p2_anc = load_two_pop_dist(sims_filepath_prefix, "p2_conditional_ancestor.txt")
sims_p2_anc_mean = sims_p2_anc %>% group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
dtwf_p2_anc_mean = dtwf_p2_anc %>% 
  filter(ancestral_freq < 1) %>%
  group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

ggplot() + 
  geom_line(data=dtwf_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
            size=line_size, show.legend = FALSE) +
  geom_point(data=sims_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
             size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestor", y = "expected descendant frequency",
       title = expression("A) Bottleneck (0.1N"['e']*")"))
ggsave("plots/two_pop_forward_A.pdf", height=56, width=56, units="mm")

### Growth (B)
sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.001/")
sims_p2_anc = load_two_pop_dist(sims_filepath_prefix, "p2_conditional_ancestor.txt")
sims_p2_anc_mean = sims_p2_anc %>% group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.001/")
dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
dtwf_p2_anc_mean = dtwf_p2_anc %>% 
  filter(ancestral_freq < 1) %>%
  group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

ggplot() + 
  geom_line(data=dtwf_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
            size=line_size, show.legend = FALSE) +
  geom_point(data=sims_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
             size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestor", y = "expected descendant frequency",
       title = expression("B) Exponential growth (r = 0.1%)"))
ggsave("plots/two_pop_forward_B.pdf", height=56, width=56, units="mm")

### Bottleneck + growth (C)
sims_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0.001/")
sims_p2_anc = load_two_pop_dist(sims_filepath_prefix, "p2_conditional_ancestor.txt")
sims_p2_anc_mean = sims_p2_anc %>% group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0.001/")
dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
dtwf_p2_anc_mean = dtwf_p2_anc %>% 
  filter(ancestral_freq < 1) %>%
  group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

ggplot() + 
  geom_line(data=dtwf_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
            size=line_size, show.legend = FALSE) +
  geom_point(data=sims_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
             size=0.75, show.legend = FALSE, shape=1) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestor", y = "expected descendant frequency",
       title = expression("C) Bottleneck + growth (0.1N"['e']*", r=0.1%)"))
ggsave("plots/two_pop_forward_C.pdf", height=56, width=56, units="mm")

### Bottleneck (D)
dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne3e3_growth0/")
dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
dtwf_p2_anc_mean = dtwf_p2_anc %>% 
  filter(ancestral_freq < 1) %>%
  group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

ggplot() + 
  geom_line(data=dtwf_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
            size=line_size, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestor", y = "expected descendant frequency",
       title = expression("D) Bottleneck (0.3N"['e']*")"))
ggsave("plots/two_pop_forward_D.pdf", height=56, width=56, units="mm")

### Growth (E)
dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.0005/")
dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
dtwf_p2_anc_mean = dtwf_p2_anc %>% 
  filter(ancestral_freq < 1) %>%
  group_by(ancestral_freq, s) %>%
  summarize(mean_p2 = sum(as.numeric(p2freq) * prob)) %>%
  mutate(ancestral_freq = as.numeric(ancestral_freq))

ggplot() + 
  geom_line(data=dtwf_p2_anc_mean, mapping=aes(ancestral_freq, mean_p2, color=s),
            size=line_size, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestor", y = "expected descendant frequency",
       title = expression("E) Exponential growth (r = 0.05%)"))
ggsave("plots/two_pop_forward_E.pdf", height=56, width=56, units="mm")

########################### Alternative demographies ###########################
### Weaker bottleneck (A)
dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne3e3_growth0/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
dtwf_anc_p2_mean = dtwf_anc_p2 %>% 
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>%
  filter(p2freq < 1)

ggplot() + 
  geom_line(data=dtwf_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
            size=line_size, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in descendant", y = "expected frequency in ancestor",
       title = expression("A) Bottleneck (0.3N"['e']*")"))
ggsave("plots/alt_demo_A.pdf", height=85, width=85, units="mm")

### Weaker growth backward (B)
dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.0005/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
dtwf_anc_p2_mean = dtwf_anc_p2 %>% 
  group_by(p2freq, s) %>%
  summarize(mean_anc = sum(as.numeric(ancestral_freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>%
  filter(p2freq < 1)

ggplot() + 
  geom_line(data=dtwf_anc_p2_mean, mapping=aes(p2freq, mean_anc, color=s),
            size=line_size, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in descendant", y = "expected frequency in ancestor",
       title = expression("B) Exponential growth (r = 0.05%)"))
ggsave("plots/alt_demo_B.pdf", height=85, width=85, units="mm")

################################ Demography CFS ################################
dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
dtwf_p1_p2 = load_two_pop_dist(dtwf_filepath_prefix, "p1_conditional_p2.txt")

dtwf_p1_p2_mean = dtwf_p1_p2 %>% group_by(p2freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>% 
  filter(p2freq < 1)
ggplot(dtwf_p1_p2_mean) + 
  geom_line(mapping=aes(p2freq, mean_p1, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "population 2 frequency", y = "expected population 1 frequency", 
       title = expression("A) Bottleneck (0.1N"['e']*")"))
ggsave("plots/demo_cfs_A.pdf", height=56, width=56, units="mm")

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.001/")
dtwf_p1_p2 = load_two_pop_dist(dtwf_filepath_prefix, "p1_conditional_p2.txt")

dtwf_p1_p2_mean = dtwf_p1_p2 %>% group_by(p2freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>% 
  filter(p2freq < 1)
ggplot(dtwf_p1_p2_mean) + 
  geom_line(mapping=aes(p2freq, mean_p1, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "population 2 frequency", y = "expected population 1 frequency", 
       title = expression("B) Exponential growth (r = 0.1%)"))
ggsave("plots/demo_cfs_B.pdf", height=56, width=56, units="mm")

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0.001/")
dtwf_p1_p2 = load_two_pop_dist(dtwf_filepath_prefix, "p1_conditional_p2.txt")

dtwf_p1_p2_mean = dtwf_p1_p2 %>% group_by(p2freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>% 
  filter(p2freq < 1)
ggplot(dtwf_p1_p2_mean) + 
  geom_line(mapping=aes(p2freq, mean_p1, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "population 2 frequency", y = "expected population 1 frequency", 
       title = expression("C) Bottleneck + growth (0.1N"['e']*", r=0.1%)"))
ggsave("plots/demo_cfs_C.pdf", height=56, width=56, units="mm")

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne3e3_growth0/")
dtwf_p1_p2 = load_two_pop_dist(dtwf_filepath_prefix, "p1_conditional_p2.txt")

dtwf_p1_p2_mean = dtwf_p1_p2 %>% group_by(p2freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>% 
  filter(p2freq < 1)
ggplot(dtwf_p1_p2_mean) + 
  geom_line(mapping=aes(p2freq, mean_p1, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "population 2 frequency", y = "expected population 1 frequency", 
       title = expression("D) Bottleneck (0.3N"['e']*")"))
ggsave("plots/demo_cfs_D.pdf", height=56, width=56, units="mm")

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.0005/")
dtwf_p1_p2 = load_two_pop_dist(dtwf_filepath_prefix, "p1_conditional_p2.txt")

dtwf_p1_p2_mean = dtwf_p1_p2 %>% group_by(p2freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>% 
  filter(p2freq < 1)
ggplot(dtwf_p1_p2_mean) + 
  geom_line(mapping=aes(p2freq, mean_p1, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "population 2 frequency", y = "expected population 1 frequency", 
       title = expression("E) Exponential growth (r = 0.05%)"))
ggsave("plots/demo_cfs_E.pdf", height=56, width=56, units="mm")

######################## Equilibrium CFS, heterozygosity #######################
dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0/")
dtwf_p1_p2 = load_two_pop_dist(dtwf_filepath_prefix, "p1_conditional_p2.txt")

dtwf_p1_p2_mean = dtwf_p1_p2 %>% group_by(p2freq, s) %>%
  summarize(mean_p1 = sum(as.numeric(p1freq) * prob)) %>%
  mutate(p2freq = as.numeric(p2freq)) %>% 
  filter(p2freq < 1)
ggplot(dtwf_p1_p2_mean) + 
  geom_line(mapping=aes(p2freq, mean_p1, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "population 2 frequency", y = "expected population 1 frequency", 
       title = "A) Conditional frequency spectrum under equilibrium demography")
ggsave("plots/equilibrium_A.pdf", height=85, width=85, units="mm")

dtwf_p1_p2_port = dtwf_p1_p2 %>% 
  group_by(s, p2freq) %>% 
  mutate(p1freq = as.numeric(p1freq)) %>%
  summarize(expected_het = sum(prob * 2 * p1freq * (1 - p1freq))) %>% 
  mutate(p2freq = as.numeric(p2freq)) %>%
  filter(p2freq > 0, p2freq < 1) %>%
  mutate(p2_het = 2 * p2freq * (1 - p2freq))
ggplot(dtwf_p1_p2_port) + 
  geom_line(mapping=aes(p2freq, p2_het), 
            linetype="dotted", show.legend = FALSE) +
  geom_line(mapping=aes(p2freq, expected_het, color=s),
            show.legend = FALSE) + 
  selection_color + universal_theme + 
  labs(x = "population 2 frequency", y = "expected population 1 heterozygosity", 
       title = "B) Expected heterozygosity under equilibrium demography")
ggsave("plots/equilibrium_B.pdf", height=85, width=85, units="mm")

#################### Backward intuition across demographies ####################
dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
backward_dtwf_ex = dtwf_anc_p2 %>% 
  mutate(ancestral_freq = as.numeric(ancestral_freq),
         p2freq = as.numeric(p2freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 100)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=backward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s, linetype=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  scale_linetype_manual(values=c("0" = "solid",
                                 "positive, hs=+5e-4" = "35",
                                 "stabilizing, hs=-5e-4" = "solid",
                                 "negative, hs=-5e-4" = "1331")) +
  labs(x = "frequency in\nancestor", 
       y = "probability",
       title = "Posterior (backward transition)")
ggsave("plots/backward_explained_A.pdf", height=45, width=45, units="mm")

dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
forward_dtwf_ex = dtwf_p2_anc %>% 
  mutate(p2freq = as.numeric(p2freq),
         ancestral_freq = as.numeric(ancestral_freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 50)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=forward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  labs(x = "frequency in\nancestor", 
       y = "likelihood",
       title = "Likelihood (forward transition)")
ggsave("plots/backward_explained_B.pdf", height=45, width=45, units="mm")

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
backward_dtwf_ex = dtwf_anc_p2 %>% 
  mutate(ancestral_freq = as.numeric(ancestral_freq),
         p2freq = as.numeric(p2freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 100)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=backward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  labs(x = "frequency in\nancestor", 
       y = "probability",
       title = "Posterior (backward transition)")
ggsave("plots/backward_explained_C.pdf", height=45, width=45, units="mm")

dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
forward_dtwf_ex = dtwf_p2_anc %>% 
  mutate(p2freq = as.numeric(p2freq),
         ancestral_freq = as.numeric(ancestral_freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 50)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=forward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  labs(x = "frequency in\nancestor", 
       y = "likelihood",
       title = "Likelihood (forward transition)") 
ggsave("plots/backward_explained_D.pdf", height=45, width=45, units="mm")

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0.001/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
backward_dtwf_ex = dtwf_anc_p2 %>% 
  mutate(ancestral_freq = as.numeric(ancestral_freq),
         p2freq = as.numeric(p2freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 100)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=backward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  labs(x = "frequency in\nancestor", 
       y = "probability",
       title = "Posterior (backward transition)") 
ggsave("plots/backward_explained_E.pdf", height=45, width=45, units="mm")

dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
forward_dtwf_ex = dtwf_p2_anc %>% 
  mutate(p2freq = as.numeric(p2freq),
         ancestral_freq = as.numeric(ancestral_freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 50)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=forward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  labs(x = "frequency in\nancestor", 
       y = "likelihood",
       title = "Likelihood (forward transition)")
ggsave("plots/backward_explained_F.pdf", height=45, width=45, units="mm")

dtwf_filepath_prefix = paste0("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/dtwf-generation2e3_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.001/")
dtwf_anc_p2 = load_two_pop_dist(dtwf_filepath_prefix, "ancestor_conditional_p2.txt")
backward_dtwf_ex = dtwf_anc_p2 %>% 
  mutate(ancestral_freq = as.numeric(ancestral_freq),
         p2freq = as.numeric(p2freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 100)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=backward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  labs(x = "frequency in\nancestor", 
       y = "probability",
       title = "Posterior (backward transition)") 
ggsave("plots/backward_explained_G.pdf", height=45, width=45, units="mm")

dtwf_p2_anc = load_two_pop_dist(dtwf_filepath_prefix, "p2_conditional_ancestor.txt")
forward_dtwf_ex = dtwf_p2_anc %>% 
  mutate(p2freq = as.numeric(p2freq),
         ancestral_freq = as.numeric(ancestral_freq)) %>%
  filter(p2freq > 0.49, p2freq < 0.51) %>% 
  mutate(bin = cut_interval(ancestral_freq, n = 50)) %>% # this accounts for uneven bin sizes
  group_by(bin, s) %>% 
  summarize(ancestral_freq = min(ancestral_freq), prob = sum(prob))
ggplot() + 
  geom_line(data=forward_dtwf_ex, 
            mapping=aes(ancestral_freq, prob, color=s),
            show.legend = FALSE) +
  universal_theme + selection_color +
  labs(x = "frequency in\nancestor", 
       y = "likelihood",
       title = "Likelihood (forward transition)") 
ggsave("plots/backward_explained_H.pdf", height=45, width=45, units="mm")

ancestral = bind_rows(read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h5e6_s-1.0e-10_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "stabilizing, hs=-5e-4", count = seq(1, 47442)),
                      read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s-1.0e-3_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "negative, hs=-5e-4", count = seq(1, 47442)),
                      read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s0.0_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "0", count = seq(1, 47442)),
                      read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/distributions/ancestral/jouganous_wo_migration_h0.5_s+1.0e-3_ancestral_sfs_count_probs.txt", col_names = c("prob")) %>% mutate(s = "positive, hs=+5e-4", count = seq(1, 47442)))
ancestral = ancestral %>% mutate(freq = count / 47442)
ggplot(ancestral, aes(freq, prob, color=s)) + 
  geom_line(show.legend=FALSE) + 
  universal_theme + selection_color + 
  scale_y_log10(limits=c(1e-27, 1)) +
  labs(x = "frequency in ancestor", y = "probability",
       title = "Prior")
ggsave("plots/backward_explained_I.pdf", height=45, width=45, units="mm")
