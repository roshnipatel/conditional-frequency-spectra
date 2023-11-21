source("~/sherlock/oak/stabilizing_selection/scripts/plotting/shared.R")

################ Load simulation-based 3d matrix: CHB, YRI | CEU ###############
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/CHBYRI_CEU.jouganous_wo_migration.likelihood_output/"
prob_df = read_tsv(paste0(filepath_prefix, "probs_h0.5_s0.0.txt")) %>% 
  mutate(s = "0") %>% mutate(prob = replace(prob, prob == 0, 1e-100))
selection = list("negative, hs=-5e-4" = "probs_h0.5_s-1.0e-3.txt",
                 "negative, hs=-2.5e-5" = "probs_h0.5_s-5.0e-4.txt",
                 "negative, hs=-5e-5" = "probs_h0.5_s-1.0e-4.txt",
                 "stabilizing, hs=-5e-4" = "probs_h5e6_s-1.0e-10.txt",
                 "stabilizing, hs=-2.5e-5" = "probs_h5e6_s-5.0e-11.txt",
                 "stabilizing, hs=-5e-5" = "probs_h5e6_s-1.0e-11.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]))
  df = df %>% mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
    mutate(s = coeff)
  prob_df = bind_rows(prob_df, df)
}
prob_df = prob_df %>% mutate(CEUfreq = round(CEUfreq, 2), 
                             CHBfreq = round(CHBfreq, 2),
                             YRIfreq = round(YRIfreq, 2))

############################# P(loss) analysis #################################
loss_df = inner_join(prob_df %>% filter(s == 0) %>% group_by(CEUfreq) %>% 
                       filter((YRIfreq < 0.05) | (YRIfreq > 0.95)) %>% 
                       summarize(YRI = sum(prob)),
                     prob_df %>% filter(s == 0) %>% group_by(CEUfreq) %>% 
                       filter((CHBfreq < 0.05) | (CHBfreq > 0.95)) %>% 
                       summarize(CHB = sum(prob))) %>%
  pivot_longer(-CEUfreq, values_to="prob")
ggplot(loss_df, aes(CEUfreq, prob, color=name)) + geom_point(size=0.5) + 
  universal_theme
ggsave("~/Documents/selection_2024/plots/yri_vs_chb_loss.pdf", height=70, width=70, units="mm")

loss_df = prob_df %>% group_by(s, CEUfreq) %>% 
  filter((YRIfreq < 0.05) | (YRIfreq > 0.95)) %>% 
  summarize(probability = sum(prob)) %>%
  filter(!grepl("-5e-5", s))
ggplot(loss_df, aes(CEUfreq, probability, color=s)) + geom_point(size=0.5) + 
  selection_color + universal_theme + labs(x = "CEU frequency")
ggsave("~/Documents/selection_2024/plots/yri_loss.pdf", height=70, width=70, units="mm")

loss_df = prob_df %>% group_by(s, CEUfreq) %>% 
  filter((CHBfreq < 0.05) | (CHBfreq > 0.95)) %>% 
  summarize(probability = sum(prob)) %>%
  filter(!grepl("-5e-5", s))
ggplot(loss_df, aes(CEUfreq, probability, color=s)) + geom_point(size=0.5) + 
  selection_color + universal_theme + labs(x = "CEU frequency")
ggsave("~/Documents/selection_2024/plots/chb_loss.pdf", height=70, width=70, units="mm")

############################ Median YRI | ancestral ############################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/jouganous_wo_migration_"
jfd_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_YRI_ancestral.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(across(everything(), cumsum)) %>%
  mutate(YRIfreq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-YRIfreq, names_to="ancestral_freq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_YRI_ancestral.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_YRI_ancestral.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_YRI_ancestral.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_YRI_ancestral.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_YRI_ancestral.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_YRI_ancestral.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(across(everything(), cumsum)) %>%
    mutate(YRIfreq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-YRIfreq, names_to="ancestral_freq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
    mutate(s = coeff)
  jfd_df = bind_rows(jfd_df, df)
}
YRI_median_df = jfd_df %>% group_by(ancestral_freq, s) %>% 
  filter(prob > 0.5) %>% 
  summarize(median_YRI = (min(YRIfreq))) %>%
  filter(!grepl("e-5", s)) %>% ungroup %>% mutate(ancestral_freq = as.numeric(ancestral_freq))
ggplot(YRI_median_df, aes(ancestral_freq, median_YRI, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestral population", y = "median YRI frequency")
ggsave("~/Documents/selection_2024/plots/YRI_ancestral_median.pdf", width=50, height=50, units="mm")

ggplot(YRI_median_df %>% filter(!grepl("stab", s)), aes(ancestral_freq, median_YRI, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestral population", y = "median YRI frequency")
ggsave("~/Documents/selection_2024/plots/YRI_ancestral_median_wo_stab.pdf", width=50, height=50, units="mm")

ggplot(YRI_median_df %>% filter(s == "0"), aes(ancestral_freq, median_YRI, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestral population", y = "median YRI frequency")
ggsave("~/Documents/selection_2024/plots/YRI_ancestral_median_wo_sel.pdf", width=50, height=50, units="mm")

############################ Median CEU | ancestral ############################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/demography_"
jfd_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_CEU_ancestral.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(across(everything(), cumsum)) %>%
  mutate(CEUfreq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-CEUfreq, names_to="ancestral_freq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_CEU_ancestral.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_CEU_ancestral.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_CEU_ancestral.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_CEU_ancestral.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_CEU_ancestral.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_CEU_ancestral.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(across(everything(), cumsum)) %>%
    mutate(CEUfreq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-CEUfreq, names_to="ancestral_freq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
    mutate(s = coeff)
  jfd_df = bind_rows(jfd_df, df)
}
CEU_median_df = jfd_df %>% group_by(ancestral_freq, s) %>% 
  filter(prob > 0.5) %>% 
  summarize(median_CEU = (min(CEUfreq))) %>%
  filter(!grepl("e-5", s)) %>% ungroup %>% mutate(ancestral_freq = as.numeric(ancestral_freq))
ggplot(CEU_median_df, aes(ancestral_freq, median_CEU, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "frequency in ancestral population", y = "median CEU frequency")
ggsave("~/Documents/selection_2024/plots/CEU_ancestral_median.pdf", width=50, height=50, units="mm")

##################### Median analysis of ancestor | CEU ########################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/jouganous_wo_migration_"
jfd_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_ancestral_YRI.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(across(everything(), cumsum)) %>%
  mutate(ancestral_freq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-ancestral_freq, names_to="YRIfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_ancestral_YRI.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_ancestral_YRI.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_ancestral_YRI.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_ancestral_YRI.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_ancestral_YRI.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_ancestral_YRI.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(across(everything(), cumsum)) %>%
    mutate(ancestral_freq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-ancestral_freq, names_to="YRIfreq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
    mutate(s = coeff)
  jfd_df = bind_rows(jfd_df, df)
}
shared_gg = c(selection_color, universal_theme, 
              list(geom_line(mapping=aes(ancestral_freq, prob, color=s), 
                             size=0.8, show.legend = FALSE),
                   labs(x = "ancestral frequency", y = "cumulative probability")))
plotting_df = jfd_df %>% group_by(YRIfreq, s) %>% 
  filter(prob > 0.5) %>% 
  summarize(median_ancestral = (min(ancestral_freq))) %>%
  filter(!grepl("e-5", s)) %>% ungroup %>% mutate(YRIfreq = as.numeric(YRIfreq))
ggplot(plotting_df, aes(YRIfreq, median_ancestral, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "YRI frequency", y = "median frequency in ancestral population")

jfd_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_ooA_CEU.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(across(everything(), cumsum)) %>%
  mutate(ooAfreq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-ooAfreq, names_to="CEUfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_ooA_CEU.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_ooA_CEU.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_ooA_CEU.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_ooA_CEU.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_ooA_CEU.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_ooA_CEU.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(across(everything(), cumsum)) %>%
    mutate(ooAfreq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-ooAfreq, names_to="CEUfreq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
    mutate(s = coeff)
  jfd_df = bind_rows(jfd_df, df)
}
shared_gg = c(selection_color, universal_theme, 
              list(geom_line(mapping=aes(ooAfreq, prob, color=s), 
                             size=0.8, show.legend = FALSE),
                   labs(x = "ooA frequency", y = "cumulative probability")))
plotting_df = jfd_df %>% group_by(CEUfreq, s) %>% 
  filter(prob > 0.5) %>% 
  summarize(median_ooA = (min(ooAfreq))) %>%
  filter(!grepl("e-5", s)) %>% ungroup %>% mutate(CEUfreq = as.numeric(CEUfreq))
ggplot(plotting_df, aes(CEUfreq, median_ooA, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "CEU frequency", y = "median frequency in ooA population")

####################### Mean analysis of ancestor | CEU ########################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/no_dem_"
jfd_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_ancestral_CEU.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(ancestral_freq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-ancestral_freq, names_to="CEUfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_ancestral_CEU.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_ancestral_CEU.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_ancestral_CEU.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_ancestral_CEU.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_ancestral_CEU.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_ancestral_CEU.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(ancestral_freq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-ancestral_freq, names_to="CEUfreq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
    mutate(s = coeff)
  jfd_df = bind_rows(jfd_df, df)
}
mean_df = jfd_df %>% group_by(CEUfreq, s) %>% 
  summarize(mean_ancestral = sum(ancestral_freq * prob))
shared_gg = c(selection_color, universal_theme, 
              list(geom_line(mapping=aes(ancestral_freq, prob, color=s), 
                             size=0.8, show.legend = FALSE),
                   labs(x = "ancestral frequency", y = "cumulative probability")))
plotting_df = mean_df %>% ungroup %>% filter(!grepl("e-5", s)) %>% 
  mutate(CEUfreq = as.numeric(CEUfreq))
a = ggplot(plotting_df, aes(CEUfreq, mean_ancestral, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "CEU frequency", y = "mean frequency in ancestral population")

jfd_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_ooA_CEU.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(ooAfreq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-ooAfreq, names_to="CEUfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_ooA_CEU.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_ooA_CEU.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_ooA_CEU.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_ooA_CEU.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_ooA_CEU.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_ooA_CEU.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(ooAfreq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-ooAfreq, names_to="CEUfreq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
    mutate(s = coeff)
  jfd_df = bind_rows(jfd_df, df)
}
mean_df = jfd_df %>% group_by(CEUfreq, s) %>% 
  summarize(mean_ooA = sum(ooAfreq * prob))
shared_gg = c(selection_color, universal_theme, 
              list(geom_line(mapping=aes(ooAfreq, prob, color=s), 
                             size=0.8, show.legend = FALSE),
                   labs(x = "ooA frequency", y = "cumulative probability")))
plotting_df = mean_df %>% ungroup %>% filter(!grepl("e-5", s)) %>% 
  mutate(CEUfreq = as.numeric(CEUfreq))
b = ggplot(plotting_df, aes(CEUfreq, mean_ooA, color=s)) + 
  geom_point(size=0.5, show.legend = FALSE) +
  universal_theme + selection_color + ylim(0, 1) +
  labs(x = "CEU frequency", y = "mean frequency in ooA population")

ggarrange(a, b)

###################### Plotting simulation KL divergence #######################
# For CEU, YRI | anc
KL_divergence = prob_df %>% 
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(ancestral_freq) %>% 
  mutate(across(-c(CEUfreq, YRIfreq, `0`), log10)) %>% 
  mutate_at(vars(-CEUfreq, -YRIfreq, -ancestral_freq, -`0`), ~ log10(`0`) - .) %>% 
  mutate_at(vars(-CEUfreq, -YRIfreq, -ancestral_freq, -`0`), ~ `0` * .) %>% 
  summarise_at(vars(-CEUfreq, -YRIfreq, -`0`), sum) %>%
  pivot_longer(-ancestral_freq, names_to="selection", values_to="KL")

ggplot(KL_divergence, aes(ancestral_freq, KL, color=selection)) + 
  geom_point() + theme_pubr() + labs(title = "CEU, YRI | anc") + selection_color +
  scale_y_log10()

# For CHB, YRI | CEU
sum_df = prob_df %>% pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% summarize_at(vars(-CHBfreq, -YRIfreq), sum) %>% 
  pivot_longer(-CEUfreq)
ggplot(sum_df, aes(CEUfreq, value, color=name)) + 
  geom_point(size=0.5) + universal_theme + labs(title = "CHB, YRI | CEU") + selection_color

KL_divergence = prob_df %>% 
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate(across(-c(CHBfreq, YRIfreq, `0`), log10)) %>% 
  mutate_at(vars(-CHBfreq, -YRIfreq, -CEUfreq, -`0`), ~ log10(`0`) - .) %>% 
  mutate_at(vars(-CHBfreq, -YRIfreq, -CEUfreq, -`0`), ~ `0` * .) %>% 
  summarise_at(vars(-CHBfreq, -YRIfreq, -`0`), sum) %>%
  pivot_longer(-CEUfreq, names_to="selection", values_to="KL")

ggplot(KL_divergence %>% filter(grepl("negative", selection)), aes(CEUfreq, KL, color=selection)) + 
  geom_point(size=0.5) + labs(title = "CHB, YRI | CEU") + 
  selection_color + universal_theme + labs(x = "CEU frequency")
ggsave("~/Documents/selection_2024/plots/sims_KL.pdf", height=85, width=85, units="mm")

# Stabilizing vs negative
KL_divergence = prob_df %>% 
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate(`stabilizing, hs=-5e-4` = 
           `negative, hs=-5e-4` * (log10(`negative, hs=-5e-4`) - log10(`stabilizing, hs=-5e-4`)),
         `stabilizing, hs=-2.5e-5` = 
           `negative, hs=-2.5e-5` * (log10(`negative, hs=-2.5e-5`) - log10(`stabilizing, hs=-2.5e-5`)),
         `stabilizing, hs=-5e-5` = 
           `negative, hs=-5e-5` * (log10(`negative, hs=-5e-5`) - log10(`stabilizing, hs=-5e-5`))) %>% 
  summarise_at(vars(starts_with("stabilizing")), sum) %>%
  pivot_longer(-CEUfreq, names_to="selection", values_to="KL")

ggplot(KL_divergence, aes(CEUfreq, KL, color=selection)) + 
  geom_point(size=0.5) + labs(title = "CHB, YRI | CEU") + 
  selection_color + universal_theme + labs(x = "CEU frequency")
ggsave("~/Documents/selection_2024/plots/sims_KL_stab_v_neg.pdf", height=85, width=85, units="mm")

chb_df = prob_df %>% group_by(CEUfreq, CHBfreq, s) %>% 
  summarize(prob = sum(prob))

KL_divergence = chb_df %>% 
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate(across(-c(CHBfreq, `0`), log10)) %>% 
  mutate_at(vars(-CHBfreq, -CEUfreq, -`0`), ~ log10(`0`) - .) %>% 
  mutate_at(vars(-CHBfreq, -CEUfreq, -`0`), ~ `0` * .) %>% 
  summarise_at(vars(-CHBfreq, -`0`), sum) %>%
  pivot_longer(-CEUfreq, names_to="selection", values_to="KL")

ggplot() + 
  geom_point(data=KL_divergence, mapping=aes(CEUfreq, KL, color=selection)) + 
  theme_pubr() + labs(title = "CHB | CEU") + selection_color + 
  geom_point(data=chb_KL, mapping=aes(median, KL), color="red")

yri_df = prob_df %>% group_by(CEUfreq, YRIfreq, s) %>% 
  summarize(prob = sum(prob))

KL_divergence = yri_df %>% 
  pivot_wider(names_from = s, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate(across(-c(YRIfreq, `0`), log10)) %>% 
  mutate_at(vars(-YRIfreq, -CEUfreq, -`0`), ~ log10(`0`) - .) %>% 
  mutate_at(vars(-YRIfreq, -CEUfreq, -`0`), ~ `0` * .) %>% 
  summarise_at(vars(-YRIfreq, -`0`), sum) %>%
  pivot_longer(-CEUfreq, names_to="selection", values_to="KL")

ggplot() + 
  geom_point(data=KL_divergence, mapping=aes(CEUfreq, KL, color=selection)) + 
  theme_pubr() + labs(title = "YRI | CEU") + selection_color + 
  geom_point(data=yri_KL, mapping=aes(median, KL), color="red")

########################## Using tree.py to plot PDFs ##########################
chb_df = prob_df %>% group_by(CEUfreq, CHBfreq, s) %>% 
  summarize(prob = sum(prob))
chb_df = chb_df %>% filter(!grepl("e-5", s))

shared_gg = c(list(geom_line(mapping=aes(CHBfreq, prob, color=s), show.legend = FALSE),
                   theme_pubr(),
                   theme(plot.title = element_text(size = 6),
                         axis.title = element_text(size = 6),
                         axis.text = element_text(size = 6),
                         legend.title = element_text(size = 6),
                         legend.text = element_text(size = 6),
                         plot.title.position = "plot"),
                   labs(x = "CHB frequency", y = "probability")), selection_color)

a = ggplot(chb_df %>% filter(CEUfreq == 0.01)) + labs(title="CEU frequency = 0.01") + shared_gg
a
ggsave("~/Documents/selection_2024/plots/chb_sims_01.pdf", height=50, width=50, units="mm")

b = ggplot(chb_df %>% filter(CEUfreq == 0.1)) + labs(title="CEU frequency = 0.1") + shared_gg

c = ggplot(chb_df %>% filter(CEUfreq == 0.2)) + labs(title="CEU frequency = 0.2") + shared_gg

d = ggplot(chb_df %>% filter(CEUfreq == 0.3)) + labs(title="CEU frequency = 0.3") + shared_gg

e = ggplot(chb_df %>% filter(CEUfreq == 0.4)) + labs(title="CEU frequency = 0.4") + shared_gg

f = ggplot(chb_df %>% filter(CEUfreq == 0.5)) + labs(title="CEU frequency = 0.5") + shared_gg

g = ggplot(chb_df %>% filter(CEUfreq == 0.6)) + labs(title="CEU frequency = 0.6") + shared_gg

h = ggplot(chb_df %>% filter(CEUfreq == 0.7)) + labs(title="CEU frequency = 0.7") + shared_gg

i = ggplot(chb_df %>% filter(CEUfreq == 0.8)) + labs(title="CEU frequency = 0.8") + shared_gg
i
ggsave("~/Documents/selection_2024/plots/chb_sims_8.pdf", height=50, width=50, units="mm")

j = ggplot(chb_df %>% filter(CEUfreq == 0.9)) + labs(title="CEU frequency = 0.9") + shared_gg

ggarrange(a, b, c, d, e, f, g, h, i, j, ncol=5, nrow=2, common.legend = TRUE, align="hv")
ggsave("~/Documents/selection_2024/plots/chb_sims.pdf", height=70, width=170, units="mm")

yri_df = prob_df %>% group_by(CEUfreq, YRIfreq, s) %>% 
  summarize(prob = sum(prob))
yri_df = yri_df %>% filter(!grepl("e-5", s)) 

shared_gg = c(list(geom_line(mapping=aes(YRIfreq, prob, color=s), show.legend=FALSE),
                   theme_pubr(),
                   theme(plot.title = element_text(size = 6),
                         axis.title = element_text(size = 6),
                         axis.text = element_text(size = 6),
                         legend.title = element_text(size = 6),
                         legend.text = element_text(size = 6),
                         plot.title.position = "plot"),
                   labs(x = "YRI frequency", y = "probability")), selection_color)

a = ggplot(yri_df %>% filter(CEUfreq == 0.01)) + labs(title="CEU frequency = 0.01") + shared_gg
a
ggsave("~/Documents/selection_2024/plots/yri_sims_01.pdf", height=50, width=50, units="mm")

b = ggplot(yri_df %>% filter(CEUfreq == 0.1)) + labs(title="CEU frequency = 0.1") + shared_gg

c = ggplot(yri_df %>% filter(CEUfreq == 0.2)) + labs(title="CEU frequency = 0.2") + shared_gg

d = ggplot(yri_df %>% filter(CEUfreq == 0.3)) + labs(title="CEU frequency = 0.3") + shared_gg

e = ggplot(yri_df %>% filter(CEUfreq == 0.4)) + labs(title="CEU frequency = 0.4") + shared_gg

f = ggplot(yri_df %>% filter(CEUfreq == 0.5)) + labs(title="CEU frequency = 0.5") + shared_gg

g = ggplot(yri_df %>% filter(CEUfreq == 0.6)) + labs(title="CEU frequency = 0.6") + shared_gg

h = ggplot(yri_df %>% filter(CEUfreq == 0.7)) + labs(title="CEU frequency = 0.7") + shared_gg
h
ggsave("~/Documents/selection_2024/plots/yri_sims_7.pdf", height=50, width=50, units="mm")

i = ggplot(yri_df %>% filter(CEUfreq == 0.8)) + labs(title="CEU frequency = 0.8") + shared_gg
i
ggsave("~/Documents/selection_2024/plots/yri_sims_8.pdf", height=50, width=50, units="mm")

j = ggplot(yri_df %>% filter(CEUfreq == 0.9)) + labs(title="CEU frequency = 0.9") + shared_gg
j
ggsave("~/Documents/selection_2024/plots/yri_sims_9.pdf", height=50, width=50, units="mm")

ggarrange(a, b, c, d, e, f, g, h, i, j, ncol=5, nrow=2, common.legend = TRUE, align="hv")
ggsave("~/Documents/selection_2024/plots/yri_sims.pdf", height=70, width=170, units="mm")

####################### Plotting control vs simulation #########################
# First load prob_df, the simulation-based 3D matrix: CHB, YRI | CEU
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/"
control_chb = read_tsv(paste0(filepath_prefix, "JFD/control.conditionUKB_WBfreq,B.popCHB.txt"), 
                       col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(CHBfreq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-CHBfreq, names_to="CEUfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(type="empirical") %>% mutate(CEUfreq = as.double(CEUfreq)) %>%
  mutate(CEUfreq = round(CEUfreq, 2), CHBfreq = round(CHBfreq, 2))
chb_df = prob_df %>% group_by(CEUfreq, CHBfreq, s) %>% 
  summarize(prob = sum(prob)) %>% filter(s == "0") %>% 
  mutate(type="simulated") %>% transmute(CEUfreq, CHBfreq, type, prob)
control_chb = bind_rows(control_chb, chb_df)
control_yri = read_tsv(paste0(filepath_prefix, "JFD/control.conditionUKB_WBfreq,B.popYRI.txt"), 
                       col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(YRIfreq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-YRIfreq, names_to="CEUfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(type="empirical") %>% mutate(CEUfreq = as.double(CEUfreq)) %>%
  mutate(CEUfreq = round(CEUfreq, 2), YRIfreq = round(YRIfreq, 2))
yri_df = prob_df %>% group_by(CEUfreq, YRIfreq, s) %>% 
  summarize(prob = sum(prob)) %>% filter(s == "0") %>% 
  mutate(type="simulated") %>% transmute(CEUfreq, YRIfreq, type, prob)
control_yri = bind_rows(control_yri, yri_df)

KL_CHB = control_chb %>% 
  pivot_wider(names_from = type, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate(across(-c(CHBfreq, simulated), log10)) %>% 
  mutate_at(vars(-CHBfreq, -CEUfreq, -simulated), ~ log10(simulated) - .) %>% 
  mutate_at(vars(-CHBfreq, -CEUfreq, -simulated), ~ simulated * .) %>% 
  summarise_at(vars(-CHBfreq, -simulated), sum) %>%
  pivot_longer(-CEUfreq, names_to="type", values_to="KL")

a = ggplot(KL_CHB, aes(CEUfreq, KL)) + 
  geom_point() + theme_pubr() + labs(title = "CHB | CEU") +
  ylim(0, 1)

KL_YRI = control_yri %>% 
  pivot_wider(names_from = type, values_from = prob) %>% 
  group_by(CEUfreq) %>% 
  mutate(across(-c(YRIfreq, simulated), log10)) %>% 
  mutate_at(vars(-YRIfreq, -CEUfreq, -simulated), ~ log10(simulated) - .) %>% 
  mutate_at(vars(-YRIfreq, -CEUfreq, -simulated), ~ simulated * .) %>% 
  summarise_at(vars(-YRIfreq, -simulated), sum) %>%
  pivot_longer(-CEUfreq, names_to="type", values_to="KL")

b = ggplot(KL_YRI, aes(CEUfreq, KL)) + 
  geom_point() + theme_pubr() + labs(title = "YRI | CEU") +
  ylim(0, 1)

ggarrange(a, b)

shared_gg = list(geom_line(mapping = aes(CHBfreq, prob, color=type), size=0.8),
                 theme_pubr())

a = ggplot(control_chb %>% filter(CEUfreq == 0.01)) + shared_gg + 
  labs(title="CEUfreq = 0.01")

b = ggplot(control_chb %>% filter(CEUfreq == 0.1)) + shared_gg + 
  labs(title="CEUfreq = 0.1")

c = ggplot(control_chb %>% filter(CEUfreq == 0.2)) + shared_gg + 
  labs(title="CEUfreq = 0.2")

d = ggplot(control_chb %>% filter(CEUfreq == 0.3)) + shared_gg + 
  labs(title="CEUfreq = 0.3")

e = ggplot(control_chb %>% filter(CEUfreq == 0.4)) + shared_gg + 
  labs(title="CEUfreq = 0.4")

f = ggplot(control_chb %>% filter(CEUfreq == 0.5)) + shared_gg + 
  labs(title="CEUfreq = 0.5")

g = ggplot(control_chb %>% filter(CEUfreq == 0.6)) + shared_gg + 
  labs(title="CEUfreq = 0.6")

h = ggplot(control_chb %>% filter(CEUfreq > 0.69, CEUfreq < 0.71)) + shared_gg + 
  labs(title="CEUfreq = 0.7")

i = ggplot(control_chb %>% filter(CEUfreq == 0.8)) + shared_gg + 
  labs(title="CEUfreq = 0.8")

j = ggplot(control_chb %>% filter(CEUfreq == 0.9)) + shared_gg + 
  labs(title="CEUfreq = 0.9")

ggarrange(a, b, c, d, e, f, g, h, i, j, ncol=5, nrow=2, common.legend = TRUE, align="hv")

shared_gg = list(geom_line(mapping = aes(YRIfreq, prob, color=type), size=0.8),
                 theme_pubr())

a = ggplot(control_yri %>% filter(CEUfreq == 0.01)) + shared_gg + 
  labs(title="CEUfreq = 0.01")

b = ggplot(control_yri %>% filter(CEUfreq == 0.1)) + shared_gg + 
  labs(title="CEUfreq = 0.1")

c = ggplot(control_yri %>% filter(CEUfreq == 0.2)) + shared_gg + 
  labs(title="CEUfreq = 0.2")

d = ggplot(control_yri %>% filter(CEUfreq == 0.3)) + shared_gg + 
  labs(title="CEUfreq = 0.3")

e = ggplot(control_yri %>% filter(CEUfreq == 0.4)) + shared_gg + 
  labs(title="CEUfreq = 0.4")

f = ggplot(control_yri %>% filter(CEUfreq == 0.5)) + shared_gg + 
  labs(title="CEUfreq = 0.5")

g = ggplot(control_yri %>% filter(CEUfreq == 0.6)) + shared_gg + 
  labs(title="CEUfreq = 0.6")

h = ggplot(control_yri %>% filter(CEUfreq > 0.69, CEUfreq < 0.71)) + shared_gg + 
  labs(title="CEUfreq = 0.7")

i = ggplot(control_yri %>% filter(CEUfreq == 0.8)) + shared_gg + 
  labs(title="CEUfreq = 0.8")

j = ggplot(control_yri %>% filter(CEUfreq == 0.9)) + shared_gg + 
  labs(title="CEUfreq = 0.9")

ggarrange(a, b, c, d, e, f, g, h, i, j, ncol=5, nrow=2, common.legend = TRUE, align="hv")

########################### Using jfd.py to plot CDFs ##########################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/demography_"
prob_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_ancestral_CEU.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(across(everything(), cumsum)) %>%
  mutate(ancestral_freq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-ancestral_freq, names_to="CEUfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_ancestral_CEU.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_ancestral_CEU.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_ancestral_CEU.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_ancestral_CEU.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_ancestral_CEU.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_ancestral_CEU.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(across(everything(), cumsum)) %>%
    mutate(ancestral_freq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-ancestral_freq, names_to="CEUfreq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
    mutate(s = coeff)
  prob_df = bind_rows(prob_df, df)
}
shared_gg = c(selection_color, universal_theme, 
              list(geom_line(mapping=aes(ancestral_freq, prob, color=s), 
                             size=0.8, show.legend = FALSE),
                   labs(x = "ancestral frequency", y = "cumulative probability")))
plotting_df = prob_df %>% filter(!grepl("e-5", s))

a = ggplot(plotting_df %>% filter(CEUfreq == 0.01)) + shared_gg +
  labs(title="CEU frequency = 0.01")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_01.pdf", width=35, height=35, units="mm")
b = ggplot(plotting_df %>% filter(CEUfreq == 0.1)) + shared_gg +
  labs(title="CEU frequency = 0.1")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_02.pdf", width=35, height=35, units="mm")
c = ggplot(plotting_df %>% filter(CEUfreq == 0.2)) + shared_gg +
  labs(title="CEU frequency = 0.2")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_03.pdf", width=35, height=35, units="mm")
d = ggplot(plotting_df %>% filter(CEUfreq == 0.3)) + shared_gg +
  labs(title="CEU frequency = 0.3")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_04.pdf", width=35, height=35, units="mm")
e = ggplot(plotting_df %>% filter(CEUfreq == 0.4)) + shared_gg +
  labs(title="CEU frequency = 0.4")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_05.pdf", width=35, height=35, units="mm")
f = ggplot(plotting_df %>% filter(CEUfreq == 0.5)) + shared_gg +
  labs(title="CEU frequency = 0.5")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_06.pdf", width=35, height=35, units="mm")
g = ggplot(plotting_df %>% filter(CEUfreq == 0.6)) + shared_gg +
  labs(title="CEU frequency = 0.6")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_07.pdf", width=35, height=35, units="mm")
h = ggplot(plotting_df %>% filter(CEUfreq == 0.7)) + shared_gg +
  labs(title="CEU frequency = 0.7")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_08.pdf", width=35, height=35, units="mm")
i = ggplot(plotting_df %>% filter(CEUfreq == 0.8)) + shared_gg +
  labs(title="CEU frequency = 0.8")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_09.pdf", width=35, height=35, units="mm")
j = ggplot(plotting_df %>% filter(CEUfreq == 0.9)) + shared_gg +
  labs(title="CEU frequency = 0.9")
ggsave("~/Documents/selection_2024/plots/ancestral_CEU_10.pdf", width=35, height=35, units="mm")

ggarrange(a, b, c, d, e, f, g, h, i, j, ncol=5, nrow=2, align="hv")

prob_df = read_tsv(paste0(filepath_prefix, "h0.5_s0.0_CHB_ooA.txt"), 
                   col_names = as.character(seq(0, 1, 0.01))) %>% 
  mutate(across(everything(), cumsum)) %>%
  mutate(CHBfreq = seq(0, 1, 0.01)) %>% 
  pivot_longer(-CHBfreq, names_to="ooAfreq", values_to="prob") %>%
  mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
  mutate(s = "0")
selection = list("negative, hs=-5e-4" = "h0.5_s-1.0e-3_CHB_ooA.txt",
                 "negative, hs=-2.5e-5" = "h0.5_s-5.0e-4_CHB_ooA.txt",
                 "negative, hs=-5e-5" = "h0.5_s-1.0e-4_CHB_ooA.txt",
                 "stabilizing, hs=-5e-4" = "h5e6_s-1.0e-10_CHB_ooA.txt",
                 "stabilizing, hs=-2.5e-5" = "h5e6_s-5.0e-11_CHB_ooA.txt",
                 "stabilizing, hs=-5e-5" = "h5e6_s-1.0e-11_CHB_ooA.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]), 
                col_names = as.character(seq(0, 1, 0.01))) %>%
    mutate(across(everything(), cumsum)) %>%
    mutate(CHBfreq = seq(0, 1, 0.01)) %>% 
    pivot_longer(-CHBfreq, names_to="ooAfreq", values_to="prob") %>%
    mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
    mutate(s = coeff)
  prob_df = bind_rows(prob_df, df)
}
shared_gg = c(selection_color, universal_theme, 
              list(geom_line(mapping=aes(CHBfreq, prob, color=s), 
                             size=0.8, show.legend = FALSE),
                   labs(x = "CHB frequency", y = "cumulative probability")))
plotting_df = prob_df %>% filter(!grepl("e-5", s))

a = ggplot(plotting_df %>% filter(ooAfreq == 0.01)) + shared_gg +
  labs(title="ooA frequency = 0.01")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_01.pdf", width=35, height=35, units="mm")
b = ggplot(plotting_df %>% filter(ooAfreq == 0.1)) + shared_gg +
  labs(title="ooA frequency = 0.1")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_02.pdf", width=35, height=35, units="mm")
c = ggplot(plotting_df %>% filter(ooAfreq == 0.2)) + shared_gg +
  labs(title="ooA frequency = 0.2")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_03.pdf", width=35, height=35, units="mm")
d = ggplot(plotting_df %>% filter(ooAfreq == 0.3)) + shared_gg +
  labs(title="ooA frequency = 0.3")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_04.pdf", width=35, height=35, units="mm")
e = ggplot(plotting_df %>% filter(ooAfreq == 0.4)) + shared_gg +
  labs(title="ooA frequency = 0.4")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_05.pdf", width=35, height=35, units="mm")
f = ggplot(plotting_df %>% filter(ooAfreq == 0.5)) + shared_gg +
  labs(title="ooA frequency = 0.5")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_06.pdf", width=35, height=35, units="mm")
g = ggplot(plotting_df %>% filter(ooAfreq == 0.6)) + shared_gg +
  labs(title="ooA frequency = 0.6")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_07.pdf", width=35, height=35, units="mm")
h = ggplot(plotting_df %>% filter(ooAfreq == 0.7)) + shared_gg +
  labs(title="ooA frequency = 0.7")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_08.pdf", width=35, height=35, units="mm")
i = ggplot(plotting_df %>% filter(ooAfreq == 0.8)) + shared_gg +
  labs(title="ooA frequency = 0.8")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_09.pdf", width=35, height=35, units="mm")
j = ggplot(plotting_df %>% filter(ooAfreq == 0.9)) + shared_gg +
  labs(title="ooA frequency = 0.9")
ggsave("~/Documents/selection_2024/plots/CHB_ooA_10.pdf", width=35, height=35, units="mm")

ggarrange(a, b, c, d, e, f, g, h, i, j, ncol=5, nrow=2, align="hv")

######################### Using tree.py to plot CDFs ###########################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/CHB_CEU.likelihood_output/"
prob_df = read_tsv(paste0(filepath_prefix, "probs_h0.5_s0.0.txt")) %>% 
  mutate(s = "0")  %>% pivot_wider(names_from=CEUfreq, values_from=prob) %>% 
  mutate(across(-c(CHBfreq, s), cumsum)) %>% 
  pivot_longer(cols=-c(CHBfreq, s), names_to="CEUfreq", values_to="cumulative_prob")
selection = list("negative, hs=-5e-4" = "probs_h0.5_s-0.001.txt",
                 "stabilizing, hs=-5e-4" = "probs_h5000000_s-1.0e-10.txt")
for (coeff in names(selection)) {
  df = read_tsv(paste0(filepath_prefix, selection[coeff][[1]]))
  df = df %>% mutate(s = coeff) %>% pivot_wider(names_from=CEUfreq, values_from=prob) %>% 
    mutate(across(-c(CHBfreq, s), cumsum)) %>% 
    pivot_longer(cols=-c(CHBfreq, s), names_to="CEUfreq", values_to="cumulative_prob")
  prob_df = bind_rows(prob_df, df)
}
shared_gg = list(geom_line(mapping=aes(CHBfreq, cumulative_prob, color=s), size=0.8),
                 theme_pubr(), 
                 scale_color_manual(values=c("0" = "#FFBC42",
                                             "negative, hs=-5e-4" = "#395B50",
                                             "negative, hs=-2.5e-5" = "#5A7684",
                                             "negative, hs=-5e-5" = "#92AFD7",
                                             "stabilizing, hs=-5e-4" = "#870058",
                                             "stabilizing, hs=-2.5e-5" = "#E3879E")))
a = ggplot(prob_df %>% filter(CEUfreq == 0.01)) + labs(title="CEUfreq = 0.01") + shared_gg

b = ggplot(prob_df %>% filter(CEUfreq == 0.1)) + labs(title="CEUfreq = 0.1") + shared_gg

c = ggplot(prob_df %>% filter(CEUfreq == 0.2)) + labs(title="CEUfreq = 0.2") + shared_gg

d = ggplot(prob_df %>% filter(CEUfreq == 0.3)) + labs(title="CEUfreq = 0.3") + shared_gg

e = ggplot(prob_df %>% filter(CEUfreq == 0.4)) + labs(title="CEUfreq = 0.4") + shared_gg

f = ggplot(prob_df %>% filter(CEUfreq == 0.5)) + labs(title="CEUfreq = 0.5") + shared_gg

g = ggplot(prob_df %>% filter(CEUfreq == 0.6)) + labs(title="CEUfreq = 0.6") + shared_gg

h = ggplot(prob_df %>% filter(CEUfreq == 0.7)) + labs(title="CEUfreq = 0.7") + shared_gg

i = ggplot(prob_df %>% filter(CEUfreq == 0.8)) + labs(title="CEUfreq = 0.8") + shared_gg

j = ggplot(prob_df %>% filter(CEUfreq == 0.9)) + labs(title="CEUfreq = 0.9") + shared_gg

ggarrange(a, b, c, d, e, f, g, h, i, j, ncol=5, nrow=2, common.legend = TRUE, align="hv")
