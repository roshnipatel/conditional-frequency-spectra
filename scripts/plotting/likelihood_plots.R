library(tidyverse)
library(ggpubr)
setwd("~/Documents/stabilizing_selection/plots/")

################### plotting CEU conditional likelihoods #######################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/"
test_dir = "simulated_h0.5_s0.0/"
neg_selection = list("-1e-2" = "probs_h0.5_s-1.0e-2.txt",
                     "-1e-3" = "probs_h0.5_s-0.001.txt",
                     "-1e-4" = "probs_h0.5_s-0.0001.txt",
                     "-1e-5" = "probs_h0.5_s-1.0e-5.txt")
prob_df = read_tsv(paste0(filepath_prefix, test_dir, "probs_h0.5_s0.0.txt"))
n = nrow(prob_df)
prob_df = prob_df %>% mutate(idx = seq(n), "0" = prob) %>% select(-prob)
for (coeff in names(neg_selection)) {
  df = read_tsv(paste0(filepath_prefix, test_dir, neg_selection[coeff][[1]]))
  df = df %>% mutate(idx = seq(n), "{coeff}" := prob) %>% select(-prob)
  prob_df = inner_join(prob_df, df)
}
prob_df = prob_df %>% mutate(CEUbin = round(CEUfreq, 2))

# might want to make this plot with the Big Simulations
grouped_df = prob_df %>% group_by(CEUbin) %>% 
  summarize(YRI = mean(YRIfreq),
            CHB = mean(CHBfreq)) %>%
  ungroup %>% pivot_longer(-CEUbin, names_to="population", values_to="frequency")
pop_palette = c("#870058", "#FF4000")
ggplot(grouped_df, aes(CEUbin, frequency, color=population)) + geom_point() + 
  theme_pubr() + scale_color_manual(values=pop_palette) +
  labs(x = "CEU frequency", y = "mean (expected) frequency")

# not sure that this shows anything useful - also might want to normalize by sample size
# or make this plot with the Big Simulations
grouped_df = prob_df %>% group_by(CEUbin) %>% 
  summarize(YRI = var(YRIfreq),
            CHB = var(CHBfreq)) %>%
  ungroup %>% pivot_longer(-CEUbin, names_to="population", values_to="frequency")
pop_palette = c("#870058", "#FF4000")
ggplot(grouped_df, aes(CEUbin, frequency, color=population)) + geom_point() + 
  theme_pubr() + scale_color_manual(values=pop_palette) +
  labs(x = "CEU frequency", y = "variance of frequency")

grouped_df = prob_df %>% select(-contains("freq"), -idx) %>% 
  mutate(across(-CEUbin, log)) %>% group_by(CEUbin) %>% 
  summarize(across(everything(), mean)) %>% 
  pivot_longer(cols=!CEUbin, names_to="s", values_to="prob")
# selection_palette = c("#395B50", "#5A7684", "#92AFD7", "#C5D1EB", "#FFBC42")
selection_palette = c("#395B50", "#5A7684", "#92AFD7", "#C5D1EB")
ggplot(grouped_df %>% filter(s != "0"), aes(CEUbin, prob, color=s)) + 
  geom_point() + labs(y = "mean log probability", x = "CEU frequency") +
  scale_color_manual(values=selection_palette) + theme_pubr() +
  geom_point(data = grouped_df %>% filter(s == "0"), mapping = aes(CEUbin, prob), color="#FFBC42")

##################### debugging unconditional likelihoods ######################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/uniform_binning/CEUCHBYRI.likelihood_output/"
test_dir = "simulated_h0.5_s0.0/"
test_dir = "simulated_h0.5_s-1e-4/"
test_dir = "simulated_h0.5_s-0.001/"
neg_selection = list("-1e-2" = "probs_h0.5_s-1.0e-2.txt",
                     "-1e-3" = "probs_h0.5_s-0.001.txt",
                     "-1e-4" = "probs_h0.5_s-0.0001.txt",
                     "-1e-5" = "probs_h0.5_s-1.0e-5.txt")
prob_df = read_tsv(paste0(filepath_prefix, test_dir, "probs_h0.5_s0.0.txt"))
n = nrow(prob_df)
prob_df = prob_df %>% mutate(idx = seq(n), "0" = prob) %>% select(-prob)
for (coeff in names(neg_selection)) {
  df = read_tsv(paste0(filepath_prefix, test_dir, neg_selection[coeff][[1]]))
  df = df %>% mutate(idx = seq(n), "{coeff}" := prob) %>% select(-prob)
  prob_df = inner_join(prob_df, df)
}
ancestral_probs = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/tmp/ancestral.txt",
                           col_names = c("prob")) %>% mutate(ancestral_freq = seq(0, 100) / 100)

# sanity checking E_p[log P(x)] > E_p[log Q(x)]
log_df = prob_df %>% select(-idx) %>% mutate(across(-contains("freq"), log))
grouped_df = log_df %>% group_by(ancestral_freq) %>% 
  summarize(across(everything(), mean))
ggplot(grouped_df) + geom_point(aes(ancestral_freq, `0`), color="black") +
  geom_point(aes(ancestral_freq, `-1e-3`), color="red") + theme_pubr() +
  labs(x = "ancestral frequency", y = "log-likelihood")
grouped_df = inner_join(grouped_df, ancestral_probs)
colSums(grouped_df %>% select(-contains("freq"), -prob) * grouped_df$prob)

# estimating distribution of selection coefficients
mle_df = prob_df %>% pivot_longer(cols=-c(contains("freq"), idx), names_to="s", values_to="prob") %>% 
  group_by(idx) %>% filter(prob == max(prob))
ggplot(mle_df, aes(s)) + geom_bar() + theme_pubr()

# likelihood ratio test
data = read_tsv(paste0(filepath_prefix, "simulated_h0.5_s0.0/simulated_ancestral_likelihoods_total.txt"))
mean_bootstraps = data %>% summarize(across(everything(), mean))
1 - pchisq((-2 * (mean_bootstraps - mean_bootstraps[[1]])) %>% as.numeric(), df=0)
data = data %>% mutate(idx = seq(100)) %>% select(idx, starts_with("h0.5"))
data = data %>% pivot_longer(cols=starts_with("h"), "selection")
tib = data %>% group_by(idx) %>% slice(which.max(value)) %>% ungroup %>% 
  group_by(selection) %>% summarise(n = n())

############### debugging ancestral conditional likelihoods ####################
filepath_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/uniform_binning/CEUCHBYRI_anc.likelihood_output/"
test_dir = "simulated_h0.5_s0.0/"
test_dir = "simulated_h0.5_s-1e-4/"
test_dir = "simulated_h0.5_s-0.001/"
neg_selection = list("-1e-2" = "probs_h0.5_s-1.0e-2.txt",
                     "-1e-3" = "probs_h0.5_s-0.001.txt",
                     "-1e-4" = "probs_h0.5_s-0.0001.txt",
                     "-1e-5" = "probs_h0.5_s-1.0e-5.txt")
prob_df = read_tsv(paste0(filepath_prefix, test_dir, "probs_h0.5_s0.0.txt"))
n = nrow(prob_df)
prob_df = prob_df %>% mutate(prob = replace(prob, prob == 0, 1e-100)) %>%
  mutate(idx = seq(n), "0" = prob) %>% select(-prob)
for (coeff in names(neg_selection)) {
  df = read_tsv(paste0(filepath_prefix, test_dir, neg_selection[coeff][[1]]))
  df = df %>% mutate(prob = replace(prob, prob == 0, 1e-100)) %>% 
    mutate(idx = seq(n), "{coeff}" := prob) %>% select(-prob)
  prob_df = inner_join(prob_df, df)
}

log_df = prob_df %>% select(-idx) %>% mutate(across(-contains("freq"), log))
grouped_df = log_df %>% group_by(ancestral_freq) %>% 
  summarize(across(everything(), mean))
pivoted_df = grouped_df %>% pivot_longer(cols=-contains("freq"), names_to="s", values_to="prob")
pivoted_df %>% group_by(ancestral_freq, ooAfreq, YRIfreq, CEUfreq, CHBfreq) %>% 
  slice(which.max(prob)) %>% ungroup %>% group_by(s) %>% count

selection_palette = c("#395B50", "#5A7684", "#92AFD7", "#FFBC42")
ggplot(pivoted_df %>% filter(s != '-1e-5'), aes(ancestral_freq, prob, color=s)) + 
  geom_point() + labs(y = "log-likelihood", x = "ancestral frequency") +
  scale_color_manual(values=selection_palette) + theme_pubr()
ggplot(grouped_df) + geom_point(aes(ancestral_freq, `0`), color="black") +
  geom_point(aes(ancestral_freq, `-1e-3`), color="red") + theme_pubr() +
  labs(x = "ancestral frequency", y = "log-likelihood")

########################## debugging P(x_A) ####################################
# Plot simulated neutral distribution for P(x_A)
dist = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/simulations/one_pop/h0.5_s0.0/all_muts.txt")
ggplot() + geom_density(data=dist, mapping=aes(freq), bins=100, alpha=0.5) + 
  geom_point(data=tibble(x = seq(0.01, 1, 0.01)) %>% mutate(y = 2 * 1.29e-8 * 1/x * 1e7), mapping=aes(x, y)) +
  theme_pubr()

# Plot probabilities for simulated neutral test set
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h0.5_s0.0.txt")
n = nrow(neut)
neut = neut %>% mutate(neut = prob, idx = seq(n)) %>% select(-prob)
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h0.5_s-0.001.txt") 
neg = neg %>% mutate(neg = prob, idx = seq(n)) %>% select(-prob)
neut_tib = inner_join(neut, neg)

# Plot probabilities for simulated negative test set
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/probs_h0.5_s0.0.txt")
n = nrow(neut)
neut = neut %>% mutate(neut = prob, idx = seq(n)) %>% select(-prob)
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/probs_h0.5_s-0.001.txt")
neg = neg %>% mutate(neg = prob, idx = seq(n)) %>% select(-prob)
neg_tib = inner_join(neut, neg)

# Plot probabilities for simulated weakly negative test set
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.0005/probs_h0.5_s0.0.txt")
n = nrow(neut)
neut = neut %>% mutate(neut = prob, idx = seq(n)) %>% select(-prob)
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.0005/probs_h0.5_s-0.001.txt")
neg = neg %>% mutate(neg = prob, idx = seq(n)) %>% select(-prob)
weak_neg_tib = inner_join(neut, neg)

# Compare YRI and ancestral frequency in simulated data (just a sanity check of sorts)
ggplot(neut_tib, aes(ancestral_freq, YRIfreq)) + geom_point(alpha=0.5) + theme_pubr()
ggplot(weak_neg_tib, aes(ancestral_freq, YRIfreq)) + geom_point(alpha=0.5) + theme_pubr()
ggplot(neg_tib, aes(ancestral_freq, YRIfreq)) + geom_point(alpha=0.5) + theme_pubr()

# Compare probabilities...?
ggplot(neut_tib) + geom_point(aes(ancestral_freq, neut), alpha=0.2) + geom_point(aes(ancestral_freq, neg), color="red", alpha=0.2) + theme_pubr()
ggplot(neg_tib) + geom_point(aes(ancestral_freq, neut), alpha=0.2) + geom_point(aes(ancestral_freq, neg), color="red", alpha=0.2) + theme_pubr()

ggplot(neut_tib) + geom_point(aes(YRIfreq, neut), alpha=0.2) + geom_point(aes(YRIfreq, neg), color="red", alpha=0.2) + theme_pubr()
ggplot(neg_tib) + geom_point(aes(YRIfreq, neut), alpha=0.2) + geom_point(aes(YRIfreq, neg), color="red", alpha=0.2) + theme_pubr()

###################### Likelihoods w Eur conditioning ##########################
directory_prefix = "~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/"

plot_binned_likelihoods = function(fp, plot_title) {
  data = read_tsv(fp)
  # data = data %>% mutate(idx = seq(1100))
  data = data %>% mutate(idx = seq(1100)) %>% select(freq_bin, count, idx, starts_with("h0.5"))
  data = data %>% pivot_longer(cols=starts_with("h"), "selection")
  tib = data %>% group_by(idx) %>% slice(which.max(value)) %>% ungroup %>% 
    group_by(freq_bin, selection) %>% summarise(n = n(), n_samp = mean(count)) %>% 
    ungroup %>% filter(freq_bin != 1)
  p = ggplot(tib, aes(x = freq_bin, y = n, fill=selection)) + 
    geom_bar(position = "stack", stat = "identity") + theme_pubr() +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
    scale_fill_manual(values=c("h0.5_s-0.001" = "#353A47", 
                               "h0.5_s-0.0005" = "#0C7489", 
                               "h0.5_s0.0" = "#CBDFBD")) +
                                # "h0.5_s0.0" = "#CBDFBD", 
                                # "h5000000_s-5.0e-11" = "#F19C79", 
                                # "h5000000_s-1.0e-10" = "#A44A3F")) +
    labs(x = "CEU frequency bin", y = "replicates", title = plot_title)
  return(p)
}

plot_total_likelihoods = function(fp, plot_title) {
  data = read_tsv(fp)
  # data = data %>% mutate(idx = seq(100))
  data = data %>% mutate(idx = seq(100)) %>% select(idx, starts_with("h0.5"))
  data = data %>% pivot_longer(cols=-c(idx, contains("freq")), "selection")
  tib = data %>% group_by(idx) %>% slice(which.max(value)) %>% ungroup %>% 
    group_by(selection) %>% summarise(n = n()) %>% ungroup %>%
    mutate(dummy = "hi")
  p = ggplot(tib, aes(x = dummy, y = n, fill=selection)) + 
    labs(y = "replicates", title=plot_title) + 
    geom_bar(position = "stack", stat = "identity") + theme_pubr() +
    scale_fill_manual(values=c("-0.001" = "#353A47", 
                               "-0.0005" = "#0C7489",
                               "0.0" = "#CBDFBD")) +
                               # "h0.5_s0.0" = "#CBDFBD", 
                               # "h5000000_s-5.0e-11" = "#F19C79", 
                               # "h5000000_s-1.0e-10" = "#A44A3F")) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  return(p)
}

#### Plots without stabilizing selection
neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/uniform_likelihoods_total.txt"), "neutral")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/uniform_likelihoods_total.txt"), "negative selection")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/uniform_likelihoods_total.txt"), "weakly negative")

ggarrange(neut_total, neg_total, neg_weak_total, common.legend = TRUE, align="hv", nrow = 1, ncol = 3)

neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/simulated_ancestral_likelihoods_total.txt"), "neutral")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/simulated_ancestral_likelihoods_total.txt"), "negative selection")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/simulated_ancestral_likelihoods_total.txt"), "weakly negative")

ggarrange(neut_total, neg_total, neg_weak_total, common.legend = TRUE, align="hv", nrow = 1, ncol = 3)

neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/control_ancestral_likelihoods_total.txt"), "neutral")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/control_ancestral_likelihoods_total.txt"), "negative selection")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/control_ancestral_likelihoods_total.txt"), "weakly negative")

ggarrange(neut_total, neg_total, neg_weak_total, common.legend = TRUE, align="hv", nrow = 1, ncol = 3)

neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/control_european_likelihoods_total.txt"), "neutral")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/control_european_likelihoods_total.txt"), "negative selection")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/control_european_likelihoods_total.txt"), "weakly negative")

ggarrange(neut_total, neg_total, neg_weak_total, common.legend = TRUE, align="hv", nrow = 1, ncol = 3)

#### Empirical data
control_bin_1 = plot_binned_likelihoods(paste0(directory_prefix, "empirical_control1/likelihoods_binned.txt"), "empirical control 1")
control_bin_2 = plot_binned_likelihoods(paste0(directory_prefix, "empirical_control2/likelihoods_binned.txt"), "empirical control 2")
control_bin_3 = plot_binned_likelihoods(paste0(directory_prefix, "empirical_control3/likelihoods_binned.txt"), "empirical control 3")
control_bin_4 = plot_binned_likelihoods(paste0(directory_prefix, "empirical_control4/likelihoods_binned.txt"), "empirical control 4")

control_total_1 = plot_total_likelihoods(paste0(directory_prefix, "empirical_control1/likelihoods_total.txt"), "empirical control 1")
control_total_2 = plot_total_likelihoods(paste0(directory_prefix, "empirical_control2/likelihoods_total.txt"), "empirical control 2")
control_total_3 = plot_total_likelihoods(paste0(directory_prefix, "empirical_control3/likelihoods_total.txt"), "empirical control 3")
control_total_4 = plot_total_likelihoods(paste0(directory_prefix, "empirical_control4/likelihoods_total.txt"), "empirical control 4")

gwas_bin = plot_binned_likelihoods(paste0(directory_prefix, "empirical_gwas/likelihoods_binned.txt"), "empirical gwas SNPs")
gwas_total = plot_total_likelihoods(paste0(directory_prefix, "empirical_gwas/likelihoods_total.txt"), "empirical gwas SNPs")

ggarrange(control_bin_1, control_total_1, control_bin_2, control_total_2, control_bin_3, control_total_3, control_bin_4, control_total_4, gwas_bin, gwas_total, 
          common.legend = TRUE, align="hv", nrow = 5, ncol = 2, widths = c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2))
ggsave("empirical_likelihoods.pdf", width = 8, height = 10)

#### Uniform
neut_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/uniform_likelihoods_binned.txt"), "neutral")
neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/uniform_likelihoods_total.txt"), "neutral")

neg_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/uniform_likelihoods_binned.txt"), "negative selection")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/uniform_likelihoods_total.txt"), "negative selection")

stab_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/uniform_likelihoods_binned.txt"), "stabilizing selection")
stab_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/uniform_likelihoods_total.txt"), "stabilizing selection")

neg_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/uniform_likelihoods_binned.txt"), "weakly negative")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/uniform_likelihoods_total.txt"), "weakly negative")

stab_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/uniform_likelihoods_binned.txt"), "weakly stabilizing")
stab_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/uniform_likelihoods_total.txt"), "weakly stabilizing")

ggarrange(neut_bin, neut_total + theme(plot.title = element_blank()), neg_bin, neg_total + theme(plot.title = element_blank()), neg_weak_bin, neg_weak_total + theme(plot.title = element_blank()), 
          common.legend = TRUE, align="hv", nrow = 3, ncol = 2, widths = c(1, 0.2, 1, 0.2, 1, 0.2))
ggarrange(neut_bin, neut_total, neg_bin, neg_total, stab_bin, stab_total, neg_weak_bin, neg_weak_total, stab_weak_bin, stab_weak_total, 
          common.legend = TRUE, align="hv", nrow = 5, ncol = 2, widths = c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2))
ggsave("uniform_likelihoods.pdf", width = 8, height = 10)

#### simulated_ancestral
neut_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/simulated_ancestral_likelihoods_binned.txt"), "neutral")
neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/simulated_ancestral_likelihoods_total.txt"), "neutral")

neg_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/simulated_ancestral_likelihoods_binned.txt"), "negative selection")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/simulated_ancestral_likelihoods_total.txt"), "negative selection")

stab_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/simulated_ancestral_likelihoods_binned.txt"), "stabilizing selection")
stab_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/simulated_ancestral_likelihoods_total.txt"), "stabilizing selection")

neg_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/simulated_ancestral_likelihoods_binned.txt"), "weakly negative")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/simulated_ancestral_likelihoods_total.txt"), "weakly negative")

stab_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/simulated_ancestral_likelihoods_binned.txt"), "weakly stabilizing")
stab_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/simulated_ancestral_likelihoods_total.txt"), "weakly stabilizing")

ggarrange(neut_bin, neut_total, neg_bin, neg_total, stab_bin, stab_total, neg_weak_bin, neg_weak_total, stab_weak_bin, stab_weak_total, 
          common.legend = TRUE, align="hv", nrow = 5, ncol = 2, widths = c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2))
ggsave("simulated_ancestral_likelihoods.pdf", width = 8, height = 10)

#### control_ancestral
neut_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/control_ancestral_likelihoods_binned.txt"), "neutral")
neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/control_ancestral_likelihoods_total.txt"), "neutral")

neg_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/control_ancestral_likelihoods_binned.txt"), "negative selection")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/control_ancestral_likelihoods_total.txt"), "negative selection")

stab_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/control_ancestral_likelihoods_binned.txt"), "stabilizing selection")
stab_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/control_ancestral_likelihoods_total.txt"), "stabilizing selection")

neg_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/control_ancestral_likelihoods_binned.txt"), "weakly negative")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/control_ancestral_likelihoods_total.txt"), "weakly negative")

stab_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/control_ancestral_likelihoods_binned.txt"), "weakly stabilizing")
stab_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/control_ancestral_likelihoods_total.txt"), "weakly stabilizing")

ggarrange(neut_bin, neut_total, neg_bin, neg_total, stab_bin, stab_total, neg_weak_bin, neg_weak_total, stab_weak_bin, stab_weak_total, 
          common.legend = TRUE, align="hv", nrow = 5, ncol = 2, widths = c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2))
ggsave("control_ancestral_likelihoods.pdf", width = 8, height = 10)

#### control_european
neut_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/control_european_likelihoods_binned.txt"), "neutral")
neut_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s0.0/control_european_likelihoods_total.txt"), "neutral")

neg_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/control_european_likelihoods_binned.txt"), "negative selection")
neg_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.001/control_european_likelihoods_total.txt"), "negative selection")

stab_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/control_european_likelihoods_binned.txt"), "stabilizing selection")
stab_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-1.0e-10/control_european_likelihoods_total.txt"), "stabilizing selection")

neg_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/control_european_likelihoods_binned.txt"), "weakly negative")
neg_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h0.5_s-0.0005/control_european_likelihoods_total.txt"), "weakly negative")

stab_weak_bin = plot_binned_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/uniform_likelihoods_binned.txt"), "weakly stabilizing")
stab_weak_total = plot_total_likelihoods(paste0(directory_prefix, "simulated_h5000000_s-5.0e-11/uniform_likelihoods_total.txt"), "weakly stabilizing")

ggarrange(neut_bin, neut_total, neg_bin, neg_total, stab_bin, stab_total, neg_weak_bin, neg_weak_total, stab_weak_bin, stab_weak_total, 
          common.legend = TRUE, align="hv", nrow = 5, ncol = 2, widths = c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2))
ggsave("control_european_likelihoods.pdf", width = 8, height = 10)


############################ Signal noise plots ################################
stab = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/h5000000_s-1.0e-10_YRI_CEU_meanvar.txt") %>%
  mutate(freq = seq(0, 1, 0.01))
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/h0.5_s-0.001_YRI_CEU_meanvar.txt") %>%
  mutate(freq = seq(0, 1, 0.01))
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/JFD/h0.5_s0.0_YRI_CEU_meanvar.txt") %>%
  mutate(freq = seq(0, 1, 0.01))

a = ggplot() + geom_line(data=stab, mapping=aes(freq, means), color="#A44A3F") + 
  geom_line(data=neg, mapping=aes(freq, means), color="#353A47") + 
  geom_line(data=neut, mapping=aes(freq, means), color="#CBDFBD") + 
  theme_pubr() + labs(x = "CEU frequency", y = "mean YRI frequency")

b = ggplot() + geom_line(data=stab, mapping=aes(freq, variances), color="#A44A3F") + 
  geom_line(data=neg, mapping=aes(freq, variances), color="#353A47") + 
  geom_line(data=neut, mapping=aes(freq, variances), color="#CBDFBD") + 
  theme_pubr() + labs(x = "CEU frequency", y = "YRI frequency variance")

ggarrange(a, b)
ggsave("signal_noise.pdf", width = 8, height = 4)

###################### Likelihoods ###########################
data = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/empirical_gwas/likelihoods.txt")
data = data %>% mutate(idx = seq(1, 100))
data = data %>% pivot_longer(-idx, "params") 
best = data %>% group_by(idx) %>% filter(value == max(value))
best %>% group_by(params) %>% count
ggplot(data, aes(params, value)) + geom_violin() + theme_pubr()

neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/likelihoods.txt")
neut = neut %>% mutate(idx = seq(1, 100))
neut = neut %>% pivot_longer(-idx, "params") 
best = neut %>% group_by(idx) %>% filter(value == max(value))
best %>% group_by(params) %>% count
ggplot(neut, aes(params, value)) + geom_violin() + theme_pubr()

neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/likelihoods.txt")
neg = neg %>% mutate(idx = seq(1, 100))
neg = neg %>% pivot_longer(-idx, "params") 
best = neg %>% group_by(idx) %>% filter(value == max(value))
best %>% group_by(params) %>% count
ggplot(neg, aes(params, value)) + geom_violin() + theme_pubr()

stab = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h5000000_s-1.0e-10/likelihoods.txt")
stab = stab %>% mutate(idx = seq(1, 100))
stab = stab %>% pivot_longer(-idx, "params") 
best = stab %>% group_by(idx) %>% filter(value == max(value))
best %>% group_by(params) %>% count
ggplot(stab, aes(params, value)) + geom_violin() + theme_pubr()

neg_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.0005/likelihoods.txt")
neg_weak = neg_weak %>% mutate(idx = seq(1, 100))
neg_weak = neg_weak %>% pivot_longer(-idx, "params") 
best = neg_weak %>% group_by(idx) %>% filter(value == max(value))
best %>% group_by(params) %>% count
ggplot(neg_weak, aes(params, value)) + geom_violin() + theme_pubr()

stab_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h5000000_s-5.0e-11/likelihoods.txt")
stab_weak = stab_weak %>% mutate(idx = seq(1, 100))
stab_weak = stab_weak %>% pivot_longer(-idx, "params") 
best = stab_weak %>% group_by(idx) %>% filter(value == max(value))
best %>% group_by(params) %>% count
ggplot(stab_weak, aes(params, value)) + geom_violin() + theme_pubr()

############################## Heat maps #######################################
tib = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h0.5_s0.0.txt")
ggplot(tib, aes(ancestral_freq, CEUfreq, fill=prob)) + geom_tile() + theme_pubr()

################## Neutral sims ######################
stab_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h5000000_s-5.0e-11.txt")
neg_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h0.5_s-0.0005.txt")
stab = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h5000000_s-1.0e-10.txt")
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h0.5_s-0.001.txt")
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s0.0/probs_h0.5_s0.0.txt")

sum(log(stab_weak$prob))
sum(log(neg_weak$prob))
sum(log(neg$prob))
sum(log(neut$prob))
sum(log(stab$prob))

################## Negative sims ######################
stab_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/probs_h5000000_s-5.0e-11.txt")
neg_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/probs_h0.5_s-0.0005.txt")
stab = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/probs_h5000000_s-1.0e-10.txt")
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/probs_h0.5_s-0.001.txt")
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h0.5_s-0.001/probs_h0.5_s0.0.txt")

sum(log(stab_weak$prob))
sum(log(neg_weak$prob))
sum(log(neg$prob))
sum(log(neut$prob))
sum(log(stab$prob))

################## Stabilizing sims ######################
stab_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h5000000_s-1.0e-10/probs_h5000000_s-5.0e-11.txt")
neg_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h5000000_s-1.0e-10/probs_h0.5_s-0.0005.txt")
stab = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h5000000_s-1.0e-10/probs_h5000000_s-1.0e-10.txt")
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h5000000_s-1.0e-10/probs_h0.5_s-0.001.txt")
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/simulated_h5000000_s-1.0e-10/probs_h0.5_s0.0.txt")

sum(log(stab_weak$prob))
sum(log(neg_weak$prob))
sum(log(neg$prob))
sum(log(neut$prob))
sum(log(stab$prob))

################## GWAS data ######################
stab_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/empirical_gwas/output_h5000000_s-5.0e-11.txt")
neg_weak = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/empirical_gwas/output_h0.5_s-0.0005.txt")
stab = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/empirical_gwas/output_h5000000_s-1.0e-10.txt")
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/empirical_gwas/output_h0.5_s-0.001.txt")
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/empirical_gwas/output_h0.5_s0.0.txt")
tib = bind_rows(tibble(idx = seq(1, 18357), prob = stab$prob, type="stabilizing"),
                tibble(idx = seq(1, 18357), prob = neg$prob, type="negative"),
                tibble(idx = seq(1, 18357), prob = neut$prob, type="neutral"))
ggplot(tib, aes(type, prob)) + geom_violin() + theme_pubr()

sum(log(stab_weak$prob))
sum(log(neg_weak$prob))
sum(log(neg$prob))
sum(log(neut$prob))
sum(log(stab$prob))

################## controls data ######################
stab = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/control/output_h5000000_s-1.0e-10.txt")
neg = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/control/output_h0.5_s-0.001.txt")
neut = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/data/likelihood_output/control/output_h0.5_s0.0.txt")
tib = bind_rows(tibble(idx = seq(1, 9900), prob = stab$prob, type="stabilizing"),
                tibble(idx = seq(1, 9900), prob = neg$prob, type="negative"),
                tibble(idx = seq(1, 9900), prob = neut$prob, type="neutral"))
ggplot(tib, aes(type, prob)) + geom_violin() + theme_pubr()

wide = tib %>% spread(type, prob)
sum(log(wide$negative))
sum(log(wide$neutral))
sum(log(wide$stabilizing))

########################## one pop sims #####################################
onepop_neutral = read_tsv("~/sherlock/oak/stabilizing_selection/data/simulations/one_pop/h0.5_s0.0/all_muts.txt")
onepop_negative = read_tsv("~/sherlock/oak/stabilizing_selection/data/simulations/one_pop/h0.5_s-0.001/all_muts.txt")
onepop_stabilizing = read_tsv("~/sherlock/oak/stabilizing_selection/data/simulations/one_pop/h5000000_s-1.0e-10/all_muts.txt")

a = ggplot() + geom_histogram(data=onepop_neutral %>% sample_n(200000), mapping=aes(freq), alpha = 0.5, binwidth= 0.01) + 
  geom_histogram(data=onepop_negative %>% sample_n(200000), mapping=aes(freq), fill="red", alpha = 0.5, binwidth= 0.01) + 
  geom_histogram(data=onepop_stabilizing %>% sample_n(200000), mapping=aes(freq), fill="blue", alpha = 0.5, binwidth= 0.01) +
  theme_pubr()

########################## ooA sims #####################################
neutral = read_tsv("~/sherlock/oak/stabilizing_selection/data/simulations/ooA/h0.5_s0.0/all_muts.txt")
negative = read_tsv("~/sherlock/oak/stabilizing_selection/data/simulations/ooA/h0.5_s-0.001/all_muts.txt")
stabilizing = read_tsv("~/sherlock/oak/stabilizing_selection/data/simulations/ooA/h5000000_s-1.0e-10/all_muts.txt")

b = ggplot() + geom_histogram(data=neutral %>% sample_n(700000), mapping=aes(ooAfreq), alpha = 0.5, binwidth= 0.01) + 
  geom_histogram(data=negative %>% sample_n(700000), mapping=aes(ooAfreq), fill="red", alpha = 0.5, binwidth= 0.01) + 
  geom_histogram(data=stabilizing %>% sample_n(700000), mapping=aes(ooAfreq), fill="blue", alpha = 0.5, binwidth= 0.01) +
  theme_pubr() + xlab("freq")

ggarrange(a, b)

ggplot() + geom_histogram(data=neutral %>% sample_n(700000) %>% filter(ooAfreq < 0.05), mapping=aes(YRIfreq), alpha = 0.5, binwidth= 0.01) + 
  geom_histogram(data=negative %>% sample_n(700000) %>% filter(ooAfreq < 0.05), mapping=aes(YRIfreq), fill="red", alpha = 0.5, binwidth= 0.01) + 
  geom_histogram(data=stabilizing %>% sample_n(700000) %>% filter(ooAfreq < 0.05), mapping=aes(YRIfreq), fill="blue", alpha = 0.5, binwidth= 0.01) +
  theme_pubr()

################## checking stationary dists ##########################
analytical = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/tmp/neutral_probs.txt", col_names=c("prob")) %>% mutate(idx = seq(1, 98))
sim = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/tmp/sims.txt", col_names=c("prob")) %>% mutate(idx = seq(0, 100))
minor = read_tsv("~/sherlock/scratch/stabilizing_selection/data_analysis/tmp/neutral_minor.txt", col_names=c("prob")) %>% mutate(idx = seq(1, 49))

minor = minor %>% mutate(adj_prob = prob / sum(minor$prob))
sim = sim %>% mutate(flip = 100-idx) %>% mutate(minor = pmin(idx, flip))
grouped_sim = sim %>% group_by(minor) %>% summarize(prob = sum(prob))

ggplot() + geom_point(data=analytical, mapping = aes(idx, prob)) + geom_point(data=sim, mapping = aes(idx, prob), color="red") + theme_pubr()
ggplot() + geom_point(data=minor, mapping = aes(idx, adj_prob)) + geom_point(data=grouped_sim, mapping = aes(minor, prob), color="red") + theme_pubr()
ggplot() + geom_point(data=analytical, mapping = aes(idx, prob)) + theme_pubr()
