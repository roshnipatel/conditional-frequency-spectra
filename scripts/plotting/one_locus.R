source("~/sherlock/oak/stabilizing_selection/scripts/plotting/shared.R")
setwd("~/Documents/selection_Genetics_2024/figures/")

N_REPS = 20

update_p = function(p, h, s) {
  q = 1 - p
  mean_fitness = 1 + 2 * p * q * h * s + p * p * s
  diff = (p * q * s * (p * (1 - h) + q * h)) / mean_fitness
  return(p + diff)
}

simulate_p = function(p, h, s, N) {
  p = update_p(p, h, s)
  count = rbinom(N, 2, p)
  new_p = sum(count) / (2 * N)
  return(new_p)
}

simulate_n_trajectories = function(n_rep, p, h, s, N, t) {
  curr_p = p
  reps = matrix(0, t, n_rep)
  for (j in 1:n_rep) {
    freqs = rep(0, t)
    p = curr_p
    for (i in 1:t) {
      freqs[i] = p
      p = simulate_p(p, h, s, N)
    }
    reps[,j] = freqs
  }
  return(reps)
}

################################### Fig 1A #####################################
p = 0.1
N = 10000 * 2
t = 2000
neutral = simulate_n_trajectories(N_REPS, p, 0.5, 0, N, t)
neutral_tib = as_tibble(neutral) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("neutral_", name),
         s = "0")

del = simulate_n_trajectories(N_REPS, p, 0.5, -0.001, N, t)
del_tib = as_tibble(del) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("neg_", name),
         s = "negative, hs=-5e-4")

underdom = simulate_n_trajectories(N_REPS, p, 5e6, -1e-10, N, t)
underdom_tib = as_tibble(underdom) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("stab_", name),
         s = "stabilizing, hs=-5e-4")

pos = simulate_n_trajectories(N_REPS, p, 0.5, 0.001, N, t)
pos_tib = as_tibble(pos) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("pos_", name),
         s = "positive, hs=+5e-4")

lowest_tib = bind_rows(neutral_tib, del_tib, underdom_tib, pos_tib) %>% 
  filter(time %% 20 == 0)
mean_lowest = lowest_tib %>% group_by(time, s) %>% 
  summarize(mean_frequency = mean(frequency))

ggplot() + 
  geom_line(data=lowest_tib, mapping=aes(time, frequency, group=name, color=s),
            show.legend = FALSE, alpha = 0.05) +
  geom_line(data=mean_lowest, mapping=aes(time, mean_frequency, group=s, color=s),
            size=line_size, show.legend = FALSE, linetype="dashed") +
  labs(x="generation", y="frequency", title="A) Trajectories beginning at 0.1") + ylim(0, 1) +
  universal_theme + selection_color
ggsave("plots/1A_alt.pdf", height=42, width=42, units="mm")

################################### Fig 1B #####################################
p = 0.5
N = 10000 * 2
t = 2000
neutral = simulate_n_trajectories(N_REPS, p, 0.5, 0, N, t)
neutral_tib = as_tibble(neutral) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("neutral_", name),
         s = "0")

del = simulate_n_trajectories(N_REPS, p, 0.5, -0.001, N, t)
del_tib = as_tibble(del) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("neg_", name),
         s = "negative, hs=-5e-4")

underdom = simulate_n_trajectories(N_REPS, p, 5e6, -1e-10, N, t)
underdom_tib = as_tibble(underdom) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("stab_", name),
         s = "stabilizing, hs=-5e-4")

pos = simulate_n_trajectories(N_REPS, p, 0.5, 0.001, N, t)
pos_tib = as_tibble(pos) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("pos_", name),
         s = "positive, hs=+5e-4")

middle_tib = bind_rows(neutral_tib, del_tib, underdom_tib, pos_tib) %>% 
  filter(time %% 20 == 0)
mean_middle = middle_tib %>% group_by(time, s) %>% 
  summarize(mean_frequency = mean(frequency))

ggplot() + 
  geom_line(data=middle_tib, mapping=aes(time, frequency, group=name, color=s),
            show.legend = FALSE, alpha = 0.05) +
  geom_line(data=mean_middle, mapping=aes(time, mean_frequency, group=s, color=s),
            size=line_size, show.legend = FALSE, linetype="dashed") +
  labs(x="generation", y="frequency", title="B) Trajectories beginning at 0.5") + ylim(0, 1) +
  universal_theme + selection_color
ggsave("plots/1B.pdf", height=42, width=42, units="mm")

################################### Fig 1C #####################################
p = 0.8
N = 10000 * 2
t = 2000
neutral = simulate_n_trajectories(N_REPS, p, 0.5, 0, N, t)
neutral_tib = as_tibble(neutral) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("neutral_", name),
         s = "0")

del = simulate_n_trajectories(N_REPS, p, 0.5, -0.001, N, t)
del_tib = as_tibble(del) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("neg_", name),
         s = "negative, hs=-5e-4")

underdom = simulate_n_trajectories(N_REPS, p, 5e6, -1e-10, N, t)
underdom_tib = as_tibble(underdom) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("stab_", name),
         s = "stabilizing, hs=-5e-4")

pos = simulate_n_trajectories(N_REPS, p, 0.5, 0.001, N, t)
pos_tib = as_tibble(pos) %>% 
  mutate(time = seq(1, t)) %>% 
  pivot_longer(-time, values_to = "frequency") %>%
  mutate(name = paste0("pos_", name),
         s = "positive, hs=+5e-4")

high_tib = bind_rows(neutral_tib, del_tib, underdom_tib, pos_tib) %>% 
  filter(time %% 20 == 0)
mean_high = high_tib %>% group_by(time, s) %>% 
  summarize(mean_frequency = mean(frequency))

ggplot() + 
  geom_line(data=high_tib, mapping=aes(time, frequency, group=name, color=s),
            show.legend = FALSE, alpha = 0.1) +
  geom_line(data=mean_high, mapping=aes(time, mean_frequency, group=s, color=s),
            size=line_size, show.legend = FALSE, linetype="dashed") +
  labs(x="generation", y="frequency", title="C) Trajectories beginning at 0.8") + ylim(0, 1) +
  universal_theme + selection_color
ggsave("plots/1C.pdf", height=42, width=42, units="mm")
