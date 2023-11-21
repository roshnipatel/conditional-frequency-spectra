library(tidyverse)
library(ggpubr)

# Color palettes
selection_color = list(scale_color_manual(values=c("0" = "#0D0630",
                                                   "positive, hs=+5e-4" = "#D90368",
                                                   "positive, hs=+2.5e-5" = "#E6579A",
                                                   "positive, hs=+5e-5" = "#F2ABCD",
                                                   "stabilizing, hs=-5e-4" = "#06A77D",
                                                   "stabilizing, hs=-2.5e-5" = "#59C4A8",
                                                   "stabilizing, hs=-5e-5" = "#ACE2D4",
                                                   "negative, hs=-5e-4" = "#3772FF",
                                                   "negative, hs=-2.5e-5" = "#7AA1FF",
                                                   "negative, hs=-5e-5" = "#BCD0FF")))

pop_color = list(scale_color_manual(values=c("gwas" = "#FBB13C",
                                             "matched" = "#0D0630")),
                 scale_fill_manual(values=c("gwas" = "#FBB13C",
                                             "matched" = "#0D0630")))

# Font sizes
paper_theme = list(theme_pubr(),
                   theme(plot.title = element_text(size = 8),
                         axis.title = element_text(size = 8),
                         axis.text = element_text(size = 6),
                         legend.title = element_text(size = 6),
                         legend.text = element_text(size = 6),
                         axis.line = element_line(size=0.5),
                         axis.ticks = element_line(size=0.5),
                         plot.title.position = "plot"))
digital_theme = list(theme_pubr(),
                     theme(plot.title.position = "plot"))
universal_theme = paper_theme

# Figure sizes
point_size = 0.5
line_size = 0.5

# Functions to load probability distributions
parse_header = function(filepath) {
  con = file(filepath, "r")
  header = readLines(con, n=3)
  close(con)
  axis0 = word(header[2], 2)
  axis0 = substr(axis0, 1, nchar(axis0) - 1)
  bins0 = as.numeric(str_split(word(header[2], 3),",")[[1]])
  axis1 = word(header[3], 2)
  axis1 = substr(axis1, 1, nchar(axis1) - 1)
  bins1 = as.numeric(str_split(word(header[3], 3),",")[[1]])
  return(list(axis0=axis0, axis1=axis1, bins0=bins0, bins1=bins1))
}

prob_df_helper = function(filepath) {
  axes = parse_header(filepath)
  df = read_delim(filepath, col_names = as.character(axes$bins1), 
                  skip=3, delim=" ") %>% 
    mutate(axis0 = as.character(axes$bins0)) %>% 
    pivot_longer(-axis0, names_to="axis1", values_to="prob")
  df = df %>% rename_at(vars(c("axis0", "axis1")), ~ c(axes$axis0, axes$axis1))
  return(df)
}

load_two_pop_dist = function(file_prefix, file_suffix) {
  df = prob_df_helper(paste0(file_prefix, "h0.5_s0.0_", file_suffix)) %>% mutate(s = "0")
  selection = list("stabilizing, hs=-5e-4" = paste0("h5e6_s-1.0e-10_", file_suffix),
                   "negative, hs=-5e-4" = paste0("h0.5_s-1.0e-3_", file_suffix),
                   "positive, hs=+5e-4" = paste0("h0.5_s+1.0e-3_", file_suffix))
  for (coeff in names(selection)) {
    temp = prob_df_helper(paste0(file_prefix, selection[coeff][[1]])) %>%
      mutate(s = coeff)
    df = bind_rows(df, temp)
  }
  return(df)
}

load_ooa_dist = function(file_prefix, file_suffix) {
  df = prob_df_helper(paste0(file_prefix, "h0.5_s0.0_", file_suffix)) %>% mutate(s = "0")
  selection = list("negative, hs=-2.5e-5" = paste0("h0.5_s-5.0e-4_", file_suffix),
                   "negative, hs=-5e-5" = paste0("h0.5_s-1.0e-4_", file_suffix),
                   "stabilizing, hs=-5e-4" = paste0("h5e6_s-1.0e-10_", file_suffix),
                   "stabilizing, hs=-2.5e-5" = paste0("h5e6_s-5.0e-11_", file_suffix),
                   "stabilizing, hs=-5e-5" = paste0("h5e6_s-1.0e-11_", file_suffix),
                   "negative, hs=-5e-4" = paste0("h0.5_s-1.0e-3_", file_suffix),
                   "positive, hs=+5e-4" = paste0("h0.5_s+1.0e-3_", file_suffix),
                   "positive, hs=+2.5e-5" = paste0("h0.5_s+5.0e-4_", file_suffix),
                   "positive, hs=+5e-5" = paste0("h0.5_s+1.0e-4_", file_suffix))
  for (coeff in names(selection)) {
    temp = prob_df_helper(paste0(file_prefix, selection[coeff][[1]])) %>%
      mutate(s = coeff)
    df = bind_rows(df, temp)
  }
  return(df)
}

marginal_helper = function(filepath, skip=2) {
  con = file(filepath, "r")
  header = readLines(con, n=2)
  close(con)
  bins = as.numeric(str_split(word(header[2], 3),",")[[1]])
  df = read_delim(filepath, col_names = c("prob"), 
                  skip=skip, delim=" ") %>% 
    mutate(freq = bins)
  return(df)
}

load_ooa_marginal = function(file_prefix, file_suffix) {
  df = marginal_helper(paste0(file_prefix, "h0.5_s0.0_", file_suffix)) %>%
    mutate(s = "0")
  selection = list("negative, hs=-2.5e-5" = paste0("h0.5_s-5.0e-4_", file_suffix),
                   "negative, hs=-5e-5" = paste0("h0.5_s-1.0e-4_", file_suffix),
                   "stabilizing, hs=-5e-4" = paste0("h5e6_s-1.0e-10_", file_suffix),
                   "stabilizing, hs=-2.5e-5" = paste0("h5e6_s-5.0e-11_", file_suffix),
                   "stabilizing, hs=-5e-5" = paste0("h5e6_s-1.0e-11_", file_suffix),
                   "negative, hs=-5e-4" = paste0("h0.5_s-1.0e-3_", file_suffix),
                   "positive, hs=+5e-4" = paste0("h0.5_s+1.0e-3_", file_suffix),
                   "positive, hs=+2.5e-5" = paste0("h0.5_s+5.0e-4_", file_suffix),
                   "positive, hs=+5e-5" = paste0("h0.5_s+1.0e-4_", file_suffix))
  for (coeff in names(selection)) {
    temp = marginal_helper(paste0(file_prefix, selection[coeff][[1]])) %>%
      mutate(s = coeff)
    df = bind_rows(df, temp)
  }
  return(df)
}

load_two_pop_marginal = function(file_prefix, file_suffix) {
  df = marginal_helper(paste0(file_prefix, "h0.5_s0.0_", file_suffix)) %>%
    mutate(s = "0")
  selection = list("stabilizing, hs=-5e-4" = paste0("h5e6_s-1.0e-10_", file_suffix),
                   "negative, hs=-5e-4" = paste0("h0.5_s-1.0e-3_", file_suffix),
                   "positive, hs=+5e-4" = paste0("h0.5_s+1.0e-3_", file_suffix))
  for (coeff in names(selection)) {
    temp = marginal_helper(paste0(file_prefix, selection[coeff][[1]])) %>%
      mutate(s = coeff)
    df = bind_rows(df, temp)
  }
  return(df)
}