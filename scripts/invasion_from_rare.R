# ---------------------------
# Purpose of script: Analyse invasion-from-rare assays
#
# What this script does:
# 1. wrangles raw dataset
# 2. tests whether relative invader fitness >1
# 3. makes Figure 2
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-03-28
#
# Copyright (c) Daniel Padfield, 2024
#
# ---------------------------
#
# Notes: None
#
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}
# load packages
librarian::shelf(cowplot, tidyverse)

## ---------------------------

# load data in ####
d <- read.csv('data/invasion_from_rare_raw.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names()

color_scheme <- c(AA = '#e41a1c', OD = '#377eb8', PC = '#4daf4a', SR = '#984ea3', VG = '#ff7f00')

labs <- c(expression(italic("Achromobacter")~sp.), expression(italic("Ochrobactrum")~sp.), expression(italic("Pseudomonas")~sp.), expression(italic("Stenotrophomonas")~sp.), expression(italic("Variovorax")~sp.))

#------------#
# wrangle ####
#------------#

# counts at end were plated from frozen, 30 µl at 10^-4. Frozen at half volume
d <- mutate_at(d, vars(ends_with('count')), function(x){x*10^4*2*(1000/30)})
glimpse(d)

# calculate fitness
d <- mutate(d, inv_m = log(inv_count / inv_t0),
            res_m = log(res_count / res_t0),
            rel_fit = inv_m / res_m,
            res_div_text = as.character(res_div))

#-----------------------------------------#
# test for ability to invade from rare ####
#-----------------------------------------#

# the test for coexistence and stability is the ability to invade from rare.
# consequently we are testing each inv_sp vs each res_sp combination against 1.

d <- mutate(d, inv_sp2 = case_when(inv_sp == 'AA' ~ "Achromobacter sp.",
                                   inv_sp == 'OD' ~ "Ochrobactrum sp.",
                                   inv_sp == 'PC' ~ "Pseudomonas sp.",
                                   inv_sp == 'SR' ~ "Stenotrophomonas sp.",
                                   inv_sp == 'VG' ~ "Variovorax sp."))

d_mods <- d %>%
  nest(-c(inv_sp, res_sp, res_div, inv_sp2)) %>%
  mutate(., mod = map(data, ~ t.test(.x$rel_fit, mu = 1)),
         tidied = map(mod, broom::tidy))

d_mods_summary <- d_mods %>%
  select(-c(mod, data)) %>%
  unnest(tidied) %>%
  mutate(padjust = p.adjust(p.value, "fdr"))

not_signif <- filter(d_mods_summary, padjust > 0.05)

#----------------#
# make a plot ####
#----------------#

# try and make a pie chart for each circle - follow this tutorial: https://stackoverflow.com/questions/43984614/rggplot2geom-points-how-to-swap-points-with-pie-charts 

head(d_mods_summary)

d_pie <- mutate(d_mods_summary, x = res_sp) %>%
  separate_rows(., res_sp, sep = '\\.') %>%
  mutate(value = 1 / res_div)

d_merge <- d_pie %>%
  select(inv_sp, x, res_div) %>%
  distinct() %>%
  group_by(inv_sp) %>%
  arrange(., desc(res_div), .by_group = TRUE) %>%
  mutate(x1 = 1:n())

d_pie <- d_pie %>%
  merge(., d_merge, by = c('inv_sp', 'x', 'res_div')) %>%
  group_by(x1, estimate, res_div, inv_sp2, inv_sp)

d_grobs <- d_pie %>% 
  do(subplots = ggplot(., aes(1, value, fill = res_sp)) + 
       geom_col(position = "fill", alpha = 1, colour = NA) + 
       coord_polar(theta = "y") +
       scale_fill_manual(values = color_scheme) +
       theme_void() + 
       guides(fill = 'none')) %>% 
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), 
                                           x = x1-0.4, y = estimate-0.4, 
                                           xmax = x1+0.4, ymax = estimate+0.4))) 

d_grobs

# try and create plot

# data
d <- merge(rename(d, x = res_sp), d_merge, by = c('inv_sp', 'x', 'res_div'))
d_mods_summary <- merge(rename(d_mods_summary, x = res_sp), d_merge, by = c('inv_sp', 'x', 'res_div'))

# plot for Achromobacter
aa_plot <- filter(d_grobs, inv_sp == 'AA') %>%
  {ggplot(data = ., aes(x1, estimate)) +
      geom_rect(aes(xmin = 5.5, xmax = 11.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_linerange(aes(x = x1, ymin = conf.low, ymax = conf.high), filter(d_mods_summary, inv_sp == 'AA')) +
      ylim(c(0, 6)) +
      scale_x_continuous(breaks = c(1, 3.5, 8.5, 13.5), labels = c(4,3,2,1)) +
      .$subgrobs +
      geom_point(aes(x1, rel_fit), filter(d, inv_sp == 'AA'), size = 0.5) +
      coord_flip() +
      theme_bw(base_size = 16) +
      theme(panel.grid.minor.y = element_blank(),
            axis.title.x = element_blank()) +
      labs(title = expression('(a)'~italic('Achromobacter')~sp.),
           x = 'resident community diversity',
           y = 'relative invader growth rate') + 
      geom_hline(aes(yintercept = 1), linetype = 2)}

# plot for Ochrobactrum
od_plot <- filter(d_grobs, inv_sp == 'OD') %>%
  {ggplot(data = ., aes(x1, estimate)) +
      geom_rect(aes(xmin = 5.5, xmax = 11.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_linerange(aes(x = x1, ymin = conf.low, ymax = conf.high), filter(d_mods_summary, inv_sp == 'OD')) +
      ylim(c(0, 6)) +
      scale_x_continuous(breaks = c(1, 3.5, 8.5, 13.5), labels = c(4,3,2,1)) +
      .$subgrobs +
      geom_point(aes(x1, rel_fit), filter(d, inv_sp == 'OD'), size = 0.5) +
      coord_flip() +
      theme_bw(base_size = 16) +
      theme(panel.grid.minor.y = element_blank(),
            axis.title.x = element_blank()) +
      labs(title = expression('(b)'~italic('Ochrobactrum')~sp.),
           x = 'resident community diversity',
           y = 'relative invader growth rate') + 
      geom_hline(aes(yintercept = 1), linetype = 2)}

# plot for Pseudomonas
pc_plot <- filter(d_grobs, inv_sp == 'PC') %>%
  {ggplot(data = ., aes(x1, estimate)) +
      geom_rect(aes(xmin = 5.5, xmax = 11.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_linerange(aes(x = x1, ymin = conf.low, ymax = conf.high), filter(d_mods_summary, inv_sp == 'PC')) +
      ylim(c(0, 6)) +
      scale_x_continuous(breaks = c(1, 3.5, 8.5, 13.5), labels = c(4,3,2,1))+
      .$subgrobs +
      geom_point(aes(x1, rel_fit), filter(d, inv_sp == 'PC'), size = 0.5) +
      coord_flip() +
      theme_bw(base_size = 16) +
      theme(panel.grid.minor.y = element_blank(),
            axis.title.x = element_blank()) +
      labs(title = expression('(c)'~italic('Pseudomonas')~sp.),
           x = 'resident community diversity',
           y = 'relative invader growth rate') + 
      geom_hline(aes(yintercept = 1), linetype = 2)}

# plot for stenotrophomonas
sr_plot <- filter(d_grobs, inv_sp == 'SR') %>%
  {ggplot(data = ., aes(x1, estimate)) +
      geom_rect(aes(xmin = 5.5, xmax = 11.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_linerange(aes(x = x1, ymin = conf.low, ymax = conf.high), filter(d_mods_summary, inv_sp == 'SR')) +
      ylim(c(0, 6)) +
      scale_x_continuous(breaks = c(1, 3.5, 8.5, 13.5), labels = c(4,3,2,1)) +
      .$subgrobs +
      geom_point(aes(x1, rel_fit), filter(d, inv_sp == 'SR'), size = 0.5) +
      coord_flip() +
      theme_bw(base_size = 16) +
      theme(panel.grid.minor.y = element_blank()) +
      labs(title = expression('(d)'~italic('Stenotrophomonas')~sp.),
           x = 'resident community diversity',
           y = 'relative invader growth rate') + 
      geom_hline(aes(yintercept = 1), linetype = 2)}

# plot for variovorax
vg_plot <- filter(d_grobs, inv_sp == 'VG') %>%
  {ggplot(data = ., aes(x1, estimate)) +
      geom_rect(aes(xmin = 5.5, xmax = 11.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf),
                fill = "grey", alpha = 0.03) +
      geom_linerange(aes(x = x1, ymin = conf.low, ymax = conf.high), filter(d_mods_summary, inv_sp == 'VG')) +
      ylim(c(0, 6)) +
      scale_x_continuous(breaks = c(1, 3.5, 8.5, 13.5), labels = c(4,3,2,1)) +
      .$subgrobs +
      geom_point(aes(x1, rel_fit), filter(d, inv_sp == 'VG'), size = 0.5) +
      coord_flip() +
      theme_bw(base_size = 16) +
      theme(panel.grid.minor.y = element_blank()) +
      labs(title = expression('(e)'~italic('Variovorax')~sp.),
           x = 'resident community diversity',
           y = 'relative invader growth rate') + 
      geom_hline(aes(yintercept = 1), linetype = 2)}

# make legend
plot_leg <- ggplot(d, aes(x1, rel_fit, col = inv_sp)) +
  geom_point() +
  scale_color_manual('Resident species', values = color_scheme, labels = labs) +
  theme_bw(base_size = 16) +
  theme(legend.text.align = 0)

plot_leg <- get_legend(plot_leg)

# combine plots
p <- plot_grid(aa_plot, od_plot, pc_plot, sr_plot, vg_plot, plot_leg, ncol = 2, align = 'hv', axis = 'b')

ggsave('plots/figure_2.png', p, width = 10, height = 13)

