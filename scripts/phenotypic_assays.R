#---------------------------#
# fresh analyses of data ####
#---------------------------#

# ideas to look at

# re whether growth rate predicts long-term abundance

#---------------------#
# load in packages ####
#---------------------#

librarian::shelf(tidyverse, patchwork, smatr, cowplot, MuMIn, padpadpadpad/MicrobioUoE, palettetown)

#---------------------#
# load in datasets ####
#---------------------#

# load in supernatant data
d_supernatant <- read.csv('data/processed/supernatant_data.csv')

# load in interaction data
d_interaction <- read.csv('data/processed/interaction_data.csv')

# load in invasion from rare data
d_invasion <- read.csv('data/processed/invasion_from_rare_data.csv')

# load in long term abundance data
d_abundance <- read.csv('data/processed/long_term_abundance_data.csv')

# combine datasets
d <- left_join(d_supernatant, d_interaction) %>%
  left_join(., d_invasion) %>%
  left_join(., d_abundance) %>%
  mutate(diversity = nchar(other_spp) + 1)

# load in carrying capacity dataset
d_carrying_capacity <- read.csv('data/processed/long_term_carrying_capacity.csv')

#-------------------------------------------#
# look at some interesting relationships ####
#-------------------------------------------#

# make colour scheme
color_scheme <- c(a = '#e41a1c', o = '#377eb8', p = '#4daf4a', s = '#984ea3', v = '#ff7f00')
labs <- c(expression(italic("Achromobacter")~sp.), expression(italic("Ochrobactrum")~sp.), expression(italic("Pseudomonas")~sp.), expression(italic("Stenotrophomonas")~sp.), expression(italic("Variovorax")~sp.))
labs2 <- c(expression(italic("Achromobacter")~sp.), expression(italic("Ochrobactrum")~sp.), expression(italic("Pseudomonas")~sp.), expression(italic("Stenotrophomonas")~sp.), expression(italic("Variovorax")~sp.))

# check whether our short-term invasion assays are just measuring long-term abundance ####
d_abundance_check <- mutate(d, prop = inv_count/long_term_abundance)

p1 <- ggplot(d_abundance_check, aes(prop)) +
  geom_vline(aes(xintercept = 1)) +
  geom_histogram(col = 'black', fill = 'grey', bins = 15) +
  theme_bw(base_size = 10) +
  labs(x = 'invasion N/equilibrium N')

p2 <- ggplot(d_abundance_check, aes(long_term_abundance, inv_count)) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  stat_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_point(aes(fill = focal_spp), size = 3, shape = 21) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1e6, 5e9)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1e6, 5e9)) +
  theme_bw() +
  labs(x = expression('equilibrium abundance (cells'~mL^-1~')'),
       y = expression('invader final abundance (cells'~mL^-1~')')) +
  scale_fill_manual('Focal species:', values = color_scheme, labels = labs) +
  theme(legend.text.align = 0)

p3 <- ggdraw(p2) +
  draw_plot(p1, .1, .58, .37, .37)

p3

# so invader final abundance is lower than equilibrium abundance in all instances
d_abundance_check %>% filter(long_term_abundance < inv_count) 
# one instance where invasion abundance > long term abundance
d_abundance_check %>% filter(inv_count/long_term_abundance > 0.8) 
# two more where this is > 0.8
# so our measure of relative invader fitness is not just an estimate of carrying capacity. GOOD!

d_coculture_sum <- filter(d, diversity == 2) %>%
  group_by(., diversity, focal_spp) %>%
  summarise(se = sd(rel_fit)/sqrt(n()),
            ave_relfit = mean(rel_fit), 
            .groups = 'drop')

# make plot of pairwise interactions
p_coculture <- ggplot(filter(d, diversity == 2)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_linerange(aes(focal_spp, ymin = ave_relfit - se, ymax = ave_relfit + se), data = d_coculture_sum) +
  geom_point(aes(focal_spp, ave_relfit, col = focal_spp), size = 5, d_coculture_sum) +
  geom_point(aes(focal_spp, rel_fit, fill = focal_spp), shape = 21, col = 'white', show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(y = 'Relative growth in coculture',
       x = 'Focal species',
       title = '(b)') +
  theme(legend.text.align = 0,
        axis.text.x = element_text(size = 8)) +
  scale_color_manual('Focal species:', values = color_scheme, labels = labs) +
  scale_fill_manual('Focal species:', values = color_scheme, labels = labs) +
  ylim(c(0, 1.6)) +
  scale_x_discrete(labels = labs2, guide = guide_axis(n.dodge = 2))

# plot of pairwise interactions from the supernatant
d_super_sum <- filter(d_supernatant, nchar(other_spp) == 1) %>%
  group_by(focal_spp) %>%
  summarise(., se = sd(inter_effect)/sqrt(n()),
            inter_effect = mean(inter_effect),
            intra_effect = unique(intra_effect)) %>%
  ungroup()

d_super_sum <- mutate(d_super_sum, lower_conf = inter_effect - 1.96*se,
                      reduction_intra = (1 - intra_effect) *100,
                      reduction_lower = (1 - lower_conf) *100,
                      intra_lower = reduction_intra / reduction_lower)


p_supernatant <- ggplot(filter(d_supernatant, nchar(other_spp) == 1), aes(x = focal_spp, y = inter_effect)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_linerange(aes(ymin = inter_effect - 1.96 *se, ymax = inter_effect + 1.96*se), data = d_super_sum) +
  geom_point(aes(col = focal_spp), size = 5, data = d_super_sum) +
  geom_point(aes(fill = focal_spp), shape = 21, col = 'white', show.legend = FALSE) + 
  geom_point(aes(y = intra_effect), shape = 21, fill = 'white', filter(d_supernatant, nchar(other_spp) == 1), size = 5) +
  labs(x = 'Focal species',
       y = 'Relative growth in supernatant',
       title = '(a)') +
  theme_bw(base_size = 12) +
  theme(legend.text.align = 0,
        axis.text.x = element_text(size = 8)) +
  scale_color_manual('Focal species:', values = color_scheme, labels = labs) +
  scale_fill_manual('Focal species:', values = color_scheme, labels = labs) +
  ylim(c(0, 1.6)) +
  scale_x_discrete(labels = labs2, guide = guide_axis(n.dodge = 2))

p_supernatant + p_coculture + plot_layout(guides = 'collect')

ggsave('figures/new/supernatant.png', p_supernatant, width = 7, height = 4)

fit_sma <- smatr::sma(rel_fit ~ inter_effect, filter(d, n_spp == 2), slope.test = 1, robust = TRUE)
preds <- data.frame(expand.grid(inter_effect = seq(min(filter(d, n_spp == 2)$inter_effect, na.rm = T), max(filter(d, n_spp == 2)$inter_effect, na.rm = T), length.out = 40), stringsAsFactors = FALSE)) %>%
  mutate(., preds = coef(fit_sma)[1] + coef(fit_sma)[2]*inter_effect)

# make plot comparing supernatant to coculture
p_coculture_supernatant <- ggplot(filter(d, n_spp == 2), aes(inter_effect, rel_fit)) +
  geom_line(aes(inter_effect, preds), preds) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_vline(aes(xintercept = 1), linetype = 2) +
  geom_point(aes(col = focal_spp, fill = focal_spp), size = 3, show.legend = FALSE) +
  labs(x = 'Relative growth in supernatant',
       y = 'Relative growth in coculture',
       title = '(c)') +
  theme_bw(base_size = 12)  +
  theme(legend.text.align = 0) +
  scale_fill_manual('Focal species:', values = color_scheme, labels = labs) +
  scale_color_manual('Focal species:', values = color_scheme, labels = labs) +
  ylim(c(0.3, 1.6)) +
  xlim(c(0.3, 1.4))

# quantify pairwise interactions as either negative or positive
d <- mutate(d, type_supernatant = ifelse(inter_effect < 1, 'negative', 'positive'),
            type_long_term = ifelse(rel_fit < 1, 'negative', 'positive'),
            short_same_as_long = ifelse(type_supernatant == type_long_term, 'yes', 'no')) 

filter(d, diversity == 2) %>%
  select(focal_spp, other_spp, rel_fit, inter_effect, short_same_as_long) %>%
  View()

# make data frame for bar plot
d_plot <- select(d, n_spp, focal_spp, type_long_term, type_supernatant) %>%
  pivot_longer(cols = starts_with('type'), names_to = 'type', values_to = 'sign', names_prefix = 'type_')

# make bar plot
p2 <- ggplot(filter(d, n_spp == 2), aes(short_same_as_long)) +
  geom_bar(fill = 'grey', col = 'black', linewidth = 0.15) +
  theme_classic(base_size = 8) +
  labs(y = 'number',
       x = 'Same interaction\nsign')
# inset plot
p3 <- p_coculture_supernatant + inset_element(p2, left = 0.05, bottom = 0.65, right = 0.35, top = 0.95)

p_supernatant + p_coculture + p3 + plot_layout(guides = 'collect')

ggsave('figures/new/interactions.png', last_plot(), height =4, width = 13)

#--------------------------------------#
# look at higher order interactions ####
#--------------------------------------#

# filter for only when there are more than two species
# classify indirect interaction as either synergistic or buffering
d_many_species <- filter(d, diversity > 2) %>%
  mutate(., hoi_class = case_when(expected_fitness_tot <= 1 & tot_hoi <= 0 ~ 'synergistic',
                                  expected_fitness_tot <= 1 & tot_hoi > 0 ~ 'buffering',
                                  expected_fitness_tot > 1 & tot_hoi <= 0 ~ 'buffering',
                                  expected_fitness_tot > 1 & tot_hoi > 0 ~ 'synergistic'))

# scale indirect interactions that are currently just the sum of observed - expected
d_many_species <- 
  mutate(d_many_species, 
         tot_hoi2 = case_when(expected_fitness_tot <= 1 & tot_hoi <= 0 ~ tot_hoi*-1,
                              expected_fitness_tot <= 1 & tot_hoi > 0 ~ tot_hoi * -1,
                              expected_fitness_tot > 1 & tot_hoi <= 0 ~ tot_hoi,
                              expected_fitness_tot > 1 & tot_hoi > 0 ~ tot_hoi))

select(d_many_species, expected_fitness_tot, rel_fit, tot_hoi, tot_hoi2, hoi_class) %>%
  View()

# calculate summary stats of total indirect effects
d_many_sum <- group_by(d_many_species, diversity, focal_spp) %>%
  summarise(se = sd(tot_hoi2, na.rm = TRUE)/sqrt(n()),
    tot_hoi2 = mean(tot_hoi2),
    .groups = 'drop')

# make plot of higher order interactions
p_hoi <- ggplot(d_many_species, aes(diversity, tot_hoi2, fill = focal_spp, col = focal_spp)) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  geom_linerange(aes(ymin = tot_hoi2 - 1.96*se, ymax = tot_hoi2 + 1.96*se), position = position_dodge(width = 0.5), col = 'black', d_many_sum) +
  geom_point(position = position_dodge(width = 0.5), data = d_many_sum, size = 3) +
  geom_point(position = position_dodge(width = 0.5), shape = 21, col = 'white', show.legend = FALSE, data = filter(d_many_species, diversity <= 4)) +
  labs(x = 'Diversity',
       y = 'Total indirect interaction',
       title = '(a)') +
  theme_bw(base_size = 12) +
  theme(legend.text.align = 0,
        legend.position = c(0.79, 0.8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.75, 'lines')) +
  scale_fill_manual('Focal species:', values = color_scheme, labels = labs) +
  scale_color_manual('Focal species:', values = color_scheme, labels = labs) +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  guides(col = guide_legend(override.aes = list(size=2))) +
  geom_text(aes(x = 2.5, y = 0, label = 'synergistic'), col = 'black', angle = 90, hjust = -0.2, stat = "unique") +
  geom_text(aes(x = 2.5, y = 0, label = 'buffering'), col = 'black', angle = 90, hjust = 1.2, stat = "unique")

p_hoi

# plot expected fitness across hoi class, with positive expected or negative expectation separate
d_many_species <- mutate(d_many_species, expected_interaction = ifelse(expected_fitness_tot <= 1, 'negative', 'positive'))

# create summary stats for expected fitness
d_many_sum2 <- group_by(d_many_species, hoi_class, expected_interaction) %>%
  summarise(se = sd(expected_fitness_tot, na.rm = TRUE)/sqrt(n()),
    expected_fitness_tot = mean(expected_fitness_tot),
    .groups = 'drop')

p_hoi2 <- ggplot(d_many_species, aes(expected_interaction, expected_fitness_tot, col = hoi_class, fill = hoi_class)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_linerange(aes(ymin = expected_fitness_tot - 1.96*se, ymax = expected_fitness_tot + 1.96*se), position = position_dodge(width = 0.3), col = 'black', d_many_sum2) +
  geom_point(position = position_dodge(width = 0.3), data = d_many_sum2, size = 5) +
  geom_point(position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.05), shape = 21, col = 'white', show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(x = 'Expected interaction sign',
       y = 'Expected fitness from\npairwise coculture',
       title = '(b)') +
  ylim(c(0, 1.5)) +
  scale_color_poke('Indirect interaction', pokemon = 'moltres') +
  scale_fill_poke('Indirect interaction', pokemon = 'moltres') +
  theme(legend.position = c(0.8, 0.15),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.75, 'lines')) +
  guides(col = guide_legend(override.aes = list(size=2)))

p_hoi2

p_hoi + p_hoi2

ggsave('figures/new/interactions_ii.png', last_plot(), height =4.5, width = 10)

mod1 <- lm(expected_fitness_tot ~ hoi_class, filter(d_many_species, expected_interaction == 'negative'))
mod2 <- lm(expected_fitness_tot ~ 1 ,filter(d_many_species, expected_interaction == 'negative'))
anova(mod1, mod2)
mod2 <- lm(expected_fitness_tot ~ 1, filter(d_many_species,  expected_fitness_tot <=1))

# have a look at number of buffering and synergisms
group_by(d_many_species, hoi_class, diversity) %>%
  summarise(n = n(),
            .groups = 'drop') %>%
  group_by(diversity) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  arrange(diversity)

lm(abs(tot_hoi) ~ hoi_class, d_many_species) %>%
  summary()


# plot expected fitness across hoi class
p_hoi4 <- ggplot(filter(d_many_species, expected_fitness_tot > 1), aes(hoi_class, expected_fitness_tot, col = hoi_class, fill = hoi_class)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(position = position_jitter(width = 0.2), shape = 21, fill = 'white') +
  theme_bw() +
  labs(x = 'Interaction type',
       y = 'Expected fitness from pairwise coculture') +
  scale_color_poke(pokemon = 'blastoise')

p_hoi4

# proportion of synergism / buffering

ggplot(d_many_species, aes(expected_fitness_tot, tot_hoi, col = hoi_class)) +
  geom_point() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  labs(x = 'expected fitness',
       y = 'Total higher order interaction') +
  theme_bw(base_size = 12) +
  theme(legend.text.align = 0) +
  facet_wrap(~diversity)

# plot total expected fitness change and then the number of buffering and synergistic interactions
d_many_species <- mutate(d_many_species, absolute_prediction = abs(1 - expected_fitness_tot))

ggplot(filter(d_many_species, expected_fitness_tot <=1), aes(hoi_class, expected_fitness_tot)) +
  geom_boxplot() +
  labs(x = 'Higher order interaction class',
       y = 'Expected fitness') +
  theme_bw(base_size = 12) +
  theme(legend.text.align = 0) +
  facet_wrap(~diversity)

mod1 <- lm(expected_fitness_tot ~ hoi_class, filter(d_many_species,  expected_fitness_tot <=1))
mod2 <- lm(expected_fitness_tot ~ 1, filter(d_many_species,  expected_fitness_tot <=1))
mod3 <- lm(expected_fitness_tot ~ hoi_class, filter(d_many_species, expected_fitness_tot <=1 & diversity == 3))
mod4 <- lm(expected_fitness_tot ~ 1, filter(d_many_species, expected_fitness_tot <=1 & diversity == 3))
mod5 <- lm(expected_fitness_tot ~ hoi_class, filter(d_many_species, expected_fitness_tot <=1 & diversity == 4))
mod6 <- lm(expected_fitness_tot ~ 1, filter(d_many_species, expected_fitness_tot <=1 & diversity == 4))

anova(mod1, mod2)
anova(mod3, mod4)
anova(mod5, mod6)