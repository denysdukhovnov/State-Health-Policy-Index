## This is the main code for the manuscript titled:
## State Health Policy Index: A Summary Indicator for Subnational Mortality Analysis
## 
## 
## -------------------

# Install and load the packages required for processing the script 
# install.packages("tidyverse")
# install.packages("GGally")
# install.packages("ggpubr")
# install.packages("scales")
# install.packages("rstatix")
library(tidyverse)
library(GGally)
library(ggpubr)
library(scales)
library(rstatix)
options(scipen = 99999, max.print = 500)

# Source auxiliary functions for computing durability kernel weighting functions
source("weights_kernel_selection.R")
source("apply_positivity_weights.R")

# Lookup table for state codes, abbreviations, and names
state.fips.lookup <- tibble(
  state.codes = c("AL","AK","AZ","AR","CA","CO","CT","DE","DC","FL",
                  "GA","HI","ID","IL","IN","IA","KS","KY","LA","ME",
                  "MD","MA","MI","MN","MS","MO","MT","NE","NV","NH",
                  "NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI",
                  "SC","SD","TN","TX","UT","VT","VA","WA","WV","WI",
                  "WY"),
  
  state.names = c("Alabama","Alaska","Arizona","Arkansas","California",
                  "Colorado","Connecticut","Delaware","District of Columbia",
                  "Florida","Georgia","Hawaii","Idaho","Illinois","Indiana",
                  "Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland",
                  "Massachusetts","Michigan","Minnesota","Mississippi",
                  "Missouri","Montana","Nebraska","Nevada","New Hampshire",
                  "New Jersey","New Mexico","New York","North Carolina",
                  "North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania",
                  "Rhode Island","South Carolina","South Dakota","Tennessee",
                  "Texas","Utah","Vermont","Virginia","Washington",
                  "West Virginia","Wisconsin","Wyoming"),
  
  state.fips = c("01","02","04","05","06","08","09","10","11","12","13","15",
                 "16","17","18","19","20","21","22","23","24","25","26","27",
                 "28","29","30","31","32","33","34","35","36","37","38","39",
                 "40","41","42","44","45","46","47","48","49","50","51","53",
                 "54","55","56")
)


#######################

# Read in the main data set comprising compiled and processed policies from all policy databases
policy.inputs <- readRDS("PolicyData_input.rds")

# Create kernel to distribute the health positivity weights over a window period and functional form (gamma, poisson, negative binomial, cauchy, or sigmoid)
kern.res <- kernel.matix.fun(window.size = length(unique(policy.inputs$year)), n.years = length(unique(policy.inputs$year)), kernel.type = "cauchy")
kernel.mat <- kern.res[[length(kern.res)]]

# Apply positivity weights and calculate the I score (may take a bit of time and may issue warnings due to NA values, which should be OK to ignore)
policy.inputs.wgt <- policy.inputs %>% 
  
  # Apply weights to the health policy positivity scores based on a kernel
  apply.positivity.weights.fun(data = ., policy.cols = 5:(ncol(.) - 1), kernel.matrix = kernel.mat) %>% 
  
  # Count the number of policies in a given state-year-category-class
  mutate(n.policies = rowSums(!is.na(.[5:ncol(.)]))) %>% 

  # Calculate the health/mortality positivity score of policies by class and category (note the number of columns to be summed over (#:ncol(.) - #))
  mutate(p = rowSums(.[5:(ncol(.) - 1)], na.rm = T)) %>% 
  
  # Compute the number of policies in a policy class in each year
  group_by(state_name, year, policy_class) %>% 
  mutate(n.policies.class = sum(n.policies, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Compute sum of policies in a year and state
  group_by(state_name, year) %>% 
  mutate(n.policies.class.sum = sum(n.policies.class, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Compute sum of policies in a year and state
  group_by(state_name, year, policy_category) %>% 
  mutate(n.policies.class.category.sum = sum(n.policies.class, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Maximum policy positivity score in any state in a given year in a category
  group_by(state_name, year, policy_category) %>% 
  mutate(sum.p.g = sum(p, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(year, policy_category) %>% 
  mutate(max.p.g = max(sum.p.g, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Ratio of policy health positivity each thematic group
  mutate(Rg = sum.p.g / max.p.g) %>% 
  
  # Sum of ratio of policy health positivity for all thematic groups
  group_by(state_name, year) %>% 
  mutate(Rg.x.n.policies.class.sum = sum(Rg, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Sum of ratio of policy health positivity for each thematic group
  group_by(state_name, year, policy_category) %>% 
  mutate(Rg.category.x.n.policies.class.sum = sum(Rg, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Weighted policy score for all thematic groups and policy classes
  mutate(I.score = Rg.x.n.policies.class.sum / n.policies.class.sum) %>% 
  
  # Compute thematic group-specific I-score
  mutate(Ig.score = Rg.category.x.n.policies.class.sum / n.policies.class.category.sum)
  

# Additionally calculate the kernel-unweighted I score 
policy.inputs.unwgt <- policy.inputs %>% select(-state) %>% 
  
  # Count the number of policies in a given state-year-category-class
  mutate(n.policies = rowSums(!is.na(.[5:ncol(.)]))) %>% 

  # Calculate the health/mortality positivity score of policies by class and category (note the number of columns to be summed over (#:ncol(.) - #))
  mutate(p = rowSums(.[5:(ncol(.) - 1)], na.rm = T)) %>% 
  
  # Compute the number of policies in a policy class in each year
  group_by(state_name, year, policy_class) %>% 
  mutate(n.policies.class = sum(n.policies, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Compute sum of policies in a year and state
  group_by(state_name, year) %>% 
  mutate(n.policies.class.sum = sum(n.policies.class, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Compute sum of policies in a year and state
  group_by(state_name, year, policy_category) %>% 
  mutate(n.policies.class.category.sum = sum(n.policies.class, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Maximum policy positivity score in any state in a given year in a category
  group_by(state_name, year, policy_category) %>% 
  mutate(sum.p.g = sum(p, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(year, policy_category) %>% 
  mutate(max.p.g = max(sum.p.g, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Ratio of policy positivity each thematic group
  mutate(Rg = sum.p.g / max.p.g) %>% 
  
  # Sum of ratio of policy positivity for all thematic groups
  group_by(state_name, year) %>% 
  mutate(Rg.x.n.policies.class.sum = sum(Rg, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Sum of ratio of policy positivity for each thematic group
  group_by(state_name, year, policy_category) %>% 
  mutate(Rg.category.x.n.policies.class.sum = sum(Rg, na.rm = T)) %>% 
  ungroup() %>% 
  
  # Weighted policy score for all thematic groups and policy classes
  mutate(I.score = Rg.x.n.policies.class.sum / n.policies.class.sum) %>% 
  
  # Compute thematic group-specific I-score
  mutate(Ig.score = Rg.category.x.n.policies.class.sum / n.policies.class.category.sum)


# Calculate the weighted I-score's rank and its normalized version
I.score.wgt <- policy.inputs.wgt %>% 
  group_by(state_name, year) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(state_name, year, I.score) %>% 
  {. ->> I.score.full.wgt} %>% 
  
  # Rank states by the policy score
  arrange(year, desc(I.score)) %>% 
  group_by(year) %>% 
  mutate(I.rank = row_number(),
         I.score.Z = scale(I.score)) %>% 
  ungroup()

# Repeat the same as above for the unweighted I score
I.score.unwgt <- policy.inputs.unwgt %>% 
  group_by(state_name, year) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(state_name, year, I.score) %>% 
  {. ->> I.score.full.unwgt} %>% 

  # Rank states by the policy score
  arrange(year, desc(I.score)) %>% 
  group_by(year) %>% 
  mutate(I.rank = row_number(),
         I.score.Z = scale(I.score)) %>% 
  ungroup()

# Prepare a plot data frame
I.score.raw.plt.df <- I.score.wgt %>% mutate(Type = "Weighted") %>% 
  bind_rows(I.score.unwgt %>% mutate(Type = "Unweighted"))

# Plot the raw form of weighted and unweighted SHP score
I.score.plt.raw <- ggplot(I.score.raw.plt.df, aes(x = year, y = I.score, color = state_name)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~Type) +
  labs(x = "Year", y = "SHP score (raw)") +
  scale_x_continuous(breaks = c(seq(1970, 2010, 10), 2019), labels = c(seq(1970, 2010, 10), 2019), 
                     limits = c(1970, 2019), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  guides(color = "none") +
  geom_vline(xintercept = 1990, color = "black", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"),
        strip.text = element_text(face = "bold", size = 14),
        panel.spacing.x = unit(1, "cm"),
        axis.title.x = element_text(size = 13, vjust = -0.5),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.text = element_text(size = 12),
        panel.grid.minor.x = element_blank())

I.score.plt.raw

# Plot the normalized weighted and unweighted I score
I.score.plt <- ggplot(I.score.raw.plt.df, aes(x = year, y = I.score.Z, color = state_name)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~Type) +
  labs(x = "Year", y = "SHP score (normalized)") +
  scale_x_continuous(breaks = c(seq(1970, 2010, 10), 2019), labels = c(seq(1970, 2010, 10), 2019), 
                     limits = c(1970, 2019), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  guides(color = "none") +
  geom_vline(xintercept = 1990, color = "black", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"),
        strip.text = element_text(face = "bold", size = 14),
        panel.spacing.x = unit(1, "cm"),
        axis.title.x = element_text(size = 13, vjust = -0.5),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.text = element_text(size = 12),
        panel.grid.minor.x = element_blank())

I.score.plt

# Ranked version
rank.in.2019.wgt <- I.score.wgt %>% 
  filter(year == 2019) %>% 
  select(state_name) %>% 
  c() %>% unlist()

rank.in.2019.unwgt <- I.score.unwgt %>% 
  filter(year == 2019) %>% 
  select(state_name) %>% 
  c() %>% unlist()

# Arrange the labels for the plot to be in the descending order of rank based on the final year
I.score.wgt.plt.df <- I.score.wgt %>% 
  mutate(state_name = factor(state_name, levels = rev(rank.in.2019.wgt)))

I.score.unwgt.plt.df <- I.score.unwgt %>% 
  mutate(state_name = factor(state_name, levels = rev(rank.in.2019.unwgt)))

# Plot the state rank
I.rank.plt <- ggplot(I.score.wgt.plt.df %>%
                       mutate(Type = "Weighted") %>% 
                       bind_rows(I.score.unwgt.plt.df %>% mutate(Type = "Unweighted")) %>% 
                       filter(year %in% c(1990, 2000, 2010, 2019)), 
                     aes(x = I.rank, y = state_name, color = as.factor(year), shape = as.factor(year))) + 
  geom_point(size = 2) +
  facet_wrap(~Type, scales = "free") +
  labs(x = "Health Policy Rank", y = "") +
  guides(color = guide_legend(title = "Year"), shape = guide_legend(title = "Year")) +
  scale_x_continuous(limits = c(1, 50), breaks = c(1, seq(5, 50, 5)), expand = c(0.01,0.01)) +
  scale_shape_manual(values = c(1, 3, 4, 5)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"),
        strip.text = element_text(face = "bold", size = 14),
        panel.spacing.x = unit(0.5, "cm"),
        axis.title.x = element_text(size = 13, vjust = -0.5),
        axis.text = element_text(size = 11),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.key.height = unit(0.5, "in"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"))

I.rank.plt


# Compute the weighted 10-year lagged I score and its rank, then save the dataset
saveRDS(I.score.wgt %>%
          arrange(state_name, year) %>% 
          group_by(state_name) %>% 
          mutate(I.score.10yr.lag = dplyr::lag(I.score),
                 I.rank.10yr.lag = dplyr::lag(I.rank)) %>%
          ungroup() %>% 
          group_by(year) %>% 
          mutate(I.score.10yr.lag.Z = scale(I.score.10yr.lag)) %>% 
          ungroup() %>% 
          filter(year >= 1990), 
        file = "I_score_wgt.rds")

# Repeat for the unweighted score
saveRDS(I.score.unwgt %>%
          arrange(state_name, year) %>% 
          group_by(state_name) %>% 
          mutate(I.score.10yr.lag = dplyr::lag(I.score),
                 I.rank.10yr.lag = dplyr::lag(I.rank)) %>%
          ungroup() %>% 
          group_by(year) %>% 
          mutate(I.score.10yr.lag.Z = scale(I.score.10yr.lag)) %>% 
          ungroup() %>% 
          filter(year >= 1990), 
        file = "I_score_unwgt.rds")

# Save the plots that were plotted earlier
ggsave(plot = I.score.plt.raw, filename = "health_policy_score_RAW_wgt_and_unwgt_plot_by_state_1970_2019.jpg", 
       width = 6.5, height = 3, units = "in", scale = 2, dpi = 300)

ggsave(plot = I.score.plt, filename = "health_policy_score_wgt_and_unwgt_plot_by_state_1970_2019.jpg", 
       width = 6.5, height = 3, units = "in", scale = 2, dpi = 300)

ggsave(plot = I.rank.plt, filename = "health_policy_rank_wgt_and_unwgt_plot_by_state_1990_2019.jpg", 
       width = 5, height = 6, units = "in", scale = 1.75, dpi = 300)


############### Repeat to plot thematic category-specific graphs ###############
Ig.score.wgt <- policy.inputs.wgt %>% 
  group_by(state_name, year, policy_category) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(state_name, year, policy_category, I.score, Ig.score) %>% 
  {. ->> Ig.score.wgt.full} %>% 

  # Rank states by the policy score
  arrange(year, desc(Ig.score)) %>% 
  group_by(year, policy_category) %>% 
  mutate(Ig.rank = row_number(),
         Ig.score.Z = scale(Ig.score)) %>% 
  ungroup()

Ig.score.unwgt <- policy.inputs.unwgt %>% 
  group_by(state_name, year, policy_category) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(state_name, year, policy_category, I.score, Ig.score) %>% 
  {. ->> Ig.score.unwgt.full} %>% 
  
  # Rank states by the policy score
  arrange(year, desc(Ig.score)) %>% 
  group_by(year, policy_category) %>% 
  mutate(Ig.rank = row_number(),
         Ig.score.Z = scale(Ig.score)) %>% 
  ungroup()

rank.g.wgt.in.2019 <- Ig.score.wgt %>% 
  filter(year == 2019) %>% 
  select(state_name, policy_category, Ig.rank) %>% 
  arrange(policy_category, Ig.rank)

rank.g.unwgt.in.2019 <- Ig.score.unwgt %>% 
  filter(year == 2019) %>% 
  select(state_name, policy_category, Ig.rank) %>% 
  arrange(policy_category, Ig.rank)

# Arrange the labels for the plot to be in the descending order of rank based on the final year
Ig.score.plt.df <- Ig.score.wgt %>% 
  mutate(Type = "Weighted") %>% 
  bind_rows(Ig.score.unwgt %>% 
              mutate(Type = "Unweighted"))


# Plot the I score for each thematic policy group 
Ig.score.plt.unwgt <- ggplot(Ig.score.plt.df %>% filter(Type == "Unweighted"), aes(x = year, y = Ig.score.Z, color = state_name)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, colour = "black") +
  labs(x = "Year", y = "SHP score (normalized)") +
  scale_x_continuous(breaks = c(seq(1970, 2010, 10), 2019), labels = c(seq(1970, 2010, 10), 2019), 
                     limits = c(1970, 2019), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-10, 10, 2), labels = seq(-10, 10, 2), 
                     expand = c(0.001, 0.001)) +
  facet_wrap(~ policy_category, ncol = 4, nrow = 2) +
  guides(color = "none") +
  geom_vline(xintercept = 1990, color = "black", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.ticks = element_line(color = "grey80"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank(), 
        panel.spacing.x = unit(1, "cm"),
        panel.border = element_rect(color = "grey80", fill = "transparent"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.placement = "inside")

Ig.score.plt.wgt <- ggplot(Ig.score.plt.df %>% filter(Type == "Weighted"), aes(x = year, y = Ig.score.Z, color = state_name)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, colour = "black") +
  labs(x = "Year", y = "SHP score (normalized)") +
  scale_x_continuous(breaks = c(seq(1970, 2010, 10), 2019), labels = c(seq(1970, 2010, 10), 2019), 
                     limits = c(1970, 2019), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-10, 10, 2), labels = seq(-10, 10, 2), 
                     expand = c(0.001, 0.001)) +
  facet_wrap(~ policy_category, ncol = 4, nrow = 2) +
  guides(color = "none") +
  geom_vline(xintercept = 1990, color = "black", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.5),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.ticks = element_line(color = "grey80"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank(), 
        panel.spacing.x = unit(1, "cm"),
        panel.border = element_rect(color = "grey80", fill = "transparent"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.placement = "inside")

Ig.score.plt.unwgt
Ig.score.plt.wgt

Ig.score.plt <- ggarrange(Ig.score.plt.unwgt, Ig.score.plt.wgt, nrow = 2,
          label.x = c(-0.025, -0.025), label.y = c(1,1), 
          labels = c("Unweighted", "Weighted"))

# Save the thematic kernel-weighted and unweighted I scores
ggsave(plot = Ig.score.plt, filename = "health_policy_score_wgt_and_unwgt_plot_by_thematic_group_1970_2019.jpg",
       dpi = 300, scale = 1.65, width = 6.5, height = 5, units = "in")

# Plot the equivalent raw/untransformed I scores for each thematic policy group 
Ig.score.plt.unwgt.raw <- ggplot(Ig.score.plt.df %>% 
                                   filter(Type == "Unweighted"), 
                                 aes(x = year, y = Ig.score, color = state_name)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, colour = "black") +
  labs(x = "Year", y = "SHP score (raw)") +
  scale_x_continuous(breaks = c(seq(1970, 2010, 10), 2019), labels = c(seq(1970, 2010, 10), 2019), 
                     limits = c(1970, 2019), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-1, 1, .2), labels = seq(-1, 1, .2), 
                     expand = c(0.001, 0.001)) +
  facet_wrap(~ policy_category, ncol = 4, nrow = 2) +
  guides(color = "none") +
  geom_vline(xintercept = 1990, color = "black", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.ticks = element_line(color = "grey80"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank(), 
        panel.spacing.x = unit(1, "cm"),
        panel.border = element_rect(color = "grey80", fill = "transparent"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.placement = "inside")

Ig.score.plt.wgt.raw <- ggplot(Ig.score.plt.df %>% filter(Type == "Weighted"), aes(x = year, y = Ig.score, color = state_name)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, colour = "black") +
  labs(x = "Year", y = "SHP score (raw)") +
  scale_x_continuous(breaks = c(seq(1970, 2010, 10), 2019), labels = c(seq(1970, 2010, 10), 2019), 
                     limits = c(1970, 2019), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-1, 1, .2), labels = seq(-1, 1, .2), 
                     expand = c(0.001, 0.001)) +
  facet_wrap(~ policy_category, ncol = 4, nrow = 2) +
  guides(color = "none") +
  geom_vline(xintercept = 1990, color = "black", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.5),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.ticks = element_line(color = "grey80"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank(), 
        panel.spacing.x = unit(1, "cm"),
        panel.border = element_rect(color = "grey80", fill = "transparent"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.placement = "inside")

Ig.score.plt.unwgt.raw
Ig.score.plt.wgt.raw

# Arrange the plots for the non-transformed/non-normalized thematic I scores and save the plot
Ig.score.plt.raw <- ggarrange(Ig.score.plt.unwgt.raw, Ig.score.plt.wgt.raw, nrow = 2,
                          label.x = c(-0.025, -0.025), label.y = c(1,1), 
                          labels = c("Unweighted", "Weighted"))

ggsave(plot = Ig.score.plt.raw, filename = "health_policy_score_RAW_wgt_and_unwgt_plot_by_thematic_group_1970_2019.jpg",
       dpi = 300, scale = 1.65, width = 6.5, height = 5, units = "in")



# Compute the time-lagged thematic group policy scores and ranks and save the dataset
saveRDS(Ig.score.wgt %>%
          group_by(state_name, policy_category) %>% 
          mutate(Ig.score.10yr.lag = dplyr::lag(Ig.score),
                 Ig.rank.10yr.lag = dplyr::lag(Ig.rank)) %>%
          ungroup() %>%
          group_by(year, policy_category) %>% 
          mutate(Ig.score.10yr.lag.Z = dplyr::lag(Ig.score.10yr.lag)) %>% 
          ungroup() %>% 
          select(-I.score) %>% 
          filter(year >= 1990), 
        file = "Ig_score_wgt.rds")

saveRDS(Ig.score.unwgt %>%
          group_by(state_name, policy_category) %>% 
          mutate(Ig.score.10yr.lag = dplyr::lag(Ig.score),
                 Ig.rank.10yr.lag = dplyr::lag(Ig.rank)) %>%
          ungroup() %>%
          group_by(year, policy_category) %>% 
          mutate(Ig.score.10yr.lag.Z = dplyr::lag(Ig.score.10yr.lag)) %>% 
          ungroup() %>% 
          select(-I.score) %>% 
          filter(year >= 1990), 
        file = "Ig_score_unwgt.rds")

#############################



########### Correlation plots of policy with mortality

mortality.state <- read.table("state_mortality.csv", sep = ",", 
                              stringsAsFactors = FALSE, header = TRUE, 
                              colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
mortality.county <- read.table("county_mortality.csv", sep = ",", 
                               stringsAsFactors = FALSE, header = TRUE,
                               colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))

I.score.wgt.rds <- readRDS("I_score_w_wgts.rds") %>% 
  group_by(year) %>% 
  mutate(I.score.quintile = factor(cut(I.score, labels = 1:5, breaks = quantile(I.score, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T)),
         I.score.10yr.lag.quintile = factor(cut(I.score.10yr.lag, labels = 1:5, breaks = quantile(I.score.10yr.lag, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T))) %>% 
  ungroup()

I.score.unwgt.rds <- readRDS("I_score_unwgt.rds") %>% 
  group_by(year) %>% 
  mutate(I.score.quintile = factor(cut(I.score, labels = 1:5, breaks = quantile(I.score, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T)),
         I.score.10yr.lag.quintile = factor(cut(I.score.10yr.lag, labels = 1:5, breaks = quantile(I.score.10yr.lag, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T))) %>% 
  ungroup()

# Thematic policy category scores
Ig.score.wgt.rds <- readRDS("Ig_score_wgt.rds") %>% 
  group_by(year) %>% 
  mutate(Ig.score.quintile = factor(cut(Ig.score, labels = 1:5, breaks = quantile(Ig.score, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T)),
         Ig.score.10yr.lag.quintile = factor(cut(Ig.score.10yr.lag, labels = 1:5, breaks = quantile(Ig.score.10yr.lag, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T))) %>% 
  ungroup()

Ig.score.unwgt.rds <- readRDS("Ig_score_unwgt.rds") %>% 
  group_by(year) %>% 
  mutate(Ig.score.quintile = factor(cut(Ig.score, labels = 1:5, breaks = quantile(Ig.score, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T)),
         Ig.score.10yr.lag.quintile = factor(cut(Ig.score.10yr.lag, labels = 1:5, breaks = quantile(Ig.score.10yr.lag, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)), right = T, include.lowest = T))) %>% 
  ungroup()

# Combine the mortality and the I score data for states, while computing log(ASDR)
policy.wgt.mort.state <- mortality.state %>% 
  left_join(state.fips.lookup, by = c("STATE" = "state.fips")) %>% 
  mutate(across(c(mx.standard.b, mx.standard.f, mx.standard.m), ~log(.x*1), .names = "log.{.col}")) %>% 
  left_join(I.score.wgt.rds, by = c("Year" = "year", "state.names" = "state_name")) %>% 
  mutate(Year = as.factor(Year))

# Combine the mortality and the I score data for county groups, while computing log(ASDR)
policy.wgt.mort.county <- mortality.county %>%
  mutate(state.fips = substr(FIPS.grp, 1, 2)) %>% 
  left_join(state.fips.lookup, by = "state.fips") %>% 
  mutate(across(c(mx.standard.b, mx.standard.f, mx.standard.m), ~log(.x*1), .names = "log.{.col}")) %>%
  left_join(I.score.wgt.rds, by = c("Year" = "year", "state.names" = "state_name")) %>%
  mutate(STATE = substr(FIPS.grp, 1, 2)) %>%
  mutate(Year = as.factor(Year)) 
  
# Check the state sample column name indices to include in the correlation plots below
for (i in 1:length(names(policy.wgt.mort.state))) print(c(i, trimws(names(policy.wgt.mort.state)[i])))

# Run and plot Pearson's correlations on state combination of ASDR and the kernel-weighted I score forms by year
corr.plt.state.mortality <- ggpairs(policy.wgt.mort.state,
                          columns = c(3:5, 11:12, 14:15),
                          aes(color = as.factor(Year), alpha = 0.5),
                          upper = list("cor"),
                          lower = list(continuous = "smooth")) +
  theme_minimal() +
  theme(plot.background = element_rect(color = "transparent", fill = "white"),
        panel.border = element_rect(color = "grey70", fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "grey70"),
        strip.text = element_text(size = 13))

corr.plt.state.mortality

# Save the state correlation plot
ggsave(plot = corr.plt.state.mortality, file = "correlation_plot_state_mortality.png", 
       height = 6, width = 6, units = "in", dpi = 150, scale = 1.7)


# Check the county sample column name indices to include in the correlation plots below
for (i in 1:length(names(policy.wgt.mort.county))) print(c(i, trimws(names(policy.wgt.mort.county)[i])))

# Run and plot Pearson's correlations on county group combination of ASDR and the kernel-weighted I score forms by year
corr.plt.county.mortality <- ggpairs(policy.wgt.mort.county,
                                    columns = c(3:5, 12:13, 15:16),
                                    aes(color = as.factor(Year), alpha = 0.5),
                                    upper = list("cor"),
                                    lower = list(continuous = "smooth")) +
  theme_minimal() +
  theme(plot.background = element_rect(color = "transparent", fill = "white"),
        panel.border = element_rect(color = "grey70", fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "grey70"),
        strip.text = element_text(size = 13))

corr.plt.county.mortality

# Save the county group correlation plot
ggsave(plot = corr.plt.county.mortality, file = "correlation_plot_county_mortality.png", 
       height = 6, width = 6, units = "in", dpi = 150, scale = 1.7)


### Repeating the above plots but for the weighted forms of thematic policy components (raw ones are not shown, but can be computed)

# Merging and processing of thematic components for state sample
policy.mort.wgt.state.g <- Ig.score.wgt.rds %>% 
  filter(year %in% c(1990, 2000, 2010, 2019)) %>%
  mutate(year = as.numeric(year)) %>% 
  left_join(policy.wgt.mort.state %>%
              mutate(Year = as.numeric(as.character(Year))),
            by = c("state_name" = "state.names", "year" = "Year"), 
            multiple = "all") %>%
  pivot_wider(id_cols = c("state_name", "year"), names_from = "policy_category", values_from = "Ig.score") %>%
  rename_with(.fn = ~ paste0("cat.", .x), .cols = 3:10) %>%
  mutate(Year = as.factor(year)) %>%
  left_join(policy.wgt.mort.state %>%   
              mutate(Year = as.numeric(as.character(Year))) %>% 
              select(state.names, Year, mx.standard.b, mx.standard.f, mx.standard.m),
            by = c("state_name" = "state.names", "year" = "Year")) %>% 
  select(state_name, Year, mx.standard.b, mx.standard.f, mx.standard.m, cat.1, cat.2, cat.3, cat.4, cat.5, cat.6, cat.7, cat.8)

# Check the column name indices to include in the correlation plots below
names(policy.mort.wgt.state.g)

# Run and plot Pearson's correlations on state combination of ASDR and the kernel-weighted thematic Ig score categories by year
corr.plt.state.mx.Ig.wgt <- ggpairs(policy.mort.wgt.state.g,
                                        columns = c(3:13),
                                        aes(color = as.factor(Year), alpha = 0.5),
                                        upper = list("cor"),
                                        lower = list(continuous = "smooth")) +
  theme_minimal() +
  theme(plot.background = element_rect(color = "transparent", fill = "white"),
        panel.border = element_rect(color = "grey70", fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "grey70"),
        strip.text = element_text(size = 13))

corr.plt.state.mx.Ig.wgt

# Save the thematic group state-level correlation plot
ggsave(plot = corr.plt.state.mx.Ig.wgt, file = "correlation_plot_state_mortality_Ig_score.png",
       height = 6, width = 6, units = "in", dpi = 150, scale = 2.3)


# Merging and processing of thematic components for county group sample
policy.mort.wgt.county.g <- Ig.score.wgt.rds %>% 
  filter(year %in% c(1990, 2000, 2010, 2019)) %>%
  mutate(year = as.numeric(year)) %>% 
  right_join(policy.wgt.mort.county %>%
              mutate(Year = as.numeric(as.character(Year))),
            by = c("state_name" = "state.names", "year" = "Year"), 
            multiple = "all", relationship = "many-to-many") %>% 
  group_by(state_name, year, policy_category) %>% 
  slice(1) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c("state_name", "year"), names_from = "policy_category", values_from = "Ig.score") %>%
  rename_with(.fn = ~ paste0("cat.", .x), .cols = 3:10) %>%
  mutate(Year = as.factor(year)) %>%
  left_join(policy.wgt.mort.county %>%   
              mutate(Year = as.numeric(as.character(Year))) %>% 
              select(FIPS.grp, state.names, Year, mx.standard.b, mx.standard.f, mx.standard.m),
            by = c("state_name" = "state.names", "year" = "Year")) %>% 
  select(FIPS.grp, state_name, Year, mx.standard.b, mx.standard.f, mx.standard.m, cat.1, cat.2, cat.3, cat.4, cat.5, cat.6, cat.7, cat.8)

# Check the column name indices to include in the correlation plots below
names(policy.mort.wgt.county.g)

# Run and plot Pearson's correlations on county group combination of ASDR and the kernel-weighted thematic Ig score categories by year
corr.plt.county.mx.Ig.wgt <- ggpairs(policy.mort.wgt.county.g,
                                         columns = c(4:6,7:14),
                                         aes(color = as.factor(Year), alpha = 0.5),
                                         upper = list("cor"),
                                         lower = list(continuous = "smooth")) +
  theme_minimal() +
  theme(plot.background = element_rect(color = "transparent", fill = "white"),
        panel.border = element_rect(color = "grey70", fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "grey70"),
        strip.text = element_text(size = 13))

corr.plt.county.mx.Ig.wgt

ggsave(plot = corr.plt.county.mx.Ig.wgt, file = "correlation_plot_county_mortality_Ig_score.png",
       height = 6, width = 6, units = "in", dpi = 150, scale = 2.3)




#######################
### State AOV analysis
#######################

# Checking normality via QQ plots 
mx.std.b.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.state, x = "mx.standard.b", add = "qqline", conf.int = T, facet.by = "Year", xlab = "")
mx.std.f.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.state, x = "mx.standard.f", add = "qqline", conf.int = T, facet.by = "Year", ylab = "")
mx.std.m.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.state, x = "mx.standard.m", add = "qqline", conf.int = T, facet.by = "Year", ylab = "", xlab = "")
qqplots.mx.std <- ggpubr::ggarrange(mx.std.b.plt, mx.std.f.plt, mx.std.m.plt, nrow = 1, ncol = 3, labels = c("Both", "F", "M"))
ggsave(qqplots.mx.std, file = "QQplots_state_mx_std_by_sex.png", width = 6, height = 3, units = 'in', dpi = 150, scale = 2.2)

# Shapiro-Wilk test (if p > 0.05 = normally distributed)
policy.wgt.mort.state %>% 
  group_by(Year) %>% 
  rstatix::shapiro_test(mx.standard.b)

# Two-way ANOVAs with quintile of SHP and lagged SHP
policy.wgt.mort.state %>%
  group_by(I.score.quintile, Year) %>% 
  rstatix::shapiro_test(mx.standard.b)

policy.wgt.mort.state %>%
  group_by(I.score.10yr.lag.quintile, Year) %>% 
  rstatix::shapiro_test(mx.standard.b)

# Check for outliers (one-way)
mortality.state %>%
  group_by(Year) %>% 
  rstatix::identify_outliers(mx.standard.b) %>% 
  ungroup() %>% 
  select(STATE, Year, mx.standard.b, is.outlier, is.extreme)

# Check for outliers (two-way)
policy.wgt.mort.state %>%
  group_by(I.score.quintile, Year) %>% 
  rstatix::identify_outliers(mx.standard.b) %>% 
  ungroup() %>% 
  select(STATE, Year, I.score.quintile, mx.standard.b, is.outlier, is.extreme)

# Checking normality via QQ plots (one-way)
log.mx.std.b.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.state, x = "log.mx.standard.b", add = "qqline", conf.int = T, facet.by = "Year", xlab = "")
log.mx.std.f.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.state, x = "log.mx.standard.f", add = "qqline", conf.int = T, facet.by = "Year", ylab = "")
log.mx.std.m.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.state, x = "log.mx.standard.m", add = "qqline", conf.int = T, facet.by = "Year", ylab = "", xlab = "")
qqplots.log.mx.std <- ggpubr::ggarrange(log.mx.std.b.plt, log.mx.std.f.plt, log.mx.std.m.plt, nrow = 1, ncol = 3, labels = c("Both", "F", "M"))
ggsave(qqplots.log.mx.std, file = "QQplots_state_log_mx_std_by_sex.png", width = 6, height = 2, units = 'in', dpi = 150, scale = 2.2)

# Shapiro-Wilk test for log-transformed mx (one-way) (if p > 0.05 = normally distributed)
policy.wgt.mort.state %>%
  group_by(Year) %>% 
  rstatix::shapiro_test(log.mx.standard.b)
policy.wgt.mort.state %>%
  group_by(Year) %>% 
  rstatix::shapiro_test(log.mx.standard.f)
policy.wgt.mort.state %>%
  group_by(Year) %>% 
  rstatix::shapiro_test(log.mx.standard.m)

# Shapro-Wilk test for log-transformed mx (two-way)
policy.wgt.mort.state %>%
  group_by(I.score.quintile, Year) %>% 
  rstatix::shapiro_test(log.mx.standard.b)


# Check for outliers with log-transformed mx (one-way)
policy.wgt.mort.state %>%
  group_by(Year) %>% 
  rstatix::identify_outliers(log.mx.standard.b) %>% 
  ungroup() %>% 
  select(STATE, Year, log.mx.standard.b, is.outlier, is.extreme)
policy.wgt.mort.state %>%
  group_by(Year) %>% 
  rstatix::identify_outliers(log.mx.standard.f) %>% 
  ungroup() %>% 
  select(STATE, Year, log.mx.standard.f, is.outlier, is.extreme)
policy.wgt.mort.state %>%
  group_by(Year) %>% 
  rstatix::identify_outliers(log.mx.standard.m) %>% 
  ungroup() %>% 
  select(STATE, Year, log.mx.standard.m, is.outlier, is.extreme)


# Check for outliers with log-transformed mx (two-way, ith I score)
policy.wgt.mort.state %>%
  group_by(I.score.10yr.lag.quintile, Year) %>% 
  rstatix::identify_outliers(log.mx.standard.b) %>% 
  ungroup() %>% 
  select(STATE, Year, I.score.quintile, log.mx.standard.b, is.outlier, is.extreme)


# Run ANOVA (one-way)
rstatix::anova_test(lm(log.mx.standard.b ~ Year, data = policy.wgt.mort.state), type = 3)
rstatix::anova_test(lm(log.mx.standard.f ~ Year, data = policy.wgt.mort.state), type = 3)
rstatix::anova_test(lm(log.mx.standard.m ~ Year, data = policy.wgt.mort.state), type = 3)

# Post-hoc t-tests (one-way)
policy.wgt.mort.state %>%
  rstatix::pairwise_t_test(log.mx.standard.b ~ Year, paired = TRUE, p.adjust.method = "bonferroni", detailed = TRUE)
policy.wgt.mort.state %>%
  rstatix::pairwise_t_test(log.mx.standard.f ~ Year, paired = TRUE, p.adjust.method = "bonferroni", detailed = TRUE)
policy.wgt.mort.state %>%
  rstatix::pairwise_t_test(log.mx.standard.m ~ Year, paired = TRUE, p.adjust.method = "bonferroni", detailed = TRUE)


# Run ANOVA (two-way)
rstatix::anova_test(lm(log.mx.standard.b ~ Year + I.score.quintile + Year:I.score.quintile, data = policy.wgt.mort.state), type = 3)
rstatix::anova_test(lm(log.mx.standard.f ~ Year + I.score.quintile + Year:I.score.quintile, data = policy.wgt.mort.state), type = 3)
rstatix::anova_test(lm(log.mx.standard.m ~ Year + I.score.quintile + Year:I.score.quintile, data = policy.wgt.mort.state), type = 3)


# Post-hoc t-tests (two-way)
print(policy.wgt.mort.state %>%
        group_by(Year) %>% 
        rstatix::pairwise_t_test(log.mx.standard.b ~ I.score.quintile, paired = TRUE, p.adjust.method = "bonferroni"),
      n = 100)
print(policy.wgt.mort.state %>%
        group_by(Year) %>% 
        rstatix::pairwise_t_test(log.mx.standard.f ~ I.score.quintile, paired = TRUE, p.adjust.method = "bonferroni"),
      n = 100)
print(policy.wgt.mort.state %>%
        group_by(Year) %>% 
        rstatix::pairwise_t_test(log.mx.standard.m ~ I.score.quintile, paired = TRUE, p.adjust.method = "bonferroni"),
      n = 100)


#######################
### County AOV analysis
#######################

mx.std.b.county.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.county, x = "mx.standard.b", add = "qqline", conf.int = T, facet.by = "Year", xlab = "")
mx.std.f.county.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.county, x = "mx.standard.f", add = "qqline", conf.int = T, facet.by = "Year", ylab = "")
mx.std.m.county.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.county, x = "mx.standard.m", add = "qqline", conf.int = T, facet.by = "Year", ylab = "", xlab = "")
qqplots.mx.std.county <- ggpubr::ggarrange(mx.std.b.county.plt, mx.std.f.county.plt, mx.std.m.county.plt, nrow = 1, ncol = 3, labels = c("Both", "F", "M"))
ggsave(qqplots.mx.std.county, file = "QQplots_county_mx_std_by_sex.png", width = 6, height = 2, units = 'in', dpi = 150, scale = 2.2)

# Shapiro-Wilk test for normality (one-way) NOT APPROPRIATE FOR LARGE SAMPLE SIZES OVER 50!!! WILL REJECT NULL
policy.wgt.mort.county %>%
  group_by(Year) %>% 
  rstatix::shapiro_test(mx.standard.b)

# Check for outliers (one-way)
print(policy.wgt.mort.county %>%
  group_by(Year) %>% 
  rstatix::identify_outliers(mx.standard.b) %>% 
  ungroup() %>% 
  select(FIPS.grp, Year, mx.standard.b, is.outlier, is.extreme),
  n = 1000)

log.mx.std.b.county.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.county, x = "log.mx.standard.b", add = "qqline", conf.int = T, facet.by = "Year", xlab = "")
log.mx.std.f.county.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.county, x = "log.mx.standard.f", add = "qqline", conf.int = T, facet.by = "Year", ylab = "")
log.mx.std.m.county.plt <- ggpubr::ggqqplot(data = policy.wgt.mort.county, x = "log.mx.standard.m", add = "qqline", conf.int = T, facet.by = "Year", ylab = "", xlab = "")
qqplots.log.mx.std.county <- ggpubr::ggarrange(log.mx.std.b.county.plt, log.mx.std.f.county.plt, log.mx.std.m.county.plt, nrow = 1, ncol = 3, labels = c("Both", "F", "M"))
ggsave(qqplots.log.mx.std.county, file = "QQplots_county_log_mx_std_by_sex.png", width = 6, height = 2, units = 'in', dpi = 150, scale = 2.2)

policy.wgt.mort.county %>%
  group_by(Year) %>% 
  rstatix::shapiro_test(log.mx.standard.b)

# Check for outliers on log-transformed mx (one-way)
print(policy.wgt.mort.county %>%
        mutate(log.mx.standard.b = log(mx.standard.b)) %>% 
        group_by(Year) %>% 
        rstatix::identify_outliers(log.mx.standard.b) %>% 
        ungroup() %>% 
        select(FIPS.grp, Year, log.mx.standard.b, is.outlier, is.extreme),
        n = 1000)
print(policy.wgt.mort.county %>%
        mutate(log.mx.standard.f = log(mx.standard.f)) %>% 
        group_by(Year) %>% 
        rstatix::identify_outliers(log.mx.standard.f) %>% 
        ungroup() %>% 
        select(FIPS.grp, Year, log.mx.standard.f, is.outlier, is.extreme),
      n = 1000)
print(policy.wgt.mort.county %>%
        mutate(log.mx.standard.m = log(mx.standard.m)) %>% 
        group_by(Year) %>% 
        rstatix::identify_outliers(log.mx.standard.m) %>% 
        ungroup() %>% 
        select(FIPS.grp, Year, log.mx.standard.m, is.outlier, is.extreme),
      n = 1000)


# Run ANOVA (one-way)
rstatix::anova_test(lm(log.mx.standard.b ~ Year, data = policy.wgt.mort.county), type = 3)
rstatix::anova_test(lm(log.mx.standard.f ~ Year, data = policy.wgt.mort.county), type = 3)
rstatix::anova_test(lm(log.mx.standard.m ~ Year, data = policy.wgt.mort.county), type = 3)


# Post-hoc t-tests (one-way)
policy.wgt.mort.county %>%
  rstatix::pairwise_t_test(log.mx.standard.b ~ Year, paired = TRUE, p.adjust.method = "bonferroni", detailed = TRUE)
policy.wgt.mort.county %>%
  rstatix::pairwise_t_test(log.mx.standard.f ~ Year, paired = TRUE, p.adjust.method = "bonferroni", detailed = TRUE)
policy.wgt.mort.county %>%
  rstatix::pairwise_t_test(log.mx.standard.m ~ Year, paired = TRUE, p.adjust.method = "bonferroni", detailed = TRUE)

# Run ANOVA (two-way, with SHP measures interacting with Year)
rstatix::anova_test(lm(log.mx.standard.b ~ Year + I.score.quintile + Year:I.score.quintile, data = policy.wgt.mort.county), type = 3)
rstatix::anova_test(lm(log.mx.standard.f ~ Year + I.score.quintile + Year:I.score.quintile, data = policy.wgt.mort.county), type = 3)
rstatix::anova_test(lm(log.mx.standard.m ~ Year + I.score.quintile + Year:I.score.quintile, data = policy.wgt.mort.county), type = 3)


