# Final script used for Flocks in Flux analysis - temporal modelling and figure production

library(lmerTest)
library(lme4)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(broom)
library(purrr)
library(scales)
library(gratia)
library(mgcv)
library(showtext)
library(ggh4x)

setwd("C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/Radar Data")

# Module fragmentation
nodes_final <- read.csv("nodes_final.csv")

# Safety filter to remove unused years
affected_years <- c(1995:2004)
nodes_final <- nodes_final %>% filter(!year %in% affected_years)
n_distinct(nodes_final$year)

roostyear <- nodes_final %>% filter(year == 2023)
hist(roostyear$mean_size)

# Filter true modules (non-zero modules & modules with >1 roost)
communities <- nodes_final %>% group_by(year) %>% 
  mutate(n_birds = sum(mean_size, na.rm = T), 
         n_communities = n_distinct(module[module != 0]),
         mean_roost_size = mean(mean_size),
         average_annual_degree = mean(annual_degree),
         na.rm = T) %>% 
  ungroup() %>%
  group_by(year, module) %>%
  summarise(
    mean_roost_size = first(mean_roost_size),
    n_birds = first(n_birds),
    gini = first(gini_coeff),
    n_module_nodes = n_distinct(CLUSTER_ID),
    module_size = sum(mean_size), 
    module_strength = sum(strength), 
    annual_network_density = first(density),
    average_annual_degree = first(average_annual_degree),
    clustering_coeff = first(clustering_coeff)) %>%
  group_by(module) %>%
  arrange(year) %>%
  mutate(delta_roosts = n_module_nodes - lag(n_module_nodes), 
         years_active = n_distinct(year),
         year = as.numeric(year)) %>%
  ungroup()
  
# module-level stats
annual_summary <- communities %>%
  group_by(year) %>% 
  mutate(n_communities = n_distinct(module)) %>%
  reframe(
    mean_roost_size = first(mean_roost_size),
    gini = mean(gini[module != 0 | n_module_nodes > 2], na.rm = T),
    n_communities = first(n_communities),
    avg_community_size = mean(module_size[module != 0 | n_module_nodes > 2], na.rm = T), 
    n_birds   = first(n_birds),
    avg_roosts_per_comm = mean(n_module_nodes[module != 0 | n_module_nodes > 2]),
    n_roosts = sum(n_module_nodes),
    prop_n_island = sum(n_module_nodes[module == 0]) / sum(n_roosts),
    prop_n_dyads = sum(n_module_nodes[n_module_nodes == 2]) / sum(n_roosts),
    density = first(annual_network_density),
    average_annual_degree = first(average_annual_degree),
    clustering_coeff =first(clustering_coeff),
    year = as.character(year)) %>% 
  ungroup() %>%
  arrange(year) %>% 
  group_by(year) %>% 
  slice_head(n = 1) %>% 
  mutate(year = as.numeric(year)) %>% ungroup()

print(1- annual_summary$clustering_coeff[12]/annual_summary$clustering_coeff[1])

summary(annual_summary)
nodes_final <- nodes_final %>% left_join(annual_summary, by = "year") %>% mutate(CLUSTER_ID = as.character(CLUSTER_ID))

s_an_sum <- annual_summary %>%
  mutate(across(where(is.numeric), scale)) %>% mutate(year = as.numeric(year), n_birds = as.numeric(n_birds))

s_an_sum <- annual_summary %>% mutate(year = scale(as.numeric(year)),
                                          n_birds = scale(n_birds))

########################################################
# WHAT ARE THE UNIVARIATE TRENDS OVER TIME?

metric_colors <- c(
  "Number of Birds" = "#79af56",
  "Number of Roosts" = "#79af56",
  "Number of\nCommunities" = "#79af56",
  "Mean Roost Size" = "#6b95c8",
  "Mean Community Size" = "#6b95c8",
  "Mean Roosts\nper Community" = "#6b95c8",
  "Mean Gini Coefficient" = "#c56d68",
  "Mean Roost Degree" = "#c56d68",
  "Clustering Coefficient" = "#c56d68")

the_plot <- annual_summary %>%
  dplyr::select(year,
         `Number of Birds` = n_birds,
         `Number of Roosts` = n_roosts,
         `Number of\nCommunities` = n_communities,
         `Mean Roost Size` = mean_roost_size,
         `Mean Community Size` = avg_community_size,
         `Mean Roosts\nper Community` = avg_roosts_per_comm,
         `Mean Gini Coefficient` = gini,
         `Mean Roost Degree` = average_annual_degree,
         `Clustering Coefficient` = clustering_coeff) %>%
  pivot_longer(cols = c(-year), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c(
    "Number of Birds",
    "Number of Roosts",
    "Number of\nCommunities",
    "Mean Roost Size",
    "Mean Community Size",
    "Mean Roosts\nper Community",
    "Mean Gini Coefficient",
    "Mean Roost Degree",
    "Clustering Coefficient")))

showtext_auto()
plot <- ggplot(the_plot, aes(x = as.numeric(year), y = value, color = metric)) +
  geom_smooth(aes(fill = metric), method = "gam", se = TRUE, linewidth = .75, alpha = 0.2) +
  geom_point(size = .75, alpha = 1) +
  scale_color_manual(values = metric_colors) +
  scale_fill_manual(values = metric_colors) +
  guides(fill = "none", color = "none") +
  labs(x = "Year", y = "Values") +
  facet_wrap2(vars(metric), scales = "free_y", remove_labels = "x", ncol = 3, nrow = 3) +
  facetted_pos_scales(
    y = list(
      metric == "Number of Birds" ~ scale_y_continuous(
        labels = scales::label_number(scale = 1e-6, suffix = "M")),
      metric == "Number of Roosts" ~ scale_y_continuous(),
      metric == "Number of\nCommunities" ~ scale_y_continuous(),
      metric == "Mean Roost Size" ~ scale_y_continuous(
        labels = scales::label_number(scale = 1e-2, suffix = "k")),
      metric == "Mean Community Size" ~ scale_y_continuous(
        labels = scales::label_number(scale = 1e-3, suffix = "k")),
      metric == "Mean Roosts\nper Community" ~ scale_y_continuous(),
      metric == "Mean Gini Coefficient" ~ scale_y_continuous(),
      metric == "Mean Roost Degree" ~ scale_y_continuous(),
      metric == "Clustering Coefficient" ~ scale_y_continuous(
        labels = scales::label_number(suffix = "%")))) +
geom_text(
    data = data.frame(
      metric = unique(the_plot$metric),
      label = paste0(letters[1:length(unique(the_plot$metric))])),
    aes(x = -Inf, y = Inf, label = label),
    hjust = 1,
    vjust = -1.5,
    family = "Kigela",
    fontface = "bold",
    inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 7) +
  theme(text = element_text(family = "Kigela"),
    legend.position = "none",
    strip.text = element_text(lineheight = 1.75),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "#e2e3e8", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_blank()
    )
plot

ggsave("Figure_4.jpeg", plot)

the_plot2 <- annual_summary %>% 
   dplyr::select(n_birds,
                 `Mean Community Size` = avg_community_size,
                 `Number of Communities` = n_communities,
                 `Number of Roosts` = n_roosts,
                 `Mean Roost Degree` = average_annual_degree,
                 `Clustering Coefficient` = clustering_coeff,
                 `Mean Roosts per Community` = avg_roosts_per_comm,
                 `Mean Gini Coefficient` = gini,
                 `Mean Roost Size` = mean_roost_size) %>%
   pivot_longer(cols = c(-n_birds), names_to = "metric", values_to = "value")
 
 
 plot2 <- ggplot(the_plot2, aes(x = n_birds, y = value, color = metric)) +
   geom_smooth(aes(fill = metric), method = "gam", se = TRUE, 
               linewidth = 2, alpha = 0.2) +
   geom_point(size = 2, alpha = 0.7) +
   scale_color_manual(values = metric_colors) +
   scale_fill_manual(values = metric_colors) +
   guides(fill = "none", color = "none") +
   scale_x_continuous(labels =  label_number(scale_cut = cut_short_scale())) +
   labs(x = "Population Size", y = "Values") +
   facet_wrap(vars(metric), scale = "free_y") +
   theme_minimal()

 plot2
 
 # Birds ~ year (baseline trend)
 AIC(gam(n_birds ~ s(year), data = annual_summary),#
     lm(n_birds ~ year, data = annual_summary))
 
 m_birds <- gam(n_birds ~ s(year), data = annual_summary)
 summary(m_birds)
 
 # Total Number of Roosts
 AIC(
   gam(n_roosts ~ te(year, n_birds), data = annual_summary),
   gam(n_roosts ~ s(year) + s(n_birds), data = annual_summary),#
   gam(n_roosts ~ s(year), data = annual_summary), 
   gam(n_roosts ~ s(n_birds), data = annual_summary),
   lm(n_roosts ~ year * n_birds, data = annual_summary),
   lm(n_roosts ~ year + n_birds, data = annual_summary),
   lm(n_roosts ~ year, data = annual_summary),
   lm(n_roosts ~ n_birds, data = annual_summary))
 
 m_n_roosts <- gam(n_roosts ~ te(year, n_birds), data = annual_summary)
 summary(m_n_roosts)
 vis.gam(m_n_roosts,theta=35, phi=30)
 
 # Number of Communities
 AIC(
   gam(n_communities ~ te(year, n_birds), data = annual_summary), #
   gam(n_communities ~ s(year), data = annual_summary),
   gam(n_communities ~ s(n_birds), data = annual_summary),
   lm(n_communities ~ year * n_birds, data = annual_summary),
   lm(n_communities ~ year, data = annual_summary),
   lm(n_communities ~ n_birds, data = annual_summary))
 
 m_n_comms <- gam(n_communities ~ te(year, n_birds),data = annual_summary)
 summary(m_n_comms)
 vis.gam(m_n_comms,theta=35, phi=30)
 
 # Average Roost Size
 AIC(
   gam(mean_roost_size ~ te(year, n_birds), data = annual_summary), #
   gam(mean_roost_size ~ s(year), data = annual_summary),
   gam(mean_roost_size ~ s(n_birds), data = annual_summary),
   lm(mean_roost_size ~ year * n_birds, data = annual_summary),
   lm(mean_roost_size ~ year, data = annual_summary),
   lm(mean_roost_size ~ n_birds, data = annual_summary))
 
 m_mean_roost <- gam(mean_roost_size ~ te(year, n_birds),data = annual_summary)
 summary(m_mean_roost)
 vis.gam(m_mean_roost,theta=35, phi=30)

# Mean Community Size
 AIC(
   gam(avg_community_size ~ te(year, n_birds), data = annual_summary), #
   gam(avg_community_size ~ s(year), data = annual_summary),
   gam(avg_community_size ~ s(n_birds), data = annual_summary),
   lm(avg_community_size ~ year * n_birds, data = annual_summary),
   lm(avg_community_size ~ year, data = annual_summary),
   lm(avg_community_size ~ n_birds, data = annual_summary))
 
 m_comm_size <- gam(avg_community_size ~ te(year, n_birds),data = annual_summary)
 summary(m_comm_size)
 vis.gam(m_comm_size,theta=35, phi=30)

#Roost loss continues even when population size stabilizes or rebounds.
 
# Average Roosts per Community
 AIC(
   gam(avg_roosts_per_comm ~ te(year, n_birds), data = annual_summary), #
   gam(avg_roosts_per_comm ~ s(year), data = annual_summary),
   gam(avg_roosts_per_comm ~ s(n_birds), data = annual_summary),
   lm(avg_roosts_per_comm ~ year * n_birds, data = annual_summary),
   lm(avg_roosts_per_comm ~ year, data = annual_summary),
   lm(avg_roosts_per_comm ~ n_birds, data = annual_summary))
 
 m_roosts_per_comm <- gam(avg_roosts_per_comm ~ te(year, n_birds), data = annual_summary)
 summary(m_roosts_per_comm)
 vis.gam(m_roosts_per_comm,theta=35, phi=30)

 # Gini Coefficient
 AIC(
   gam(gini ~ te(year, n_birds), data = annual_summary), 
   gam(gini ~ s(year), data = annual_summary),
   gam(gini ~ s(n_birds), data = annual_summary), #
   lm(gini ~ year * n_birds, data = annual_summary),
   lm(gini ~ year, data = annual_summary),
   lm(gini ~ n_birds, data = annual_summary))
 

 m_gini <- gam(gini ~ s(n_birds), data = annual_summary)
 summary(m_gini)

 # Average degree
 AIC(
   gam(average_annual_degree ~ te(year, n_birds), data = annual_summary),
   gam(average_annual_degree ~ s(year), data = annual_summary),#
   gam(average_annual_degree ~ s(n_birds), data = annual_summary),
   lm(average_annual_degree ~ year * n_birds, data = annual_summary),
   lm(average_annual_degree ~ year, data = annual_summary),
   lm(average_annual_degree ~ n_birds, data = annual_summary))
 
 m_degree <- gam(average_annual_degree ~ s(year), data = annual_summary)
 summary(m_degree)

 # Clustering Coefficient
 AIC(
   gam(clustering_coeff ~ te(year, n_birds), data = annual_summary),#
   gam(clustering_coeff ~ s(year), data = annual_summary),
   gam(clustering_coeff ~ s(n_birds), data = annual_summary),
   lm(clustering_coeff ~ year * n_birds, data = annual_summary), 
   lm(clustering_coeff ~ year, data = annual_summary),
   lm(clustering_coeff ~ n_birds, data = annual_summary))
 
 summary(
   lm(clustering_coeff ~ n_birds, data = annual_summary))
 
 m_cl_coef <- gam(clustering_coeff ~ te(year, n_birds,k=4), data = annual_summary)
 summary(m_cl_coef)
 vis.gam(m_cl_coef)
 