library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(lubridate)
library(tidyr)
source("./simulations/stack_functions.R")

bma_bps <- readRDS("./flu_analysis/bma_bps_crps.rds")
sgp_eq <- readRDS("./flu_analysis/sgp_eq.rds") %>% 
  filter(week < 30)
all_flu <- read.csv("../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv")

# sgp_eq %>% 
#   mutate(rel_stack = stack_crps/eq_crps) %>% 
#   group_by(location) %>% 
#   summarise(mcrps = mean(rel_stack)) %>% 
#   ggplot() +
#   geom_point(aes(x = mcrps, y = location))


all_crps <- sgp_eq %>% 
  mutate(season_week = week) %>% 
  select(location, season_week, forecast_date, stack_crps, eq_crps) %>% 
  left_join(bma_bps %>% 
              select(location, season_week, bma_crps, bps_crps),
            by = c("location", "season_week"))



loc_mean_crps <- all_crps %>% 
  # mutate(bma_crps = bma_crps/stack_crps,
  #        bps_crps = bps_crps/stack_crps,
  #        eq_crps = eq_crps/stack_crps,
  #        stack_crps = stack_crps/stack_crps) %>% 
  pivot_longer(4:7, names_to = "method", values_to = "crps") %>% 
  group_by(location, method) %>% 
  summarise(mcrps = mean(crps)) %>% 
  mutate(ind = 1) %>% 
  arrange(location, mcrps) %>% 
  group_by(location) %>% 
  mutate(rank = cumsum(ind)) %>% 
  mutate(base_crps = mcrps[method == "stack_crps"]) %>% 
  ungroup() %>% 
  mutate(mcrps = mcrps/base_crps) %>% 
  left_join(all_flu %>% 
              select(location, location_name) %>% 
              unique(), by = "location")

sgp_order <- loc_mean_crps %>% 
  group_by(location_name) %>%
  mutate(mcrps = mcrps/min(mcrps)) %>%
  filter(method == "stack_crps") %>% 
  arrange(mcrps)


loc_mean_crps %>% 
  # filter(method != "bma_crps") %>%
  # filter(!(method %in% c("bma_crps", "eq_crps"))) %>%
  mutate(method = ifelse(method == "bma_crps", "BMA", 
                         ifelse(method == "bps_crps", "AVS", 
                                ifelse(method == "eq_crps", "EQW",
                                       ifelse(method == "stack_crps", "SGP", 
                                              NA))))) %>% 
  mutate(location_name = factor(location_name, 
                                levels = sgp_order$location_name)) %>%
  ggplot() +
  geom_point(aes(y = location_name, x = mcrps, fill = method, shape = method),
             size = 2.5ß) +
  scale_shape_manual(values=c(21:24)) +
  xlab("CRPS") +
  ylab("Region") +
  labs(fill = "Method", shape = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size = 14),
        axis.title=element_text(size=20),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = c(.15,.25)
        # ,legend.position = "none"
        )


loc_scores <- loc_mean_crps %>% 
  # filter(method != "bma_crps") %>%
  # filter(!(method %in% c("bma_crps", "eq_crps"))) %>%
  mutate(method = ifelse(method == "bma_crps", "BMA", 
                         ifelse(method == "bps_crps", "AVS", 
                                ifelse(method == "eq_crps", "EQW",
                                       ifelse(method == "stack_crps", "SGP", 
                                              NA))))) %>% 
  mutate(location_name = factor(location_name, 
                                levels = sgp_order$location_name))

loc_scores %>% 
  ggplot() +
  geom_boxplot(aes(y = mcrps, x = method))


loc_scores %>% 
  # filter(location_name != "New Hampshire") %>% 
  arrange(location_name, mcrps) %>% 
  mutate(ind = 1) %>% 
  group_by(location_name) %>% 
  mutate(rank = cumsum(ind)) %>% 
  group_by(method) %>% 
  reframe(table(rank))


unique_loc <- unique(all_crps$location)



all_crps %>% 
  pivot_longer(4:7, names_to = "method", values_to = "crps") %>% 
  mutate(method = ifelse(method == "bma_crps", "BMA", 
                         ifelse(method == "bps_crps", "AVS", 
                                ifelse(method == "eq_crps", "EQW",
                                       ifelse(method == "stack_crps", "SGP", 
                                              NA))))) %>%
  # filter(method != "bma_crps") %>% 
  group_by(season_week, method) %>% 
  # filter(location == unique_loc[i]) %>%
  summarise(mcrps = mean(crps)) %>% 
  group_by(season_week) %>% 
  # mutate(base_crps = mcrps[method == "SGP"]) %>% 
  # mutate(mcrps = mcrps/base_crps) %>% 
  ggplot() +
  geom_line(aes(x = season_week, y = mcrps, colour = method,
                linetype = method), size = .7) +
  xlab("Week")+
  ylab("CRPS") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size = 14),
        axis.title=element_text(size=20),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.position = c(.8,.7))


all_crps %>% 
  pivot_longer(4:7, names_to = "method", values_to = "crps") %>% 
  mutate(method = ifelse(method == "bma_crps", "BMA", 
                         ifelse(method == "bps_crps", "AVS", 
                                ifelse(method == "eq_crps", "EQW",
                                       ifelse(method == "stack_crps", "SGP", 
                                              NA))))) %>%
  # filter(method != "bma_crps") %>% 
  group_by(season_week, method) %>% 
  # filter(location == unique_loc[i]) %>%
  summarise(mcrps = mean(crps)) %>% 
  mutate(base_crps = mcrps[method == "SGP"]) %>% 
  mutate(mcrps = mcrps/base_crps) %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = mcrps))


wk_scores <- all_crps %>% 
  filter(location != "33") %>% 
  pivot_longer(4:7, names_to = "method", values_to = "crps") %>% 
  mutate(method = ifelse(method == "bma_crps", "BMA", 
                         ifelse(method == "bps_crps", "AVS", 
                                ifelse(method == "eq_crps", "EQW",
                                       ifelse(method == "stack_crps", "SGP", 
                                              NA))))) %>%
  # filter(method != "bma_crps") %>% 
  group_by(season_week, method) %>% 
  # filter(location == unique_loc[i]) %>%
  summarise(mcrps = mean(crps))
  
wk_scores %>% 
  arrange(season_week, mcrps) %>% 
  mutate(ind = 1) %>% 
  group_by(season_week) %>% 
  mutate(rank = cumsum(ind)) %>% #filter(rank == 1, method == "SGP")
  group_by(method) %>% 
  reframe(table(rank))


all_crps %>% 
  # filter(location != "33") %>% 
  pivot_longer(4:7, names_to = "method", values_to = "crps") %>% 
  mutate(method = ifelse(method == "bma_crps", "BMA", 
                         ifelse(method == "bps_crps", "AVS", 
                                ifelse(method == "eq_crps", "EQW",
                                       ifelse(method == "stack_crps", "SGP", 
                                              NA))))) %>% 
  group_by(method) %>% 
  summarise(mean(crps))



all_crps_long <- all_crps %>% 
  pivot_longer(4:7, names_to = "method", values_to = "crps") %>% 
  mutate(method = ifelse(method == "bma_crps", "BMA", 
                         ifelse(method == "bps_crps", "AVS", 
                                ifelse(method == "eq_crps", "EQW",
                                       ifelse(method == "stack_crps", "SGP", 
                                              NA)))))


all_crps_long %>% 
  arrange(location, season_week, crps) %>% 
  filter(season_week != 1) %>% 
  mutate(ind = 1) %>% 
  group_by(location, season_week) %>% 
  mutate(rank = cumsum(ind)) %>% filter(rank == 1) %>% group_by(method) %>% summarise(n())
  group_by(method) %>% 
  reframe(table(rank))

  
  
  
  
  
