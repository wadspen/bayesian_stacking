library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(lubridate)
library(tidyr)
source("./simulations/stack_functions.R")

bma_bps <- readRDS("./flu_analysis/bma_bps_crps.rds")
sgp_eq <- readRDS("./flu_analysis/sgp_eq.rds")
all_flu <- read.csv("../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv")

sgp_eq %>% 
  mutate(rel_stack = stack_crps/eq_crps) %>% 
  group_by(location) %>% 
  summarise(mcrps = mean(rel_stack)) %>% 
  ggplot() +
  geom_point(aes(x = mcrps, y = location))


all_crps <- sgp_eq %>% 
  mutate(season_week = week) %>% 
  select(location, season_week, forecast_date, stack_crps, eq_crps) %>% 
  left_join(bma_bps %>% 
              select(location, season_week, bma_crps, bps_crps),
            by = c("location", "season_week"))



loc_mean_crps <- all_crps %>% 
  pivot_longer(4:7, names_to = "method", values_to = "crps") %>% 
  group_by(location, method) %>% 
  summarise(mcrps = mean(crps)) %>% 
  left_join(all_flu %>% 
              select(location, location_name) %>% 
              unique(), by = "location")

sgp_order <- loc_mean_crps %>% 
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
             size = 4) +
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
        legend.position = c(.13,.8))



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
  ggplot() +
  geom_line(aes(x = season_week, y = mcrps, colour = method,
                linetype = method), size = 1.3) +
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
  mutate(rank = cumsum(ind)) %>% 
  group_by(method) %>% 
  summarise(table(rank))

  
  
  
  
  
