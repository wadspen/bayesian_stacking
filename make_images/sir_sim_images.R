library(ggplot2)
library(dplyr)
setwd(paste0(here::here(), "/make_images"))
sim_scores <- read.csv("../simulations/sir_res/sir_res.csv") %>%
select(-X)
# sim_scores$time <- rep(rep(1:38, 4), 1000)
source("../simulations/stack_functions.R")



sim_uwd1s <- sim_scores %>% 
  group_by(method, seq, rep) %>% 
  summarise(uwd1 = unit_wass_dist(ecdf(pit)))


pit_box <-sim_uwd1s %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = uwd1), size = 1) +
  ylab("UWD1") +
  xlab("") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 18),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))



pit_hist <- sim_scores %>% 
  # filter(replicate == 6) %>% 
  ggplot() +
  geom_histogram(aes(x = pit, y = ..density..)) +
  facet_wrap(~method, scales = "free") +
  ylab("") +
  xlab("PIT") +
  scale_x_continuous(breaks = seq(0,1, by = .5)) +
  theme_bw() +
  theme(axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(size = 11),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))



logsp <- sim_scores %>%
  filter(logs != Inf) %>% 
  filter(time <= 30 & time > 1) %>%
  # filter(methods == "BMA") %>%
  group_by(method, time) %>%
  summarise(mlogs = mean(logs, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = time, y = mlogs, colour = method,
                linetype = method), size = 1.1) +
  scale_colour_manual(name = "Method",
                      labels = c("AVS", "BMA", "EQW", "SGP", "SGP50")
                      ,values = c("grey80", "grey60", "grey40",
                                  "grey20", "grey0")) +
  scale_linetype_manual(name= "Method",
                        values=c("twodash", "dotdash", "longdash",
                                 "solid", "dotted"),
                        labels=c("AVS", "BMA", "EQW", "SGP", "SGP50")) +
  ylab("LogS") +
  xlab("") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 15),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
        ,legend.position = c(.85, .72)
  )

crpsp <- sim_scores %>%
  filter(logs != Inf) %>% 
  filter(time <= 30 & time > 1) %>%
  # filter(methods == "BMA") %>%
  group_by(method, time) %>%
  summarise(mcrps = mean(crps, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = time, y = mcrps, colour = method,
  linetype = method), size = 1.1) +
  scale_colour_manual(name = "Method",
                      labels = c("AVS", "BMA", "EQW", "SGP", "SGP50")
                      ,values = c("grey80", "grey60", "grey40", 
                                  "grey20", "grey0")) +
  scale_linetype_manual(name= "Method",
                        values=c("twodash", "dotdash", "longdash", 
                                 "solid", "dotted"),
                        labels=c("AVS", "BMA", "EQW", "SGP", "SGP50")) +
  ylab("CRPS") +
  xlab("Time") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 15),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
        ,legend.position = "none"
       )



cowplot::plot_grid(logsp, pit_hist, crpsp, pit_box, nrow = 2, 
                   ncol = 2, rel_heights = c(1,1))




crpsp <- sim_scores %>%
  filter(logs != Inf) %>% 
  filter(time <= 30 & time > 1) %>%
  # filter(methods == "BMA") %>%
  group_by(method, time) %>%
  summarise(mcrps = mean(crps, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = time, y = mcrps, colour = method,
                linetype = method), size = 1.1) +
  scale_colour_manual(name = "Method",
                      labels = c("AVS", "BMA", "EQW", "SGP", "SGP50")
                      ,values = c("#E69F00", "#56B4E9", "#009E73", 
                                  "#D55E00", "#CC79A7")) +
  scale_linetype_manual(name= "Method",
                        values=c("twodash", "dotdash", "longdash", 
                                 "solid", "dotted"),
                        labels=c("AVS", "BMA", "EQW", "SGP", "SGP50")) +
  ylab("CRPS") +
  xlab("Time") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 15),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
        ,legend.position = c(.85, .8)
  )

cowplot::plot_grid(crpsp, pit_box, nrow = 1, 
                   ncol = 2, rel_heights = c(1,1))


################################################
#################median plots###################
################################################


logsp <- sim_scores %>%
  filter(logs != Inf) %>% 
  filter(time <= 30 & time > 1) %>%
  # filter(methods == "BMA") %>%
  group_by(method, time) %>%
  summarise(mlogs = median(logs, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = time, y = mlogs, colour = method,
                linetype = method), size = 1.1) +
  scale_colour_manual(name = "Method",
                      labels = c("AVS", "BMA", "EQW", "SGP", "SGP50")
                      ,values = c("grey80", "grey60", "grey40",
                                  "grey20", "grey0")) +
  scale_linetype_manual(name= "Method",
                        values=c("twodash", "dotdash", "longdash",
                                 "solid", "dotted"),
                        labels=c("AVS", "BMA", "EQW", "SGP", "SGP50")) +
  ylab("LogS") +
  xlab("") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 15),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
        ,legend.position = c(.87, .75)
  )

crpsp <- sim_scores %>%
  filter(logs != Inf) %>% 
  filter(time <= 30 & time > 1) %>%
  # filter(methods == "BMA") %>%
  group_by(method, time) %>%
  summarise(mcrps = median(crps, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = time, y = mcrps, colour = method,
                linetype = method), size = 1.1) +
  scale_colour_manual(name = "Method",
                      labels = c("AVS", "BMA", "EQW", "SGP", "SGP50")
                      ,values = c("grey80", "grey60", "grey40", 
                                  "grey20", "grey0")) +
  scale_linetype_manual(name= "Method",
                        values=c("twodash", "dotdash", "longdash", 
                                 "solid", "dotted"),
                        labels=c("AVS", "BMA", "EQW", "SGP", "SGP50")) +
  ylab("CRPS") +
  xlab("Time") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 15),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
        ,legend.position = "none"
  )



cowplot::plot_grid(logsp, crpsp, nrow = 1, 
                   ncol = 2, rel_heights = c(1,1))

sim_scores %>% 
  filter(logs != Inf) %>% 
  filter(time > 1, time <= 30) %>% 
  group_by(method) %>% 
  summarise(mean(crps))


sim_scores %>%
  filter(time > 1) %>% # & time <= 30) %>%
  arrange(seq, replicate, time, logs) %>%
  mutate(ind = 1) %>%
  group_by(seq, rep, time) %>%
  mutate(rank = cumsum(ind)) %>%
  group_by(method) %>%
  summarise(median(rank))




