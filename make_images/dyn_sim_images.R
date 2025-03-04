library(dplyr)
library(ggplot2)
library(tidyr)

setwd(paste0(here::here(), "/make_images"))

dyn <- readRDS("../simulations/dynamic_stacks_all_weeks_fin.rds") %>% 
  filter(time < 50)

dyn <- dyn %>% 
  mutate(method = ifelse(method == "avs", "AVS",
                         ifelse(method == "bma", "BMA",
                                ifelse(method == "eqw", "EQW", 
                                       ifelse(method == "sgp", "SGP", NA)))))

dyn_uwd1s <- dyn %>% 
  group_by(method, replicate) %>% 
  summarise(uwd1 = unit_wass_dist(ecdf(pit)))


pit_box <- dyn_uwd1s %>% 
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

pit_hist <- dyn %>% 
  # filter(replicate == 6) %>% 
  ggplot() +
  geom_histogram(aes(x = pit, y = ..density..)) +
  facet_wrap(~method, scales = "free") +
  ylab("") +
  xlab("PIT") +
  theme_bw() +
  theme(axis.text.y=element_text(size=11),
        axis.text.x=element_text(size = 11),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))


# dyn %>% 
#   ggplot() +
#   geom_boxplot(aes(x = method, y = mlogs))


logsp <- dyn %>% 
  # filter(time > 1) %>% 
  group_by(method, time) %>% 
  summarise(mcrpss = mean(mcrps),
            mlogss = mean(mlogs),
            selogs = sd(mlogs)/500) %>% 
  ggplot() + 
  geom_line(aes(x = time, y = mlogss, 
                colour = method, linetype = method), size = 1.1) +
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
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.position = c(.91,.78))



crpsp <- dyn %>% 
  # filter(time > 1) %>% 
  group_by(method, time) %>% 
  summarise(mcrpss = mean(mcrps),
            mlogss = mean(mlogs),
            selogs = sd(mlogs)/500,
            secrps = sd(mcrps)/500) %>% 
  ggplot() + 
  geom_line(aes(x = time, y = mcrpss, 
                colour = method, linetype = method), size = 1.1) +
  scale_colour_manual(name = "Method",
                      labels = c("AVS", "BMA", "EQW", "SGP", "SGP50")
                      ,values = c("grey80", "grey60", "grey40",
                                  "grey20", "grey0")) +
  scale_linetype_manual(name= "Method",
                        values=c("twodash", "dotdash", "longdash",
                                 "solid", "dotted"),
                        labels=c("AVS", "BMA", "EQW", "SGP", "SGP50")) +
  ylab("CRPS") +
  xlab("t") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 15),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.position = "none")

cowplot::plot_grid(logsp, pit_hist, crpsp, pit_box, nrow = 2, 
                   ncol = 2, rel_heights = c(1,1))


# dyn %>% 
#   group_by(time, replicate) %>% 
#   mutate(mincrps = min(mlogs)) %>% 
#   mutate(rel_crps = mlogs/mincrps) %>% 
#   # group_by(method) %>% 
#   # summarise(m = mean(rel_crps))
#   ggplot() +
#   geom_boxplot(aes(x = method, y = rel_crps)) #+
#   facet_wrap(~method)



library(kableExtra)
dyn_sum <- dyn %>% 
  group_by(method) %>% 
  summarise(mcrps = mean(mcrps),
            mlogs = mean(mlogs),
            uwd1 = unit_wass_dist(ecdf(pit)))

kbl(dyn_sum)









