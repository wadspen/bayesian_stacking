library(dplyr)
library(ggplot2)
library(tidyr)



dyn <- readRDS("./simulations/dynamic_stacks_all_weeks.rds")

dyn <- dyn %>% 
  mutate(method = ifelse(method == "avs", "AVS",
                         ifelse(method == "bma", "BMA",
                                ifelse(method == "eqw", "EQW", 
                                       ifelse(method == "sgp", "SGP", NA)))))

dyn_uwd1s <- dyn %>% 
  group_by(method, replicate) %>% 
  summarise(uwd1 = unit_wass_dist(ecdf(pit)))


dyn_uwd1s %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = uwd1), size = 1) +
  ylab("UWD1") +
  xlab("Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=17),
        axis.text.x=element_text(size = 21),
        axis.title=element_text(size=23),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))

dyn %>% 
  # filter(replicate == 6) %>% 
  ggplot() +
  geom_histogram(aes(x = pit, y = ..density..)) +
  facet_wrap(~method, scales = "free") +
  ylab("") +
  xlab("PIT") +
  theme_bw() +
  theme(axis.text.y=element_text(size=17),
        axis.text.x=element_text(size = 15),
        axis.title=element_text(size=23),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 19),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))


dyn %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = mlogs))


dyn %>% 
  # filter(time > 1) %>% 
  group_by(method, time) %>% 
  summarise(mcrpss = mean(mcrps),
            mlogss = mean(mlogs),
            selogs = sd(mlogs)/500) %>% 
  ggplot() + 
  geom_line(aes(x = time, y = mlogss, 
                colour = method, linetype = method), size = 1.1) +
  ylab("LogS") +
  xlab("t") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 19),
        axis.title=element_text(size=23),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = c(.91,.78))



dyn %>% 
  # filter(time > 1) %>% 
  group_by(method, time) %>% 
  summarise(mcrpss = mean(mcrps),
            mlogss = mean(mlogs),
            selogs = sd(mlogs)/500,
            secrps = sd(mcrps)/500) %>% 
  ggplot() + 
  geom_line(aes(x = time, y = mcrpss, 
                colour = method, linetype = method), size = 1.1) +
  ylab("CRPS") +
  xlab("t") +
  labs(colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 19),
        axis.title=element_text(size=23),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.position = "none")


dyn %>% 
  group_by(time, replicate) %>% 
  mutate(mincrps = min(mlogs)) %>% 
  mutate(rel_crps = mlogs/mincrps) %>% 
  # group_by(method) %>% 
  # summarise(m = mean(rel_crps))
  ggplot() +
  geom_boxplot(aes(x = method, y = rel_crps)) #+
  facet_wrap(~method)



library(kableExtra)
dyn_sum <- dyn %>% 
  group_by(method) %>% 
  summarise(mcrps = mean(mcrps),
            mlogs = mean(mlogs),
            uwd1 = unit_wass_dist(ecdf(pit)))

kbl(dyn_sum)









