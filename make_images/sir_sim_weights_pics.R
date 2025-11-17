
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)


sim_res <- read.csv("./sir_wts/seq_1.csv")
sim_res <- readRDS("./sir_wts/all_wt_res.rds")


m_wts <- sim_res %>% 
  group_by(model, week, alpha) %>% 
  summarise(
    m = mean(wt),
    upp = quantile(wt, probs = 0.975),
    low = quantile(wt, probs = 0.025)
  )


m_wts %>% 
  filter(week < 31) %>%
  mutate(SGP = ifelse(alpha == 1, "SGP", "SGP50")) %>%
  # filter(SGP == "SGP50") %>% 
  ggplot() +
  geom_segment(aes(x = low, xend = upp,
                   y = as.numeric(factor(model)) + 
                     ifelse(SGP == "SGP", 0.1, -0.1), 
                   colour = SGP), size = 1.7) +
  scale_colour_manual(name = "Method",
                      labels = c("SGP", "SGP50")
                      ,values = c("#E69F00", "#56B4E9")) +
  xlab("") +
  facet_wrap(~week) +
  theme_bw() +
  theme(legend.position = c(.925, .13),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15))


m_wts %>% 
  mutate(model = ifelse(model == "sir", "MB", 
                        ifelse(model == "asg", "HB",
                               ifelse(model == "arima", "AM", 
                                      "RW")))) %>% 
  filter(week < 31) %>%
  mutate(SGP = ifelse(alpha == 1, "SGP", "SGP50")) %>%
  ggplot() +
  geom_errorbarh(aes(y = model, xmin = low, 
                    xmax = upp, color = SGP),
                 size = 1.6,
                 height = 0,
                position = position_dodge(width = .6)) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0.25, 0.75)
  ) +
  scale_colour_manual(name = "SGP",
                      labels = c("SGP", "SGP50")
                      ,values = c("#E69F00", "#56B4E9")) +
  ylab("") +
  xlab("") +
  facet_wrap(~(week - 1)) +
  theme_bw() +
  theme(legend.position = c(.915, .07),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 12))



sim_res %>% 
  mutate(model = ifelse(model == "sir", "MB", 
                        ifelse(model == "asg", "HB",
                               ifelse(model == "arima", "AM", 
                                      "RW")))) %>% 
  mutate(SGP = ifelse(alpha == 1, "SGP", "SGP50")) %>%
  # filter(alpha == 50) %>% 
  # filter(week == 30) %>% 
  filter(week %in% c(5, 15, 25, 35)) %>% 
  ggplot() +
  geom_histogram(aes(x = wt, fill = SGP)) +
  scale_fill_manual(name = "SGP",
                      labels = c("SGP", "SGP50")
                      ,values = c("#E69F00", "#56B4E9")) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0.25, 0.75)
  ) +
  geom_vline(aes(xintercept = 0.25)) +
  facet_grid(model~(week - 1), scales = "free_y") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(
        # legend.position = c(.915, .07),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 14))



srep <- sample(60, 1)

#let srep = 24
sim_res %>% 
  filter(step == 4, rep == srep) %>% 
  mutate(model = ifelse(model == "sir", "MB", 
                        ifelse(model == "asg", "HB",
                               ifelse(model == "arima", "AM", 
                                      "RW")))) %>% 
  filter(week < 31) %>%
  mutate(SGP = ifelse(alpha == 1, "SGP", "SGP50")) %>%
  ggplot() +
  geom_errorbarh(aes(y = model, xmin = low, 
                     xmax = upp, color = SGP),
                 size = 1.6,
                 height = 0,
                 position = position_dodge(width = .6)) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0.25, 0.75)
  ) +
  scale_colour_manual(name = "SGP",
                      labels = c("SGP", "SGP50")
                      ,values = c("#E69F00", "#56B4E9")) +
  ylab("") +
  xlab("") +
  facet_wrap(~(week - 1)) +
  theme_bw() +
  theme(legend.position = c(.915, .07),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 11))











