library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

iid <- readRDS("./simulations/iid_stacks_res.rds")

# iid %>% 
#   # filter(N > 20) %>% 
#   filter(method != "bma") %>% 
#   group_by(method, N) %>% 
#   summarise(mmcrps = median(mcrps),
#             muwd1 = median(uwd1)) %>% 
#   ggplot() +
#   geom_line(aes(x = N, y = muwd1, colour = method))


iid %>% 
  # filter(mcrps < 1.5) %>%
  mutate(method = ifelse(method == "bma", "BMA",
                         ifelse(method == "avs", "AVS",
                                ifelse(method == "msgp", "SGP", NA)))) %>% 
  mutate(method = factor(method, levels = c("BMA", "AVS", "SGP"))) %>% 
  # filter(N %in% c(10, 50, 200, 400)) %>% 
  mutate(N = factor(N)) %>% 
  ggplot() +
  geom_boxplot(aes(x = N, y = mcrps, colour = method), size = .8) +
  xlab("n") +
  ylab("CRPS") +
  labs(colour = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 19),
        axis.title=element_text(size=23),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = c(.81,.78))


iid %>% 
  # filter(mcrps < 1.5) %>%
  mutate(method = ifelse(method == "bma", "BMA",
                         ifelse(method == "avs", "AVS",
                                ifelse(method == "msgp", "SGP", NA)))) %>% 
  mutate(method = factor(method, levels = c("BMA", "AVS", "SGP"))) %>% 
  # filter(N %in% c(10, 50, 200, 400)) %>% 
  mutate(N = factor(N)) %>% 
  ggplot() +
  geom_boxplot(aes(x = N, y = mlogs, colour = method), size = .8) +
  xlab("n") +
  ylab("LogS") +
  labs(colour = "Method") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 19),
        axis.title=element_text(size=23),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = c(.81,.78))








  

x <- seq(-1, 12, length.out = 1001)
ys <- c()
for (i in 1:length(x)) {
  ys[i] <- sum(comps[6,]*dnorm(x[i], mus, sigmas))
}
lines(ys ~ x, col = "orange")

comps <- iid %>% 
  select(contains("comp"))


rep_samp <- sample(unique(iid$replicate), 1)

iid_rep <- iid %>% 
  filter(replicate == rep_samp)
Ns <- unique(iid_rep$N)
all_comp_lines <- data.frame()
for (i in 1:length(Ns)) {
  comp_samp <- iid_rep %>% 
    filter(N == Ns[i]) %>% 
    select(contains("comp"))
  
  comp_lines <- apply(comp_samp, MARGIN = 1, FUN = dmixnorm, x = x, 
        mus = mus, sigmas = sigmas) %>% as.data.frame() %>% 
    pivot_longer(1:3, names_to = "method", values_to = "y") %>% 
    mutate(method = str_replace_all(method, "[:digit:]", "")) %>% 
    arrange(method)
  
  comp_lines$x <- rep(x, 3)
  comp_lines$N <- Ns[i]
  all_comp_lines <- rbind(all_comp_lines, comp_lines)
}

true_dist <- data.frame(x = x, y = dmixnorm(x, tmus, tsigmas, tws))
all_comp_lines %>% 
  mutate(method = ifelse(method == "wavs", "AVS", 
                         ifelse(method == "wmsgp", "SGP",
                                ifelse(method == "wpmp", "BMA", NA)))) %>% 
  mutate(method = factor(method, levels = c("BMA", "AVS", "SGP"))) %>% 
  ggplot() +
  geom_line(data = true_dist, aes(x = x, y = y), 
            colour = "grey", size = 1.1) +
  geom_line(aes(x = x, y = y), size = 1.1) + 
  facet_grid(method~N) +
  xlab("x") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12),
        axis.title=element_text(size=21),
        strip.text.y = element_text(size = 14,),
        strip.text.x = element_text(size = 12),
        legend.position = "none")

