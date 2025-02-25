library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(distr)

iid <- readRDS("./simulations/iid_stacks_res.rds")

# iid %>% 
#   # filter(N > 20) %>% 
#   filter(method != "bma") %>% 
#   group_by(method, N) %>% 
#   summarise(mmcrps = median(mcrps),
#             muwd1 = median(uwd1)) %>% 
#   ggplot() +
#   geom_line(aes(x = N, y = muwd1, colour = method))


crpsp <- iid %>% 
  # filter(mcrps < 1.5) %>%
  mutate(method = ifelse(method == "bma", "BMA",
                         ifelse(method == "avs", "AVS",
                                ifelse(method == "msgp", "SGP", 
                                       ifelse(method == "eqw", 
                                              "EQW", NA))))) %>% 
  mutate(method = factor(method, levels = c("BMA", "AVS", "SGP", "EQW"))) %>% 
  # filter(N %in% c(10, 50, 200, 400)) %>% 
  mutate(N = factor(N)) %>% 
  ggplot() +
  geom_boxplot(aes(x = N, y = mcrps, colour = method), size = .8) +
  scale_colour_manual(name = "Model",
                      labels = c("BMA", "AVS", "EQW", "SGP")
                      ,values = c("grey60", "grey80", "grey40", 
                                  "grey20")) +
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
        legend.position = "none")


logsp <- iid %>% 
  # filter(mcrps < 1.5) %>%
  mutate(method = ifelse(method == "bma", "BMA",
                         ifelse(method == "avs", "AVS",
                                ifelse(method == "msgp", "SGP", 
                                       ifelse(method == "eqw", 
                                              "EQW", NA))))) %>% 
  mutate(method = factor(method, levels = c("BMA", "AVS", "SGP", "EQW"))) %>% 
  # filter(N %in% c(10, 50, 200, 400)) %>% 
  mutate(N = factor(N)) %>% 
  ggplot() +
  geom_boxplot(aes(x = N, y = mlogs, colour = method), size = .8) +
  scale_colour_manual(name = "Model",
                      labels = c("BMA", "AVS", "EQW", "SGP")
                      ,values = c("grey60", "grey80", "grey40", 
                                  "grey20")) +
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



cowplot::plot_grid(logsp, crpsp, ncol = 2)



mix <- UnivarMixingDistribution(Norm(3, 1),
                         Norm(6.5, 1),
                         mixCoeff = c(.65, .35))


tmus <- c(3, 6.5)
tsigmas <- c(1, 1)
tws <- c(.65, .35)
dmixnorm <- function(x, mus, sigmas, ws) {
  y <- 0
  for (i in 1:length(mus)) {
    y <- y + dnorm(x, mus[i] + sigmas[i])*ws[i]
  }
  return(y)
}  

x <- seq(-1, 12, length.out = 1001)
data.frame(x) %>% 
  mutate

xs <- seq(-3, 14, length.out = 1001)
ys <- apply(matrix(mus, nrow = length(mus)), MARGIN = 1, FUN = dnorm, x = xs)
colnames(ys) <- paste0("y", 1:length(mus))

true_comps <- data.frame(x = xs, ys) %>% 
  mutate(y_true = dmixnorm(x, tmus, tsigmas, tws)) %>% 
  ggplot() +
  geom_line(aes(x = x, y = y1), colour = "grey80", size = 1.1) + 
  geom_line(aes(x = x, y = y2), colour = "grey80", size = 1.1) +
  geom_line(aes(x = x, y = y3), colour = "grey80", size = 1.1) +
  geom_line(aes(x = x, y = y4), colour = "grey80", size = 1.1) +
  geom_line(aes(x = x, y = y5), colour = "grey80", size = 1.1) +
  geom_line(aes(x = x, y = y6), colour = "grey80", size = 1.1) +
  geom_line(aes(x = x, y = y_true), size = 1.3) +
  ylab("p(x)") +
  xlab("x") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12),
        axis.title=element_text(size=21),
        strip.text.y = element_text(size = 14,),
        strip.text.x = element_text(size = 12),
        legend.position = "none")


comps <- iid %>% 
  dplyr::select(contains("comp"))


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
    pivot_longer(1:4, names_to = "method", values_to = "y") %>% 
    mutate(method = str_replace_all(method, "[:digit:]", "")) %>% 
    arrange(method)
  
  comp_lines$x <- rep(x, 4)
  comp_lines$N <- Ns[i]
  all_comp_lines <- rbind(all_comp_lines, comp_lines)
}

true_dist <- data.frame(x = x, y = dmixnorm(x, tmus, tsigmas, tws))
fit_lines <- all_comp_lines %>% 
  mutate(method = ifelse(method == "wavs", "AVS", 
                         ifelse(method == "wmsgp", "SGP",
                                ifelse(method == "wpmp", "BMA", 
                                       ifelse(method == "weq", "EQW", 
                                              NA))))) %>% 
  mutate(method = factor(method, levels = c("BMA", "AVS", "EQW", "SGP"))) %>% 
  ggplot() +
  geom_line(data = true_dist, aes(x = x, y = y), 
            colour = "grey", size = .8) +
  geom_line(aes(x = x, y = y), size = .8) + 
  facet_grid(method~N, scales = "free") +
  xlab("x") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12),
        axis.title=element_text(size=21),
        strip.text.y = element_text(size = 14,),
        strip.text.x = element_text(size = 12),
        legend.position = "none")


cowplot::plot_grid(true_comps, fit_lines, ncol = 2, 
                   rel_widths = c(.6,1),
                   rel_heights = c(.6,1))





