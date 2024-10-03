fit_sum <- fits_list[[k]]$summary()
forcs <- unique(comp_forcs$model[comp_forcs$location == "US"])

all_fits <- data.frame()
for (i in 1:length(fits_list)) {
  fitsi <- fits_list[[i]]$summary(variables = "omegas") %>% 
    select(variable, q5, q95) 
  fitsi$model <- forcs
  fitsi$date <- sub_dates[i]
  
  fitsi <- fitsi %>% 
    pivot_longer(2:3, values_to = "quantile")
  
  all_fits <- rbind(all_fits, fitsi)
}   


formatter <- function(...){
  function(x) round(x, 1)#format(x, ..., scientific = T, digit = 2)
}

all_fits %>% 
  ggplot() +
  geom_line(aes(x = quantile, y = model), size = 1) +
  facet_wrap(~date) +
  scale_x_continuous(labels = formatter(nsmall = 2)) +
  ylab("") +
  xlab("Weight") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 17),
        strip.text = element_text(size = 12))


fits_list[[29]]$summary(variables = "omegas") %>% 
  select(variable, q5, q95) %>% 
  pivot_longer(2:3, values_to = "quantile") %>% 
  ggplot() +
  geom_line(aes(x = quantile, y = variable), size = 3) +
  theme_bw()

fits_list[[3]]$draws(variables = "omegas", format = "df") %>% 
  pivot_longer(1:11, names_to = "comp", values_to = "weight") %>% 
  select(comp, weight) %>% 
  ggplot() +
  geom_histogram(aes(x = weight)) + 
  facet_wrap(~comp) +
  theme_bw()