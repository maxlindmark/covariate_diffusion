library(ggplot2); theme_set(theme_light())
library(dplyr)
library(ggtext)
library(patchwork)

# Read data
res <- read.csv(paste0(here::here(), "/results/2024-04-04_identity_LNP/Results_cz.csv")) %>%
  mutate(deltaAIC_sign = ifelse(deltaAIC > 0,
                                "Without diffusion favored",
                                "Diffusion favored"),
         deltaAIC_sign_sig = ifelse(abs(deltaAIC) > 2,
                                    "Clear",
                                    "Not clear")) %>%
  mutate(X = ifelse(X == "a_yfs",
                    paste0("<span style=\"color: ", "tomato3", "\">", X, "</span>"),
                    X))

p1 <- ggplot(res, aes(corr_depth, reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(x = "Correlation between depth\nand diffused depth",
       y = "Species")

p2 <- ggplot(res, aes(deltaAIC, reorder(X, desc(deltaAIC)), color = deltaAIC_sign, shape = deltaAIC_sign_sig)) +
  geom_point(fill = NA) +
  labs(x = "ΔAIC", y = "Species") +
  scale_shape_manual(values = c(19, 21), name = "ΔAIC difference") +
  scale_color_manual(values = c("steelblue", "tomato3"), name = "ΔAIC sign")

# Quite some range in ranges :)
p3 <- ggplot(res, aes(log(range), reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(y = "Species")

p1 + p2 + p3 + plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom",
        axis.text.y = ggtext::element_markdown())
