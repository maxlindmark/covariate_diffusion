library(ggsidekick)
library(ggtext)
library(viridis)
library(rnaturalearth)
library(ggplot2)

map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", continent = "North America")

pnw <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -125, ymin = 25, xmax = -90, ymax = 55))))

map_plot_bb <- ggplot(pnw) +
  geom_sf(linewidth = 0.4, color = "gray40") +
  theme_sleek() +
  coord_sf(expand = FALSE) +
  xlim(-125, -101.5) +
  ylim(30.5, 49.9) +
  labs(x = "Longitude", y = "Latitude")

mp_bb_s <- map_plot_bb +
  guides(fill = guide_colorbar(position = "inside", title.position = "top", title.hjust = 0.5)) +
  theme(legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.direction = "horizontal",
        legend.key.spacing = unit(0, "mm"),
        strip.text = element_markdown(),
        legend.position.inside = c(0.16, 0.019)) +
  facet_wrap(~abb_name, ncol = 1) +
  scale_fill_viridis(option = "mako") +
  NULL

mp_bb_fc <- map_plot_bb +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  theme(
    legend.direction = "horizontal",
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.position = "bottom",
    legend.key.width = unit(0.7, "cm"),
    legend.key.height = unit(0.4, "cm"),
    axis.text.x = element_text(angle = 90),
    strip.text = element_markdown()
  ) +
  facet_wrap(~abb_name) +
  scale_fill_gradient2() +
  NULL
