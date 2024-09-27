library(ggsidekick)
library(ggtext)
library(viridis)
library(rnaturalearth)
library(ggplot2)

map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", continent = "North America")

alaska <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -179, ymin = 50, xmax = -150, ymax = 70))))

map_plot_ebs <- ggplot(alaska) +
  geom_sf(linewidth = 0.4, color = "gray40") +
  theme_sleek() +
  coord_sf(expand = FALSE) +
  xlim(-178, -155) +
  ylim(53, 66) +
  labs(x = "Longitude", y = "Latitude")

mp_ebs_s <- map_plot_ebs +
  guides(fill = guide_colorbar(position = "inside", title.position = "top", title.hjust = 0.5),
         color = "none") +
  theme(legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.direction = "horizontal",
        strip.text = element_markdown(),
        legend.position.inside = c(0.17, 0.03)) +
  facet_wrap(~species, ncol = 1) +
  NULL

mp_ebs_fc <- map_plot_ebs +
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
  facet_wrap(~species) + # FIXME: use the abb_name once I get the latin names
  scale_fill_gradient2() +
  NULL

# FIXME: In the EBS data, I think because the grid is so high-res, I see cell boundaries from geom_sf, even with color = null.
# The only way to remove this is to set color also to a variable (the fill). I wonder if this distorts the pattern, because the line has a width...
