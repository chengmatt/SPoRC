# Purpose: To create a map of sablefish management areas used in the spatial collapsibility script
# Creator: Matthew LH. Cheng
# UAF-CFOS


# Setup -------------------------------------------------------------------

library(here)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyverse)

# Read in maps here
west <- ne_states(c("United States of America", "Russia", "Canada"), returnclass = "sf")
west <- st_shift_longitude(west) # shift ongitude for plotting

# Read in stat areas
nmfs_areas <- read_sf(dsn = here("dev", "sporc_manuscript_demonstrations", "NMFS_Stat_Areas", "Sablefish_Longline_Area"), layer = "Sablefish_Longline_Area")
nmfs_areas <- nmfs_areas %>% mutate(GEN_NAME = ifelse(NAME %in% c("East Yakutat / Southeast Alaska", "West Yakutat"), "Eastern Gulf of Alaska", "A")) %>%
  mutate(         NAME = case_when(
    NAME == "Aleutian Islands" ~ "AI",
    NAME == "Bering Sea" ~ "BS",
    NAME == "Western Gulf of Alaska" ~ "WGOA",
    NAME == "Central Gulf of Alaska" ~ "CGOA",
    NAME == "West Yakutat" ~ "EGOA",
    NAME == "East Yakutat / Southeast Alaska" ~ "EGOA"
  ), NAME = factor(NAME, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  group_by(NAME) %>%
  summarise(geometry = st_union(geometry))

# Coerce longline areas
nmfs_areas <- st_make_valid(nmfs_areas) # make valid so that vertices aren't duplicated
nmfs_areas <- nmfs_areas %>% st_transform(4326) # transform to crs 4326
nmfs_areas <- st_shift_longitude(nmfs_areas) # shift longitude for plotting

# get centroids of the geometry for plotting
centroids <- nmfs_areas %>%
  group_by(NAME) %>%
  summarise(geometry = st_centroid(geometry)) %>%
  ungroup()

colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"
)

# make map
map_plot <- ggplot() +
  geom_sf(data = nmfs_areas, alpha = 0.55, aes(fill = NAME)) +
  geom_sf(data = west, color = "black", linewidth = 0.1) +
  coord_sf(ylim = c(45.2, 72.5), xlim = c(165, 230), expand = FALSE) +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.925, 0.125),
        legend.background = element_blank()) +
  labs(
    x = expression(Longitude~(degree)),
    y = expression(Latitude~(degree)),
    fill = "Region"
  )

ggsave(
  here("dev", "sporc_manuscript_demonstrations", "figs", "map.png"),
  map_plot, width = 13, height = 10,
)

