# Map of the core region: Cape Rock-jumper reporting rate per pentad.
# Pentad positions come from the pentad code (DDMM_DDMM = latitude S, longitude E),
# so no coordinates file is needed. Run from the project root.
#
# install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata"))

library(ggplot2)
library(sf)
library(rnaturalearth)

raw <- readRDS("data/clean/cleaned_multispp_occ_data.RDS")   # adjust path if different

# ---- 1. pentad summaries ----

pent <- aggregate(cbind(det = raw$Spp.caro != "-") ~ Pentad, data = raw,
                  FUN = function(x) c(n = length(x), k = sum(x)))
pent <- data.frame(Pentad = pent$Pentad, n = pent$det[, "n"], k = pent$det[, "k"])
pent$rate <- pent$k / pent$n

# ---- 2. pentad code -> centre coordinates ----
# The code gives the north-west corner. The centre is 2.5 arc-minutes south and east.
# (A shift of this size is invisible at the scale of the map.)

code_to_deg <- function(x) as.numeric(substr(x, 1, 2)) + (as.numeric(substr(x, 3, 4)) + 2.5) / 60
pent$lat <- -code_to_deg(sub("_.*", "", pent$Pentad))
pent$lon <-  code_to_deg(sub(".*_", "", pent$Pentad))
half <- 2.5 / 60                                               # half a pentad in degrees

pent$rate_plot <- ifelse(pent$k == 0, NA, pent$rate)

# ---- 3. coastline and range labels ----

za <- ne_countries(scale = "medium", country = "South Africa", returnclass = "sf")

ranges <- data.frame(
  name  = c("Cederberg", "Paarl and\nStellenbosch", "Langeberg", "Swartberg"),
  lon   = c(19.85, 18.52, 21.2, 22.4),     # where the label sits
  lat   = c(-32.00, -33.70, -34.2, -33.5),
  tlon  = c(19.18, 18.95, 20.9, 22.2),     # the range it points at (approximate)
  tlat  = c(-32.45, -33.85, -33.95, -33.35),
  sx    = c(-0.12, 0.30, 0, 0),            # where the line leaves the label
  sy    = c(-0.08, -0.12, 0, 0)
)   # approximate positions, adjust to taste

xlim <- range(pent$lon) + c(-0.6, 0.35)
ylim <- range(pent$lat) + c(-0.3, 0.3)

# keep only labels that fall inside the plotted area
ranges <- ranges[ranges$lon > xlim[1] & ranges$lon < xlim[2] &
                   ranges$lat > ylim[1] & ranges$lat < ylim[2], ]

# ---- 4. plot ----

sea  <- "#D9E6EE"
land <- "#F6F3EC"
none <- "#E6E0D2"     # pentads with no Cape Rock-jumper records

p <- ggplot() +
  geom_sf(data = za, fill = land, colour = "grey55", linewidth = 0.3) +
  geom_tile(data = pent[is.na(pent$rate_plot), ], aes(lon, lat),
            width = 2 * half, height = 2 * half, fill = none,
            colour = "white", linewidth = 0.1) +
  geom_tile(data = pent[!is.na(pent$rate_plot), ], aes(lon, lat, fill = rate_plot),
            width = 2 * half, height = 2 * half, colour = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.05, end = 0.85,
                       name = "Share of checklists with the Cape Rock-jumper",
                       labels = scales::percent_format(accuracy = 1),
                       trans = "sqrt", breaks = c(0.05, 0.25, 0.5, 0.75)) +
  geom_segment(data = ranges, aes(x = lon + sx, y = lat + sy, xend = tlon, yend = tlat),
               colour = "grey25", linewidth = 0.3) +
  geom_text(data = ranges, aes(lon, lat, label = name), size = 3.2,
            fontface = "italic", colour = "grey10", lineheight = 0.9) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top", barwidth = 12, barheight = 0.5)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = sea, colour = NA),
        panel.grid = element_line(colour = scales::alpha("white", 0.5), linewidth = 0.25),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        plot.margin = margin(4, 8, 4, 4))

# scale bar if the package is installed
if (requireNamespace("ggspatial", quietly = TRUE)) {
  p <- p + ggspatial::annotation_scale(location = "br", width_hint = 0.25,
                                       bar_cols = c("grey20", "white"))
}

dir.create("figs", showWarnings = FALSE)
ggsave("figs/core-region-map.png", p, width = 5.2, height = 7.4, dpi = 300)