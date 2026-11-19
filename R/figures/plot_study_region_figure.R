# Two-panel study-region figure.
#   A  how the region was chosen (all Western Cape pentads, bounding box, retained pentads)
#   B  Cape Rock-jumper reporting rate in the retained pentads (from plot_core_region_map.R)
# Run from the project root.
#
# install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "patchwork", "ABAP"))

library(ggplot2)
library(patchwork)

# Builds panel B (object `p`), plus `pent`, `half`, `za`, `code_to_deg`, `land`, `sea`.
source("R/figures/plot_core_region_map.R")

# ---- all Western Cape pentads ----
# EDIT THIS BLOCK to load the pentads the same way you did for the old Figure 3.1.
# The only thing needed is a data frame with a column of pentad codes like "3410_1900".
all_pent <- tryCatch(
  ABAP::getRegionPentads(.region_type = "province", .region = "Western Cape"),
  error = function(e) stop("Could not load all pentads with ABAP. Load them the way the old ",
                           "Figure 3.1 script did and call the code column `Pentad`.")
)
names(all_pent)[grep("pentad", names(all_pent), ignore.case = TRUE)[1]] <- "Pentad"

all_pent$lat <- -code_to_deg(sub("_.*", "", all_pent$Pentad))
all_pent$lon <-  code_to_deg(sub(".*_", "", all_pent$Pentad))
all_pent$kept <- all_pent$Pentad %in% pent$Pentad

# bounding box implied by the retained pentads
box <- data.frame(xmin = min(pent$lon) - half, xmax = max(pent$lon) + half,
                  ymin = min(pent$lat) - half, ymax = max(pent$lat) + half)

pA <- ggplot() +
  geom_sf(data = za, fill = land, colour = "grey55", linewidth = 0.3) +
  geom_tile(data = all_pent[!all_pent$kept, ], aes(lon, lat),
            width = 2 * half, height = 2 * half, fill = "#DDD8CB",
            colour = "white", linewidth = 0.05) +
  geom_tile(data = all_pent[all_pent$kept, ], aes(lon, lat),
            width = 2 * half, height = 2 * half, fill = "#3B7EA1",
            colour = "white", linewidth = 0.05) +
  geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "firebrick", linewidth = 0.6) +
  coord_sf(xlim = range(all_pent$lon) + c(-0.4, 0.4),
           ylim = range(all_pent$lat) + c(-0.4, 0.4), expand = FALSE) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = sea, colour = NA),
        panel.grid = element_line(colour = scales::alpha("white", 0.5), linewidth = 0.25))

pB <- p + labs(tag = "B") +
  guides(fill = guide_colourbar(title.position = "top", barwidth = 8, barheight = 0.5)) +
  theme(legend.title = element_text(size = 8))

combined <- pA + pB + plot_layout(widths = c(1.6, 1))

ggsave("figs/study-region.png", combined, width = 10, height = 6.4, dpi = 300)
ggsave("figs/study-region-selection.png", pA, width = 6.5, height = 5, dpi = 300)
