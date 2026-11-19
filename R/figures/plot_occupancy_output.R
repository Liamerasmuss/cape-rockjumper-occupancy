# Model output figures for Model B (original-prior fit).
#   1. Map of the posterior probability that each pentad is occupied.
#   2. Probability of occupancy after K visits with no detection, with and without the indicator.
# Both use only psi, q0 and q1 (the identified quantities), from the saved draws.
# Run from the project root. Keep plot_core_region_map.R next to this file.

library(ggplot2)

source("R/figures/plot_core_region_map.R")     # gives raw, pent, za, half, code_to_deg, land, sea, xlim, ylim

draws <- readRDS("output/realfit_modelB_original_draws.rds")
psi <- draws$psi; q0 <- draws$q0; q1 <- draws$q1

# ---- 1. posterior occupancy probability per pentad ----

rocky  <- (raw$Spp.sert != "-") | (raw$Spp.grwo != "-")
fynbos <- (raw$Spp.obsu != "-") | (raw$Spp.casi != "-")
raw$synth <- as.integer(rocky & fynbos)

K0  <- tapply(raw$synth == 0, raw$Pentad, sum)
K1  <- tapply(raw$synth == 1, raw$Pentad, sum)
det <- tapply(raw$Spp.caro != "-", raw$Pentad, any)
ids <- names(K0)

occ_prob <- function(k0, k1, detected) {
  if (detected) return(rep(1, length(psi)))
  l1 <- (1 - q0)^k0 * (1 - q1)^k1          # probability of no detections if occupied
  psi * l1 / (psi * l1 + (1 - psi))
}

occ <- sapply(ids, function(id) occ_prob(K0[id], K1[id], det[id]))   # draws x pentads

occ_df <- data.frame(Pentad = ids,
                     mean = colMeans(occ),
                     lo = apply(occ, 2, quantile, 0.05),
                     hi = apply(occ, 2, quantile, 0.95))
occ_df$lat <- -code_to_deg(sub("_.*", "", occ_df$Pentad))
occ_df$lon <-  code_to_deg(sub(".*_", "", occ_df$Pentad))

cat("Mean of pentad occupancy probabilities:", round(mean(occ_df$mean), 3),
    " (compare psi =", round(mean(psi), 3), ")\n")

occ_df$detected <- as.logical(det[occ_df$Pentad])

pm <- ggplot() +
  geom_sf(data = za, fill = land, colour = "grey55", linewidth = 0.3) +
  geom_tile(data = occ_df[!occ_df$detected, ], aes(lon, lat, fill = mean),
            width = 2 * half, height = 2 * half, colour = "white", linewidth = 0.1) +
  geom_tile(data = occ_df[occ_df$detected, ], aes(lon, lat),
            width = 2 * half, height = 2 * half, fill = "#7A1F3D",
            colour = "white", linewidth = 0.1) +
  scale_fill_gradientn(colours = c("#F6F3EC", "#F2D49B", "#E9965C", "#C8553D"),
                       limits = c(0, 0.3), oob = scales::squish,
                       name = "Probability occupied (pentads with no records)",
                       labels = scales::percent_format(accuracy = 1)) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top", barwidth = 12, barheight = 0.5)) +
  labs(x = NULL, y = NULL,
       caption = "Dark red squares are pentads where the species was recorded.") +
  theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = sea, colour = NA),
        panel.grid = element_line(colour = scales::alpha("white", 0.5), linewidth = 0.25),
        legend.position = "bottom", legend.title = element_text(size = 9),
        plot.caption = element_text(size = 8, colour = "grey30", hjust = 0))

dir.create("figs", showWarnings = FALSE)
ggsave("figs/occupancy-map.png", pm, width = 5.2, height = 7.4, dpi = 300)

# ---- 2. occupancy probability after K visits with no detection ----

K <- 0:60
curve <- function(ind_present) {
  t(sapply(K, function(k) {
    l1 <- if (ind_present) (1 - q1)^k else (1 - q0)^k
    pr <- psi * l1 / (psi * l1 + (1 - psi))
    c(mean = mean(pr), lo = quantile(pr, 0.05), hi = quantile(pr, 0.95))
  }))
}
c0 <- as.data.frame(curve(FALSE)); c0$K <- K; c0$group <- "Indicator absent on every visit"
c1 <- as.data.frame(curve(TRUE));  c1$K <- K; c1$group <- "Indicator present on every visit"
names(c0)[2:3] <- names(c1)[2:3] <- c("lo", "hi")
cv <- rbind(c0, c1)

pc <- ggplot(cv, aes(K, mean, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = mean(psi), linetype = 3, colour = "grey40") +
  annotate("text", x = 59, y = mean(psi) + 0.025, label = "overall occupancy", hjust = 1,
           size = 3, colour = "grey30") +
  scale_colour_manual(values = c("#C8553D", "#3B7EA1"), name = NULL) +
  scale_fill_manual(values = c("#C8553D", "#3B7EA1"), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.3)) +
  labs(x = "Number of visits with no Cape Rock-jumper recorded",
       y = "Probability the pentad is occupied") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave("figs/occupancy-vs-visits.png", pc, width = 6.5, height = 4.2, dpi = 300)