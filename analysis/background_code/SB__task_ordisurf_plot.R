

library(ggplot2)
library(vegan)
library(grid)

data("dune")
data("dune.env")
ord <- metaMDS(dune)
dune.sf <- ordisurf(ord ~ dune.env$A1, plot = FALSE, scaling = 3)

# get out
species.scores <- as.data.frame(scores(ord, "species"))
species.scores$species <- rownames(species.scores)
names(species.scores)[c(1, 2)] <- c("x", "y")
species.scores$z <- NA

# get surf out
extract.xyz <- function(obj) {   # obj <- dune.sf
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz) }
contour.vals <- extract.xyz(obj = dune.sf)

# identify stuff
labelz <- data.frame(x = c(-0.85, -0.8, -0.45, -0.15, 0.15, 0.5, 0.85), y = c(0.05, 
                            1.1, 1.1, 1.1, 1, 0.75, 0.65), z = NA, 
                     labels = c("3.5", "4.0", "4.5", "5.0", "5.5", "6.0", "6.5"))

# plot stuff
ggplot(data = contour.vals, aes(x, y, z = z)) +
  theme_minimal() +
  stat_contour(aes(colour = ..level..)) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1.5)) +
  geom_text(data = species.scores,
            aes(x = x, y = y, label = species),
            colour = "red") + 
  coord_equal() +
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(panel.border = element_rect(fill = NA),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none") +
  geom_text(data = labelz, 
            aes(x = x, y = y, label = labels), 
            angle = -80, 
            size = 6)

# # could try
# library(directlabels)
# direct.label(p)

