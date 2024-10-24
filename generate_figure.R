library(ggplot2)
library(grid)

png(filename = "/Users/adina/Desktop/figure.png", width = 6, height = 5, units = "in", res = 600, type = "cairo")

data <- read.csv('/Users/adina/Desktop/input.csv')

data$migrant <- factor(data$migrant, levels = unique(data$migrant))
data$type <- factor(data$type, levels = c("add", "dim_epi", "rein_epi", "pla"))
group_labels <- data.frame(
  group = 1:9,
  label = c(" A. Diminishing Epistasis", "", "", "B. Reinforcing Epistasis", "", "", "C. Phenotypic Plasticity", "", "")
)

ggplot(data) +
  geom_bar(stat = "identity", 
           aes(x = migrant, y = HL, fill = type, color = type),
           position = position_dodge(width = 0.7), 
           width = 0.6) +
  facet_wrap(~ group, nrow = 3, scales = "free_x") +
  scale_color_manual(values = c("add" = "blue", "dim_epi" = "orange", "rein_epi" = "indianred", "pla" = "red"), guide = "none") +
  scale_fill_manual(values = c("add" = "white", "dim_epi" = "orange", "rein_epi" = "indianred", "pla" = "red"), labels = c("Additivity (No Epistasis)", "Diminishing Epistasis", "Reinforcing Epistasis", "Plasticity")) +
  guides(fill = guide_legend(override.aes = list(fill = c("add" = "white", "dim_epi" = "orange", "rein_epi" = "indianred", "pla" = "red"), color = c("add" = "blue", "dim_epi" = "orange", "rein_epi" = "indianred", "pla" = "red")))) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    labels = c("0.1", "1", "10", "100", "1000", "10000")
  ) +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 0, size = 9), 
    axis.text.y = element_text(size = 9), 
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9, margin = margin(r = 8)),
    panel.grid = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.35, "cm"),
    legend.position = c(0.48, 1.12),
    legend.direction = "horizontal",
    plot.margin = margin(1.5, 0.3, 0.3, 0.3, "cm"), #top, right, bottom, left
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "\nMigrants",
    y = expression(italic("H")[L] ~ "Model / Drift"),
    fill = ""
  ) +
  geom_text(data = group_labels, aes(label = label),
            x = Inf, y = Inf, hjust = 1.12, vjust = 1, size = 3,
            fontface = "bold")

dev.off()