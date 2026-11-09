# Redraws fig15_wdbc_subsets from the saved set-level results (no rerun of the simulation).
# Run from the folder that contains real_data_wdbc_subsets/ and figures/.
library(ggplot2)
res <- read.csv("real_data_wdbc_subsets/set_level_results.csv")
TAU <- 0.5; NTR <- c(25, 50, 100)
f <- rbind(data.frame(p = res$p, n = res$n, comparison = "Common covariance\n(PE - LDA)", diff = res$d_PE_LDA),
           data.frame(p = res$p, n = res$n, comparison = "Class-specific\n(PE - QDA)", diff = res$d_PE_QDA))
f$p <- factor(paste("p =", f$p)); f$n <- factor(paste("n =", f$n), levels = paste("n =", NTR))
g <- ggplot(f, aes(diff, comparison)) +
  geom_vline(xintercept = 0, colour = "#52514e", linewidth = 0.3) +
  geom_vline(xintercept = c(-TAU, TAU), colour = "#52514e", linewidth = 0.3, linetype = "dashed") +
  geom_point(colour = "#2a78d6", alpha = 0.6, size = 1.3,
             position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  facet_grid(p ~ n) +
  labs(x = "Difference in test error (pp; negative: PE better)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.spacing.x = unit(1.4, "lines"),
        strip.text = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11), axis.title.x = element_text(size = 12))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig15_wdbc_subsets.pdf", g, width = 8, height = 4.8)
ggsave("figures/fig15_wdbc_subsets.png", g, width = 8, height = 4.8, dpi = 300)
cat("fig15 redrawn\n")
