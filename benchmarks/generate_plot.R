#!/usr/bin/env Rscript
# Generate performance comparison bar chart in ggplot2 theme_bw style

library(ggplot2)

df <- data.frame(
    Implementation = factor(c(
        "Python m3_light\n(1 core)",
        "C++17\n(1 core)",
        "C++17\n(4 cores)",
        "C++17\n(8 cores)",
        "C++17\n(12 cores)"
    ), levels = c(
        "Python m3_light\n(1 core)",
        "C++17\n(1 core)",
        "C++17\n(4 cores)",
        "C++17\n(8 cores)",
        "C++17\n(12 cores)"
    )),
    Time_min = c(35, 17, 5, 2.5, 2),
    Speedup = c("1x", "2x", "7x", "14x", "18x"),
    Type = c("Python", "C++", "C++", "C++", "C++")
)

p <- ggplot(df, aes(x = Implementation, y = Time_min, fill = Type)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.3) +
    geom_text(aes(label = Speedup), vjust = -0.5, size = 3.5, fontface = "bold") +
    geom_text(aes(label = paste0(Time_min, " min")), vjust = 1.5, size = 3,
              color = "white", fontface = "bold") +
    scale_fill_manual(values = c("Python" = "#e05555", "C++" = "#2166ac"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                       breaks = seq(0, 35, 5)) +
    labs(
        title = "Tetramer search: m3_light (Python) vs rnamotifs_search (C++17)",
        subtitle = "NOVA dataset, mm9 genome, 512 tetramers | Intel i7-8700, 64 GB RAM, Ubuntu 24.04",
        x = NULL, y = "Wall-clock time (minutes)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40", size = 9),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9)
    )

ggsave("benchmarks/performance_comparison.pdf",
       p, width = 8, height = 4.5)
ggsave("benchmarks/performance_comparison.svg",
       p, width = 8, height = 4.5)

cat("Benchmark plot saved to benchmarks/performance_comparison.{pdf,svg}\n")
