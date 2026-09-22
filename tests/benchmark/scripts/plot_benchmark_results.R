#!/usr/bin/env Rscript
# ============================================================
# plot_benchmark_results.R
#
# Regenerates every benchmark figure from the three committed
# workbooks under tests/benchmark/results/.
#
# Project layout assumed:
#   vcf2fasta-rust/
#   ├── docs/
#   │   └── figures/                       ← OUTPUT: PDF + PNG
#   └── tests/benchmark/
#       ├── scripts/
#       │   └── plot_benchmark_results.R   ← THIS SCRIPT
#       └── results/                       ← INPUT: three .xlsx
#           ├── results_benchmark.xlsx
#           ├── results_scaling.xlsx
#           └── results_unphased.xlsx
#
# Usage (from anywhere):
#   Rscript path/to/tests/benchmark/scripts/plot_benchmark_results.R
#
# Override default paths:
#   RESULTS_DIR=/some/other/results \
#   FIGURES_DIR=/tmp/figures \
#     Rscript plot_benchmark_results.R
#
# Required packages (install once):
#   install.packages(c("readxl","dplyr","tidyr","stringr",
#                      "ggplot2","patchwork","scales","forcats"))
# ============================================================


# ============================================================
# 0) Resolve paths relative to this script
# ============================================================
# We want the script to work whether it is invoked as
#   `Rscript tests/benchmark/scripts/plot_benchmark_results.R`
# or sourced from an interactive R session, or run from a
# different working directory.

get_script_dir <- function() {
  # 1. Rscript invocation: --file=<path>
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  matches <- grep(file_arg, args)
  if (length(matches) > 0) {
    path <- sub(file_arg, "", args[matches[1]])
    return(dirname(normalizePath(path, mustWork = FALSE)))
  }
  # 2. source() invocation: sys.frame(1)$ofile
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile, mustWork = FALSE)))
  }
  # 3. Fallback — assume the caller is in scripts/
  warning("Could not determine script directory; falling back to getwd().")
  getwd()
}

script_dir   <- get_script_dir()
# script_dir = <project>/tests/benchmark/scripts
# project_root = <project>
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."),
                              mustWork = FALSE)

# Defaults, in project-relative terms:
#   input  -> tests/benchmark/results/
#   output -> docs/figures/
default_results_dir <- normalizePath(file.path(script_dir, "..", "..", "results"),
                                     mustWork = FALSE)
default_figures_dir <- normalizePath(file.path(project_root, "docs", "figures"),
                                     mustWork = FALSE)

results_dir <- normalizePath(
  Sys.getenv("RESULTS_DIR", unset = default_results_dir),
  mustWork = FALSE
)
figures_dir <- normalizePath(
  Sys.getenv("FIGURES_DIR", unset = default_figures_dir),
  mustWork = FALSE
)

if (!dir.exists(results_dir)) {
  stop("RESULTS_DIR does not exist: ", results_dir)
}
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}

message("Script directory:  ", script_dir)
message("Project root:      ", project_root)
message("Results directory: ", results_dir)
message("Figures directory: ", figures_dir)


# ============================================================
# 1) Packages
# ============================================================
required_pkgs <- c("readxl","dplyr","tidyr","stringr",
                   "ggplot2","patchwork","scales","forcats")
missing_pkgs <- required_pkgs[!vapply(required_pkgs,
                                      requireNamespace,
                                      logical(1),
                                      quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them with:\n",
    '  install.packages(c(',
    paste0('"', missing_pkgs, '"', collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(forcats)
})


# ============================================================
# 2) Global theme + colors
# ============================================================
theme_set(
  theme_classic(base_size = 10) +
    theme(
      text             = element_text(family = "sans"),
      axis.title       = element_text(size = 10),
      axis.text        = element_text(size = 9),
      legend.text      = element_text(size = 9),
      legend.title     = element_blank(),
      legend.position  = "top",
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey30"),
      axis.line        = element_line(linewidth = 0.4),
      axis.ticks       = element_line(linewidth = 0.4),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = 9)
    )
)

# Same colors used in every figure
method_colors <- c(
  "vcflib" = "#888888",
  "CPU 1T" = "#1f77b4",
  "CPU 2T" = "#2ca02c",
  "GPU 1T" = "#d62728",
  "GPU 2T" = "#9467bd"
)


# ============================================================
# 3) Helpers: parse "MMm SS.sss" / "mean ± sd" strings
# ============================================================

# "00m 22.21s" -> 22.21  (also handles "1h 02m 03.45s" and "0.50s")
parse_one <- function(x) {
  x <- as.character(x)
  h <- as.numeric(str_match(x, "(\\d+)\\s*h")[, 2])       # optional hours
  m <- as.numeric(str_match(x, "(\\d+)\\s*m(?!s)")[, 2])   # optional minutes
  s <- as.numeric(str_match(x, "([\\d.]+)\\s*s")[, 2])     # seconds
  h[is.na(h)] <- 0
  m[is.na(m)] <- 0
  h * 3600 + m * 60 + s
}

# "00m 22.21s ± 00m 00.77s" -> tibble(mean_s, sd_s)
parse_ms <- function(x) {
  parts <- str_split_fixed(x, "±", 2)
  tibble(
    mean_s = parse_one(parts[, 1]),
    sd_s   = parse_one(parts[, 2])
  )
}

# Apply to one column, returning <prefix>_mean / <prefix>_sd
split_ms <- function(df, col, prefix) {
  out <- parse_ms(df[[col]])
  names(out) <- c(paste0(prefix, "_mean"), paste0(prefix, "_sd"))
  bind_cols(df, out)
}

# ----- Power-of-10 tick labels (with unit "s") -----
fmt_log10_s <- function(x) {
  parse(text = paste0("10^", round(log10(x)), "~s"))
}

# ----- Plain "Nx" labels for speedup -----
fmt_plain_x <- function(x) paste0(x, "x")

# Shared break positions
time_breaks  <- c(0.1, 1, 10, 100, 1000)
speed_breaks <- c(1, 10, 100, 1000)

# ----- Safe lower error bound -----
# If SD > mean, mean - sd goes negative and log-scale plots break.
# Cap at mean/5 so high-variance points stay visible without fake spikes.
safe_lower <- function(mean, sd) pmax(mean - sd, mean / 5, 0.01)

# Quick diagnostic: print min / max / ratio per figure
range_check <- function(df, col, label) {
  r <- range(df[[col]], na.rm = TRUE)
  cat(sprintf("%-22s  min=%8.3f s  max=%8.2f s  ratio=%8.1f\n",
              label, r[1], r[2], r[2] / r[1]))
}

# Save a figure as both PDF and PNG, in FIGURES_DIR
save_fig <- function(plot, name, width, height, dpi = 300) {
  pdf_path <- file.path(figures_dir, paste0(name, ".pdf"))
  png_path <- file.path(figures_dir, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height)
  ggsave(png_path, plot, width = width, height = height, dpi = dpi)
  message("Wrote ", png_path)
}

# Check that an input workbook exists, with a friendly message
require_workbook <- function(file_name) {
  path <- file.path(results_dir, file_name)
  if (!file.exists(path)) {
    stop("Required workbook not found: ", path,
         "\nExpected it under RESULTS_DIR = ", results_dir,
         "\nSee docs/BENCHMARKING.md §3 for the workbook layout.")
  }
  path
}


# ============================================================
# 4) BENCHMARK (Platinum, two samples per chromosome)
# ============================================================
bench_xlsx <- require_workbook("results_benchmark.xlsx")
bench_raw <- read_excel(bench_xlsx, sheet = "Combined")

bench <- bench_raw %>%
  split_ms("vcflib",      "vcflib") %>%
  split_ms("cpu_t1_time", "cpu_t1") %>%
  split_ms("cpu_t2_time", "cpu_t2") %>%
  split_ms("gpu_t1_time", "gpu_t1") %>%
  split_ms("gpu_t2_time", "gpu_t2") %>%
  mutate(
    Chromosome = factor(Chromosome,
                        levels = c(paste0("chr", 1:22), "chrX", "chrY")),
    Sample     = factor(Sample)
  ) %>%
  arrange(Chromosome, Sample)

glimpse(bench)

# Long format for ggplot
bench_long <- bench %>%
  select(Chromosome, Sample, Variants, Ploidy,
         vcflib_mean, vcflib_sd,
         cpu_t1_mean, cpu_t1_sd,
         cpu_t2_mean, cpu_t2_sd,
         gpu_t1_mean, gpu_t1_sd,
         gpu_t2_mean, gpu_t2_sd) %>%
  pivot_longer(
    cols = -c(Chromosome, Sample, Variants, Ploidy),
    names_to = c("method", ".value"),
    names_pattern = "(.*)_(mean|sd)"
  ) %>%
  mutate(
    method = recode(method,
                    vcflib = "vcflib",
                    cpu_t1 = "CPU 1T",
                    cpu_t2 = "CPU 2T",
                    gpu_t1 = "GPU 1T",
                    gpu_t2 = "GPU 2T"),
    method = factor(method,
                    levels = c("vcflib","CPU 1T","CPU 2T","GPU 1T","GPU 2T"))
  )

range_check(bench_long, "mean", "Fig1 (benchmark)")

# Diagnostic: how many cells have SD > mean?
bench_long %>%
  filter(!is.na(mean), !is.na(sd), sd > mean) %>%
  count(method) %>%
  print()

# ---------- Figure 1: runtime per chromosome, faceted by Sample ----------
p1 <- ggplot(bench_long,
             aes(x = Chromosome, y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75, colour = "black", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = safe_lower(mean, sd), ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2, linewidth = 0.35
  ) +
  facet_wrap(~ Sample, ncol = 1, scales = "free_y") +
  scale_y_log10(breaks = time_breaks, labels = fmt_log10_s) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "Chromosome",
    y = "Runtime (log scale)",
    fill = NULL,
    title = "Platinum — runtime across chromosomes",
    subtitle = "Bars: mean of n = 4 runs; error bars: SD"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save_fig(p1, "fig1_runtime_chromosomes", width = 13, height = 7)

# ---------- Figure 1-alt: chromosome on facet, Sample on X ----------
p1_alt <- ggplot(bench_long,
                 aes(x = Sample, y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75, colour = "black", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = safe_lower(mean, sd), ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2, linewidth = 0.35
  ) +
  facet_wrap(~ Chromosome, ncol = 6, scales = "free_y") +
  scale_y_log10(breaks = time_breaks, labels = fmt_log10_s) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = NULL,
    y = "Runtime (log scale)",
    fill = NULL,
    title = "Platinum — runtime per chromosome and sample",
    subtitle = "Bars: mean of n = 4 runs; error bars: SD"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save_fig(p1_alt, "fig1_alt", width = 13, height = 8)

# ---------- Figure 3: speedup vs vcflib, faceted by Sample ----------
bench_speedup <- bench %>%
  mutate(
    sp_cpu1 = vcflib_mean / cpu_t1_mean,
    sp_cpu2 = vcflib_mean / cpu_t2_mean,
    sp_gpu1 = vcflib_mean / gpu_t1_mean,
    sp_gpu2 = vcflib_mean / gpu_t2_mean
  ) %>%
  select(Chromosome, Sample, sp_cpu1, sp_cpu2, sp_gpu1, sp_gpu2) %>%
  pivot_longer(-c(Chromosome, Sample),
               names_to = "method", values_to = "speedup") %>%
  mutate(
    method = recode(method,
                    sp_cpu1 = "CPU 1T",
                    sp_cpu2 = "CPU 2T",
                    sp_gpu1 = "GPU 1T",
                    sp_gpu2 = "GPU 2T"),
    method = factor(method,
                    levels = c("CPU 1T","CPU 2T","GPU 1T","GPU 2T"))
  )

p3 <- ggplot(bench_speedup,
             aes(x = Chromosome, y = speedup, fill = method)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ Sample, ncol = 1, scales = "free_y") +
  scale_y_log10(breaks = speed_breaks, labels = fmt_plain_x) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "Chromosome",
    y = "Speedup vs vcflib (log scale)",
    fill = NULL,
    title = "Speedup across chromosomes",
    subtitle = "Dashed line = parity with vcflib"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save_fig(p3, "fig3_speedup", width = 13, height = 7)


# ============================================================
# 5) SCALING (1000 Genomes, samples vs runtime)
# ============================================================
scal_xlsx <- require_workbook("results_scaling.xlsx")
scal_raw <- read_excel(scal_xlsx, sheet = "Combined")

scal <- scal_raw %>%
  split_ms("cpu_t1_time",       "cpu_t1") %>%
  split_ms("cpu_t2_time",       "cpu_t2") %>%
  split_ms("gpu_dev0_t1_time",  "gpu_t1") %>%
  split_ms("gpu_dev0_t2_time",  "gpu_t2")

# Average across chromosomes at each sample count.
# Pooled SD across groups: sqrt(mean(sd^2)).
scal_summary <- scal %>%
  group_by(Samples) %>%
  summarise(
    cpu_t1_mean = mean(cpu_t1_mean),
    cpu_t2_mean = mean(cpu_t2_mean),
    gpu_t1_mean = mean(gpu_t1_mean),
    gpu_t2_mean = mean(gpu_t2_mean),
    cpu_t1_sd   = sqrt(mean(cpu_t1_sd^2)),
    cpu_t2_sd   = sqrt(mean(cpu_t2_sd^2)),
    gpu_t1_sd   = sqrt(mean(gpu_t1_sd^2)),
    gpu_t2_sd   = sqrt(mean(gpu_t2_sd^2)),
    .groups = "drop"
  ) %>%
  arrange(Samples)

scal_long <- scal_summary %>%
  pivot_longer(
    cols = -Samples,
    names_to = c("method", ".value"),
    names_pattern = "(.*)_(mean|sd)"
  ) %>%
  mutate(
    method = recode(method,
                    cpu_t1 = "CPU 1T",
                    cpu_t2 = "CPU 2T",
                    gpu_t1 = "GPU 1T",
                    gpu_t2 = "GPU 2T"),
    method = factor(method,
                    levels = c("CPU 1T","CPU 2T","GPU 1T","GPU 2T"))
  )

range_check(scal_long, "mean", "Fig2 (scaling)")

# ---------- Figure 2: scaling with sample count (log-log) ----------
p2 <- ggplot(scal_long,
             aes(x = Samples, y = mean,
                 colour = method, shape = method)) +
  geom_errorbar(aes(ymin = safe_lower(mean, sd), ymax = mean + sd),
                width = 0.02, linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(breaks = sort(unique(scal_long$Samples))) +
  scale_y_log10(breaks = time_breaks, labels = fmt_log10_s) +
  scale_colour_manual(values = method_colors) +
  scale_shape_manual(values = c(
    "CPU 1T" = 16, "CPU 2T" = 15, "GPU 1T" = 17, "GPU 2T" = 18
  )) +
  labs(
    x = "Number of samples",
    y = "Runtime",
    colour = NULL, shape = NULL,
    title = "Runtime scaling with sample count on chromosome 22",
    subtitle = "1000 Genomes Project; mean over runs, error bars: pooled SD (n = 4)"
  )

save_fig(p2, "fig2_scaling", width = 6.5, height = 5)


# ============================================================
# 6) UNPHASED workloads
# ============================================================
unph_xlsx <- require_workbook("results_unphased.xlsx")
unph_raw <- read_excel(unph_xlsx, sheet = "Combined")

unph <- unph_raw %>%
  split_ms("cpu_t1_time",       "cpu_t1") %>%
  split_ms("cpu_t2_time",       "cpu_t2") %>%
  split_ms("gpu_dev0_t1_time",  "gpu_t1") %>%
  split_ms("gpu_dev0_t2_time",  "gpu_t2") %>%
  mutate(
    Dataset = factor(Dataset,
                     levels = c("tiny","small","medium","larger","ploidy_edge")),
    x_label = paste0(Dataset, "\n(S=", Samples, ", C=", Chromosomes, ")"),
    work    = Samples * Chromosomes
  ) %>%
  arrange(Dataset)

unph_long <- unph %>%
  select(Dataset, x_label, Samples, Chromosomes, work,
         cpu_t1_mean, cpu_t1_sd,
         cpu_t2_mean, cpu_t2_sd,
         gpu_t1_mean, gpu_t1_sd,
         gpu_t2_mean, gpu_t2_sd) %>%
  pivot_longer(
    cols = -c(Dataset, x_label, Samples, Chromosomes, work),
    names_to = c("method", ".value"),
    names_pattern = "(.*)_(mean|sd)"
  ) %>%
  mutate(
    method = recode(method,
                    cpu_t1 = "CPU 1T",
                    cpu_t2 = "CPU 2T",
                    gpu_t1 = "GPU 1T",
                    gpu_t2 = "GPU 2T"),
    method = factor(method,
                    levels = c("CPU 1T","CPU 2T","GPU 1T","GPU 2T"))
  )

range_check(unph_long, "mean", "Fig4 (unphased)")

# ---------- Figure 4: unphased workloads (log) ----------
p4 <- ggplot(unph_long,
             aes(x = x_label, y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75, colour = "black", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = safe_lower(mean, sd), ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2, linewidth = 0.35
  ) +
  scale_y_log10(breaks = c(1, 10, 100, 1000), labels = fmt_log10_s) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "Dataset",
    y = "Runtime (log scale)",
    fill = NULL,
    title = "Runtime on synthetic unphased workloads of varying size",
    subtitle = "Bars: mean; error bars: SD. X labels: (S = samples, C = chromosomes)."
  ) +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.9))

save_fig(p4, "fig4_unphased", width = 9, height = 4.5)

# ---------- Figure 4-alt (optional): workload size on X ----------
p4_work <- ggplot(unph_long,
                  aes(x = work, y = mean,
                      colour = method, shape = method)) +
  geom_errorbar(aes(ymin = safe_lower(mean, sd), ymax = mean + sd),
                width = 0.02, linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_x_log10() +
  scale_y_log10(breaks = c(1, 10, 100, 1000), labels = fmt_log10_s) +
  scale_colour_manual(values = method_colors) +
  scale_shape_manual(values = c(
    "CPU 1T" = 16, "CPU 2T" = 15, "GPU 1T" = 17, "GPU 2T" = 18
  )) +
  labs(
    x = "Workload (Samples x Chromosomes, log scale)",
    y = "Runtime (log scale)",
    colour = NULL, shape = NULL,
    title = "Unphased workloads vs workload size",
    subtitle = "Mean; error bars: SD"
  )

save_fig(p4_work, "fig4_workload", width = 6.5, height = 5)


# ============================================================
# 7) Combined multi-panel figure
#    A: Fig1  B: Fig3  C: Fig2  D: Fig4
# ============================================================
combined <- (p1 / p3) | (p2 / p4) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(2, 1))

save_fig(combined, "fig_combined", width = 16, height = 10)


# ============================================================
# 8) Preflight checks
# ============================================================

cat("\n--- Preflight -------------------------------------------------\n")

# a) any NA in parsed means/SDs?
na_counts <- bind_rows(
  bench_long %>% filter(is.na(mean) | is.na(sd)) %>%
    count(method) %>% mutate(figure = "Fig1 benchmark"),
  scal_long  %>% filter(is.na(mean) | is.na(sd)) %>%
    count(method) %>% mutate(figure = "Fig2 scaling"),
  unph_long  %>% filter(is.na(mean) | is.na(sd)) %>%
    count(method) %>% mutate(figure = "Fig4 unphased")
)
if (nrow(na_counts) == 0) {
  cat("NA in parsed means/SDs: none\n")
} else {
  cat("NA in parsed means/SDs:\n")
  print(na_counts)
}

# b) how many cells have SD > mean (cause of the spike problem)?
sd_gt <- bind_rows(
  bench_long %>%
    filter(!is.na(mean), !is.na(sd), sd > mean) %>%
    count(method, name = "n_sd_gt_mean") %>% mutate(figure = "Fig1 benchmark"),
  unph_long %>%
    filter(!is.na(mean), !is.na(sd), sd > mean) %>%
    count(method, name = "n_sd_gt_mean") %>% mutate(figure = "Fig4 unphased")
)
if (nrow(sd_gt) == 0) {
  cat("Cells where SD > mean: none\n")
} else {
  cat("Cells where SD > mean (expected for GPU on tiny contigs):\n")
  print(sd_gt)
}

# c) sample count per chromosome (should be 2 for Platinum)
n_samples_per_chrom <- bench %>% count(Chromosome, Sample) %>% count(Chromosome)
cat("Samples per chromosome (Platinum):\n")
print(n_samples_per_chrom)

# d) parse sanity checks
cat("Parse sanity checks:\n")
print(parse_ms("00m 22.21s ± 00m 00.77s"))
print(parse_ms("01m 38.99s ± 00m 06.30s"))
print(parse_ms("00m 17.45s ± 00m 22.55s"))
print(parse_ms("15m 00.00s ± 00m 12.40s"))
print(parse_ms("0.50s ± 0.02s"))

cat("---------------------------------------------------------------\n")
cat("All figures written to: ", figures_dir, "\n", sep = "")
cat("Done.\n")