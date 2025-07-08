intensity<-read.delim("FragPipe Input/E733_DIA_Report_MinImp_pool_removed.tsv")
library(tidyverse)

df_long <- intensity %>%
  pivot_longer(
    cols = 7:ncol(.),  # intensity columns
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    # Extract IP type (IgG or GR)
    IP = case_when(
      str_detect(Sample, "_GR_") | str_detect(Sample, "_GR_IP") ~ "GR",
      str_detect(Sample, "_IgG_") | str_detect(Sample, "_IgG_IP") ~ "IgG",
      TRUE ~ "Other"
    ),
    # Extract cell line before the IP tag
    CellLine = Sample %>%
      str_remove("_GR_.*|_IgG_.*|_GR_IP_.*|_IgG_IP_.*") %>%
      str_replace_all("\\.+", "")  # remove any stray dots from names
  ) %>%
  filter(IP != "Other")  # keep only IgG and GR samples

df_long$log10Intensity <- log10(df_long$Intensity)
  
df_long <- df_long %>%mutate(
            CellLine = str_replace(CellLine, "^X231", "MDA-MB-231")
                )
                  


p<-ggplot(df_long, aes(x = interaction(CellLine, IP), y = log10Intensity, fill = IP)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c("IgG" = "#1f77b4", "GR" = "#ff7f0e")) +
  labs(
    x = "Cell Line and IP Type",
    y = expression(log[10]*"(Protein Intensity)"),
    title = "Protein Intensities by Cell Line and IP Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4))

#Unused
#ggsave(p, file="Rebuttal Figure X - Protein Intesities By Cell line.png",dpi=300, width=170, height=120, unit="mm")

df_long_pdx <- df_long %>%
  mutate(
    # Assign CellLine: if underscore exists, extract before it, else keep as is
    CellLine = if_else(str_detect(CellLine, "_"), str_remove(CellLine, "_.*$"), CellLine),
    # Assign PDX label for those 3 cell lines
    CellLine = if_else(CellLine %in% c("BB3RC31", "BB7", "HBC34"), "PDX", CellLine),
    # Extract replicate number from sample name where possible
    Replicate = str_extract(Sample, "_\\d+$") %>% str_remove("_")
  ) %>%
  # Now assign replicate numbers explicitly for those PDX samples based on original cell line
  mutate(
    Replicate = case_when(
      CellLine == "PDX" & str_detect(Sample, "BB3RC31") ~ "1",
      CellLine == "PDX" & str_detect(Sample, "BB7") ~ "2",
      CellLine == "PDX" & str_detect(Sample, "HBC34") ~ "3",
      TRUE ~ Replicate
    )
  )

df_summary <- df_long_pdx %>%
  group_by(CellLine, IP, Replicate) %>%
  summarise(
    AvgIntensity = mean(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(log10AvgIntensity = log10(AvgIntensity + 1e-6))

p2 <- ggplot(df_summary, aes(x = CellLine, y = log10AvgIntensity, color = IP)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.6, outlier.shape = NA) +  # boxplot without outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2, alpha = 0.8) +  # overlay points
  scale_color_manual(values = c("IgG" = "#1f77b4", "GR" = "#ff7f0e")) +
  labs(
    x = "Cell Line",
    y = expression(log[10]*"(Average Protein Intensity)"),
    title = "Average Protein Intensity per Replicate by Cell Line and IP"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p2, file="Rebuttal Figure 1 - Average Protein Intesities By Cell line.png",dpi=300, width=170, height=120, unit="mm")
