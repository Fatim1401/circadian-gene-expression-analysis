## Load data -------------------------------------------------------------------
library(dplyr)
library(ggplot2)
rm(list = ls())
link.list_leaf <- read.table("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/link_list_dynGENIE3.txt", header = TRUE)
link.list_seed <- read.table("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/link_list_dynGENIE3_seed.txt", header = TRUE)
GO_description <- read.csv("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/GO_description.csv",header = TRUE)
GO_annotation_v2 <- read.csv2("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/Gene_annotation_v2.csv", header =TRUE)


################################################################################
Go_abiotic_L <- "GO:0009628" 
## filter the list of results dyngenie3 ----------------------------------------
link.list_leaf_0.025 <- link.list_leaf %>%
  filter(weight >= 0.025) 

## filter the interest GO in annotation table ----------------------------------
Go_annot_abiotic_L <- GO_annotation_v2 %>%
  filter(GO_annotation_v2$GO %in% Go_abiotic_L) 

## filter the target gene in the abiotic annotation table ----------------------
GO_network_abiotic_Target_L <- link.list_leaf_0.025 %>%
  filter(link.list_leaf_0.025$target.gene %in% Go_annot_abiotic_L$GID)

filter_reg_tar <- GO_network_abiotic_Target_L %>% 
  filter(GO_network_abiotic_Target_L$regulatory.gene %in% Go_annot_abiotic_L$GID)
length(unique(filter_reg_tar$regulatory.gene))
length(unique(filter_reg_tar$target.gene))
count.regulator <- data.frame(table(filter_reg_tar$regulatory.gene)) 

 # export data ------------------------------

write.table(filter_reg_tar,
            "E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/filter_reg_tar.txt",
            sep = "\t",
            row.names = FALSE)              

 write.table(count.regulator,
            "E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/count.regulator1.txt",
            sep = "\t",
            row.names = FALSE) 
 
library(openxlsx)
write.xlsx(count.regulator,
           "E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/count.regulator1.xlsx",
           row.names = FALSE)

################################################################################
  
## Ad phase information ---------------------------
# Merge

phase_table <- read.table("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/Leaf_Metacycle_results.txt", header = TRUE)
merged_data_TF <- filter_reg_tar %>% 
  left_join(phase_table, by = c("regulatory.gene" = "CycID"))

merged_data_TF <- merged_data_TF %>%
  mutate(
    JTK_adjphase = case_when(
      JTK_adjphase == 24 ~ 0,
      JTK_adjphase == 26 ~ 2,
      TRUE ~ JTK_adjphase
    )
  )
## unique TF -------------------------------------------------------------------
TF_unique <- merged_data_TF %>%
  distinct(regulatory.gene,JTK_adjphase)

## Count --------------------------
count.regulator_phase_TF <- data.frame(table(TF_unique$JTK_adjphase))
colnames(count.regulator_phase_TF) <- c("phase", "n")

# Plot
cairo_pdf("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/figures/Phase_TF_leaf.pdf",
          width = 9, height = 9)

ggplot(count.regulator_phase_TF, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_TF$n, na.rm = TRUE), " Facteurs de transcription"),
    x = "Phase (ZT)",
    y = "Nombre de gènes régulateurs"
  ) 

dev.off()

################################################################################

merged_data_target <- filter_reg_tar %>% 
  left_join(phase_table, by = c("target.gene" = "CycID"))

merged_data_target <- merged_data_target %>%
  mutate(
    JTK_adjphase = case_when(
      JTK_adjphase == 24 ~ 0,
      JTK_adjphase == 26 ~ 2,
      TRUE ~ JTK_adjphase
    )
  )

## target unique ---------------------------------------------------------------
Target_unique <- merged_data_target %>%
  distinct(target.gene,JTK_adjphase)

## Count --------------------------
count.regulator_phase_target <- data.frame(table(Target_unique$JTK_adjphase))
colnames(count.regulator_phase_target) <- c("phase", "n")

# Plot
cairo_pdf("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/figures/Phase_Target_leaf.pdf",
          width = 9, height = 9)

ggplot(count.regulator_phase_target, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_target$n, na.rm = TRUE), " gènes cibles"),
    x = "Phase (ZT)",
    y = "Nombre de gènes cibles" 
  ) 

dev.off() 

################################################################################

################################################################################


## graines  --------------------------------------------------------------------

## Define GO_term --------------------------------------------------------------
Go_abiotic_S <- "GO:0009628" 

## filter the list of results dyngenie3 ----------------------------------------
link.list_seed_0.025 <- link.list_seed %>%
  filter(weight >= 0.025)

## filter the interest GO in annotation table ----------------------------------
Go_annot_abiotic_S <- GO_annotation_v2 %>%
  filter(GO_annotation_v2$GO %in% Go_abiotic_S) 

## filter the target gene in the abiotic annotation table ----------------------
GO_network_abiotic_Target_S <- link.list_seed_0.025 %>%
  filter(link.list_seed_0.025$target.gene %in% Go_annot_abiotic_S$GID)
filter_reg_tar_s <- GO_network_abiotic_Target_S %>% 
  filter(GO_network_abiotic_Target_S$regulatory.gene %in% Go_annot_abiotic_S$GID)

length(unique(filter_reg_tar_s$regulatory.gene))
length(unique(filter_reg_tar_s$target.gene))
count.regulator_s <- data.frame(table(filter_reg_tar_s$regulatory.gene)) 

## export data --------------------------

write.table(filter_reg_tar_s,
            "E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/filter_reg_tar_s.txt",
            row.names = FALSE,
            sep = "\t")

write.table(count.regulator_s,
            "E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/count.regulator_s.txt",
            sep = "\t",
            row.names = FALSE) 

library(openxlsx)
write.xlsx(count.regulator_s,
           "E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/count.regulator_s.xlsx",
           row.names = FALSE)

################################################################################
## Ad phase information ---------------------------
# Merge

phase_table <- read.table("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/Leaf_Metacycle_results.txt", header = TRUE)
merged_data_TF_s <- filter_reg_tar_s %>% 
  left_join(phase_table, by = c("regulatory.gene" = "CycID"))

merged_data_TF_s <- merged_data_TF_s %>%
  mutate(
    JTK_adjphase = case_when(
      JTK_adjphase == 24 ~ 0,
      JTK_adjphase == 26 ~ 2,
      TRUE ~ JTK_adjphase
    )
  )
## unique TF -------------------------------------------------------------------
TF_unique_s <- merged_data_TF_s %>%
  distinct(regulatory.gene,JTK_adjphase)

## Count --------------------------
count.regulator_phase_TF_s <- data.frame(table(TF_unique_s$JTK_adjphase))
colnames(count.regulator_phase_TF_s) <- c("phase", "n")

# Plot
cairo_pdf("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/figures/Phase_TF_seed.pdf",
          width = 9, height = 9)

ggplot(count.regulator_phase_TF_s, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_TF_s$n, na.rm = TRUE), " Facteurs de transcription"),
    x = "Phase (ZT)",
    y = "Nombre de gènes régulateurs"
  ) 

dev.off()

################################################################################

merged_data_target_s <- filter_reg_tar_s %>% 
  left_join(phase_table, by = c("target.gene" = "CycID"))

merged_data_target_s <- merged_data_target_s %>%
  mutate(
    JTK_adjphase = case_when(
      JTK_adjphase == 24 ~ 0,
      JTK_adjphase == 26 ~ 2,
      TRUE ~ JTK_adjphase
    )
  )

## target unique ---------------------------------------------------------------
Target_unique_s <- merged_data_target_s %>%
  distinct(target.gene,JTK_adjphase)

## Count --------------------------
count.regulator_phase_target_s<- data.frame(table(Target_unique_s$JTK_adjphase))
colnames(count.regulator_phase_target_s) <- c("phase", "n")

# Plot
cairo_pdf("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/figures/Phase_Target_seed.pdf",
          width = 9, height = 9)

ggplot(count.regulator_phase_target_s, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_target_s$n, na.rm = TRUE), " gènes cibles"),
    x = "Phase (ZT)",
    y = "Nombre de gènes cibles" 
  ) 

dev.off() 
############################################################################
# Overlaf between Top 15 tf seed and leaf 
library(readxl)

Top15_seed <- read_excel("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/TOP_ 15_TF_S1.xlsx")

Top15_leaf <- read_excel("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/data/Top_15 _TF_1.xlsx")

Top15_leaf$condition <- "Feuilles"
Top15_seed$condition <- "Graines" 
library(dplyr)

Top15_seed <- Top15_seed %>%
  select(IDPS_V2, IDPS_V1, DESCRIPTION_PS, Freq, condition)

Top15_leaf <- Top15_leaf %>%
  select(IDPS_V2, IDPS_V1, DESCRIPTION_PS, Freq,condition)
data <- bind_rows(Top15_seed, Top15_leaf)


library(ggplot2)
cairo_pdf("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/figures/top15_graine.PDF",width = 6, height =7)
 ggplot(Top15_seed, aes(x = reorder(DESCRIPTION_PS, Freq), y = Freq)) +
  
  geom_bar(stat = "identity",
           fill = "#e31a1c",
           width = 0.7) +
  
  coord_flip() +
  
  labs(
    title = "Top 15 Facteurs de trancription - Graines",
    x = "Facteurs de transcription",
    y = "Nombres de gènes régulés"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10)
  )
dev.off() 

cairo_pdf("E:/Stage_Fatoumata/Projet-de-Stage-/outputs/figures/top15_feuilles.PDF",width = 6, height =7)
ggplot(Top15_leaf, aes(x = reorder(DESCRIPTION_PS, Freq), y = Freq)) +
  
  geom_bar(stat = "identity",
           fill = "#1f78b4",
           width = 0.7) +
  
  coord_flip() +
  
  labs(
    title = "Top 15 Facteurs de transcription - Feuilles",
    x = "Facteurs de transcription",
    y = "Nombres de gènes régulés"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10)
  )
 dev.off()
 
 
 library(VennDiagram)
 library(grid)
 
 # TF communs
 overlap_leaf_seed <- intersect(Top15_seed$DESCRIPTION_PS, Top15_leaf$DESCRIPTION_PS)
 
 # Voir les TF communs
 length(overlap_leaf_seed)
 overlap_leaf_seed
 
 # Création du PDF
 cairo_pdf(
   "E:/Stage_Fatoumata/Projet-de-Stage-/outputs/figures/Venn_Graine_feuilles_top15.pdf",
   width = 6,
   height = 6
 )
 
 # Dessiner le Venn
 draw.pairwise.venn(
   area1 = length(Top15_leaf$DESCRIPTION_PS),
   area2 = length(Top15_seed$DESCRIPTION_PS),
   cross.area = length(overlap_leaf_seed),
   category = c("Feuilles", "Graines"),
   fill = c("#1f78b4", "#e31a1c"),
   alpha = 0.7,
   cex = 1.3,
   cat.cex = 1.2,
   lty = "blank"
 )
 
 # Ajouter les noms des TF communs dans l'intersection
 grid.text(
   paste(overlap_leaf_seed, collapse = "\n"),
   x = 0.5,
   y = 0.5,
   gp = gpar(fontsize = 8)
 )
 
 dev.off()
 
 
 
 
 
  
  
  
