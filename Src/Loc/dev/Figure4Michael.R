# -------------------------------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# -------------------------------------------------------------------------
library(tidyverse)
library(igraph)
library(ggraph)
library(msigdbr)

# -------------------------------------------------------------------------
# 2. CONFIGURATION & LOADING
# -------------------------------------------------------------------------
# Replace with your actual file
raw_data <- read_csv("ComplexDietCentralAnalysisGeneSets.csv")
raw_data <- readRDS("OUtput/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCombinedGOresultsWithAlternates-KeggReactome.rds")

target_column <- "Fig4"

# -------------------------------------------------------------------------
# 3. FILTER AND PROCESS DATA DYNAMICALLY
# -------------------------------------------------------------------------
node_data <- raw_data %>%
  rename(pathway = 1) %>% 
  
  # 1. Primary activity filter (Fig4)
  filter(.data[[target_column]] == 1) %>%
  
  # FIX: Force columns to numeric first so R reads the decimals correctly
  mutate(
    across(c(CH_stat, HI_stat, HV_stat, CH, HI, HV), as.numeric)
  ) %>%
  
  # 2. CONDITIONAL FILTERING: Only enforce the stat direction if the group is active (== 1).
  # If the group is inactive (== 0), it automatically passes (TRUE).
  filter(
    if_else(CH == 1, CH_stat < 0, TRUE),
    if_else(HI == 1, HI_stat > 0, TRUE),
    if_else(HV == 1, HV_stat > 0, TRUE)
  ) %>%
  
  # 3. Create the overlap/master categorizations based on CH, HI, HV
  mutate(
    master_group = case_when(
      CH == 1 & HI == 0 & HV == 0 ~ "CH Only",
      CH == 0 & HI == 1 & HV == 0 ~ "HI Only",
      CH == 0 & HI == 0 & HV == 1 ~ "HV Only",
      CH == 1 & HI == 1 & HV == 0 ~ "CH + HI",
      CH == 1 & HI == 0 & HV == 1 ~ "CH + HV",
      CH == 0 & HI == 1 & HV == 1 ~ "HI + HV",
      CH == 1 & HI == 1 & HV == 1 ~ "All Three",
      TRUE                        ~ "None"
    ),
    
    # Create distinct binary flags for the individual plots
    plot_CH = if_else(CH == 1, "CH Active", "Inactive"),
    plot_HI = if_else(HI == 1, "HI Active", "Inactive"),
    plot_HV = if_else(HV == 1, "HV Active", "Inactive")
  )

# -------------------------------------------------------------------------
# 4. FETCH CANONICAL GENES FROM MSigDB
# -------------------------------------------------------------------------
human_msigdb <- msigdbr(species = "Homo sapiens")

filtered_ref <- human_msigdb %>%
  filter(gs_name %in% node_data$pathway) %>%
  select(pathway = gs_name, gene = gene_symbol)

node_sizes <- filtered_ref %>%
  group_by(pathway) %>%
  summarise(total_reference_genes = n(), .groups = "drop")

node_data <- node_data %>%
  left_join(node_sizes, by = "pathway") %>%
  drop_na(total_reference_genes)

# -------------------------------------------------------------------------
# 5. CALCULATE OVERLAP (EDGES) BETWEEN PATHWAYS
# -------------------------------------------------------------------------
pathway_gene_lists <- filtered_ref %>%
  group_by(pathway) %>%
  summarise(genes = list(gene), .groups = "drop")

pairs <- combn(nrow(pathway_gene_lists), 2)
jaccard <- function(a, b) { length(intersect(a, b)) / length(union(a, b)) }

edge_data <- tibble(
  from = pathway_gene_lists$pathway[pairs[1, ]],
  to = pathway_gene_lists$pathway[pairs[2, ]],
  similarity = map2_dbl(pathway_gene_lists$genes[pairs[1, ]], 
                        pathway_gene_lists$genes[pairs[2, ]], jaccard)
) %>% 
  filter(similarity > 0.15)

# -------------------------------------------------------------------------
# 6. BUILD GRAPH, FILTER SIZE, AND **FREEZE THE LAYOUT**
# -------------------------------------------------------------------------
# Build initial graph to find clusters and filter sizes
initial_net <- graph_from_data_frame(d = edge_data, vertices = node_data, directed = FALSE)
cluster_assignments <- components(initial_net)$membership

# --- NEW: Drop clusters with fewer than X nodes ---
min_cluster_size <- 1 # Adjust this number to change your drop threshold

filtered_nodes <- node_data %>%
  mutate(cluster_id = cluster_assignments[pathway]) %>%
  # Count sizes per cluster group
  group_by(cluster_id) %>%
  mutate(cluster_size = n()) %>%
  # Filter step: drops clusters below the node threshold AND keeps your size rule
  # (Note: Remove 'total_reference_genes > 1500' check if it drops too many nodes)
  filter(cluster_size >= min_cluster_size, any(total_reference_genes > 0)) %>% 
  
  mutate(clean_name = str_remove(pathway, "^[^_]+_")) %>%
  mutate(
    is_shortest = row_number() == which.min(nchar(clean_name)),
    cluster_label = if_else(is_shortest, clean_name, "")
  ) %>% 
  ungroup()

filtered_edges <- edge_data %>%
  filter(from %in% filtered_nodes$pathway & to %in% filtered_nodes$pathway)

# Rebuild final network
net <- graph_from_data_frame(d = filtered_edges, vertices = filtered_nodes, directed = FALSE)

# Compute and save coordinates
fixed_layout <- create_layout(net, layout = "fr")


# =========================================================================
# 7. PLOTTING BLOCKS
# =========================================================================

# Shared aesthetic themes for consistency
base_theme <- theme_void() + 
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey30", size = 12),
        legend.title = element_text(face = "bold", size = 10))

# Custom color palette matching the individual plots
master_colors <- c(
  "CH Only"   = "#EE7733",  # Exact color from plot_CH
  "HI Only"   = "#377EB8",  # Exact color from plot_HI
  "HV Only"   = "#4DAF4A",  # Exact color from plot_HV
  "CH + HI"   = "#984EA3",  # Purple blend
  "CH + HV"   = "#FF7F00",  # Orange/Yellow-Green blend
  "HI + HV"   = "#00A08A",  # Teal blend
  "All Three" = "#A65628",  # Brown/Neutral mix
  "None"      = "grey80"
)

# -------------------------------------------------------------------------
# PLOT 1: MASTER NETWORK (WITH CUSTOM INHERITED COLORS)
# -------------------------------------------------------------------------
plot_master <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  geom_node_point(aes(size = total_reference_genes, color = master_group)) + 
  geom_node_text(aes(label = if_else(total_reference_genes > 1, clean_name, ""),
                     alpha = if_else(total_reference_genes > 1, 1.0, 0.0)), 
                 repel = TRUE, 
                 size = 3, 
                 fontface = "bold", 
                 max.overlaps = 100) +
  # NEW: Overriding the color scales using our manual mapping
  scale_color_manual(values = master_colors, name = "Overlap Group") +
  scale_size_continuous(range = c(10, 30)) + guides(edge_alpha = "none") +
  base_theme + labs(size = "Category Size", edge_width = "Category Overlap")

plot_master_nolab <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  geom_node_point(aes(size = total_reference_genes, color = master_group)) + 
  geom_node_text(aes(label = "")) +
  # NEW: Overriding the color scales using our manual mapping
  scale_color_manual(values = master_colors, name = "Overlap Group") +
  scale_size_continuous(range = c(10, 30)) + guides(edge_alpha = "none") +
  base_theme + labs(size = "Category Size", edge_width = "Category Overlap")

# -------------------------------------------------------------------------
# PLOT 2: CH COMPARISON ONLY
# -------------------------------------------------------------------------
plot_CH <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  geom_node_point(aes(size = total_reference_genes, color = plot_CH, alpha = plot_CH)) + 
  geom_node_text(aes(label = cluster_label, alpha = plot_CH), repel = TRUE, size = 3, fontface = "bold") +
  scale_color_manual(values = c("CH Active" = "#EE7733", "Inactive" = "grey80"), name = "CH Status") +
  scale_alpha_manual(values = c("CH Active" = 1.0, "Inactive" = 0.2), guide = "none") +
  scale_size_continuous(range = c(10, 30)) + guides(edge_alpha = "none") +
  base_theme + labs(title = "CH Comparison Network", size = "Reference Genes", edge_width = "Gene Overlap")

# -------------------------------------------------------------------------
# PLOT 3: HI COMPARISON ONLY
# -------------------------------------------------------------------------
plot_HI <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  geom_node_point(aes(size = total_reference_genes, color = plot_HI, alpha = plot_HI)) + 
  geom_node_text(aes(label = cluster_label, alpha = plot_HI), repel = TRUE, size = 3, fontface = "bold") +
  scale_color_manual(values = c("HI Active" = "#377EB8", "Inactive" = "grey40"), name = "HI Status") +
  scale_alpha_manual(values = c("HI Active" = 1.0, "Inactive" = 0.2), guide = "none") +
  scale_size_continuous(range = c(10, 30)) + guides(edge_alpha = "none") +
  base_theme + labs(title = "HI Comparison Network", size = "Reference Genes", edge_width = "Gene Overlap")

# -------------------------------------------------------------------------
# PLOT 4: HV COMPARISON ONLY
# -------------------------------------------------------------------------
plot_HV <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  geom_node_point(aes(size = total_reference_genes, color = plot_HV, alpha = plot_HV)) + 
  geom_node_text(aes(label = cluster_label, alpha = plot_HV), repel = TRUE, size = 3, fontface = "bold") +
  scale_color_manual(values = c("HV Active" = "#4DAF4A", "Inactive" = "grey40"), name = "HV Status") +
  scale_alpha_manual(values = c("HV Active" = 1.0, "Inactive" = 0.2), guide = "none") +
  scale_size_continuous(range = c(10, 30)) + guides(edge_alpha = "none") +
  base_theme + labs(title = "HV Comparison Network", size = "Reference Genes", edge_width = "Gene Overlap")

# Print them
#plot_master
#plot_master_nolab
#plot_CH
#plot_HI
#plot_HV

ggsave("Fig4.jpg", plot_master, height = 20, width = 20, dpi = 300)
ggsave("Fig4_nolab.jpg", plot_master_nolab, height = 20, width = 20, dpi = 300)
