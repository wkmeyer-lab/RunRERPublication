# =========================================================================
# 1. LOAD REQUIRED PACKAGES
# =========================================================================
library(tidyverse)
library(igraph)
library(ggraph)
library(msigdbr)
library(stringr)
library(ggforce)  # For pie chart geometry on coordinates
library(scales)   # For color_ramp gradient generation

# =========================================================================
# 2. CONFIGURATION & LOADING
# =========================================================================
raw_data <- read_csv("ComplexDietCentralAnalysisGeneSets.csv", 
                     col_types = cols(.default = "c"))

# =========================================================================
# 3. FILTER AND PROCESS DATA DYNAMICALLY
# =========================================================================
node_data <- raw_data %>%
  rename(pathway = 1) %>% 
  
  # Safely forces text stats to numeric, retaining all tiny decimals, then replaces NA with 0
  mutate(
    across(c(CH_stat, HI_stat, HV_stat, 
             CH_PadjNumSignificant, HI_PadjNumSignificant, HV_PadjNumSignificant), 
           ~ suppressWarnings(as.numeric(.))),
    across(c(CH_stat, HI_stat, HV_stat, 
             CH_PadjNumSignificant, HI_PadjNumSignificant, HV_PadjNumSignificant), 
           ~ replace_na(., 0))
  ) %>%
  
  # CONDITIONAL DIRECTION CHECK
  filter(
    if_else(CH_PadjNumSignificant > 0, CH_stat <= 0, TRUE),
    if_else(HI_PadjNumSignificant > 0, HI_stat >= 0, TRUE),
    if_else(HV_PadjNumSignificant > 0, HV_stat >= 0, TRUE)
  ) %>%
  
  # OVERALL SIGNIFICANCE FILTER
  filter(
    CH_PadjNumSignificant > 0 | HI_PadjNumSignificant > 0 | HV_PadjNumSignificant > 0
  ) %>%
  
  # EXPLICIT COLOR ASSIGNMENT FOR EVERY INTERSECTION
  rowwise() %>%
  mutate(
    is_CH = CH_PadjNumSignificant > 0,
    is_HI = HI_PadjNumSignificant > 0,
    is_HV = HV_PadjNumSignificant > 0,
    
    master_hex_color = case_when(
      is_CH & is_HI & is_HV ~ "#E7298A", # All Three: Magenta
      is_CH & is_HI         ~ "#984EA3", # CH & HI overlap: Purple
      is_CH & is_HV         ~ "#E6AB02", # CH & HV overlap: Gold
      is_HI & is_HV         ~ "#00CED1", # HI & HV overlap: Teal
      is_CH                 ~ "#EE7733", # CH only: Orange
      is_HI                 ~ "#377EB8", # HI only: Blue
      is_HV                 ~ "#4DAF4A", # HV only: Green
      TRUE                  ~ "#808080"
    ),
    
    # ABSOLUTE SHADING (Gradient logic based on the highest significance value)
    max_sig = pmax(CH_PadjNumSignificant, HI_PadjNumSignificant, HV_PadjNumSignificant),
    master_alpha = 0.35 + 0.65 * (pmin(max_sig, 100) / 100)
  ) %>%
  ungroup()

# =========================================================================
# 4. FETCH CANONICAL GENES FROM MSigDB
# =========================================================================
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

# =========================================================================
# 5. CALCULATE OVERLAP (EDGES) BETWEEN PATHWAYS
# =========================================================================
pathway_gene_lists <- filtered_ref %>%
  group_by(pathway) %>%
  summarise(genes = list(gene), .groups = "drop")

pairs <- combn(nrow(pathway_gene_lists), 2)
# Jaccard index: Identifies overlaps > 15% to build network springs
jaccard <- function(a, b) { length(intersect(a, b)) / length(union(a, b)) }

edge_data <- tibble(
  from = pathway_gene_lists$pathway[pairs[1, ]],
  to = pathway_gene_lists$pathway[pairs[2, ]],
  similarity = map2_dbl(pathway_gene_lists$genes[pairs[1, ]], 
                        pathway_gene_lists$genes[pairs[2, ]], jaccard)
) %>% 
  filter(similarity > 0.15)

# -------------------------------------------------------------------------
# 6. BUILD GRAPH, FILTER SIZE, AND FREEZE THE LAYOUT
# -------------------------------------------------------------------------
initial_net <- graph_from_data_frame(d = edge_data, vertices = node_data, directed = FALSE)
cluster_assignments <- components(initial_net)$membership

min_cluster_size <- 1 

filtered_nodes <- node_data %>%
  mutate(cluster_id = cluster_assignments[pathway]) %>%
  group_by(cluster_id) %>%
  mutate(cluster_size = n()) %>%
  filter(cluster_size >= min_cluster_size, any(total_reference_genes > 0)) %>% 
  
  mutate(
    clean_name = str_remove(pathway, "^[^_]+_"),
    cluster_label = clean_name  # <-- UPDATED: Every node gets its own label now
  ) %>% 
  ungroup()

filtered_edges <- edge_data %>%
  filter(from %in% filtered_nodes$pathway & to %in% filtered_nodes$pathway)

# Rebuild final network
net <- graph_from_data_frame(d = filtered_edges, vertices = filtered_nodes, directed = FALSE)

# Compute and save coordinates so every plot has the exact same architecture
fixed_layout <- create_layout(net, layout = "fr")

# =========================================================================
# 7. PLOTTING BLOCKS 
# =========================================================================

base_theme <- theme_void() + 
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey30", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

# Create dummy data for the Master Plot's color key legend, detailing all 7 categories
dummy_legend_data <- data.frame(
  group = c("CH Only", 
            "HI Only", 
            "HV Only", 
            "CH & HI", 
            "CH & HV", 
            "HI & HV", 
            "CH & HI & HV"),
  x = 0, y = 0
)

# -------------------------------------------------------------------------
# PLOT 1: MASTER NETWORK
# -------------------------------------------------------------------------
plot_master <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  
  # Layer 1: Filled node with mapped alpha shading (gradient)
  geom_node_point(aes(size = total_reference_genes, color = master_hex_color, alpha = master_alpha)) + 
  
  # Layer 2: Opaque stroke inheriting the EXACT color of the blended node
  geom_node_point(aes(size = total_reference_genes, color = master_hex_color), shape = 1, stroke = 1.2) +
  
  # Layer 3: Invisible points to force a base-color legend
  geom_point(data = dummy_legend_data, aes(x = x, y = y, fill = group), alpha = 0, stroke = 0) +
  
  geom_node_text(aes(label = cluster_label), repel = TRUE, size = 3, fontface = "bold") +
  
  # Scales
  scale_color_identity() + 
  scale_alpha_identity() +
  scale_size_continuous(range = c(5, 16)) +
  scale_fill_manual(
    name = "Base Color Key",
    values = c("CH Only" = "#EE7733", 
               "HI Only"   = "#377EB8", 
               "HV Only"  = "#4DAF4A", 
               "CH & HI" = "#984EA3",
               "CH & HV"   = "#E6AB02",
               "HI & HV"   = "#00CED1",
               "CH & HI & HV" = "#E7298A"),
    # Keep the legend items ordered as defined in the dummy data factor
    breaks = c("CH Only", "HI Only", "HV Only", 
               "CH & HI", "CH & HV", "HI & HV", "CH & HI & HV"),
    guide = guide_legend(override.aes = list(alpha = 1, size = 6, shape = 21, color = "black"))
  ) +
  guides(edge_alpha = "none") +
  base_theme + 
  labs(size = "Category Size", edge_width = "Category Overlap")

# Save output
ggsave("Fig4B.jpg", plot_master, height = 20, width = 20, dpi = 300)

# -------------------------------------------------------------------------
# PLOT 2: CH COMPARISON 
# -------------------------------------------------------------------------
plot_CH <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  
  geom_node_point(aes(size = total_reference_genes, color = CH_PadjNumSignificant)) + 
  
  geom_node_text(aes(label = cluster_label, alpha = if_else(CH_PadjNumSignificant > 0, 1, 0)), repel = TRUE, size = 3, fontface = "bold") +
  
  scale_color_gradient(low = "grey85", high = "#EE7733", name = "CH Alternates", limits = c(0, 100), na.value = "#EE7733", guide = "colorbar") +
  scale_alpha_identity() + scale_size_continuous(range = c(5, 16)) + guides(edge_alpha = "none") + base_theme + 
  labs(title = "CH Network by Significance Frequency")

# Save output
ggsave("Fig4CH.jpg", plot_CH, height = 20, width = 20, dpi = 300)

# -------------------------------------------------------------------------
# PLOT 3: HI COMPARISON
# -------------------------------------------------------------------------
plot_HI <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  
  geom_node_point(aes(size = total_reference_genes, color = HI_PadjNumSignificant)) + 
  
  geom_node_text(aes(label = cluster_label, alpha = if_else(HI_PadjNumSignificant > 0, 1, 0)), repel = TRUE, size = 3, fontface = "bold") +
  scale_color_gradient(low = "grey85", high = "#377EB8", name = "HI Alternates", limits = c(0, 100), na.value = "#377EB8", guide = "colorbar") +
  scale_alpha_identity() + scale_size_continuous(range = c(5, 16)) + guides(edge_alpha = "none") + base_theme + 
  labs(title = "HI Network by Significance Frequency")

# Save output
ggsave("Fig4HI.jpg", plot_HI, height = 20, width = 20, dpi = 300)

# -------------------------------------------------------------------------
# PLOT 4: HV COMPARISON
# -------------------------------------------------------------------------
plot_HV <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  
  geom_node_point(aes(size = total_reference_genes, color = HV_PadjNumSignificant)) + 
  
  geom_node_text(aes(label = cluster_label, alpha = if_else(HV_PadjNumSignificant > 0, 1, 0)), repel = TRUE, size = 3, fontface = "bold") +
  scale_color_gradient(low = "grey85", high = "#4DAF4A", name = "HV Alternates", limits = c(0, 100), na.value = "#4DAF4A", guide = "colorbar") +
  scale_alpha_identity() + scale_size_continuous(range = c(5, 16)) + guides(edge_alpha = "none") + base_theme + 
  labs(title = "HV Network by Significance Frequency")

# Save output
ggsave("FigHV.jpg", plot_HV, height = 20, width = 20, dpi = 300)

# =========================================================================
# 7. ALTERNATIVE PLOT: PIE CHART NODES WITH INDEPENDENT GRADIENTS
# =========================================================================

# Step 1: Extract coordinates from fixed_layout and scale a radius for the pies
layout_df <- as_tibble(fixed_layout)
x_range <- diff(range(layout_df$x, na.rm = TRUE))

# Radius scaling for nodes
layout_df <- layout_df %>%
  mutate(
    pie_radius = scales::rescale(sqrt(total_reference_genes), 
                                 to = c(x_range * 0.005, x_range * 0.016))
  )

# Step 2: Wrangle data into long format to calculate pie slices and colors
pie_data <- layout_df %>%
  select(name, x, y, pie_radius, 
         CH_PadjNumSignificant, HI_PadjNumSignificant, HV_PadjNumSignificant) %>%
  pivot_longer(
    cols = ends_with("Significant"),
    names_to = "comparison", 
    values_to = "sig_val"
  ) %>%
  filter(sig_val > 0) %>% 
  group_by(name) %>%
  mutate(
    # Slice math
    n_slices = n(),
    slice_angle = 2 * pi / n_slices,
    end_angle = cumsum(slice_angle),
    start_angle = end_angle - slice_angle,
    
    # Gradient math for the inside fill
    norm_val = pmin(sig_val / 50, 1),
    slice_color = case_when(
      str_detect(comparison, "CH") ~ scales::colour_ramp(c("grey85", "#EE7733"))(norm_val),
      str_detect(comparison, "HI") ~ scales::colour_ramp(c("grey85", "#377EB8"))(norm_val),
      str_detect(comparison, "HV") ~ scales::colour_ramp(c("grey85", "#4DAF4A"))(norm_val)
    ),
    
    # NEW: Full opacity base color for the split outer stroke
    stroke_color = case_when(
      str_detect(comparison, "CH") ~ "#EE7733",
      str_detect(comparison, "HI") ~ "#377EB8",
      str_detect(comparison, "HV") ~ "#4DAF4A"
    ),
    
    # Internal slice dividers (white if sliced, clear if single)
    border_color = if_else(n_slices > 1, "white", NA_character_)
  ) %>%
  ungroup()

# Step 3: Draw the plot
plot_pie_network <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = similarity, width = similarity), color = "grey40") + 
  scale_edge_width(range = c(1.2, 3.0)) +
  
  # Invisible point layer solely to generate the ggplot size legend
  geom_node_point(aes(size = total_reference_genes), color = "transparent") +
  
  # 1. Draw the pie slices (fills and internal white dividers)
  geom_arc_bar(
    data = pie_data,
    aes(x0 = x, y0 = y, r0 = 0, r = pie_radius, 
        start = start_angle, end = end_angle, fill = slice_color, color = border_color),
    linewidth = 0.4, inherit.aes = FALSE
  ) +
  
  # 2. NEW: Draw the split outer stroke matching each slice's base category color
  geom_arc(
    data = pie_data,
    aes(x0 = x, y0 = y, r = pie_radius, 
        start = start_angle, end = end_angle, color = stroke_color),
    linewidth = 1.2, inherit.aes = FALSE
  ) +
  
  geom_node_text(aes(label = cluster_label), repel = TRUE, size = 3, fontface = "bold") +
  
  # Tells ggplot to use our hex codes directly for fill/color without altering them
  scale_fill_identity() + 
  scale_color_identity() + 
  
  # Formats the size legend explicitly to match your other plots
  scale_size_continuous(range = c(5, 16), name = "Category Size") + 
  guides(edge_alpha = "none", edge_width = "none") +
  base_theme

# Save output
ggsave("Fig4B_PieNodes.jpg", plot_pie_network, height = 20, width = 20, dpi = 300)