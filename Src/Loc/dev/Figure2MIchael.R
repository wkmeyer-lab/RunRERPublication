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
raw_data <- readRDS("OUtput/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCombinedGOresultsWithAlternates-KeggReactome.rds")

# >>> SET YOUR DYNAMIC COLUMNS HERE <<<
target_column <- "IV-significantRobust"
histat_column <- "IV-stat"

raw_data$pathway = rownames(raw_data)
raw_data = raw_data[,c(195, 1:194)]

# -------------------------------------------------------------------------
# 3. FILTER AND PROCESS DATA DYNAMICALLY
# -------------------------------------------------------------------------
node_data <- raw_data %>%
  rename(pathway = 1) %>% 
  
  # 1. Dynamically keep rows active in your specified target column
  filter(.data[[target_column]] == 1) %>%
  
  # 2. Convert your specified HI_stat column to a strict binary category
  mutate(
    diet_direction = case_when(
      .data[[histat_column]] < 0  ~ "Invertivore", # FIXED: Removed 's' to match color scale
      .data[[histat_column]] > 0  ~ "Vertivore"    # FIXED: Removed 's' to match color scale
    ),
    # Keep it as a plain character string so igraph doesn't corrupt it
    diet_direction = as.character(diet_direction),
    
    # Extract prefix database category (KEGG, REACTOME, PID)
    category = str_extract(pathway, "^[^_]+")
  ) %>%
  select(pathway, diet_direction, category, everything())

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

num_pathways <- nrow(pathway_gene_lists)
pairs <- combn(num_pathways, 2)

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

edge_data <- tibble(
  from = pathway_gene_lists$pathway[pairs[1, ]],
  to = pathway_gene_lists$pathway[pairs[2, ]],
  similarity = map2_dbl(pathway_gene_lists$genes[pairs[1, ]], 
                        pathway_gene_lists$genes[pairs[2, ]], 
                        jaccard)
)

edge_data_filtered <- edge_data %>% 
  filter(similarity > 0.15)

# -------------------------------------------------------------------------
# 6. BUILD AND PLOT THE NETWORK GRAPH
# -------------------------------------------------------------------------
# 1. Build the initial network object
initial_net <- graph_from_data_frame(d = edge_data_filtered, 
                                     vertices = node_data, 
                                     directed = FALSE)

# 2. Identify distinct clusters (islands) in the graph
cluster_assignments <- components(initial_net)$membership

# 3. Figure out which clusters contain at least one node > 1000 genes
valid_clusters <- node_data %>%
  mutate(cluster_id = cluster_assignments[pathway]) %>%
  group_by(cluster_id) %>%
  summarise(has_large_node = any(total_reference_genes > 0)) %>%
  #filter(has_large_node == TRUE) %>%
  pull(cluster_id)

# 4. Filter the node dataset to keep only nodes belonging to valid clusters
 filtered_nodes <- node_data %>%
   mutate(cluster_id = cluster_assignments[pathway]) %>%
   filter(cluster_id %in% valid_clusters) %>%
   # Dynamically strip out the prefix (e.g. "KEGG_", "REACTOME_", "PID_")
   mutate(clean_name = str_remove(pathway, "^[^_]+_")) %>%

   # Group by the cluster ID to evaluate names locally within each island
   group_by(cluster_id) %>%
   mutate(
     # row_number() == which.min(...) flags only the single shortest string in the group
     is_shortest = row_number() == which.min(nchar(clean_name)),
     cluster_label = if_else(is_shortest, clean_name, "")
   ) %>%
   ungroup()

 # 5. Filter the edge dataset to match only remaining nodes
 filtered_edges <- edge_data_filtered %>%
   filter(from %in% filtered_nodes$pathway & to %in% filtered_nodes$pathway)

 # 6. Rebuild the final filtered network object
 net <- graph_from_data_frame(d = filtered_edges,
                              vertices = filtered_nodes,
                              directed = FALSE)

 p <- ggraph(net, layout = "fr") + 
   # CHANGED: Width maps back to the Jaccard similarity fraction
   geom_edge_link(aes(alpha = similarity, width = similarity), 
                  color = "grey40", show.legend = TRUE) + 
   scale_edge_width(range = c(0.5, 2.5)) +
   
   geom_node_point(aes(size = total_reference_genes, 
                       color = diet_direction)) + 
   
   # Use clean_name or cluster_label depending on which variation you are running
   geom_node_text(aes(label = clean_name), 
                  repel = TRUE, 
                  size = 2.5, 
                  fontface = "bold",
                  max.overlaps = 100) +
   
   scale_color_manual(
     values = c(
       "Invertivore" = "#0077BB", 
       "Vertivore"   = "#E7298A"
     ),
     name = "Diet"
   ) +
   
   scale_size_continuous(range = c(5, 16)) + 
   
   # Strips out the duplicate alpha legend, leaving only edge_width
   guides(edge_alpha = "none") +
   
   theme_void() +
   labs(
     size = "Total Genes in Reference",
     # CHANGED: Clear legend title reflecting the similarity metric
     edge_width = "Gene Overlap"
   ) +
   theme(
     legend.position = "right",
     plot.title = element_text(face = "bold", size = 16),
     legend.title = element_text(face = "bold", size = 10)
   ) +
   guides(
     edge_alpha = "none",                   # Still hides the duplicate edge legend
     color = guide_legend(order = 1),       # Forces Diet legend to the very top
     size = guide_legend(order = 2),        # Forces Node Size legend to the middle
     edge_width = guide_legend(order = 3)   # Forces Gene Overlap legend to the bottom
   )
p

ggsave("Output/COmplexDietCentralAnalysis/Visualization/Figure2AllLegends.jpg", plot = p, width = 20, height = 20, dpi = 600)

p <- ggraph(net, layout = "fr") + 
  # CHANGED: Width maps back to the Jaccard similarity fraction
  geom_edge_link(aes(alpha = similarity, width = similarity), 
                 color = "grey40", show.legend = TRUE) + 
  scale_edge_width(range = c(0.5, 2.5)) +
  
  geom_node_point(aes(size = total_reference_genes, 
                      color = diet_direction)) + 
  
  # Use clean_name or cluster_label depending on which variation you are running
  geom_node_text(aes(label = ""), 
                 repel = TRUE, 
                 size = 2.5, 
                 fontface = "bold",
                 max.overlaps = 100) +
  
  scale_color_manual(
    values = c(
      "Invertivore" = "#0077BB", 
      "Vertivore"   = "#E7298A"
    ),
    name = "Diet"
  ) +
  
  scale_size_continuous(range = c(5, 16)) + 
  
  # Strips out the duplicate alpha legend, leaving only edge_width
  guides(edge_alpha = "none") +
  
  theme_void() +
  labs(
    size = "Total Genes in Reference",
    # CHANGED: Clear legend title reflecting the similarity metric
    edge_width = "Gene Overlap"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(
    edge_alpha = "none",                   # Still hides the duplicate edge legend
    color = guide_legend(order = 1),       # Forces Diet legend to the very top
    size = guide_legend(order = 2),        # Forces Node Size legend to the middle
    edge_width = guide_legend(order = 3)   # Forces Gene Overlap legend to the bottom
  )
p

ggsave("Output/COmplexDietCentralAnalysis/Visualization/Figure2NoLabs.png", plot = get_last_plot(), width = 20, height = 20, dpi = 600)
