library(tidyverse)
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(cowplot)
library(ggpubr)
library(igraph)
library(phyloseq)
# Network analyses --------------------------------------------------------
source("libs/libs.R")
load("data/osd2014_sparcc_graph.Rda", verbose = TRUE)
load("data/osd2014_amp_mg_intersect.Rda")
load("data/osd2014_plancto_asvs_beta.Rda")
load("data/osd2014_dada2_phyloseq_beta.Rda")
load("data/osd2014_plancto_phyloseq_beta.Rda")
load("data/osd2014_plancto_pcoa.Rda")
load("data/com_colors.Rda")
selnodes <- V(osd2014_sparcc_graph)[asv %in% osd2014_plancto_asvs_beta$asv]
# get their network neighborhood
selegoV <- ego(osd2014_sparcc_graph, order=1, nodes = selnodes, mode = "all", mindist = 0)

# turn the returned list of igraph.vs objects into a graph
selegoG <- induced_subgraph(osd2014_sparcc_graph, unlist(selegoV))

# Get the taxonomic classification
selegoG_tax <- selegoG %>%
  as_tbl_graph() %>%
  inner_join(as(tax_table(osd2014_dada2_phyloseq_beta), "matrix") %>% as_tibble(rownames = "asv") %>% dplyr::select(asv, Phylum, Class, Order, Family, Genus)) %>% # %>% unite(tax, c(Phylum, Class, Order, Family), sep = ";")) %>%
  dplyr::select(-asv)

# Create the letter-value plot
plancto_lv <- selegoG_tax %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  dplyr::mutate(name = Phylum) %>%
  as_data_frame(what = "edges") %>%
  filter(grepl("Planctomycetes", from) | grepl("Planctomycetes", to)) %>%
  dplyr::mutate(item1 = ifelse(from == "Planctomycetes", from, to),
                item2 =  ifelse(from != "Planctomycetes", from, to)) %>%
  as_tibble() %>%
  unite(associations, c(item1, item2), sep = "|", remove = FALSE) %>%
  group_by(item2) %>%
  add_count() %>%
  ungroup() %>%
  dplyr::mutate(item2 = paste0(item2, " (", n, ")"),
                item2 = fct_reorder(item2, n, .desc = FALSE)) %>%
  ggplot(aes(item2, weight_orig)) +
  lvplot::geom_lv(width.method = "height", color = "#404040",
                  fill = "#CB566A",
                  outlier.colour = "#CB566A",
                  outlier.size = 1, outlier.stroke = 0.2,  varwidth=FALSE, size = 0.2) +
  lvplot::geom_lv(width.method = "height",
                  color = "#404040",
                  fill = "#CB566A",
                  outlier.colour = "#404040",
                  outlier.size = 1, outlier.stroke = 0.2, outlier.shape = 21,  varwidth=FALSE, size = 0.2) +
  ggpubr::rotate() +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)
  ) +
  xlab("") +
  ylab("SPARCC correlation")

# Save the letter-value plot
cowplot::save_plot(filename = "figures/osd2014_planctomycetes_lvplot.pdf", plot = plancto_lv, base_width = 7, base_height = 4)

# Create plots and stats for each Phyla associated to any Planctomycetes
plots <- lapply(V(selegoG_tax)$Phylum %>% unique(), plot_com_tax,
                G = selegoG, physeq = osd2014_plancto_phyloseq_beta,
                physeq_all = osd2014_dada2_phyloseq_beta,
                ord = osd2014_plancto_pcoa)
names(plots) <- V(selegoG_tax)$Phylum %>% unique()

# COmbine and plot
p <- ggpubr::ggarrange(plotlist = plots[1:5], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_1.pdf"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[6:10], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_2.pdf"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[11:15], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_3.pdf"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[16:20], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_4.pdf"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[21:25], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_5.pdf"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[26:29], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_5.pdf"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)

p <- ggpubr::ggarrange(plotlist = plots[1:5], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_1.png"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[6:10], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_2.png"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[11:15], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_3.png"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[16:20], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_4.png"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[21:25], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_5.png"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)
p <- ggpubr::ggarrange(plotlist = plots[26:29], common.legend = FALSE, ncol = 1, nrow = 5)
fname <- file.path("figures", paste0("osd2014_planctomycetes_com_plot_5.png"))
cowplot::save_plot(fname, p, base_height = 20, base_width = 11, ncol = 1, nrow = 1)



# Investigate Plantomycetes graph -----------------------------------------
G <- selegoG
tax <- "Planctomycetes"

# Create a subgraph for the Planctomycetes
G_sel <- G %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  dplyr::mutate(name = asv) %>%
  inner_join(as(tax_table(osd2014_plancto_phyloseq_beta), "matrix") %>% as.data.frame() %>% as_tibble(rownames = "asv")) %>%
  activate(edges) %>%
  mutate(Phylum_from = .N()$Phylum[from],
         Phylum_to = .N()$Phylum[to]) %>%
  activate(nodes)

# Preapare for collapsing the nodes
g1 <- G_sel %>%
  activate(nodes) %>%
  dplyr::mutate(comb = paste0(Phylum, com, sep = "|")) %>%
  dplyr::select(name, com, Phylum, comb) %>%
  activate(edges) %>%
  dplyr::select(-sign, -Phylum_from, -Phylum_to)

# Collapse the nodes
g.c <- contract.vertices(as.igraph(g1), as.factor(V(as.igraph(g1))$comb),
                         vertex.attr.comb=list(name= function(x) paste(x, collapse="|"),
                                               com= function(x) unique(x),
                                               comb = function(x) unique(x),
                                               Phylum = function(x) unique(x))
)

# Merge the edges (median weight)
g.c <- igraph::simplify(g.c, edge.attr.comb = list(weight = function(x) median(x),
                                                   weight_orig =  function(x) median(x),
                                                   pvalue = function(x) median(x)), remove.loops = TRUE) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  dplyr::mutate(sign = ifelse(weight_orig < 0, "NEG", "POS"))


# Let's plot the graph
sign_color <- c("#404040", "#B7144B")
names(sign_color) <- c("POS", "NEG")

class_colors <- c('#111111', '#65ADC2', '#233B43',
                  '#E84646', '#C29365', '#362C21',
                  '#316675', '#168E7F', '#109B37')

names(class_colors) <- G_sel %>% activate(nodes) %>% as_tibble() %>% .$Class %>% unique()

# Get the graph layout of the graph so it will be used for all the different plots
layout <- create_layout(G_sel, layout = 'fr')

# Plot the graph for Class
g_class <- ggraph(layout) +
  geom_edge_link0(aes(edge_width = weight, color = sign)) +
  geom_node_point(aes(fill = Class), size = 2, shape = 21, color = "black") +
  scale_edge_alpha(guide = FALSE) +
  scale_edge_width_continuous(range = c(0.2, 1.2)) +
  scale_edge_color_manual(values = sign_color) +
  scale_color_manual(values = class_colors) +
  theme_graph(base_family = 'Helvetica', base_size = 9) +
  theme(legend.position = "none",
        legend.key.width = unit(2, units = "mm"),
        legend.key.height = unit(5*length(labels), units = "mm"),
        legend.box = "vertical",
        text = element_text(size = 9)) +
  coord_fixed() +
  ggtitle("Class")

# Plot the graph with communities
g_com <- ggraph(layout) +
  geom_edge_link0(aes(edge_width = weight, color = sign))+
  geom_node_point( aes(fill = com), size = 2, shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = shapes, guide = guide_legend(direction = "vertical",
                                                           title.position = 'top',
                                                           reverse = FALSE,
                                                           label.position = "bottom")) +
  scale_edge_alpha(guide = FALSE) +
  scale_edge_width_continuous(range = c(0.2, 1.2)) +
  scale_edge_color_manual(values = sign_color) +
  theme_graph(base_family = 'Helvetica', base_size = 9) +
  theme(legend.position = "none",
        legend.key.width = unit(2, units = "mm"),
        legend.key.height = unit(5*length(labels), units = "mm"),
        #legend.justification=c(0, 1),
        legend.box = "vertical",
        text = element_text(size = 9)) +
  coord_fixed() +
  ggtitle("Louvain communities")

# Combine plots
p_g_class_com <- plot_grid(g_class, g_com, ncol = 1)

# Get proportions
physeq_prop <- transform_sample_counts(osd2014_dada2_phyloseq_beta, function(x) x/sum(x))


# Get the count data
graph_prop <- dplyr::as_data_frame(G_sel %>% activate(nodes)) %>%
  group_by(Class) %>%
  add_count() %>%
  ungroup() %>%
  inner_join(qpsmelt(physeq_prop) %>% dplyr::rename(asv = OTU))

graph_prop_agg <- graph_prop %>%
  dplyr::select(label, com, Abundance, Class, n) %>%
  dplyr::mutate(label = fct_relevel(label, rev(samples_sorted$label)))

plot_planct <- graph_prop_agg %>%
  dplyr::mutate(Class = paste0(Class, " (", n, ")")) %>%
  dplyr::select(-n) %>%
  #filter(counts >=5) %>%
  group_by(com,label, Class) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>%
  #filter(Abundance >0) %>%
  #mutate(Abundance = ifelse(Abundance == 0, 0.00001, Abundance)) %>%
  #droplevels() %>%
  ungroup() %>%
  ggplot(aes(label, Abundance, fill = com)) +
  geom_col(width = 1, color = "#404040") +
  facet_wrap(~Class, ncol = 6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  #theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        text = element_text(size = 9),
        legend.key.size = unit(0.5, "cm")) +
  ylab("Proportion") +
  xlab("OSD samples") +
  ggpubr::rotate() +
  scale_fill_manual(values = colors, name = NULL, guide = TRUE) +
  guides(fill = guide_legend(nrow = 1))

# Combina grah plots with the barplots
p_g_comb <- plot_grid(p_g_class_com, plot_grid(plot_planct, NULL, ncol = 1, rel_heights = c(0.6, 0.4)), ncol = 2, rel_widths = c(0.3, 0.7))

# Get legend
g_com_legend <- get_legend(g_com + theme_bw() + theme(legend.position = "top", text = element_text(size = 11)))

# Create graphs by community highlighting the Class of the nodes
g_com_list <- lapply(sort(G_sel %>% activate(nodes) %>% as_tibble() %>% .$com %>% unique()), map_coms, layout = layout)

# Combine plots
plot_g_com_list <- plot_grid(
  plotlist = g_com_list,
  align = 'vh',
  ncol = 3,
  nrow = 3
)

# Get legend
plot_tmp <- ggraph(G_sel) +
  geom_node_point(aes(fill = Class), shape = 21, color = "#404040", size = 3) +
  scale_fill_manual(values = class_colors)

legend <- get_legend(
  # create some space to the left of the legend
  plot_tmp +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)
plot_g_com <- plot_grid(legend, plot_g_com_list, ncol = 1, rel_heights = c(.1, 1))

# save plots
save_plot(p_g_comb, filename = "figures/osd2014_plancto_Planctomycetes_Class_com_graphs.pdf", base_height = 8, base_width = 9.5)
save_plot(plot_g_com, filename = "figures/osd2014_plancto_Planctomycetes_Class_subgraphs.pdf", base_height = 11, base_width = 11)

# save legend
save_plot(g_com_legend , filename = "figures/osd2014_plancto_Planctomycetes_Class_com_graphs_legend.pdf", base_height = 8, base_width = 11)

