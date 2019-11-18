surf <- function(X, Y) {
  ordi <- vegan::ordisurf(X, Y, plot = FALSE, bs="ds")
  ordi.grid <- ordi$grid #extracts the ordisurf object
  #str(ordi.grid) #it's a list though - cannot be plotted as is
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
  ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
  ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
  ordi.mite.na
}

generate_grid <- function(max_lat = -90, min_lat = 90, max_lon = 195, min_lon = -180, ticks = 100){
  require(maps)
  grid <- expand.grid(seq(min_lon, max_lon, length = ticks),
                      seq(min_lat, max_lat, length = ticks))
  colnames(grid) <- c("lon", "lat")
  rownames(grid) <- 1:nrow(grid)
  return(grid)
}

plot_dispersion <- function(betadspr, title = NULL){
  library(ggConvexHull)
  # extract the centroids and the site points in multivariate space.
  centroids <- data.frame(group = rownames(betadspr$centroids), data.frame(betadspr$centroids)) %>% as_tibble()
  vectors <- data.frame(group = betadspr$group, data.frame(betadspr$vectors)) %>% as_tibble()


  # to create the lines from the centroids to each point we will put it in a format that ggplot can handle
  seg.data <- centroids %>%
    as_tibble() %>%
    dplyr::select("group","PCoA1","PCoA2") %>%
    inner_join(vectors %>% as_tibble() %>% dplyr::select("group","PCoA1","PCoA2") %>% setNames(c("group", "v.PCoA1","v.PCoA2"))) %>%
    dplyr::mutate(group = fct_inorder(group))
  lvls <- levels(seg.data$group)
  seg.data <- seg.data %>%
    dplyr::mutate(group1 = group,
                  group = fct_relevel(group, after = 1)) %>%
    bind_rows(seg.data %>% dplyr::mutate(group1 = "All") ) %>%
    dplyr::mutate(group1 = fct_relevel(group1, c("All", lvls)))
  ggthemr::ggthemr_reset()

  ratio_values <- get_ratio(x = seg.data$v.PCoA1, y = seg.data$v.PCoA2)
  ggplot() +
    geom_segment(data = seg.data, aes(x= v.PCoA1, xend= PCoA1, y= v.PCoA2, yend= PCoA2, group = group), alpha=0.30) +
    geom_convexhull(data = seg.data, aes(x = v.PCoA1, y = v.PCoA2, group = group), alpha = 0, color = "black", linetype="dashed", size = 0.3) +
    geom_point(data = seg.data, aes(x = PCoA1, y = PCoA2, group = group), size = 3, colour = "black", fill = "#B7144B", shape=21) +
    geom_point(data = seg.data, aes(x= v.PCoA1, y= v.PCoA2, group = group), size = 1, colour = "black", fill = "#404040", shape=21) +
    facet_wrap(~group1) +
    theme_bw() +
    coord_fixed(ratio = ratio_values) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          text = element_text(face = "bold", size = 8)) +
    ggtitle(title)
}

get_ratio <- function(x = x, y = y, display = 4/3){
  ratio_display <- display
  ratio_values <- ((max(x) - min(x)))/((max(y) - min(y)))
  ratio_values/ratio_display
}

get_ratio_log10 <- function(x = x, y = y, display = 4/3){
  ratio_display <- display
  ratio_values <- (log10(max(x) - min(x)))/(log10(max(y) - min(y)))
  ratio_values/ratio_display
}

write_beta <- function(beta_permut = beta_permut, beta_anova = beta_anova, permanova = permanova, label = label){

  xls_fname <- file.path(paste0("results/osd2014_planctomycetes_beta_", label, ".xlsx"))

  beta_permut$tab %>%
    as.data.frame() %>%
    write.xlsx(file = xls_fname, sheetName = "BETADIV PERMUT", row.names = TRUE)

  beta_anova %>%
    as.data.frame() %>%
    write.xlsx(file = xls_fname, sheetName = "BETADIV ANOVA", row.names = TRUE, append = TRUE)

  permanova %>%
    as.data.frame() %>%
    write.xlsx(file = xls_fname, sheetName = "PERMANOVA MP", row.names = TRUE, append = TRUE)
}


# Functions for graph analysis --------------------------------------------
qpsmelt <- function(X) {
  if (taxa_are_rows(X)) {
    count_table <- as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "OTU") %>%
      gather(label, Abundance, -OTU)
  }else{
    count_table <-as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "label") %>%
      gather(OTU, Abundance, -label)
  }
  sample_table <- as(sample_data(X), "matrix") %>%
    as_tibble()
  taxa_table <- as(tax_table(X), "matrix") %>%
    as_tibble(rownames = "OTU")

  count_table %>%
    left_join(sample_table) %>%
    left_join(taxa_table)
}

plot_com_tax <- function(X, G = G, physeq = physeq, physeq_all = physeq_all, ord = ord){

  tax <- X
  cat(paste("Processing", X, "\n"))

  tax_sel <- as(tax_table(physeq_all), "matrix") %>%
    as_tibble(rownames = "asv") %>%
    dplyr::select(asv, Phylum, Class, Order, Family, Genus, asv_name) %>%
    filter(Phylum == "Planctomycetes" | Phylum == tax)


  if(tax != "Planctomycetes"){
    G_sel <- G %>%
      as_tbl_graph() %>%
      activate(nodes) %>%
      filter(asv %in% tax_sel$asv) %>%
      dplyr::mutate(name = asv) %>%
      inner_join(as(tax_table(physeq_all), "matrix") %>% as.data.frame() %>% as_tibble(rownames = "asv")) %>%
      activate(edges) %>%
      mutate(Phylum_from = .N()$Phylum[from],
             Phylum_to = .N()$Phylum[to]) %>%
      filter(Phylum_from != Phylum_to) %>%
      activate(nodes) %>%
      filter(!node_is_isolated())
  }else{
    G_sel <- G %>%
      as_tbl_graph() %>%
      activate(nodes) %>%
      filter(asv %in% tax_sel$asv) %>%
      dplyr::mutate(name = asv) %>%
      inner_join(as(tax_table(physeq_all), "matrix") %>% as.data.frame() %>% as_tibble(rownames = "asv")) %>%
      activate(edges) %>%
      mutate(Phylum_from = .N()$Phylum[from],
             Phylum_to = .N()$Phylum[to]) %>%
      activate(nodes) %>%
      filter(!node_is_isolated())
  }


  l_p <- transform_sample_counts(physeq_all, function(x) x/sum(x))


  l <- G_sel %>%
    activate(nodes) %>%
    as_tibble() %>%
    group_by(Phylum) %>%
    add_count() %>%
    ungroup() %>%
    inner_join(qpsmelt(l_p) %>% dplyr::rename(asv = OTU)) %>%
    filter(label %in% sample_names(physeq))

  samples_sorted <- ord$vectors %>%
    as_tibble(rownames = "label")  %>%
    dplyr::select(label, Axis.1, Axis.2) %>%
    dplyr::arrange(Axis.1, -Axis.2)

  ggthemr::ggthemr(layout = "scientific", palette = "fresh")

  if( tax == "Planctomycetes"){
    l_pl <- l %>%
      dplyr::select(label, com, Abundance, Phylum, n) %>%
      dplyr::mutate(label = fct_relevel(label, samples_sorted$label),
                    com = fct_relevel(com, names(colors)),
                    Phylum = paste0(Phylum, " (", n, ")")) %>%
      dplyr::select(-n) %>%
      #filter(counts >=5) %>%
      group_by(com, label, Phylum) %>%
      dplyr::summarise(Abundance = sum(Abundance)) %>%
      ungroup() %>%
      ggplot(aes(label, Abundance, fill = com)) +
      geom_bar(stat = "identity", width = 1) +
      facet_wrap(~Phylum, nrow = 2) +
      scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "top",
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm")) +
      ylab("Proportion") +
      xlab("OSD samples") +
      scale_fill_manual(values = colors, name = NULL)
  }else{
    l_pl <- l %>%
      dplyr::select(label, com, Abundance, Phylum, n) %>%
      dplyr::mutate(Phylum = fct_relevel(Phylum, c(tax, "Planctomycetes")),
                    Phylum1 = paste0(Phylum, " (", n, ")"))
    taxs <- l_pl %>% dplyr::select(Phylum, Phylum1) %>% unique() %>% arrange(as.factor(Phylum))
    l_pl <- l_pl %>%
      dplyr::mutate(Phylum1 = fct_relevel(Phylum1, taxs$Phylum1),
                    label = fct_relevel(label, samples_sorted$label),
                    com = fct_relevel(com, names(colors))) %>%
      dplyr::select(-n) %>%
      #filter(counts >=5) %>%
      group_by(com, label, Phylum1) %>%
      dplyr::summarise(Abundance = sum(Abundance)) %>%
      ungroup() %>%
      ggplot(aes(label, Abundance, fill = com)) +
      geom_bar(stat = "identity", width = 1, color = "#404040") +
      facet_wrap(~Phylum1, nrow = 2, scales = "free") +
      scale_y_continuous(labels = scales::percent) +
      #theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "top",
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm")) +
      ylab("Proportion") +
      xlab("OSD samples") +
      scale_fill_manual(values = colors, name = NULL)
  }

  g1 <- G_sel %>%
    activate(nodes) %>%
    dplyr::mutate(comb = paste0(Phylum, com, sep = "|")) %>%
    dplyr::select(name, com, Phylum, comb) %>%
    activate(edges) %>%
    dplyr::select(-sign, -Phylum_from, -Phylum_to)

  g.c <- contract.vertices(as.igraph(g1), as.factor(V(as.igraph(g1))$comb),
                           vertex.attr.comb=list(name= function(x) paste(x, collapse="|"),
                                                 com= function(x) unique(x),
                                                 comb = function(x) unique(x),
                                                 Phylum = function(x) unique(x))
  )

  g.c <- igraph::simplify(g.c, edge.attr.comb = list(weight = function(x) median(x),
                                                     weight_orig =  function(x) median(x),
                                                     pvalue = function(x) median(x)), remove.loops = TRUE) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    dplyr::mutate(sign = ifelse(weight_orig < 0, "NEG", "POS"))

  shapes <- c(21, 22)
  names(shapes) <- c("Planctomycetes", tax)
  sign_color <- c("#404040", "#B7144B")
  names(sign_color) <- c("POS", "NEG")


  p_g <- ggraph::ggraph(g.c, layout = "fr") +
    geom_edge_link0(aes(edge_width = weight, color = sign))+
    #geom_edge_fan(aes(color = sign, width = weight_orig), show.legend = TRUE) +
    geom_node_point(aes(fill = com, shape = Phylum), size = 3, color = "black") +
    scale_fill_manual(values = colors, guide = FALSE) +
    # scale_edge_color_gradientn(colours = rev(pal),
    #                            limits = c(-1,1),
    #                            breaks = c(-0.5, 0, 0.5),
    #                            guide = guide_edge_colourbar(
    #                              direction = "horizontal",
    #                              title.position = 'top',
    #                              reverse = FALSE,
    #                              label.position = "bottom"
    #                            ), name = "SPARCC correlation") +
    scale_shape_manual(values = shapes, guide = guide_legend(direction = "vertical",
                                                             title.position = 'top',
                                                             reverse = FALSE,
                                                             label.position = "bottom")) +
    scale_edge_alpha(guide = FALSE) +
    scale_edge_width_continuous(range = c(0.2, 1.2)) +
    scale_edge_color_manual(values = sign_color) +
    theme_graph(base_family = 'Helvetica', base_size = 9) +
    theme(legend.position = "right",
          legend.key.width = unit(2, units = "mm"),
          legend.key.height = unit(5*length(labels), units = "mm"),
          #legend.justification=c(0, 1),
          legend.box = "vertical") +
    coord_fixed()

  #p_g + l_pl

  # Write graph data to an excel file
  library(xlsx)

  # In each xlsx we want:
  # Plancto nodes
  # Phylum nodes
  # All connections
  ntax <- gsub("[^[:alnum:] ]", "_", tax)
  xls_fname <- file.path("results/", paste0("osd2014_graph_", ntax, ".xlsx"))
  G_sel %>%
    activate(nodes) %>%
    as_tibble() %>%
    dplyr::select(asv_name, asv, Phylum, Class, Order, Family, Genus) %>%
    arrange(Class) %>%
    as_tibble() %>%
    filter(Phylum == "Planctomycetes") %>%
    as.data.frame() %>%
    write.xlsx(file=xls_fname, sheetName = "Planctomycetes nodes", row.names = FALSE)

  if (tax != "Planctomycetes"){
    G_sel %>%
      activate(nodes) %>%
      as_tibble() %>%
      dplyr::select(asv_name, asv, Phylum, Class, Order, Family, Genus) %>%
      as_tibble() %>%
      filter(Phylum == tax) %>%
      arrange(Class) %>%
      as.data.frame() %>%
      write.xlsx(file=xls_fname, sheetName=paste(tax, "nodes"), row.names=FALSE, append = TRUE)
  }

  G_sel %>%
    dplyr::select(asv_name, asv, Phylum, Class, Order, Family, Genus) %>%
    dplyr::mutate(name = asv_name) %>%
    as_data_frame(what = "edges") %>%
    write.xlsx(file=xls_fname, sheetName=paste("Planctomycetes", tax, "graph"), row.names=FALSE, append = TRUE)

  cowplot::plot_grid(p_g, l_pl, nrow = 1, ncol = 2, align = "v", axis = "t")
}

map_coms <- function(f_com, layout = layout){

  v_com <-  G_sel %>% activate(nodes) %>%
    filter(com == f_com) %>%
    activate(nodes) %>%
    as_tibble()


  ggraph(layout) +
    geom_edge_link0(aes(edge_width = weight, color = sign))+
    #geom_edge_fan(aes(color = sign, width = weight_orig), show.legend = TRUE) +
    geom_node_point(fill = "#C0C2C2", shape = 21, color = "#404040", size = 2) +
    #geom_node_point(aes(color = Class), size = 3, alpha = 0.7) +
    geom_node_point(data = layout %>% filter(name %in% v_com$name), aes(fill = Class), size = 2, shape = 21, color = "black") +

    #scale_fill_manual(values = colors, guide = FALSE) +
    # scale_edge_color_gradientn(colours = rev(pal),
    #                            limits = c(-1,1),
    #                            breaks = c(-0.5, 0, 0.5),
    #                            guide = guide_edge_colourbar(
    #                              direction = "horizontal",
    #                              title.position = 'top',
    #                              reverse = FALSE,
    #                              label.position = "bottom"
    #                            ), name = "SPARCC correlation") +
  scale_shape_manual(values = shapes, guide = guide_legend(direction = "vertical",
                                                           title.position = 'top',
                                                           reverse = FALSE,
                                                           label.position = "bottom")) +
    scale_edge_alpha(guide = FALSE) +
    scale_edge_width_continuous(range = c(0.2, 1.2)) +
    scale_edge_color_manual(values = sign_color) +
    scale_fill_manual(values = class_colors) +
    theme_graph(base_family = 'Helvetica', base_size = 9) +
    theme(legend.position = "none",
          #legend.key.height = unit(2, units = "mm"),
          #legend.key.width = unit(10*length(labels), units = "mm"),
          #legend.justification=c(0, 1),
          legend.box = "horizontal",
          text = element_text(size = 9)) +
    ggtitle(f_com) + coord_fixed()
}

