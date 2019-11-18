library(tidyverse)
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(vegan)
library(ape)
library(raster)
library(geoR)
library(gstat)
library(rgeos)
library(wesanderson)
library(coenocliner)
library(broom)
library(pals)
library(cowplot)
library(ggpubr)
library(xlsx)

# Load necessary libraries for TINA/PINA calculation
source("libs/functions_com_sim.R")
source("libs/libs.R")

# Load SPARCC matrices
load("data/osd2104_global_TINA_PINA.Rda")
load("data/osd2014_meow_regions.Rda")

# Load phyloseq objects
#load("~/Desktop/osd2014_analyses_repo/osd2014_analysis/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata")

#load("~/Desktop/osd2014_analyses_repo/osd2014_analysis/osd2014_16S_asv/data//osd2014_basic_cdata.Rdata", verbose = TRUE)
load(file = "data/osd2014_cdata.Rda")
load(file = "data/osd2014_amp_mg_intersect.Rda")
load(file = "data/osd2014_dada2_phyloseq_alpha.Rda")
load(file = "data/osd2014_dada2_phyloseq_beta.Rda")


# Add some modified meow regions
osd2014_cdata <- osd2014_cdata %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% left_join(osd2014_meow_regions)

# Objects used for alpha and bar plotting
osd2014_plancto_asvs_alpha <- as(tax_table(osd2014_dada2_phyloseq_alpha), "matrix") %>%
  as_tibble(rownames = "asv") %>%
  filter(Phylum == "Planctomycetes")

osd2014_plancto_phyloseq_alpha <- prune_taxa(osd2014_plancto_asvs_alpha$asv, osd2014_dada2_phyloseq_alpha)
osd2014_plancto_phyloseq_alpha <- prune_samples(sample_sums(osd2014_plancto_phyloseq_alpha) > 0, osd2014_plancto_phyloseq_alpha)

# Objects used for beta and bar plotting
osd2014_plancto_asvs_beta <- as(tax_table(osd2014_dada2_phyloseq_beta), "matrix") %>%
  as_tibble(rownames = "asv") %>%
  filter(Phylum == "Planctomycetes")

osd2014_plancto_phyloseq_beta <- prune_taxa(osd2014_plancto_asvs_beta$asv, osd2014_dada2_phyloseq_beta)
osd2014_plancto_phyloseq_beta <- prune_samples(sample_sums(osd2014_plancto_phyloseq_beta) > 0, osd2014_plancto_phyloseq_beta)

# Get how many ASVs and sequence we have in each sample
osd2014_plancto_asvXsample <- as(otu_table(osd2014_plancto_phyloseq_beta), "matrix") %>%
  as_tibble(rownames = "label") %>%
  gather(key = asv, value = "counts", -label) %>%
  filter(counts > 0) %>%
  group_by(label) %>%
  dplyr::summarise(n = n(), N = sum(counts)) %>%
  ungroup()


# Plot the number of ASVs and sequence we have in each sample
kelly.colours <- kelly()

names(kelly.colours) <- c("Removed", sort(unique(osd2014_cdata$meow_province)))

osd2014_plancto_asvXsample <- osd2014_plancto_asvXsample %>%
  inner_join(osd2014_cdata %>% dplyr::select(label, meow_province)) %>%
  mutate(keep = ifelse((n > 3 | N >= 50), meow_province, "Removed"),
         keep = fct_relevel(keep, c("Removed", sort(unique(osd2014_cdata$meow_province)))))



ratio_values <- get_ratio_log10(x = osd2014_plancto_asvXsample$N, y = osd2014_plancto_asvXsample$n)

plot_osd2014_plancto_asvXsample <- ggplot(osd2014_plancto_asvXsample, aes(N, n, fill = keep)) +
  geom_point(size = 3, shape = 21, color = "black") +
  scale_y_log10() +
  scale_x_log10(labels = scales::comma) +
  xlab("Sequence counts") +
  ylab("Number of ASVs") +
  scale_fill_manual(values = kelly.colours, name = "") +
  theme_bw()


osd2014_plancto_asvXsample_legend <- cowplot::get_legend(plot_osd2014_plancto_asvXsample +
                                                           theme(panel.grid = element_blank(),
                                                                 legend.position = "right") +
                                                           guides(fill = guide_legend(ncol= 2)))

osd2014_plancto_asvXsample_legend <- as_ggplot(osd2014_plancto_asvXsample_legend)

plot_osd2014_plancto_asvXsample <- plot_osd2014_plancto_asvXsample +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  coord_fixed(ratio = ratio_values)

save_plot(filename = "figures/osd2014_plancto_asvXsample_legend.pdf", osd2014_plancto_asvXsample_legend)

# We remove samples with only not seen in more than 3 samples and at least with 50 sequences
osd2014_plancto_phyloseq_beta <- prune_samples(osd2014_plancto_asvXsample %>% filter((n > 3 | N >= 50)) %>% .$label,
                                               osd2014_plancto_phyloseq_beta)



# Start TINA/PINA inference -----------------------------------------------

# Get TINA
PARAM <- list();
PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$sample.steps <- c(0.01, 0.02, 0.05, 0.1, 0.15, seq(0.2, 0.9, by=0.1));
PARAM$use.cores <- 4;

#Set parameters
size.thresh <- 1;
pseudocount <- 10^-6;
nblocks <- 400;
use.cores <- 4;

# We have three different matrices for TINA/PINA
# S.sparcc: as calculated by TINA scripts
# S.phylo: as calculated by PINA scripts
# S.sparcc.filt: SPARCC graph filtered in a way that weak links are removed
#                but the network integrity remains intact

S.sparcc <- S.sparcc[osd2014_plancto_asvs_beta$asv, osd2014_plancto_asvs_beta$asv]
S.phylo <- S.phylo[osd2014_plancto_asvs_beta$asv, osd2014_plancto_asvs_beta$asv]
S.sparcc.filt <- S.sparcc.filt[osd2014_plancto_asvs_beta$asv, osd2014_plancto_asvs_beta$asv]

############################
#Preallocate classical indices to calculate
get.cs <- list(); t <- 1;

############################
#Preallocate corrected indices to calculate
t <- start.t <- length(get.cs) + 1;
#Jaccard index, SparCC-corrected, unweighted, normalized
#=> unweighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, unweighted";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, SparCC-corrected, weighted, normalized
#=> weighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, weighted";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;
#Jaccard index, SparCC-corrected-filtered, unweighted, normalized
#=> unweighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, unweighted filtered";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, SparCC-corrected-filtered, weighted, normalized
#=> weighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, weighted filtered";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;
#Jaccard index, phylo-corrected, unweighted, normalized
#=> unweighted PINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "PINA, unweighted";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, phylo-corrected, weighted, normalized
get.cs[[t]] <- list();
get.cs[[t]]$name <- "PINA, weighted";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;

o.t <- t(as(object = otu_table(osd2014_plancto_phyloseq_beta, taxa_are_rows = TRUE), "matrix"))
ot <- t(as(object = otu_table(osd2014_plancto_phyloseq_beta, taxa_are_rows = TRUE), "matrix"))
#ot <- ot[plancto_asvs$asv,]

############################
#Iterate through (corrected) indices and calculate pairwise community similarities
for (t in seq(start.t, length(get.cs))) {
  print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
  #Calculate all pairwise similarities
  if (get.cs[[t]]$name %in% c("TINA, unweighted", "TINA, weighted")) {
    curr.cs <- community.similarity.corr.par(ot, S = S.sparcc, distance = get.cs[[t]]$call, blocksize = 10, use.cores = PARAM$use.cores)
  } else if (get.cs[[t]]$name %in% c("PINA, unweighted", "PINA, weighted")) {
    curr.cs <- community.similarity.corr.par(ot, S = S.phylo, distance = get.cs[[t]]$call, blocksize = 1000, use.cores = PARAM$use.cores)
  } else if (get.cs[[t]]$name %in% c("TINA, unweighted filtered", "TINA, weighted filtered")) {
    curr.cs <- community.similarity.corr.par(ot, S = S.sparcc.filt, distance = get.cs[[t]]$call, blocksize = 1000, use.cores = PARAM$use.cores)
  }

  #Correct for rounding errors
  curr.cs[curr.cs < 0] <- 0;
  #Export histogram
  curr.data <- curr.cs[upper.tri(curr.cs)];
  curr.plot <- ggplot(data.frame(data=sample(curr.data, 2000)), aes(x = data)) +
    geom_density(alpha = 0.2, fill = "green") +
    xlim(0,1) +
    ggtitle(paste("Distribution of Pairwise", get.cs[[t]]$name, "Distances")) +
    ylab("Density");
  #ggsave(curr.plot, width=10, height=10, filename=paste(PARAM$folder.output, "hmp_samples.similarity_hist.", get.cs[[t]]$name, ".pdf", sep = ""), useDingbats=F);
  #Store current similarities
  get.cs[[t]]$cs <- curr.cs;
  get.cs[[t]]$plot <- curr.plot
  #Tidy up
  rm(curr.data, curr.cs);
  print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
  gc();
}

osd2014_plancto_pina_tina_results <- get.cs


# Beta diversity analyses -------------------------------------------------

# Calculate PCoA
osd2014_plancto_pcoa <- pcoa(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs))
# We need to run it with the stats version to be able to get the ordisurf
ord <- cmdscale(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs))

# Scree plot
osd2014_plancto_pcoa_eig <- osd2014_plancto_pcoa$values$Relative_eig %>%
  enframe(name = "PC", value = "prop") %>%
  mutate(PC = as.character(PC), PC = fct_reorder(.desc = TRUE, PC, prop))

osd2014_plancto_pcoa_scree <- osd2014_plancto_pcoa_eig %>%
  top_n(n = 10, wt = prop) %>%
  mutate(PC = fct_inorder(PC)) %>%
  ggplot(aes(x = PC, y = prop)) +
  geom_col(fill = "indianred", color = "black") +
  theme_light() +
  scale_y_continuous(labels = scales::percent) +
  xlab("Principal components") +
  ylab("Percent explained")

# Filter the contextual data to the samples used
osd2014_cdata_filt <- osd2014_cdata %>%
  filter(label %in% rownames(osd2014_plancto_pcoa$vectors)) %>%
  mutate(wt_cat = cut(water_temperature, breaks = 6),
         wt_cat = fct_inorder(wt_cat))
osd2014_cdata_filt_df <- osd2014_cdata_filt %>%
  as.data.frame() %>%
  column_to_rownames("label")

ordisurf_data <- surf(ord, osd2014_cdata_filt_df[rownames(ord),]$water_temperature)

# Data for plotting the PCoA
data4pcoa <- osd2014_plancto_pcoa$vectors %>%
  as_tibble(rownames = "label") %>%
  dplyr::select(label, Axis.1, Axis.2) %>%
  inner_join(osd2014_cdata %>% dplyr::select(label, start_lon, start_lat, water_temperature, longhurst_biome, meow_province)) %>%
  dplyr::rename(lon = start_lon, lat = start_lat)

# Get colors
pal <- (viridis::magma(8, alpha = 0.8)[2:7])

ratio_values <- get_ratio(x = data4pcoa$Axis.1, y = data4pcoa$Axis.2)

# Plot the PCoA
plot_pcoa <- ggplot() +
  stat_contour(data = ordisurf_data, aes(x = x, y = y, z = z, color = ..level..),  binwidth = 1) +
  geom_point(data = data4pcoa, aes(x = Axis.1, y = Axis.2, fill = meow_province), shape = 21, size = 3, color = "black") +
  theme_bw() +
  xlab(paste0("PC1 (", scales::percent(osd2014_plancto_pcoa_eig[1,]$prop), ")")) +
  ylab(paste0("PC2 (", scales::percent(osd2014_plancto_pcoa_eig[2,]$prop), ")")) +
  scale_fill_manual(name = "Longhurst biome",
                    values = kelly.colours, #c("#B7144B", "#F09735", "#1F629A", "#404040"),
                    guide = FALSE) +
  scale_color_gradientn(colours = pal,
                        name = "Water temperature (C)",
                        guide = guide_colorbar(
                          direction = "horizontal",
                          title.position = 'top',
                          reverse = FALSE,
                          label.position = "bottom"
                        ))

# Get legend
plot_pcoa_legend <- cowplot::get_legend(plot_pcoa +
                                          theme(legend.position= "top",
                                                legend.key.height = unit(2, units = "mm"),
                                                legend.key.width = unit(10*length(labels), units = "mm"),
                                                #legend.justification=c(0, 1),
                                                legend.box = "horizontal",
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank()))
# Plot legend
plot_pcoa_legend <- as_ggplot(plot_pcoa_legend)
save_plot(filename = "figures/osd2014_plot_pcoa_legend_legend.pdf", plot_pcoa_legend)


plot_pcoa <- plot_pcoa +
  theme(legend.position = "none",
        legend.key.height = unit(2, units = "mm"),
        legend.key.width = unit(10*length(labels), units = "mm"),
        #legend.justification=c(0, 1),
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1.15*(ratio_values))



# Interpolation of the PCoA PC1 -------------------------------------------

########################################
#Function to generate grids
########################################

# Get the data from the ordination
osd2014_cdata_filt <- osd2014_cdata %>%
  filter(label %in% rownames(osd2014_plancto_pcoa$vectors))

# Prepare data for the interpolation
idw_data <- osd2014_plancto_pcoa$vectors %>%
  as_tibble(rownames = "label") %>%
  dplyr::select(label, Axis.1) %>%
  inner_join(osd2014_cdata %>% dplyr::select(label, start_lon, start_lat)) %>%
  dplyr::rename(lon = start_lon, lat = start_lat) %>%
  as.data.frame() %>%
  column_to_rownames("label")


# Create the grid where we will interpolate the values
max_lon <- ifelse(round(max(osd2014_cdata_filt$start_lon) + 10) > 180, 180, round(max(osd2014_cdata_filt$start_lon) + 10))
min_lon <- ifelse(round(min(osd2014_cdata_filt$start_lon) - 10) < -180, -180, round(min(osd2014_cdata_filt$start_lon) - 10))
max_lat <- ifelse(round(max(osd2014_cdata_filt$start_lat) + 10) > 90, 90, round(max(osd2014_cdata_filt$start_lat) + 10))
min_lat <- ifelse(round(min(osd2014_cdata_filt$start_lat) - 10) < -90, -90, round(min(osd2014_cdata_filt$start_lat) - 10))
idw_grid <- as.data.frame(generate_grid(max_lat = max_lat, min_lat = min_lat,
                                        max_lon = max_lon, min_lon = min_lon,
                                        ticks = 200))

# Crop the raster to be centered in the regions of interest
data("wrld_simpl", package = "maptools")
wm <- crop(gBuffer(wrld_simpl, byid = TRUE, width = 0), extent(min_lon, max_lon, min_lat, max_lat))

# Transform df into spatial data
dat <- idw_data
coordinates(dat) <- ~lon + lat
coordinates(idw_grid) <- ~lon + lat
gridded(idw_grid) <- TRUE


# Interpolation
idw.variable <- idw(formula = Axis.1 ~ 1, locations = dat, newdata = idw_grid, maxdist = 8, idp = 2)
idw.output = as.data.frame(idw.variable)
names(idw.output)[1:3] <- c("long", "lat", "var1.pred")

# Map the predicted values
# Some colors for the rasters
pal <- wes_palette("Zissou1", 100, type = "continuous")

# Get the raster from the predictions
idw.r <- rasterFromXYZ(idw.output[, c("long", "lat", "var1.pred")])

# Convert it to a df
idw.msk.dfr <- as.data.frame(rasterToPoints(idw.r))

ggthemr::ggthemr_reset()
# Let's plot
# Adapted from https://stackoverflow.com/questions/48816773/polar-stereographic-map-in-r
plot_idw <- ggplot() +
  geom_tile(data = idw.msk.dfr, aes(x = x, y = y, color = var1.pred), size = 1) +
  geom_polygon(data = wm, aes(long, lat, group = group), fill = '#262626', color = '#3A3A3A', size = 0.1) +
  geom_point(data = idw_data %>% rownames_to_column("label") %>% inner_join(osd2014_cdata) %>% as_tibble(),
             aes(x = lon, y = lat, fill = meow_province), size = 0.5, alpha = 1, shape = 21, color = "black") +
  # Convert to polar coordinates
  #coord_map("mercator") +
  coord_equal(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) +
  # Removes Axes and labels
  xlab("") +
  ylab("") +
  # Adds labels
  #geom_text(aes(x = 180, y = seq(min_lat, max_lat - 10, by = 10), hjust = -0.2, label = paste0(seq(min_lat, max_lat - 10, by = 10), "°N"))) +
  #geom_text(aes(x = x_lines, y = min_lat - 15, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +
  # Adds axes
  #geom_hline(aes(yintercept = min_lat - 5), size = 1)  +
  #geom_segment(aes(y = min_lat - 10, yend = 90, x = x_lines, xend = x_lines), linetype = "dashed") +
  scale_color_gradientn(colours = pal,breaks=c(-0.25,0,0.25),labels=c(-0.25,0,0.25),
                        limits=c(-0.5,0.5), name = "PC1") +
  scale_fill_manual(name = "Longhurst biome",
                    values = kelly.colours, #c("#B7144B", "#F09735", "#1F629A", "#404040"),
                    guide = FALSE) +
  # Change theme to remove axes and ticks
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "#f5f5f2", color = NA))

# Plot temperature against the PC1
pc1_temp_data <- osd2014_plancto_pcoa$vectors %>%
  as_tibble(rownames = "label") %>%
  dplyr::select(label, Axis.1) %>%
  inner_join(osd2014_cdata %>% dplyr::select(label, start_lon, start_lat, water_temperature)) %>%
  dplyr::rename(lon = start_lon, lat = start_lat) %>%
  filter(water_temperature > 0) %>%
  as.data.frame() %>%
  column_to_rownames("label")

ratio_values <- get_ratio(x = pc1_temp_data$water_temperature, y =pc1_temp_data$Axis.1)

plot_pc1_temp <- pc1_temp_data %>%
  ggplot(aes(water_temperature, Axis.1)) +
  geom_point(shape = 21, size = 3, fill = "#315C71", color = "black", alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", fill = "black", alpha = 0.2) +
  ylab("PC1") +
  xlab("Water temperature (C)") +
  theme_bw() +
  theme(legend.position= "none",
        legend.key.height = unit(2, units = "mm"),
        legend.key.width = unit(10*length(labels), units = "mm"),
        #legend.justification=c(0, 1),
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio = ratio_values)


beta_panel <- cowplot::plot_grid(plot_osd2014_plancto_asvXsample, plot_pcoa, plot_idw, plot_pc1_temp,
                                 align = "hv", axis = "t",  nrow = 2, ncol = 2)

save_plot(filename = "figures/osd2014_planctomycetes_beta_panel_raw.pdf", plot = beta_panel, base_width = 11, base_height = 8)

osd2014_plancto_pcoa$vectors %>%
  as_tibble(rownames = "label") %>%
  dplyr::select(label, Axis.1) %>%
  inner_join(osd2014_cdata %>% dplyr::select(label, start_lon, start_lat, water_temperature)) %>%
  dplyr::rename(lon = start_lon, lat = start_lat) %>%
  filter(water_temperature > 0) %>%
  do(fit = lm(Axis.1 ~ water_temperature, data = .)) %>%
  glance(fit)


# Betadisper --------------------------------------------------------------
set.seed(12345)
# Be sure the contextual data has the same order as the distance matrices
osd2014_cdata_filt_df <- osd2014_cdata_filt_df[rownames(osd2014_plancto_pina_tina_results[[2]]$cs),]

# Let's test the homogeneity of multivariate dispersions
# We have quite an unbalanced setup so might be difficult to be sure that difference we observe are owing to
# between or within group
betadisper_tina_mp <- betadisper(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs), osd2014_cdata_filt_df$meow_province, type = "median")
p_betadisper_tina_mp <- permutest(betadisper_tina_mp)
anova_betadisper_tina_mp <- anova(betadisper_tina_mp)
plot_betadisper_tina_mp <- plot_dispersion(betadisper_tina_mp, title = "MEOW provinces")
save_plot(filename = "figures/osd2014_planctomycetes_dispersion_mp.pdf", plot = plot_betadisper_tina_mp, base_width = 8, base_height = 8)


betadisper_tina_mr <- betadisper(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs), osd2014_cdata_filt_df$meow_realm, type = "median")
p_betadisper_tina_mr <- permutest(betadisper_tina_mr)
anova_betadisper_tina_mr <- anova(betadisper_tina_mr)
plot_betadisper_tina_mr <- plot_dispersion(betadisper_tina_mr, title = "MEOW realm")
save_plot(filename = "figures/osd2014_planctomycetes_dispersion_mr.pdf", plot = plot_betadisper_tina_mr, base_width = 8, base_height = 8)

betadisper_tina_wt <- betadisper(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs), osd2014_cdata_filt_df$wt_cat, type = "median")
p_betadisper_tina_wt <- permutest(betadisper_tina_wt)
anova_betadisper_tina_wt <- anova(betadisper_tina_wt)
plot_betadisper_tina_wt <- plot_dispersion(betadisper_tina_wt, title = "Water temperature")
save_plot(filename = "figures/osd2014_planctomycetes_dispersion_wt.pdf", plot = plot_betadisper_tina_wt, base_width = 8, base_height = 8)

betadisper_tina_lhc <- betadisper(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs), osd2014_cdata_filt_df$longhurst_code, type = "median")
p_betadisper_tina_lhc <- permutest(betadisper_tina_lhc)
anova_betadisper_tina_lhc <- anova(betadisper_tina_lhc)
plot_betadisper_tina_lhc <- plot_dispersion(betadisper_tina_lhc, title = "Longhurst code")
save_plot(filename = "figures/osd2014_planctomycetes_dispersion_lhc.pdf", plot = plot_betadisper_tina_lhc, base_width = 8, base_height = 8)

betadisper_tina_lhb <- betadisper(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs), osd2014_cdata_filt_df$longhurst_biome, type = "median")
p_betadisper_tina_lhb <- permutest(betadisper_tina_lhb)
anova_betadisper_tina_lhb <- anova(betadisper_tina_lhb)
plot_betadisper_tina_lhb <- plot_dispersion(betadisper_tina_lhb, title = "Lonhgurst biome")
save_plot(filename = "figures/osd2014_planctomycetes_dispersion_lhb.pdf", plot = plot_betadisper_tina_lhb, base_width = 7, base_height = 7)


# Permanova ---------------------------------------------------------------

perm_tina_mp <- adonis2(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs) ~ meow_province, data = osd2014_cdata_filt_df)
perm_tina_mr <- adonis2(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs) ~ meow_realm, data = osd2014_cdata_filt_df)
perm_tina_wt <- adonis2(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs) ~ wt_cat, data = osd2014_cdata_filt_df)
perm_tina_lhc <- adonis2(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs) ~ longhurst_code, data = osd2014_cdata_filt_df)
perm_tina_lhb <- adonis2(as.dist(osd2014_plancto_pina_tina_results[[2]]$cs) ~ longhurst_biome, data = osd2014_cdata_filt_df)


# Save results betadiv and permanova in excel

write_beta(beta_permut = p_betadisper_tina_mp, beta_anova = anova_betadisper_tina_mp, permanova = perm_tina_mp, label = "mp")
write_beta(beta_permut = p_betadisper_tina_mr, beta_anova = anova_betadisper_tina_mr, permanova = perm_tina_mr, label = "mr")
write_beta(beta_permut = p_betadisper_tina_wt, beta_anova = anova_betadisper_tina_wt, permanova = perm_tina_wt, label = "wt")
write_beta(beta_permut = p_betadisper_tina_lhc, beta_anova = anova_betadisper_tina_lhc, permanova = perm_tina_lhc, label = "lhc")
write_beta(beta_permut = p_betadisper_tina_lhb, beta_anova = anova_betadisper_tina_lhb, permanova = perm_tina_lhb, label = "lhb")

