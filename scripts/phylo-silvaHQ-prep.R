library(tidyverse)

plan_silva <- read_tsv("data/planctomycetes_silva_ssu-quality.txt.gz", col_names = TRUE) %>%
  separate(path, into = c("domain",
                          "phylum",
                          "class",
                          "order",
                          "family",
                          "genus",
                          "species"), sep = ";")
plan_silva$phylum %>% table()
plan_silva %>%
  group_by(class) %>%
  count()

ggthemr::ggthemr(palette = "fresh", layout = "scientific")
plan_silva %>%
  mutate(Quality = ifelse((quality > 90 & (regionLength > 900 & regionLength < 2000)), "HQ", "notHQ")) %>%
  ggplot(aes(regionLength, quality, fill = Quality)) +
  geom_point(shape = 21, color = "black", alpha = 0.6) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(.~class) +
  theme_bw() +
  xlab("Length") +
  ylab("Quality score")

plan_silva %>%
  filter(quality > 90 & (regionLength > 900 & regionLength < 2000)) %>%
  filter(!is.na(genus), genus != "") %>%
  arrange(regionLength) %>%
  select(primaryAccession) %>%
  write_tsv(path = "data/plancto_HQ_silva.r132.tsv", col_names = FALSE)

plan_silva %>%
  filter(quality > 90 & (regionLength > 900 & regionLength < 2000)) %>%
  group_by(class) %>%
  count() %>%
  ungroup() %>%
  mutate(class = fct_rev(class)) %>%
  ggplot(aes(class, n)) +
  geom_col() +
  scale_y_log10() +
  ggpubr::rotate()

plan_silva %>%
  filter(quality > 90 & (regionLength > 900 & regionLength < 2000)) %>%
  select(domain, phylum, class, order, family, genus, species) %>%
  write_tsv(path = "data/plancto_HQ_silva.r132.alluvial.tsv", col_names = TRUE)

  group_by(genus) %>%
  count() %>%
  ungroup() %>%
  mutate(genus = fct_rev(genus)) %>%
  ggplot(aes(genus, n)) +
  geom_col() +
  scale_y_log10() +
  ggpubr::rotate()
