library(ggordiplots)

spp_PCoA <- fishFD$x.axes #extract coordinates from PCoA
spp_groups <- fishFD$spfgr %>% as.factor() #extract groups

# create basic ordiplot object to pull into advanced plot later
PCoA_ordiplot <- gg_ordiplot(ord = spp_PCoA,
                             groups = spp_groups,
                             ellipse = FALSE,
                             hull = FALSE,
                             spiders = TRUE)

PCoA_spiders <- PCoA_ordiplot$df_spiders %>% 
  rename(A1 = x, A2 = y) %>%
  as_tibble

# basic plot building block
ord_gg <- ggplot() + theme_bw(16) + 
  geom_vline(xintercept = 0, lty=3, color="darkgrey") +
  geom_hline(yintercept = 0, lty=3, color="darkgrey") +
  theme(panel.grid=element_blank())

# pal <- pnw_palette("Bay", type = "discrete")

# advanced plot with spiders, sans vectors
PCoA_plot <- ord_gg +
  geom_text(data=spp_PCoA,
            aes(x=A1, y=A2,
                colour=spp_groups,
                label= colnames(fish.list$abund)),
            #label.padding=unit(0.3,"lines"),
            fontface="bold")  +
  # geom_point(data=PCoA_list$species,
  #            aes(x=A1, y=A2,
  #                colour=spp_groups))  +
  labs(x="PCoA Axis 1",
       y="PCoA Axis 2") +
  # scale_color_manual(values = pal) +
  geom_segment(data=PCoA_spiders,
               aes(x=cntr.x, y=cntr.y,
                   xend=A1, yend=A2, color = Group),
               size=1); PCoA_plot
