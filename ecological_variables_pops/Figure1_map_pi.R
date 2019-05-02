library(data.table)
library(dplyr)
library(ggplot2)
library(raster)
library(rgdal)
library(viridis)
library(scatterpie)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(ggsn)
# import data -------------------------------------------------------------

r.her.coords <- read.csv("output/rh1_poolAF_secchi_sal_20190327.csv", header = TRUE)
r.her.coords <- r.her.coords %>% rename(AF = F261Y)
pop.order <- fread("data/pop_order_df.txt", header = TRUE) %>% rename(num_pop = V1)

r.her.coords <- left_join(r.her.coords, pop.order, by=c("pop" = "pop_order")) %>% dplyr::select(-X)


# pie charts -------------------------------------------------------

lon.ext <- c(3, 30.5)
lat.ext <- c(52, 68)

#backgroundmaps
world <- map_data('world')
wrld_map <- readOGR("data/10m-admin-0-countries/10m_admin_0_countries.shp")
df_wrld_map <- fortify(wrld_map)


adj.her.coords <- r.her.coords %>%
  group_by(name, lon) %>% 
  mutate(lon.adj = jitter(lon, amount = .5)) %>% 
  mutate(lat.adj = jitter(lat, amount = .5))# %>% 


adj.her.coords <- as.data.frame(adj.her.coords)

adj.her.coords <- adj.her.coords %>% mutate(season = ifelse(autumn == "Y", "Autumn", "Spring/Summer"))

pie <- adj.her.coords %>% dplyr::select(lon.adj, lat.adj, AF) %>% mutate(Y = 1 - AF, radius = 0.7) %>%
  rename(F = AF)

pie.list <- pie %>% 
  tidyr::gather(F261Y, value, -lon.adj, -lat.adj, -radius) %>% 
  tidyr::nest(F261Y, value) %>% 
  mutate(pie.grob = purrr::map(data,function(d) ggplotGrob(ggplot(d, aes(x = 1, y = value, fill = F261Y)) + 
                                                             geom_col(show.legend = FALSE) + 
                                                             coord_polar(theta = "y") + 
                                                             scale_fill_manual(values = c("black", "white")) +
                                                             theme_void()))) %>%
  rowwise() %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = lon.adj - radius, xmax = lon.adj + radius,
                                          ymin = lat.adj - radius, ymax = lat.adj + radius)))


p <- ggplot(world, aes(long, lat)) + xlim(lon.ext) + ylim(lat.ext) + 
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(text = element_text(size=40), 
        axis.text=element_text(colour="black"),
        axis.ticks.length=unit(.65, "cm"),
        panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.2)) 

pie.overlay <- p + 
  geom_tile(data = pie %>% tidyr::gather(F261Y, value, -lon.adj, -lat.adj, -radius),
            aes(x = lon.adj, y = lat.adj, fill = F261Y), 
            color = "black", width = 0.01, height = 0.01, 
            inherit.aes = FALSE) + 
  scale_fill_manual(values = c("black", "white")) +
  pie.list$subgrob +
  scale_color_manual(name = "Season", values = c("#d8b365","#5ab4ac")) +
  geom_point(data = adj.her.coords, aes(x = lon.adj, y = lat.adj), size = 23, fill = "white", shape = 1, stroke = 2) +
  geom_label_repel(data = adj.her.coords, aes(x = lon.adj, y = lat.adj, label = num_pop, color = season), fill = "white", 
                   box.padding = .65, point.padding = 3.5, segment.size = 0.8, size = 10, nudge_x	= 1, nudge_y = 1, label.size = 2.2, force = 1)
#pie.overlay

pdf("output/piechart_plotting/grob_pie_plots.pdf", width = 19.69333333333333, height = 24.746666666666666, useDingbats = FALSE)
pie.overlay
dev.off()

point.ref <- p +
  geom_point(data = adj.her.coords, aes(x = lon, y = lat), size = 2) 

pdf("output/piechart_plotting/point.ref.pdf", width = 19.69333333333333, height = 24.746666666666666, useDingbats = FALSE)
point.ref
dev.off()
