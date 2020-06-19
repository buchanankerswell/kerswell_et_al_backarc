# Load libraries
library(ggplot2)
library(fields)
library(viridis)
# Load data
load('backarc_heatflow_data.RData')
# Draw Plots

# Open pdf for drawing plots
pdf(file = 'global_heat_flow_plots.pdf', width = 11, height = 8.5, useDingbats = FALSE)

# Set plot parameters
# 2x3 plots per page
par(mfrow = c(2,1),
    mar = c(1, 1, 1, 1)
)

# P. Bird boundaries ----
# # WGS84 EPSG: 4326 +proj=longlat
# plot(st_geometry(bbox), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap), col = 'grey', lwd = 1, add = TRUE)
# title('P. Bird (2003, G3) Subduction Zones')
# # plot(st_geometry(pb_bound), lwd = 1.5, col = 'black', add = TRUE)
# plot(st_geometry(filter(pb_bound, Type == 'subduction')), lwd = 2, col = 'lemonchiffon', add = TRUE)
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('P. Bird (2003, G3) Subduction Zones')
plot(st_union(pb_bound_robin_buffer_sz), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(st_geometry(pb_bound_robin), lwd = 1.5, col = 'black', add = TRUE)
plot(st_geometry(pb_bound_robin_sz), lwd = 2, col = 'lemonchiffon', add = TRUE)

# USGS boundaries ----
# # WGS84 EPSG: 4326 +proj=longlat
# plot(st_geometry(bbox), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap), col = 'grey', lwd = 1, add = TRUE)
# title('USGS Convergent Boundaries')
# # plot(st_geometry(usgs_bound), lwd = 1.5, col = 'black', add = TRUE)
# plot(st_geometry(filter(usgs_bound, Type == 'convergent')), lwd = 2, col = 'lemonchiffon', add = TRUE)
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('USGS Convergent Boundaries')
plot(st_union(usgs_bound_robin_buffer_sz), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(st_geometry(usgs_bound_robin), lwd = 1.5, col = 'black', add = TRUE)
plot(st_geometry(usgs_bound_robin_sz), lwd = 2, col = 'lemonchiffon', add = TRUE)

# UTIG boundaries ----
# # WGS84 EPSG: 4326 +proj=longlat
# plot(st_geometry(bbox), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap), col = 'grey', lwd = 1, add = TRUE)
# title('UTIG "PLATES" Trenches')
# # plot(st_geometry(utig_bound), lwd = 1.5, col = 'black', add = TRUE)
# plot(st_geometry(filter(utig_bound, Type == 'TR')), lwd = 2, col = 'lemonchiffon', add = TRUE)
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('UTIG "PLATES" Trenches')
plot(st_union(utig_bound_robin_buffer_sz), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(st_geometry(utig_bound_robin), lwd = 1.5, col = 'black', add = TRUE)
plot(st_geometry(utig_bound_robin_sz), lwd = 2, col = 'lemonchiffon', add = TRUE)

# Syracuse and Abers (2006) boundaries ----
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('Syracuse & Abers (2006) Segments')
for(i in 1:length(segs)){
  plot(segs_buffer[[i]], col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
  plot(segs[[i]][[1]][1], lwd = 2, col = 'lemonchiffon', add = TRUE)
  plot(segs[[i]][[1]][2:length(segs[[i]][[1]])], lwd = 1, col = rgb(0, 0, 0, 0.5), add = TRUE)
}

# Global heat flow data (IHFC, 2018) ----
# # WGS84 EPSG: 4326 +proj=longlat
# plot(st_geometry(bbox), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap), col = 'grey', lwd = 1, add = TRUE)
# title(main = 'Global Heat Flow Data')
# # plot(st_geometry(pb_bound), lwd = 1.5, col = 'black', add = TRUE)
# plot(global_hf_filter['Heat_Flow'], cex = 0.1, pch = 20, add = TRUE)
# image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title(main = 'Global Heat Flow Data')
# plot(st_geometry(pb_bound_robin), lwd = 1.5, col = 'black', add = TRUE)
plot(global_hf_robin_filter['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE)
image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)

# P Bird global heat flow buffer robin projection ----
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('P. Bird (2003, G3) Subduction Zones')
plot(st_union(pb_bound_robin_buffer_sz), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
plot(hf_in_pb_buffer_sz['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(pb_bound_robin), lwd = 1.5, col = 'black', add = TRUE)
plot(st_geometry(pb_bound_robin_sz), lwd = 2, col = 'lemonchiffon', add = TRUE)
image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)

# # P. Bird segments 1:13 ----
# # WGS84 EPSG: NA +proj=robin
# plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
# title('P. Bird (2003, G3) Subduction Zones Segments 1-13')
# plot(st_union(filter(pb_bound_robin_buffer_sz, segment_index >= 1 & segment_index <= 13)), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(filter(hf_in_pb_buffer_sz, segment_index >= 1 & segment_index <= 13)['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(filter(pb_bound_robin_sz, segment_index >= 1 & segment_index <= 13)), lwd = 2, col = 'lemonchiffon', add = TRUE)
# image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)
# 
# # P. Bird segments 14:26 ----
# # WGS84 EPSG: NA +proj=robin
# plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
# title('P. Bird (2003, G3) Subduction Zones Segments 14-26')
# plot(st_union(filter(pb_bound_robin_buffer_sz, segment_index >= 4 & segment_index <= 26)), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(filter(hf_in_pb_buffer_sz, segment_index >= 14 & segment_index <= 26)['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(filter(pb_bound_robin_sz, segment_index >= 14 & segment_index <= 26)), lwd = 2, col = 'lemonchiffon', add = TRUE)
# image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)
# 
# # P. Bird segments 27:39 ----
# # WGS84 EPSG: NA +proj=robin
# plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
# title('P. Bird (2003, G3) Subduction Zones Segments 27-39')
# plot(st_union(filter(pb_bound_robin_buffer_sz, segment_index >= 27 & segment_index <= 39)), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(filter(hf_in_pb_buffer_sz, segment_index >= 27 & segment_index <= 39)['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(filter(pb_bound_robin_sz, segment_index >= 27 & segment_index <= 39)), lwd = 2, col = 'lemonchiffon', add = TRUE)
# image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)
# 
# # P. Bird segments 40:52 ----
# # WGS84 EPSG: NA +proj=robin
# plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
# title('P. Bird (2003, G3) Subduction Zones Segments 40-52')
# plot(st_union(filter(pb_bound_robin_buffer_sz, segment_index >= 40 & segment_index <= 52)), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(filter(hf_in_pb_buffer_sz, segment_index >= 40 & segment_index <= 52)['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(filter(pb_bound_robin_sz, segment_index >= 40 & segment_index <= 52)), lwd = 2, col = 'lemonchiffon', add = TRUE)
# image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)
# 
# # P. Bird segments 53:65 ----
# # WGS84 EPSG: NA +proj=robin
# plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
# plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
# plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
# title('P. Bird (2003, G3) Subduction Zones Segments 53-65')
# plot(st_union(filter(pb_bound_robin_buffer_sz, segment_index >= 53 & segment_index <= 65)), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
# plot(filter(hf_in_pb_buffer_sz, segment_index >= 53 & segment_index <= 65)['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(filter(pb_bound_robin_sz, segment_index >= 53 & segment_index <= 65)), lwd = 2, col = 'lemonchiffon', add = TRUE)
# image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)

# USGS global heat flow buffer robin projection ----
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('USGS Convergent Boundaries')
plot(st_union(usgs_bound_robin_buffer_sz), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
plot(hf_in_usgs_buffer_sz['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(usgs_bound_robin), lwd = 1.5, col = 'black', add = TRUE)
plot(st_geometry(usgs_bound_robin_sz), lwd = 2, col = 'lemonchiffon', add = TRUE)
image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)

# UTIG global heat flow buffer robin projection ----
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('UTIG "PLATES" Trenches')
plot(st_union(utig_bound_robin_buffer_sz), col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
plot(hf_in_utig_buffer_sz['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
# plot(st_geometry(utig_bound_robin), lwd = 1.5, col = 'black', add = TRUE)
plot(st_geometry(utig_bound_robin_sz), lwd = 2, col = 'lemonchiffon', add = TRUE)
image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)

# Syracuse and Abers (2006) global heat flow buffer robin projection ----
# WGS84 EPSG: NA +proj=robin
plot(st_geometry(bbox_robin), lwd = 1.5, col = rgb(176, 196, 222, 127.5, maxColorValue = 255))
plot(st_geometry(grats_robin), col = 'grey', lwd = 0.5, add = TRUE)
plot(st_geometry(worldMap_robin), col = 'grey', lwd = 1, add = TRUE)
title('Syracuse & Abers (2006) Segments')
for(i in 1:length(segs)){
  plot(segs_buffer[[i]], col = rgb(0.5, 0.5, 0.5, 0.2), border = 'NA', add = TRUE)
  plot(segs_hf[[i]]['Heat_Flow'], cex = 0.2, pch = 20, add = TRUE, main = bquote('Heat Flow'~(mW/m^2)))
  plot(segs[[i]][[1]][1], lwd = 2, col = 'lemonchiffon', add = TRUE)
  plot(segs[[i]][[1]][2:length(segs[[i]][[1]])], lwd = 1, col = rgb(0, 0, 0, 0.5), add = TRUE)
}
image.plot(horizontal = FALSE, legend.mar = 9.1, legend.shrink = 0.5, legend.line = 3, legend.only = TRUE, zlim = c(0,200), col = plasma(256, alpha = 0.8), legend.lab = bquote('Heat Flow'~mW/m^2), add = TRUE)

# # P. Bird ggplot ----
# pb_hf_plot <- ggplot() +
#   geom_sf(data = grats_robin, color = 'grey', size = 0.25) +
#   geom_sf(data = bbox_robin, color = 'black', size = 1, fill = 'lightsteelblue', alpha = 0.5) +
#   geom_sf(data = worldMap_robin) +
#   # geom_sf(data = pb_bound_robin, color = 'black', size = 0.5) +
#   geom_sf(data = st_union(pb_bound_robin_buffer_sz), fill = 'grey50', color = 'NA', alpha = 0.2) +
#   geom_sf(data = hf_in_pb_buffer_sz, aes(color = Heat_Flow), size = 0.2, alpha = 0.2) +
#   geom_sf(data = pb_bound_robin_sz, color = 'lemonchiffon', size = 1) +
#   scale_color_viridis(option = 'magma') +
#   labs(color = bquote('Heat'~'Flow'~(mW/m^2))) +
#   ggtitle('Global Heat Flow Data', subtitle = '(IHFC, 2018) 1000km buffer') +
#   theme_bw(base_size = 14) +
#   theme(
#     panel.border = element_blank(),
#     axis.text = element_text(face = 'plain', color = 'black'),
#     axis.ticks = element_line(size = 1, color = 'black'),
#     axis.title = element_text(size = 12, face = 'plain'),
#     legend.position = 'bottom',
#     legend.title = element_text(size = 12, face = 'plain')
#   )
# pb_hf_plot

# # USGS ggplot ----
# usgs_hf_plot <- ggplot() +
#   geom_sf(data = grats_robin, color = 'grey', size = 0.25) +
#   geom_sf(data = bbox_robin, color = 'black', size = 1, fill = 'lightsteelblue', alpha = 0.5) +
#   geom_sf(data = worldMap_robin) +
#   # geom_sf(data = usgs_bound_robin, color = 'black', size = 0.5) +
#   geom_sf(data = st_union(usgs_bound_robin_buffer_sz), fill = 'grey50', color = 'NA', alpha = 0.2) +
#   geom_sf(data = hf_in_usgs_buffer_sz, aes(color = Heat_Flow), size = 0.2, alpha = 0.2) +
#   geom_sf(data = usgs_bound_robin_sz, color = 'lemonchiffon', size = 1) +
#   scale_color_viridis(option = 'magma') +
#   labs(color = bquote('Heat'~'Flow'~(mW/m^2))) +
#   ggtitle('Global Heat Flow Data', subtitle = '(IHFC, 2018) 1000km buffer') +
#   theme_bw(base_size = 14) +
#   theme(
#     panel.border = element_blank(),
#     axis.text = element_text(face = 'plain', color = 'black'),
#     axis.ticks = element_line(size = 1, color = 'black'),
#     axis.title = element_text(size = 12, face = 'plain'),
#     legend.position = 'bottom',
#     legend.title = element_text(size = 12, face = 'plain')
#   )
# usgs_hf_plot

# # UTIG "PLATES" ggplot ----
# utig_hf_plot <- ggplot() +
#   geom_sf(data = grats_robin, color = 'grey', size = 0.25) +
#   geom_sf(data = bbox_robin, color = 'black', size = 1, fill = 'lightsteelblue', alpha = 0.5) +
#   geom_sf(data = worldMap_robin) +
#   # geom_sf(data = utig_bound_robin, color = 'black', size = 0.5) +
#   geom_sf(data = st_union(utig_bound_robin_buffer_sz), fill = 'grey50', color = 'NA', alpha = 0.2) +
#   geom_sf(data = hf_in_utig_buffer_sz, aes(color = Heat_Flow), size = 0.2, alpha = 0.2) +
#   geom_sf(data = utig_bound_robin_sz, color = 'lemonchiffon', size = 1) +
#   scale_color_viridis(option = 'magma') +
#   labs(color = bquote('Heat'~'Flow'~(mW/m^2))) +
#   ggtitle('Global Heat Flow Data', subtitle = '(IHFC, 2018) 1000km buffer') +
#   theme_bw(base_size = 14) +
#   theme(
#     panel.border = element_blank(),
#     axis.text = element_text(face = 'plain', color = 'black'),
#     axis.ticks = element_line(size = 1, color = 'black'),
#     axis.title = element_text(size = 12, face = 'plain'),
#     legend.position = 'bottom',
#     legend.title = element_text(size = 12, face = 'plain')
#   )
# utig_hf_plot

# Close pdf ----
dev.off()
