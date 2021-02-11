rm(list = ls())
source('functions.R')
load('data/hf.RData')

# Global

# Summarise global hf data
print('Global heat flow data')
d.hf %>%
  group_by(segment) %>%
  summarise(n = n())

# Summary plots
# Number of hf data points
p.no <- d.hf %>%
  group_by(segment) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(x = n, y = forcats::fct_rev(factor(segment)), fill = segment), stat = 'identity', alpha = 0.8, show.legend = F) +
  labs(x = NULL, y = NULL, title = 'Number of Heat Flow Data Points') +
  labs(y = 'Segment', x = NULL, title = bquote(Data~Points)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )

# Print
p.no

# Heat flow
p.hf <- d.hf %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = hf, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(Heat~Flow~~mWm^-2)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )

# Print
p.hf

# grad
p.grad <- d.hf %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = grad, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(grad~~mKm^-1)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )

# Print
p.grad

# con
p.con <- d.hf %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = con, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(con~~Wm^-1~K^-1)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )

# Print
p.con

# ele
p.el <- d.hf %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = ele, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(ele~~m)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )

# Print
p.el

# No. of Temps
p.nt <- d.hf %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = tempNo, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(Number~of~Temps)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )

# Print
p.nt

# Composite summary plot
(p.no + theme(plot.margin = margin(8, 0, 8, 8), axis.text.y = element_text(margin = margin(0, 4, 0, 0)))) +
  (p.hf + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) +
  (p.grad + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) +
  (p.con + theme(plot.margin = margin(8, 0, 8, 8), axis.text.y = element_text(margin = margin(0, 4, 0, 0)))) +
  (p.el + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) +
  (p.nt + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) &
  theme(panel.spacing.x = unit(0, 'pt'), axis.line.y.left = element_line(color = rgb(0, 0, 0, 0.2)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA))

# Save plot
ggsave('figs/global-summary.png', device = 'png', dpi = 330, width = 11, height = 7, bg = 'transparent')

# Summary by variable
print('Data summary')
print(d.hf %>%
        select(c(hf, grad, con, ele, tempNo, segment)) %>%
        pivot_longer(-segment, 'variable') %>%
        mutate(variable = factor(variable, levels = c('hf', 'grad', 'con', 'ele', 'tempNo'))) %>%
        group_by(segment, variable) %>%
        summarise(
          min = min(value, na.rm = T),
          max = max(value, na.rm = T),
          median = median(value, na.rm = T),
          mean = mean(value, na.rm = T)),
      n = 65)