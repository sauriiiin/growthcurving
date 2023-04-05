##### LIQ SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : saurinbp@gmail.com
##### Date   : 04/03/2023

##### INITIALIZE
shhh <- suppressPackageStartupMessages
options(warn=-1)

shhh({
library(dplyr)
library(reshape2)
library(growthcurver)
library(lattice)
library(deSolve)
library(growthrates)
library(ggplot2)})

`%notin%` <- Negate(`%in%`)

system("mkdir -p output")
system("mkdir -p output/figs")

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GATHER DATA
cat("What is the experiment ID: ")
expt_id <- readLines(con = 'stdin', n = 1)
expt_id <- expt_id[[1]]

platemap <- data.frame(expt_id = expt_id, read.csv(sprintf('input/%s_Platemap.csv',expt_id)))
platemap$smudge[is.na(platemap$smudge)] <- 'no'

data <- data.frame(expt_id = expt_id, read.csv(sprintf('input/%s_Data.csv',expt_id)))

data2 <- merge(data %>% melt(id.vars = c('expt_id','Time'), variable.name = 'well_ID', value.name = 'OD'),
               platemap, by = c('expt_id','well_ID'))
# unique(data2$condition)
cnd_ids <- c("YPDA (pH 7.5)","DMSO (pH 7.5)","1ug/ml Tunicamycin (pH 7.5)")
data2$condition <- factor(data2$condition, levels = cnd_ids)

# unique(data2$sample)
orf_ids <- c("BY4741","BY4742","DELHO::KAN (AARON PAPER)",
             "YBR KAN DEL", "YBR KAN DEL (AARON PAPER)",
             "YBR CRSPR (YARC492)", "YBR CRSPR (YARC493)",
             "YBR196C-A ATG", "YBR196C-A AAG",
             "YBR CRSPR (YARC141)", "YBR CRSPR (YARC142)",
             "YDL204W-A ATG", "YDL204W-A AAG",
             "BLANK")
data2$sample <- factor(data2$sample, levels = orf_ids)

##### GROWTH ANALYSIS
data.gc <- data.frame(expt_id = expt_id, SummarizeGrowthByPlate(data[,-1]))
data.gr <- NULL
for (e in unique(data$expt_id)) {
  if (data$Time[data$expt_id == e][2] - data$Time[data$expt_id == e][1] > 30) {
    h <- 3
  } else {
    h <- 8
  }
  for (i in 3:dim(data[data$expt_id == e,])[2]) {
    if (sum(data[data$expt_id == e,i] < 0.005) < 5) {
      fit0 <- fit_easylinear(data$Time[data$expt_id == e], data[data$expt_id == e,i], h = h, quota = 1);
      
      temp_res <- data.frame(expt_id = e, sample = colnames(data[i]), maxgr = coef(fit0)[[3]],
                             dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
      data.gr <- rbind(data.gr, temp_res)
    }
  }
}
# head(data.gc)
# head(data.gr)

res.liq <- merge(platemap, merge(data.gc, data.gr, by = c('expt_id','sample')), by.x = c('expt_id','well_ID'), by.y = c('expt_id','sample'))
res.liq <- res.liq %>% filter(smudge != 'yes', note == '')
res.liq <- res.liq[,colnames(res.liq) %notin% c('smudge','note')]
# head(res.liq)

res.liq$condition <- factor(res.liq$condition, levels =  cnd_ids)
res.liq$sample <- factor(res.liq$sample, levels = orf_ids)

##### GROWTH CURVES
# head(data2)

plot.gc <- data2 %>%
  filter(sample != 'BLANK', smudge != 'yes') %>%
  ggplot() +
  stat_summary(aes(x = Time, y = OD, col = sample),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = Time, y = OD, fill = sample),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_grid(background~condition) +
  labs(x="Minute",y="OD")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))
ggsave(sprintf("output/figs/%s_growth_curves.jpg",expt_id), plot.gc,
       height = 180, width = 300,
       units = 'mm',
       bg = 'white',
       dpi = 300)


##### GROWTH DYNAMICS
plot.auc <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = sample, y = auc_l)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  labs(x = '', y = 'AUC') +
  facet_grid(background~condition) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"))
ggsave(sprintf("output/figs/%s_area_under_curve.jpg",expt_id), plot.auc,
       height = 180, width = 300,
       units = 'mm',
       bg = 'white',
       dpi = 300)


plot.dtime <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = sample, y = t_gen)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_log10() +
  labs(x = '', y = 'Doubling Time (min.)') +
  facet_grid(background~condition) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"))
ggsave(sprintf("output/figs/%s_doubling_time.jpg",expt_id), plot.dtime,
       height = 180, width = 300, units = 'mm',
       bg = 'white',
       dpi = 300)



plot.ltime <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = sample, y = ltime)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_log10() +
  labs(x = '', y = 'Lag Time (min.)') +
  facet_grid(background~condition) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"))
ggsave(sprintf("output/figs/%s_lag_time.jpg",expt_id), plot.ltime,
       height = 180, width = 300, units = 'mm',
       bg = 'white',
       dpi = 300)


write.csv(file = sprintf('output/%s_results.csv', expt_id), res.liq, row.names = F)