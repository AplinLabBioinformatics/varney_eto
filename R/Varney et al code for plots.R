library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
library(writexl)
library(ggsignif)
library(RColorBrewer)
library(WRS2)
# library(userfriendlyscience)
library(rstatix)
library(survival)
library(ggsurvfit)
library(readxl)
library(factoextra)
library(cluster)
library(grid)
library(robustbase)

setwd(file.path(dirname(getwd()),"data"))

##### Figure 2 GSDM cleavage plots

GSDM.cl.data <- read.csv("MP.Meti.GSDM.cleave.r.csv")

draw <-function(plot, x_in = 4.24, y_in = 2.6) {
  
  grid::grid.newpage()
  
  grid::rectGrob(gp = grid::gpar(fill = "gray")) |>
    grid::grid.draw()
  
  grid::viewport(width  = ggplot2::unit(x_in, "in"), 
                 height = ggplot2::unit(y_in, "in")) |>
    grid::pushViewport()
  
  ggplot2::ggplot_build(plot) |>
    ggplot2::ggplot_gtable() |>
    grid::grid.draw()
}


# Calculate mean and standard deviation
avg_data <- GSDM.cl.data %>%
  group_by(Cells, GSDM_Frag, Condition, GSDM) %>%
  summarize(AveragePercent = mean(Percentage, na.rm = TRUE),
            StdDev = sd(Percentage, na.rm = TRUE))

avg_data$Condition <- factor(avg_data$Condition, 
                             levels = c("DMSO","ETO","IACs","WZB","6AN"))

avg_data$GSDM <- factor(avg_data$GSDM, 
                        levels = c("GSDME","GSDMD"))

avg_data$GSDM_Frag <- factor(avg_data$GSDM_Frag, 
                             levels = c("GSDME-FL","GSDME-NT","GSDMD-FL","GSDMD-43kDa","GSDMD-NT"))


# Filter the data for "GSDME-NT" and "GSDMD-NT"
filtered_data.GSDME <- GSDM.cl.data %>%
  filter(GSDM_Frag %in% c("GSDME-NT"))

filtered_data.GSDMD <- GSDM.cl.data %>%
  filter(GSDM_Frag %in% c("GSDMD-NT"))

cell_lines <- unique(GSDM.cl.data$Cells)

###

perform_anova_tukey <- function(data, cell_line) {
  cell_data <- data %>% filter(Cells == cell_line)
  
  # Perform one-way ANOVA
  anova_result <- aov(Percentage ~ Condition, data = cell_data)
  print(paste("ANOVA results for", cell_line))
  print(summary(anova_result))
  
  # Perform Tukey's HSD post hoc test
  tukey_result <- TukeyHSD(anova_result)
  print(paste("Tukey HSD results for", cell_line))
  print(tukey_result)
}


for (cell_line in cell_lines) {
  perform_anova_tukey(filtered_data.GSDME, cell_line)
}



########## Plots for GSDM cleavage

GSDM.cl.data$Condition <- factor(GSDM.cl.data$Condition, 
                                 levels = c("DMSO","ETO","IACs","WZB","6AN"))

custom_colors <- c("DMSO" = "#696969", "ETO" = "#696969", "IACs" = "#696969",
                   "WZB" = "#696969","6AN" = "#696969")


GSDME.cleave<-ggplot(GSDM.cl.data%>% filter(GSDM_Frag == "GSDME-NT"), aes(x = Condition, y = Percentage, fill = Condition)) +
  geom_boxplot(aes(group = Condition), color = "black", outlier.size = 0) +
  geom_dotplot(binaxis = "y",  
               color = "black",
               fill = 'black',
               binwidth = max(GSDM.cl.data$Percentage)*.03, stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center")+
  facet_grid(~ Cells)+
  labs(title = "Cleaved GSDME-NT Percent",
       x = NULL,
       y = "GSDME-NT % of total") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = custom_colors,
                    guide = NULL)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size=10, color = "black"),
        plot.title=element_text(size=12),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10),
        panel.grid = element_line(color = "#F7F7F7"),
        strip.background = element_rect("white", margin(0,0,0,0)))+
  scale_y_continuous(limits = c(-0.01,1.05))+
  geom_hline(yintercept = 0, color = "black", size = 1)+
  geom_line(data = tibble(x = c(1,2), y = c(.76, .76)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1,1), y = c(.25, .76)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1,4), y = c(.95, .95)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(4,4), y = c(.87, .95)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(2,2), y = c(.73, .76)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38", x = 1.8, y = 0.83),
            aes(x = x, y = y),
            label = "p = 0.012",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38", x = 2.5, y = 1.02),
            aes(x = x, y = y),
            label = "n.s.",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP46", x = 1.8, y = 0.83),
            aes(x = x, y = y),
            label = "p = 0.005",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP46", x = 2.5, y = 1.02),
            aes(x = x, y = y),
            label = "p = 0.002",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65", x = 1.9, y = 0.83),
            aes(x = x, y = y),
            label = "p = 0.0008",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65", x = 2.5, y = 1.02),
            aes(x = x, y = y),
            label = "n.s.",
            size = 2.5,
            inherit.aes = FALSE)


draw(GSDME.cleave)
ggsave("GSDME.cleave.png", plot = GSDME.cleave,
       width = 4.24, height = 2.6, units = "in", dpi = 300)


###### now for GSDMD

GSDMD.cleave<-ggplot(GSDM.cl.data%>% filter(GSDM_Frag == "GSDMD-NT"), aes(x = Condition, y = Percentage, fill = Condition)) +
  geom_boxplot(aes(group = Condition), color = "black", outlier.size = 0) +
  geom_dotplot(binaxis = "y",  
               color = "black",
               fill = 'black',
               binwidth = max(GSDM.cl.data$Percentage)*.03, stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center")+
  facet_grid(~ Cells)+
  labs(title = "Cleaved GSDMD-NT Percent",
       x = NULL,
       y = "GSDMD-NT % of total") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = custom_colors,
                    guide = NULL)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size=10, color = "black"),
        plot.title=element_text(size=12),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10),
        panel.grid = element_line(color = "#F7F7F7"),
        strip.background = element_rect("white", margin(0,0,0,0)))+
  scale_y_continuous(limits = c(-0.01,1.05))+
  geom_hline(yintercept = 0, color = "black", size = 1)

draw(GSDMD.cleave)
ggsave("GSDMD.cleave.png", plot = GSDMD.cleave,
       width = 4.24, height = 2.6, units = "in", dpi = 300)

# END Fig 2 #




# Figure 3 CPT1A activity related




# Panel A TIM


#Data for Panel B
UVM.CPT1A.surv.data <- read.csv("UVM.SCNA.met.cpt1.pscore.csv")


#Stats for figure 3B

UVM.CPT1A.surv.data %>% t_test(CPT1AinhibDN ~ BAP1.SCNA,
                               p.adjust.method = "bonferroni")

# plot for Figure 3B 

CPT1A.BAP1.BP <- ggplot(UVM.CPT1A.surv.data, aes(x = BAP1.SCNA, y = CPT1AinhibDN, fill = BAP1.SCNA)) +
  geom_boxplot(aes(group = BAP1.SCNA), color = "black") +
  geom_dotplot(binaxis = "y",  
               aes(color = BAP1.SCNA), 
               binwidth = max(UVM.CPT1A.surv.data$CPT1AinhibDN)*.06, stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center") +  # Adjust jitter width here
  labs(title = "CPT1A activity in UM (TCGA-UVM)",
       x = NULL, y = "CPT1A-inhib. down \n signature score",
       fill = "SCNA cluster\n(BAP1 Status)") +  # Add this line to set the legend title
  scale_color_manual(values = rep("black", 10),
                     guide = "none") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), labels = c("SCNA 1/2\n(BAP1 WT)", "SCNA 3/4\n(BAP1-Def)")) +
  theme_bw() +
  scale_y_continuous(limits = c(min(UVM.CPT1A.surv.data$CPT1AinhibDN), 
                                max(UVM.CPT1A.surv.data$CPT1AinhibDN)*1.3)) +
  theme(axis.text.x = element_text(size=10, color = "black"),
        plot.title = element_text(size=12, hjust = 0.25, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=8, color = "black"),
        legend.text = element_text(size=8, color = "black"),
        legend.title = element_text(size=9, hjust = 0.5, color = "black"),
        legend.position = "right",
        legend.margin = margin(-5, -5, -5, -5),
        legend.box.spacing = unit(.5, "lines"),
        legend.key.spacing.y = unit(0.5, "lines")) +
  geom_line(data = tibble(x = c(1, 2), y = c(0.48, 0.48)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE) +
  geom_line(data = tibble(x = c(1, 1), y = c(.43, 0.48)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE) +
  geom_text(data = tibble(x = 1.5, y = 0.55),
            aes(x = x, y = y),
            label = "p < 0.001",
            size = 2.5,
            inherit.aes = FALSE)



draw(CPT1A.BAP1.BP)
ggsave("CPT1A.BAP1.BP.png", plot = CPT1A.BAP1.BP,
       width = 3.15, height = 2.6, units = "in", dpi = 300)


###### Fig 3C  Survival data 


# NOTE THIS WILL ONLY WORK AFTER RUNNING CODE For km clusters (figure s2)
# NOTE THIS WILL ONLY WORK AFTER RUNNING CODE For km clusters (figure s2)
# NOTE THIS WILL ONLY WORK AFTER RUNNING CODE For km clusters (figure s2)

UM.specific.MET <- UVM.CPT1A.surv.data %>%
  filter(X1..Time.to.UM.Metastatsis != "[Unknown]")

UM.specific.MET$Time.to.UM.Met.or.last.followup <-as.numeric(UM.specific.MET$Time.to.UM.Met.or.last.followup)
UM.specific.MET$Time.to.UM.Met.or.last.followup.mo <- UM.specific.MET$Time.to.UM.Met.or.last.followup / 30.44


###### using time to MET 
surv_object <- Surv(time = UM.specific.MET$Time.to.UM.Met.or.last.followup.mo, 
                    event = UM.specific.MET$Metastatic.Disease == "Yes")

fit <- survfit(surv_object ~ OrderedCluster, data = UM.specific.MET)

logrank_test.T <- survdiff(surv_object ~ OrderedCluster, data = UM.specific.MET)

print(logrank_test.T)


custom_colors.T<- c("#FE9929", "#D7191C", "#004C99")
### group order classification 

draw <- function(plot, x_in = 4, y_in = 2.5) {
  
  grid::grid.newpage()
  
  grid::rectGrob(gp = grid::gpar(fill = "gray")) |>
    grid::grid.draw()
  
  grid::viewport(width  = ggplot2::unit(x_in, "in"), 
                 height = ggplot2::unit(y_in, "in")) |>
    grid::pushViewport()
  
  ggplot2::ggplot_build(plot) |>
    ggplot2::ggplot_gtable() |>
    grid::grid.draw()
}



CPT1A.cluster.METsurv<-ggsurvfit(fit, linewidth = 1) +
  labs(title = "CPT1A 'activity' Signature (GSVA)",
       x = "Time (Months)",
       y = "UM Specific Metastasis",
       color = "CPT1A Act.\nKM-Cluster\n(cases/events)", # Change legend title
       fill = "CPT1A Cluster \n(cases/events)") + # Change legend title
  theme_bw() +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5), # Adjust plot margins
        legend.position = "right") + # Position legend at the bottom
  annotate("label", x = Inf, y = 1, label = paste("Log-rank p-value =", 
                                                  signif(logrank_test.T$pvalue, 3)),
           hjust = 1.1, vjust = 1, alpha = 0.5, fill = "white", size = 2.5) + # Semi-transparent background
  scale_color_manual(values = custom_colors.T, labels = c("Cluster 1\n(22/3)", "Cluster 2\n(30/4)", "Cluster 3\n(18/9)")) + # Set legend labels
  scale_fill_manual(values = custom_colors.T, labels = c("Cluster 1\n(22/3)", "Cluster 2\n(30/4)", "Cluster 3\n(18/9)")) +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, color = "black",  hjust = 0.4),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(hjust = .5,size = 10, color = "black"),
        legend.position = "right",
        legend.margin = margin(-0.5,-1,-2,-5),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        legend.box.spacing = unit(.5, "lines"),
        legend.key.spacing.y = unit(0.5, "lines")) +
  scale_y_continuous(limits = c(0.25, 1))+
  add_censor_mark() 


draw(CPT1A.cluster.METsurv)
ggsave("CPT1A.cluster.METsurv.png", plot = CPT1A.cluster.METsurv,
       width = 3.8, height = 2.59, units = "in", dpi = 300)


# END FIGURE 3

#Start Figure 4

MP38data <- read.csv("MP38data.Rformat.csv")
MP46data <- read.csv("MP46data.Rformat.csv")



Treatment_order<-c("DMSO","ETO 50 µM","ETO 75 µM","ETO 100 µM")

MP46data[ MP46data == "ETO 50uM"] <-"ETO 50 µM"
MP46data[ MP46data == "ETO 75uM"] <-"ETO 75 µM"
MP46data[ MP46data == "ETO 100uM"] <-"ETO 100 µM"
MP38data[ MP38data == "ETO 50 \xb5M"] <-"ETO 50 µM"
MP38data[ MP38data == "ETO 75 \xb5M"] <-"ETO 75 µM"
MP38data[ MP38data == "ETO 100 \xb5M"] <-"ETO 100 µM"


DF.46<-MP46data
DF.38<-MP38data
### setting drawing sizes for panel of figure
draw <- function(plot, x_in = 2.83, y_in = 2.44) {
  
  grid::grid.newpage()
  
  grid::rectGrob(gp = grid::gpar(fill = "gray")) |>
    grid::grid.draw()
  
  grid::viewport(width  = ggplot2::unit(x_in, "in"), 
                 height = ggplot2::unit(y_in, "in")) |>
    grid::pushViewport()
  
  ggplot2::ggplot_build(plot) |>
    ggplot2::ggplot_gtable() |>
    grid::grid.draw()
}

#Making compiled line plots

avg_values.38 <- aggregate(red.per ~ Hour + Treat, data = DF.38, 
                           FUN = function(x) c(mean = mean(x), sd = sd(x)))

avg_values.46 <- aggregate(red.per ~ Hour + Treat, data = DF.46, 
                           FUN = function(x) c(mean = mean(x), sd = sd(x)))



#PLot for Figure 4B

MP38.ETO.lineplot<-ggplot(avg_values.38, 
                          aes(x = Hour, y = red.per[, "mean"]*100, 
                              color = Treat, group = Treat)) +
  geom_line(linetype = "solid", linewidth = 1) +
  geom_point(size = 1) +
  labs(title = "PI uptake in MP38",
       x = "Hours",
       y = "PI+ Area / Phase Area") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "lightcoral", "red", "darkred"), 
                     breaks = Treatment_order,
                     guide = guide_legend(title = NULL)) +
  guides(fill = "none", shape = "none") +
  geom_errorbar(aes(ymin = red.per[, "mean"]*100 - red.per[, "sd"]*100, 
                    ymax = red.per[, "mean"]*100 + red.per[, "sd"]*100), 
                width = 0.5)+ 
  scale_x_continuous(breaks = c(0,24,48,72))+
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "inside",
        legend.position.inside = c(0.04,0.95),
        legend.justification = c("left", "top"),
        legend.key.spacing = unit(0.5, "lines"),
        legend.key.spacing.y = unit(-0.5, "lines"),
        legend.box.background = element_rect(color = "black", linewidth = 0.8),
        legend.margin = margin(0.5,0.5,0.5,0.5),
        legend.text = element_text(size=8, color = "black"),
        legend.title = NULL)+
  theme(plot.title = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=10, color = "black"),
        axis.text.x = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"))

draw(MP38.ETO.lineplot)

ggsave("MP38.ETO.lineplot.png", plot = MP38.ETO.lineplot,
       width = 2.83, height = 2.44, units = "in", dpi = 300)


#plot for Figure 4E

MP46.ETO.lineplot<-ggplot(avg_values.46, 
                          aes(x = Hour, y = red.per[, "mean"]*100, 
                              color = Treat, group = Treat)) +
  geom_line(linetype = "solid", linewidth = 1) +
  geom_point(size = 1) +
  labs(title = "PI uptake in MP46",
       x = "Hours",
       y = "PI+ Area / Phase Area") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "lightcoral", "red", "darkred"), 
                     breaks = Treatment_order,
                     guide = guide_legend(title = NULL)) +
  guides(fill = "none", shape = "none") +
  geom_errorbar(aes(ymin = red.per[, "mean"]*100 - red.per[, "sd"]*100, 
                    ymax = red.per[, "mean"]*100 + red.per[, "sd"]*100), 
                width = 0.5)+ 
  scale_x_continuous(breaks = c(0,24,48,72))+
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "inside",
        legend.position.inside = c(0.04,0.95),
        legend.justification = c("left", "top"),
        legend.key.spacing = unit(0.5, "lines"),
        legend.key.spacing.y = unit(-0.5, "lines"),
        legend.box.background = element_rect(color = "black", linewidth = 0.8),
        legend.margin = margin(0.5,0.5,0.5,0.5),
        legend.text = element_text(size=8, color = "black"),
        legend.title = NULL)+
  theme(plot.title = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=10, color = "black"),
        axis.text.x = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"))

draw(MP46.ETO.lineplot)
ggsave("MP46.ETO.lineplot.png", plot = MP46.ETO.lineplot,
       width = 2.83, height = 2.44, units = "in", dpi = 300)



# set up for creating AUC plots and setting sum timepoint limits

timepoint <- 72  # Define the timepoint (in hours) you want to limit


auc_38 <- MP38data %>%
  filter(Hour <= timepoint) %>%  # Filter data based on the timepoint
  group_by(Exp.N, Treat) %>%
  summarize(AUC = sum((red.per[-1] + red.per[-n()]) * diff(Hour))/2)

auc_46 <- MP46data %>%
  filter(Hour <= timepoint) %>%  # Filter data based on the timepoint
  group_by(Exp.N, Treat) %>%
  summarize(AUC = sum((red.per[-1] + red.per[-n()]) * diff(Hour))/2)


# this creates a factor column that will properly order the treatment groups
auc_38$Treat.2 <- factor(auc_38$Treat, 
                         levels = c("DMSO","ETO 50 µM","ETO 75 µM","ETO 100 µM"))
auc_46$Treat.2 <- factor(auc_46$Treat, 
                         levels = c("DMSO","ETO 50 µM","ETO 75 µM","ETO 100 µM"))

#Aves and stdev for error bars hopefully
AUC38aves <-aggregate(AUC ~ Treat.2, 
                      data = auc_38, 
                      FUN = function(x) c(mean = mean(x), 
                                          sd = sd(x)))
AUC46aves <-aggregate(AUC ~ Treat.2, 
                      data = auc_46, 
                      FUN = function(x) c(mean = mean(x), 
                                          sd = sd(x)))



# now do the anovas

aov_model.38 <- aov(AUC ~ Treat, data = auc_38)
anova_summary.38 <- summary(aov_model.38)
p_value.38 <- anova_summary.38[[1]][["Pr(>F)"]]  # Extracting the p-value
print(p_value.38)
# 1.068014e-10 (timepoint set to 72h)

aov_model.46 <- aov(AUC ~ Treat, data = auc_46)
anova_summary.46 <- summary(aov_model.46)
p_value.46 <- anova_summary.46[[1]][["Pr(>F)"]]  # Extracting the p-value
print(p_value.46)
# 7.891553e-09 (timepoint set to 72h)

# Tukey time!

tukey_test.38 <- TukeyHSD(aov_model.38)
print(tukey_test.38)

tukey_test.46 <- TukeyHSD(aov_model.46)
print(tukey_test.46)

# first set draw dims #
draw <- function(plot, x_in = 2.34, y_in = 2.44) {
  
  grid::grid.newpage()
  
  grid::rectGrob(gp = grid::gpar(fill = "gray")) |>
    grid::grid.draw()
  
  grid::viewport(width  = ggplot2::unit(x_in, "in"), 
                 height = ggplot2::unit(y_in, "in")) |>
    grid::pushViewport()
  
  ggplot2::ggplot_build(plot) |>
    ggplot2::ggplot_gtable() |>
    grid::grid.draw()
}


SF <-0.7

max.2 <-max(auc_46$AUC[auc_46$Treat.2 == "ETO 50 µM"])+ SF
max.3 <- max(auc_46$AUC[auc_46$Treat.2 == "ETO 75 µM"]) + SF
max.4 <- max(auc_46$AUC[auc_46$Treat.2 == "ETO 100 µM"]) + SF 

max.38.2 <-max(auc_38$AUC[auc_38$Treat.2 == "ETO 50 µM"])+ SF
max.38.3 <- max(auc_38$AUC[auc_38$Treat.2 == "ETO 75 µM"]) + SF
max.38.4 <- max(auc_38$AUC[auc_38$Treat.2 == "ETO 100 µM"]) + SF


##### Plot for Figure 4C



MP38.ETO.PI.boxplot<-ggplot(auc_38, aes(x = Treat.2, y = AUC, fill = Treat.2)) +
  geom_boxplot(color = "black",outlier.size = 0) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               aes(color = Treat.2), 
               binwidth = 0.3, stroke = 2) +
  labs(title = "PI uptake in MP38 (72 h)",x = NULL, y = "PI+/Phase (AUC)") +
  scale_color_manual(values = c("black", "black", "black", "black"),
                     guide = "none")+
  scale_fill_manual(values = c("blue", "lightcoral", "red", "darkred"),
                    guide = "none") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45,hjust = 1  , size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"))+
  geom_line(data = tibble(x = c(1, 2), y = c(max.38.2, max.38.2)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_line(data = tibble(x = c(1, 1), y = c(max.38.2, max.38.2*.8)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_text(data = tibble(x = 1.5, y = max.38.2+0.5),
            aes(x = x, y = y),
            label = "n.s.",
            size = 3,
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1, 3), y = c(max.38.3*.95, max.38.3*.95)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_line(data = tibble(x = c(1, 1), y = c(max.38.3*.95, max.38.3*.7)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_text(data = tibble(x = 2, y = max.38.3*1.05),
            aes(x = x, y = y),
            label = "p < 0.0001",
            size = 3,
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1, 4), y = c(max.38.4*.95, max.38.4*.95)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_line(data = tibble(x = c(1, 1), y = c(max.38.4*.95, max.38.4*.65)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_text(data = tibble(x = 2.5, y = max.38.4),
            aes(x = x, y = y),
            label = "p < 0.0001",
            size = 3,
            inherit.aes = FALSE)

draw(MP38.ETO.PI.boxplot)

ggsave("MP38.ETO.PI.boxplot.png", plot = MP38.ETO.PI.boxplot,
       width = 2.35, height = 2.44, units = "in", dpi = 300)

#### Plot for Figure 4F

MP46.ETO.PI.boxplot<-ggplot(auc_46, aes(x = Treat.2, y = AUC, fill = Treat.2)) +
  geom_boxplot(color = "black",outlier.size = 0) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               aes(color = Treat.2), 
               binwidth = 0.4, stroke = 2) +
  labs(title = "PI uptake in MP46 (72 h)",x = NULL, y = "PI+/Phase (AUC)") +
  scale_color_manual(values = c("black", "black", "black", "black"),
                     guide = "none")+
  scale_fill_manual(values = c("blue", "lightcoral", "red", "darkred"),
                    guide = "none") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45,hjust = 1  , size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"))+
  geom_line(data = tibble(x = c(1, 2), y = c(max.2, max.2)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_line(data = tibble(x = c(1, 1), y = c(max.2, max.2*.8)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_text(data = tibble(x = 1.5, y = max.2+0.5),
            aes(x = x, y = y),
            label = "n.s.",
            size = 3,
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1, 3), y = c(max.3*.95, max.3*.95)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_line(data = tibble(x = c(1, 1), y = c(max.3*.95, max.3*.7)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_text(data = tibble(x = 2, y = max.3*1.05),
            aes(x = x, y = y),
            label = "p < 0.0001",
            size = 3,
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1, 4), y = c(max.4*.95, max.4*.95)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_line(data = tibble(x = c(1, 1), y = c(max.4*.95, max.4*.65)),
            aes(x = x, y = y),
            size = 0.5,color = 'grey41',
            inherit.aes = FALSE) +
  geom_text(data = tibble(x = 2.5, y = max.4),
            aes(x = x, y = y),
            label = "p < 0.0001",
            size = 3,
            inherit.aes = FALSE)

draw(MP46.ETO.PI.boxplot)

ggsave("MP46.ETO.PI.boxplot.png", plot = MP46.ETO.PI.boxplot,
       width = 2.35, height = 2.44, units = "in", dpi = 300)


#### End Figure 4


#### Start Figure 5

MP38.siCPT1A <- read.csv("MP38.siCPT1A.PI.data.csv")
MP65.siCPT1A <- read.csv("MP65.siCPT1A.PI.data.csv")

Treatment_order<-c("siControl","siCPT1A #1","siCPT1A #2")


AUC_values.38 <- aggregate(red.per ~ Exp.N + Treat, data = MP38.siCPT1A, 
                           FUN = function(x) AUC = sum(x))

AUC_values.65 <- aggregate(red.per ~ Exp.N + Treat, data = MP65.siCPT1A, 
                           FUN = function(x) AUC = sum(x))


AUC_values.38$Cells <- rep("MP38", 9)
AUC_values.65$Cells <- rep("MP65", 9)

siCPT1A.PI.AUC<-rbind(AUC_values.38, AUC_values.65)



# MP38 AUC stats #

aov_siCPT1A.38 <- aov(red.per ~ Treat, data = AUC_values.38)
anova_summary_siCPT1A.38 <- summary(aov_siCPT1A.38)
print(anova_summary_siCPT1A.38)

tukey_test_siCPT1A.38 <- TukeyHSD(aov_siCPT1A.38)
print(tukey_test_siCPT1A.38)

# MP65 AUC stats #

aov_siCPT1A.65 <- aov(red.per ~ Treat, data = AUC_values.65)
anova_summary_siCPT1A.65 <- summary(aov_siCPT1A.65)
print(anova_summary_siCPT1A.65)

tukey_test_siCPT1A.65 <- TukeyHSD(aov_siCPT1A.65)
print(tukey_test_siCPT1A.65)

# first set draw dims #
draw <-function(plot, x_in = 6, y_in = 2.75) {
  
  grid::grid.newpage()
  
  grid::rectGrob(gp = grid::gpar(fill = "gray")) |>
    grid::grid.draw()
  
  grid::viewport(width  = ggplot2::unit(x_in, "in"), 
                 height = ggplot2::unit(y_in, "in")) |>
    grid::pushViewport()
  
  ggplot2::ggplot_build(plot) |>
    ggplot2::ggplot_gtable() |>
    grid::grid.draw()
}


#### Plot for Figure 5B


siCPT1A.PI.AUC.bp<-ggplot(siCPT1A.PI.AUC, aes(x = Treat, y = red.per, fill = Treat)) +
  geom_boxplot(aes(group = Treat), color = "black", outlier.size = 0) +
  geom_dotplot(binaxis = "y",  
               aes(color = Treat), 
               stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center")+
  facet_wrap(~ Cells, scales = "free")+  # Adjust jitter width here
  labs(title = "PI uptake after 72 h",
       x = NULL, y = "PI/Phase (AUC)") +
  scale_color_manual(values = rep("black", 6),
                     guide = "none") +
  scale_fill_manual(values = c("#AA4499", "#fE732F", "#FE9929"),
                    guide = NULL) +
  theme_bw()+
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0, 0.2)))+
  theme(plot.title = element_text(size = 12, hjust = 0.5, color = "black", margin = margin(1,0,0,0)),
        axis.text.x = element_text(size=10, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        strip.background = element_rect("white", margin(0,0,0,0)),
        strip.text = element_text(size = 10),
        plot.margin = margin(t = -1, r = 1, b = 1, l = 1),
        legend.box.spacing = unit(0.1, "lines"),
        legend.key.width = unit(1, "lines"))+
  geom_line(data = tibble(Cells = "MP38",x = c(1, 2), y = c(3.55, 3.55)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1,3), y = c(4.2, 4.2)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1,1), y = c(4.2, 3.8)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1,1), y = c(3.55, 3.0)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 1.6, y = 3.75),
            aes(x = x, y = y),
            label = "p = 0.0007",
            size = 3,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 2.1, y = 4.37),
            aes(x = x, y = y),
            label = "p = 0.0248",
            size = 3,
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65",x = c(1, 2), y = c(1.75, 1.75)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65",x = c(1,3), y = c(2.1, 2.1)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65",x = c(1,1), y = c(2.1, 1.9)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65",x = c(1,1), y = c(1.75, 1.5)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65",x = 1.55, y = 1.88),
            aes(x = x, y = y),
            label = "p = 0.004",
            size = 3,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65",x = 2, y = 2.22),
            aes(x = x, y = y),
            label = "p = 0.001",
            size = 3,
            inherit.aes = FALSE)


draw(siCPT1A.PI.AUC.bp)

ggsave("siCPT1A.PI.AUC.bp.png", plot = siCPT1A.PI.AUC.bp,
       width = 6, height = 2.75, units = "in", dpi = 300)




#read in data figure 5D
siCPT1A.HMGB1 <- read.csv("2024.09.26.siCPT1A.HMGB1.csv")


#STATS  for figure 5D
#### MP38

MP38_siCPT1A_HMGB1 <- siCPT1A.HMGB1 %>%
  filter(Cells %in% c("MP38"))

aov_siCPT1A.HMGB1.38 <- aov(HMGB1 ~ siRNA, data = MP38_siCPT1A_HMGB1)
anova_summary_siCPT1A.HMGB1.38 <- summary(aov_siCPT1A.HMGB1.38)
print(anova_summary_siCPT1A.HMGB1.38)

tukey_test_siCPT1A.HMGB1.38 <- TukeyHSD(aov_siCPT1A.HMGB1.38)
print(tukey_test_siCPT1A.HMGB1.38)

######### MP65

MP65_siCPT1A_HMGB1 <- siCPT1A.HMGB1 %>%
  filter(Cells %in% c("MP65"))

aov_siCPT1A.HMGB1.65 <- aov(HMGB1 ~ siRNA, data = MP65_siCPT1A_HMGB1)
anova_summary_siCPT1A.HMGB1.65 <- summary(aov_siCPT1A.HMGB1.65)
print(anova_summary_siCPT1A.HMGB1.65)

tukey_test_siCPT1A.HMGB1.65 <- TukeyHSD(aov_siCPT1A.HMGB1.65)
print(tukey_test_siCPT1A.HMGB1.65)



#### Plot for Figure 5D


siCPT1A.bp<-ggplot(siCPT1A.HMGB1, aes(x = siRNA, y = HMGB1, fill = siRNA)) +
  geom_boxplot(aes(group = siRNA), color = "black", outlier.size = 0) +
  geom_dotplot(binaxis = "y",  
               aes(color = siRNA), 
               binwidth = max(siCPT1A.HMGB1$HMGB1)*.04, stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center")+
  facet_grid(~ Cells, scales = "free")+  # Adjust jitter width here
  labs(title = "Secreted HMGB1",
       x = NULL, y = "HMGB1/Coomassie A.U.") +
  scale_color_manual(values = rep("black", 10),
                     guide = "none") +
  scale_fill_manual(values = c("#AA4499", "#fE732F", "#FE9929"),
                    guide = "none") +
  theme_bw()+
  scale_y_continuous(limits = c(0, max(siCPT1A.HMGB1$HMGB1)*1.3))+
  theme(plot.title = element_text(size = 12, hjust = 0.5, color = "black", margin = margin(1,0,0,0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        legend.title = element_text(size=10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.position = "right",
        legend.margin = margin(0,0,0,0),
        strip.background = element_rect("white", margin(0,0,0,0)),
        strip.text = element_text(size = 10),
        plot.margin = margin(t = -1, r = 1, b = 1, l = 1))+
  geom_line(data = tibble(Cells = "MP38", x = c(1, 2), y = c(0.4, 0.4)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1,3), y = c(0.46, 0.46)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1,1), y = c(0.2, 0.46)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 1.68, y = 0.435),
            aes(x = x, y = y),
            label = "p = 0.0003",
            size = 3,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 2.1, y = 0.497),
            aes(x = x, y = y),
            label = "p = 0.0443",
            size = 3,
            inherit.aes = FALSE)+
  
  geom_line(data = tibble(Cells = "MP65", x = c(1, 2), y = c(0.36, 0.36)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65", x = c(1,3), y = c(0.44, 0.44)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65", x = c(1,1), y = c(0.2, 0.44)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65", x = 1.62, y = 0.393),
            aes(x = x, y = y),
            label = "p = 0.008",
            size = 3,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65", x = 2, y = 0.475),
            aes(x = x, y = y),
            label = "p = 0.01",
            size = 3,
            inherit.aes = FALSE)

draw(siCPT1A.bp)
ggsave("HMGB1.siCPT1A.bp.png", plot = siCPT1A.bp,
       width = 3.3, height = 3.5, units = "in", dpi = 300)



# END figure 5


# Start Figure 6

### read in data for combined plot
auc_MPs.si<- read.csv("combindMPsidata.Rformat.Inna.final.csv")
auc_MPs.si$Treat <- as.factor(auc_MPs.si$Treat)
auc_MPs.si$siRNA <- factor(auc_MPs.si$siRNA, 
                           levels = c("Control", "GSDME#4", "GSDME#18",
                                      "GSDMD#1", "GSDMD#19"))


Treatment_siorder<-c("siCon_DMSO","siE#4_DMSO","siE#18_DMSO","siD#1_DMSO","siD#19_DMSO",
                     "siCon_ETO","siE#4_ETO","siE#18_ETO","siD#1_ETO","siD#19_ETO")

auc_MPs.si$Group <-factor(auc_MPs.si$Group, levels = Treatment_siorder)

auc_MPs.si <- auc_MPs.si %>%
  group_by(Cells) %>%
  mutate(binwidth = max(AUC) * 0.025)



#Stats for MP38
auc_38.si<- read.csv("MP38sidata.Rformat.Inna.final.csv")

auc_38.si$Treat <- as.factor(auc_38.si$Treat)
auc_38.si$siRNA <- factor(auc_38.si$siRNA, 
                          levels = c("Control", "GSDME#4", "GSDME#18",
                                     "GSDMD#1", "GSDMD#19"))


Treatment_siorder<-c("siCon_DMSO","siE#4_DMSO","siE#18_DMSO","siD#1_DMSO","siD#19_DMSO",
                     "siCon_ETO","siE#4_ETO","siE#18_ETO","siD#1_ETO","siD#19_ETO")

auc_38.si$Group <-factor(auc_38.si$Group, levels = Treatment_siorder)

# Fit the robust regression model with interaction
fit.38AUC <- lmrob(AUC ~ siRNA * Treat, data = auc_38.si)

# Obtain a summary of the model
fit_summary <- summary(fit.38AUC)
print(fit_summary)

# Extract coefficients, robust standard errors, t-statistics, and p-values
coef_table <- fit_summary$coefficients
print(coef_table)





#Stats for MP46

auc_46.si<- read.csv("MP46siGSDMs.PI.AUC.csv")

auc_46.si$Treat <- as.factor(auc_46.si$Treat)
auc_46.si$siRNA <- factor(auc_46.si$siRNA, 
                          levels = c("Control", "GSDME#4", "GSDME#18",
                                     "GSDMD#1", "GSDMD#19"))


Treatment_siorder<-c("siCon_DMSO","siE#4_DMSO","siE#18_DMSO","siD#1_DMSO","siD#19_DMSO",
                     "siCon_ETO","siE#4_ETO","siE#18_ETO","siD#1_ETO","siD#19_ETO")

auc_46.si$Group <-factor(auc_46.si$Group, levels = Treatment_siorder)
# Fit the robust regression model with interaction
fit.46AUC <- lmrob(AUC ~ siRNA * Treat, data = auc_46.si)

# Obtain a summary of the model
fit_summary.46 <- summary(fit.46AUC)
print(fit_summary.46)

# Extract coefficients, robust standard errors, t-statistics, and p-values
coef_table.46 <- fit_summary.46$coefficients
print(coef_table.46)




############### PLot for fig 6B siGSDM PI uptake combined plot #########

siGSMD.PI.combined<-ggplot(auc_MPs.si, aes(x = Treat, y = AUC, fill = siRNA)) +
  geom_boxplot(aes(group = Group), color = "black", outlier.size = 0) +
  geom_dotplot(binaxis = "y",  
               aes(color = Group),
               stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center") +
  facet_wrap(~ Cells, scales = "free_y") +
  labs(title = "PI uptake induced by ETO",
       x = NULL, y = "PI+/Phase (AUC)") +
  scale_color_manual(values = rep("black", 10),
                     guide = "none") +
  scale_fill_manual(values = c("#AA4499", "#108A8D", "#00BF99", "#fE732F", "#FE9929"),
                    guide = guide_legend(title = "siRNA")) +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2)))+ # need to fix
  theme(plot.title = element_text(size = 10, hjust = 0.5, color = "black", margin = margin(1,0,0,0)),
        axis.text.x = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        legend.title = element_text(size=8, hjust = 0.5, color = "black", face = "bold"),
        legend.text = element_text(size = 8, color = "black"),
        legend.position = "right",
        legend.margin = margin(0,0,0,0),
        strip.background = element_rect("white", margin(0,0,0,0)),
        strip.text = element_text(size = 10),
        plot.margin = margin(t = -1, r = 1, b = 1, l = 1),
        legend.box.spacing = unit(0.1, "lines"),
        legend.key.width = unit(1, "lines"))+
  
  geom_line(data = tibble(Cells = "MP38",x = c(.7, 1.7), y = c(1.5, 1.5)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1.7,1.85), y = c(1.6, 1.6)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1.7,2.), y = c(1.8, 1.8)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 1.2, y = 1.6),
            aes(x = x, y = y),
            label = "p = 2.7e-06",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 1.8, y = 1.7),
            aes(x = x, y = y),
            label = "p = 0.07",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 1.85, y = 1.9),
            aes(x = x, y = y),
            label = "p = 0.03",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP46",x = c(.7, 1.7), y = c(6.2, 6.2)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP46",x = c(1.7,1.85), y = c(6.7, 6.7)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP46",x = c(1.7,2.), y = c(7.7, 7.7)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP46",x = 1.15, y = 6.6),
            aes(x = x, y = y),
            label = "p = 7.34e-10",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP46",x = 1.75, y = 7.1),
            aes(x = x, y = y),
            label = "p = 0.004",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP46",x = 1.85, y = 8.1),
            aes(x = x, y = y),
            label = "p = 0.0006",
            size = 2.5,
            inherit.aes = FALSE)

draw(siGSMD.PI.combined)

ggsave("siGSMD.PI.combined.png", plot = siGSMD.PI.combined,
       width = 5.4, height = 2.8, units = "in", dpi = 300)


#HMGB1 data for panel 6D
### STats

### MP38
MP38.siGs.ETO.noD<- read.csv("HMGB1.MP38.siGSDMs.ETO.noD.csv")



MP38.siGs.ETO.noD$siRNA <- factor(MP38.siGs.ETO.noD$siRNA, 
                                  levels = c("Control", "siGSDME #4", "siGSDME #18"))

MP38.siGs.ETO.noD$Treat <- factor(MP38.siGs.ETO.noD$Treat, 
                                  levels = c("DMSO", "ETO"))


anova_result.38.noD <- aov(HMGB1 ~ siRNA * Treat, data = MP38.siGs.ETO.noD)
summary(anova_result.38.noD)
tukey_result.38.noD <- TukeyHSD(anova_result.38.noD)
print(tukey_result.38.noD)

# MP65

MP65.siGs.ETO.noD<- read.csv("HMGB1.MP65.siGSDMs.ETO.noD.csv")

MP65.siGs.ETO.noD$siRNA <- factor(MP65.siGs.ETO.noD$siRNA, 
                                  levels = c("Control", "siGSDME #4", "siGSDME #18"))

MP65.siGs.ETO.noD$Treat <- factor(MP65.siGs.ETO.noD$Treat, 
                                  levels = c("DMSO", "ETO"))

anova_result.65.noD <- aov(HMGB1 ~ siRNA * Treat, data = MP65.siGs.ETO.noD)
summary(anova_result.65.noD)
tukey_result.65.noD <- TukeyHSD(anova_result.65.noD)
print(tukey_result.65.noD)


### combined data for plot

siGSDME.HMGB1 <- read.csv("HMGB1.MPs.siGSDMs.ETO.noD.csv")

siGSDME.HMGB1$Treat <- as.factor(siGSDME.HMGB1$Treat)
siGSDME.HMGB1$siRNA <- factor(siGSDME.HMGB1$siRNA, 
                              levels = c("Control", "GSDME#4", "GSDME#18"))


Treatment_siorder<-c("DMSO_Control","DMSO_siGSDME #4","DMSO_siGSDME #18",
                     "ETO_Control","ETO_siGSDME #4","ETO_siGSDME #18")

siGSDME.HMGB1$Group <-factor(siGSDME.HMGB1$Group, levels = Treatment_siorder)





# first set draw dims #
draw <- function(plot, x_in = 3.5, y_in = 2.8) {
  
  grid::grid.newpage()
  
  grid::rectGrob(gp = grid::gpar(fill = "gray")) |>
    grid::grid.draw()
  
  grid::viewport(width  = ggplot2::unit(x_in, "in"), 
                 height = ggplot2::unit(y_in, "in")) |>
    grid::pushViewport()
  
  ggplot2::ggplot_build(plot) |>
    ggplot2::ggplot_gtable() |>
    grid::grid.draw()
}



siGSDME.HMGB1.bp<-ggplot(siGSDME.HMGB1, aes(x = Treat, y = HMGB1, fill = siRNA)) +
  geom_boxplot(aes(group = Group), color = "black", outlier.size = 0) +
  geom_dotplot(binaxis = "y",  
               aes(color = Group), 
               stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center")+
  facet_wrap(~ Cells, scales = "free")+  # Adjust jitter width here
  labs(title = "Secreted HMGB1",
       x = NULL, y = "HMGB1/PonceauS A.U.") +
  scale_color_manual(values = rep("black", 6),
                     guide = "none") +
  scale_fill_manual(values = c("#AA4499", "#108A8D", "#00BF99"),
                    guide = guide_legend(title = "siRNA")) +
  theme_bw()+
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0, 0.2)))+
  theme(plot.title = element_text(size = 12, hjust = 0.5, color = "black", margin = margin(1,0,0,0)),
        axis.text.x = element_text(size=10, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        legend.title = element_text(size=8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.position = "bottom",
        legend.margin = margin(-.5,0,0,0),
        strip.background = element_rect("white", margin(0,0,0,0)),
        strip.text = element_text(size = 10),
        plot.margin = margin(t = -1, r = 1, b = 1, l = 1),
        legend.box.spacing = unit(0.1, "lines"),
        legend.key.width = unit(1, "lines"))+
  
  geom_line(data = tibble(Cells = "MP38", x = c(.75, 1.75), y = c(0.235, 0.235)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1.75,2), y = c(0.265, 0.265)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP38",x = c(1.75,2.25), y = c(0.305, 0.305)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 1.25, y = 0.253),
            aes(x = x, y = y),
            label = "p = 5.93E-05",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 1.875, y = 0.285),
            aes(x = x, y = y),
            label = "p = 0.0049",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP38",x = 2, y = 0.325),
            aes(x = x, y = y),
            label = "p = 0.0006",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65", x = c(.75, 1.75), y = c(2.45, 2.45)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65", x = c(1.75,2), y = c(2.75, 2.75)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(Cells = "MP65", x = c(1.75,2.25), y = c(3.03, 3.03)),
            aes(x = x, y = y),
            size = .75, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65",x = 1.25, y = 2.635),
            aes(x = x, y = y),
            label = "p = 0.001",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65",x = 1.875, y = 2.88),
            aes(x = x, y = y),
            label = "p = 0.006",
            size = 2.5,
            inherit.aes = FALSE)+
  geom_text(data = tibble(Cells = "MP65",x = 2, y = 3.2),
            aes(x = x, y = y),
            label = "p = 0.027",
            size = 2.5,
            inherit.aes = FALSE)


draw(siGSDME.HMGB1.bp)
ggsave("siGSDME.HMGB1.bp.png", plot = siGSDME.HMGB1.bp,
       width = 3.5, height = 2.8, units = "in", dpi = 300)

#### END figure 6



#### Start figure 7


Morrison_data <- read_xlsx("Morrison data with pyroptosis score simple for Glenn.xlsx")


df <- Morrison_data[, c("treat sort", "PyroScore", "subject.id", "response")]
df$`treat sort`[df$`treat sort` == 1] <- "Pre" 
df$`treat sort`[df$`treat sort` == 2] <- "On" 
df$`treat sort` <- factor(df$`treat sort`, levels = c("Pre", "On"))
df$`response`[df$`response` == "PD"] <- "Progressive" 
df$`response`[df$`response` == "SD"] <- "Stable" 
df$`response` <- factor(df$`response`, levels = c("Progressive", "Stable"))
df <- df[!(df$subject.id %in% c("CA209064-4-64174_064")), ] # remove row without pair
df$subject.id <- as.factor(df$subject.id)

# Combine 'treat sort' and 'response' into a new variable 'group'
df$group <- factor(paste(df$`treat sort`, df$response, sep = "_"),
                   levels = c("Pre_Progressive", "On_Progressive", "Pre_Stable", "On_Stable"))


draw <- function(plot, x_in = 4, y_in = 2.4) {
  
  grid::grid.newpage()
  
  grid::rectGrob(gp = grid::gpar(fill = "gray")) |>
    grid::grid.draw()
  
  grid::viewport(width  = ggplot2::unit(x_in, "in"), 
                 height = ggplot2::unit(y_in, "in")) |>
    grid::pushViewport()
  
  ggplot2::ggplot_build(plot) |>
    ggplot2::ggplot_gtable() |>
    grid::grid.draw()
}

##### plot Figure 7B Note: stats were performed by Inna Chervoneva




PScore.ICi.plot.stats<-ggplot(df, aes(x = group, y = PyroScore, fill = response)) +
  geom_boxplot(alpha = 0.25)+
  geom_line(aes(group = subject.id, color = response), size = 1.2, alpha = 1) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.0), color = "black") +
  geom_signif(comparisons = list(c("Pre_Stable", "On_Stable")), 
              map_signif_level = TRUE, 
              annotations = "p = 0.016",
              y_position = 0.5,
              textsize = 3,
              vjust = -0.2) +
  geom_signif(comparisons = list(c("Pre_Progressive", "On_Progressive")), 
              map_signif_level = TRUE, 
              annotations = "n.s.",
              y_position = 0.0,
              textsize = 3,
              vjust = -0.2) +
  labs(title = "PScore response to new ICb treatment",
       x = "ICb treatment status",
       y = "PScore (GSVA)",
       color = "Disease Response", # Change legend title
       fill = "Disease Response") +
  theme_bw() +
  scale_x_discrete(labels = c("Pre", "On", "Pre", "On"))+
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.1, 0.13)))+
  theme(legend.title = element_text(hjust = 0.5,size=10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size=10, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        axis.title.x = element_text(size=10, color = "black"),
        plot.title = element_text(size = 12, hjust = 0.1, color = "black")) 

draw(PScore.ICi.plot.stats)


ggsave("PScore.ICi.plot.stats.png", plot = PScore.ICi.plot.stats,
       width = 4, height = 2.4, units = "in", dpi = 300)






#### Fig S2 

# Panel S2A
#STats and plots for CPT1AinhibDN by SCNA


UVM.data <- read.csv("UVM_GSVAs_r.csv")

UVM.data$SCNA <- factor(UVM.data$SCNA)

aov_GSVA.CPT1A.SCNA <- aov(CPT1AinhibDN ~ SCNA, data = UVM.data)
anova_summary_GSVA.CPT1A.SCNA <- summary(aov_GSVA.CPT1A.SCNA)
print(anova_summary_GSVA.CPT1A.SCNA)

tukey_test_GSVA.CPT1A.SCNA <- TukeyHSD(aov_GSVA.CPT1A.SCNA)
print(tukey_test_GSVA.CPT1A.SCNA)

UVM.CPT1A.surv.data$SCNA <-as.factor(UVM.CPT1A.surv.data$SCNA)

CPT1A.SNCA.BP<-ggplot(UVM.CPT1A.surv.data, aes(x = SCNA, y = CPT1AinhibDN, fill = SCNA)) +
  geom_boxplot(aes(group = SCNA), color = "black") +
  geom_dotplot(binaxis = "y",  
               aes(color = SCNA), 
               binwidth = max(UVM.CPT1A.surv.data$CPT1AinhibDN)*.06, stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center") +  # Adjust jitter width here
  labs(title = "CPT1A 'activity' in UM SCNA clusters",
       x = NULL, y = "CPT1A-inhib. down sig. (A.U.)") +
  scale_color_manual(values = rep("black", 10),
                     guide = "none") +
  scale_fill_manual(values = c("#FE9929", "#D7191C", "#6BAED6","#004C99"),
                    guide = guide_legend(title = "SCNA Cluster")) +
  theme_bw()+
  scale_y_continuous(limits = c(min(UVM.CPT1A.surv.data$CPT1AinhibDN), 
                                max(UVM.CPT1A.surv.data$CPT1AinhibDN)*1.75))+
  theme(axis.text.x = element_text(size=10, color = "black"),
        plot.title=element_text(size=12, hjust = 0.3, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=8, color = "black"),
        legend.text = element_text(size=10, color = "black"),
        legend.title = element_text(size=9,hjust = 0.5, color = "black"),
        legend.position = "right",
        legend.margin = margin(-5,-5,-5,-5))+
  geom_line(data = tibble(x = c(3, 4), y = c(0.5, 0.5)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(3, 3), y = c(0.45, 0.5)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(2, 2), y = c(0.37, 0.42)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(2,3), y = c(0.42, 0.42)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1,4), y = c(.7, 0.7)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(1,1), y = c(.65, 0.7)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(2,4), y = c(.6, 0.6)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_line(data = tibble(x = c(2,2), y = c(.55, 0.6)),
            aes(x = x, y = y),
            linewidth = .8, color = 'grey41',
            inherit.aes = FALSE)+
  geom_text(data = tibble(x = 3.5, y = 0.548+0.008),
            aes(x = x, y = y),
            label = "p = 0.006",
            size = 2.4,
            color = "black",
            inherit.aes = FALSE)+
  geom_text(data = tibble(x = 2.45, y = 0.468+0.006),
            aes(x = x, y = y),
            label = "p = 0.004",
            size = 2.4,
            color = "black",
            inherit.aes = FALSE)+
  geom_text(data = tibble(x = 2.5, y = 0.743+0.02),
            aes(x = x, y = y),
            label = "p < 0.001",
            size = 2.4,
            color = "black",
            inherit.aes = FALSE)+
  geom_text(data = tibble(x = 3, y = 0.648+0.01),
            aes(x = x, y = y),
            label = "p < 0.001",
            size = 2.4,
            color = "black",
            inherit.aes = FALSE)


draw(CPT1A.SNCA.BP)
ggsave("CPT1A.SNCA.BP.png", plot = CPT1A.SNCA.BP,
       width = 3, height = 2.51, units = "in", dpi = 300)



# Attempt to find the correct number of clusters based on CPT1A activity
# CODE For km clusters (figure s2)
# CODE For km clusters (figure s2)
# CODE For km clusters (figure s2)
# Scale the data

data_scaled <- scale(UVM.CPT1A.surv.data$CPT1AinhibDN)

# ELBOW PLOT Panel S2B
fviz_nbclust(data_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2) +  # Adjust based on the elbow plot
  labs(subtitle = "Elbow method")


# Perform K-Means Clustering
set.seed(123)
kmeans_result <- kmeans(data_scaled, centers = 3, nstart = 25)  # Adjust 'centers' based on the elbow plot
# Add the cluster labels to the original data
UVM.CPT1A.surv.data$Cluster <- as.factor(kmeans_result$cluster)


# Calculate the mean value for each cluster
cluster_means <- aggregate(data_scaled, by = list(cluster = kmeans_result$cluster), FUN = mean)
print(cluster_means)

# Sort the clusters by their mean values in increasing order
cluster_means <- cluster_means[order(cluster_means$V1, decreasing = FALSE), ]
print(cluster_means)

# Create a mapping from old cluster numbers to new cluster numbers
old_to_new <- setNames(seq_along(cluster_means$cluster), cluster_means$cluster)

# Reorder the clusters in the original data
UVM.CPT1A.surv.data$OrderedCluster <- as.factor(sapply(UVM.CPT1A.surv.data$Cluster, function(x) old_to_new[as.character(x)]))

# Verify the reordering
print(UVM.CPT1A.surv.data$OrderedCluster)
print(UVM.CPT1A.surv.data$Cluster)


UVM.CPT1A.surv.data$Cluster == UVM.CPT1A.surv.data$OrderedCluster

levels(UVM.CPT1A.surv.data$OrderedCluster) <- c("Cluster 1", "Cluster 2", "Cluster 3")

# Verify the changes
print(UVM.CPT1A.surv.data$OrderedCluster)


############### plot for S2C
Q

CPT1A.cluster.BP<-ggplot(UVM.CPT1A.surv.data, aes(x = OrderedCluster, y = CPT1AinhibDN, fill = OrderedCluster)) +
  geom_boxplot(aes(group = OrderedCluster), color = "black") +
  geom_dotplot(binaxis = "y",  
               aes(color = OrderedCluster), 
               binwidth = max(UVM.CPT1A.surv.data$CPT1AinhibDN)*.08, stroke = 2,
               position = position_dodge(0.75),
               stackdir = "center") +  # Adjust jitter width here
  labs(title = "CPT1A 'Activity' in k-mean clusters",
       x = NULL, y = "CPT1A-inhib. down sig. (A.U.)",
       color = "CPT1A Activity Cluster", # Change legend title
       fill = "CPT1A Activity Cluster") +
  scale_color_manual(values = rep("black", 3),
                     guide = "none") +
  scale_fill_manual(values = c("#FE9929", "#D7191C", "#004C99"),
                    guide = guide_legend(title = "CPT1A Act. \n KM-Cluster")) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.35)) +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(hjust = .4,size = 10, color = "black"),
        legend.position = "right",
        legend.margin = margin(-0.5,-0.5,-5,-5),
        axis.title.y = element_text(size = 9, color = "black"))

# margin(t = -1, r = 1, b = 1, l = 1)

draw(CPT1A.cluster.BP)
ggsave("CPT1A.cluster.BP.png", plot = CPT1A.cluster.BP,
       width = 3.15, height = 2.4, units = "in", dpi = 300)






### Data prep for panel S2D

UVM.CPT1A.surv.data$SCNA <- as.factor(UVM.CPT1A.surv.data$SCNA)

SNCA.counts <- UVM.CPT1A.surv.data %>%
  group_by(OrderedCluster, SCNA) %>%
  summarise(count = n()) %>%
  ungroup()

SNCA.counts <- SNCA.counts %>%
  group_by(OrderedCluster) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()


### Plot for panel S2D
SCNA.count.CPT1Aclust<-ggplot(UVM.CPT1A.surv.data, aes(x = OrderedCluster, fill = SCNA)) +
  geom_bar(position = "fill", alpha = 0.75) +
  geom_text(data = SNCA.counts, aes(x = OrderedCluster, y = prop, label = count),
            position = position_stack(vjust = 0.5), size = 2.5, color = "black") +
  labs(title = "SCNA Percentage by CPT1A Cluster",
       x = "CPT1A activity Cluster",
       y = "Percentage",
       fill = "SCNA") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#FE9929", "#D7191C", "#6BAED6","#004C99"),
                    guide = guide_legend(title = "SCNA Cluster")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, color = "black", hjust = 0.4),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(hjust = .5,size = 10, color = "black"),
        legend.position = "right",
        legend.margin = margin(-0.5,-1,-2,-5),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"))



draw(SCNA.count.CPT1Aclust)
ggsave("SCNA.count.CPT1Aclust.png", plot = SCNA.count.CPT1Aclust,
       width = 3.15, height = 3.15, units = "in", dpi = 300)


# end fig S2






