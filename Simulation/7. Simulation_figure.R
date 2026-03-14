############################################
## MGWCR
## Figure2: Bandwidth distribution
## Figure3: MSE distribution 
## Figure4: Parameter surfaces and regional estimates 
## Figure5: Analytical and simulated local standard errors 
############################################

library(ggplot2)
library(openxlsx)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(scales)
select <- dplyr::select


####################################
######## Figure 2  #################
####################################
# Bandwidth

## Simulation 1
bdw1 <- GWCR_bws_fix1 %>% left_join(MGWCR2_bws_fix1, by='r')
bdw1$r <- NULL

p <- ggplot() +
  geom_boxplot(aes(y = bdw2$gwcr, x = "GWCR"), fill = "#80796bff", outlier.colour = "grey50") +
  geom_boxplot(aes(y = bdw2$mgwcr1, x = "MGWCR_beta1"), fill = "#80796bff", outlier.colour = "grey50") +
  geom_boxplot(aes(y = bdw2$mgwcr2, x = "MGWCR_beta2"), fill = "#80796bff", outlier.colour = "grey50") +
  theme_classic() +
  labs(title = "Box plots",
       y = "Bandwidth") +
  scale_x_discrete(labels = c(expression(paste(GW, " ",Cox)),
                              expression(paste(MGWCR, " ", beta[1])),
                              expression(paste(MGWCR, " ", beta[2])))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", vjust = 1.2,size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0), size=18))
print(p)



## Simulation 2
bdw1 <- GWCR_bws_fix2 %>% left_join(MGWCR_bws_fix2, by='r')
bdw1$r <- NULL

p <- ggplot() +
  geom_boxplot(aes(y = bdw1$gwcr, x = "GW Cox"), fill = "#80796bff", outlier.colour = "grey50") +
  geom_boxplot(aes(y = bdw1$mgwcr1, x = "MGWCR_beta1"), fill = "#80796bff", outlier.colour = "grey50") +
  geom_boxplot(aes(y = bdw1$mgwcr2, x = "MGWCR_beta2"), fill = "#80796bff", outlier.colour = "grey50") +
  geom_boxplot(aes(y = bdw1$mgwcr3, x = "MGWCR_beta3"), fill = "#80796bff", outlier.colour = "grey50") +
  theme_classic() +
  labs(title = "Box plots",
       y = "Bandwidth") +
  scale_x_discrete(labels = c(expression(paste(GW, " ",Cox)),
                              expression(paste(MGWCR, " ", beta[1])),
                              expression(paste(MGWCR, " ", beta[2])),
                              expression(paste(MGWCR, " ", beta[3])))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", vjust = 1.2,size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0), size=18))
print(p)



####################################
######## Figure3  ##################
####################################
# MAB, MSE 
## Simulation 1
df.mse <- cbind(mse_glb1,mse_gwcr1,mse_mgwcr1)

df_long <- melt(df.mse)
df_long$group <- rep(c("GLB", "GWCR", "MGWCR"), each = 1000, times = 2)
df_long$beta <- rep(c(1, 2), each = 3000)

plots <- list()
for (beta_value in 1:2) {
  df_beta <- df_long %>% filter(beta == beta_value)
  
  p <- ggplot(df_beta, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("GLB" = "#374e5550", "GWCR" = "#6a659950", "MGWCR" = "#df8f4450")) +
    theme_classic() +
    labs(x = "", y = "MSE") +
    scale_x_discrete(labels = c("GLB", "GW Cox", "MGWCR")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", vjust = 1.2,size=13),
          axis.text.y = element_text(color = "black",size=13),
          legend.position = "none",
          axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0), size=15))
  
  plots[[beta_value]] <- p
}

grid.arrange(grobs = plots, ncol = 2)

png(filename="~/figure3_1.png",width=466,height=300,unit="px",bg="transparent")
grid.arrange(grobs = plots, ncol = 2)
dev.off()


## Simulation 2
df.mse <- cbind(mse_glb2,mse_gwcr2,mse_mgwcr2)

df_long <- melt(df.mse)
df_long$group <- rep(c("GLB", "GWCR", "MGWCR"), each = 1000, times = 3)
df_long$beta <- rep(c(1, 2, 3), each = 3000)


plots <- list()
for (beta_value in 1:3) {
  df_beta <- df_long %>% filter(beta == beta_value)
  
  p <- ggplot(df_beta, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("GLB" = "#374e5550", "GWCR" = "#6a659950", "MGWCR" = "#df8f4450")) +
    theme_classic() +
    labs(x = "", y = "MSE") +
    scale_x_discrete(labels = c("GLB", "GW Cox", "MGWCR")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", vjust = 1.2,size=13),
          axis.text.y = element_text(color = "black",size=13),
          legend.position = "none",
          axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0), size=15))
  
  plots[[beta_value]] <- p
}

grid.arrange(grobs = plots, ncol = 3)


png(filename="~/figure3_2.png",width=700,height=300,unit="px",bg="transparent")
grid.arrange(grobs = plots, ncol = 3)
dev.off()



#############################
#####  Figure 4  #############
#############################
data_list <- list(
  betas_true1,   
  betas_true2,  
  betas_glb1, 
  betas_glb2,  
  betas_gwcr1, 
  betas_gwcr2,  
  betas_mgwcr1, 
  betas_mgwcr2
)

mg_fun <- function(idx, var, lim){
  griddf <- as.data.frame(cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                                      lngcoords = seq(1, 25, 1))), loc = 1:625))
  values <- data_list[[idx]]
  
  griddf0 <- merge(griddf, values, by = "loc")
  plot <- ggplot(griddf0, aes(x = latcoords, y = lngcoords, fill = get(var))) +
    geom_tile() +
    scale_fill_gradientn(colors = color, limits = lim) + 
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) 
  return(plot)
}

color <- rev(brewer.pal(10, "RdBu"))

#### Simulation 1
plot1 <- mg_fun(1, "V1", c(0, 2.1)) + theme(legend.position = "none") # True1
plot2 <- mg_fun(1, "V2", c(-1, 1)) + theme(legend.position = "none")
plot3 <- mg_fun(3, "V1", c(0, 2.1)) + theme(legend.position = "none") # Glb1
plot4 <- mg_fun(3, "V2", c(-1, 1)) + theme(legend.position = "none")
plot5 <- mg_fun(5, "V1", c(0, 2.1)) + theme(legend.position = "none") # GWCR1
plot6 <- mg_fun(5, "V2", c(-1, 1)) + theme(legend.position = "none")
plot7 <- mg_fun(7, "V1", c(0, 2.1)) + theme(legend.position = "none") # MGWCR1
plot8 <- mg_fun(7, "V2", c(-1, 1)) + theme(legend.position = "none")

legend1 <- get_legend(mg_fun(1, "V1", c(0, 2.1)) + theme(legend.position = "right"))
legend2 <- get_legend(mg_fun(5, "V2", c(-1, 1)) + theme(legend.position = "right"))

legends_combined <- plot_grid(legend1, legend2, ncol = 1)
combined_plot <- plot_grid(
  plot1, plot2, plot3, plot4,
  plot5, plot6, plot7, plot8,
  ncol = 4, byrow = FALSE,
  align="hv", axis="tblr",
  rel_widths=c(1,1,1,1), rel_heights=c(1,1,1)
)
final_plot <- plot_grid(combined_plot, legends_combined, ncol = 2, rel_widths = c(4, 0.4))
print(final_plot)

#### Simulation 2
plot1 <- mg_fun(2, "V1", c(0, 2.2)) + theme(legend.position = "none") # True2
plot2 <- mg_fun(2, "V2", c(0, 2.2)) + theme(legend.position = "none")
plot3 <- mg_fun(2, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")
plot4 <- mg_fun(4, "V1", c(0, 2.2)) + theme(legend.position = "none") # Glb2
plot5 <- mg_fun(4, "V2", c(0, 2.2)) + theme(legend.position = "none")
plot6 <- mg_fun(4, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")
plot7 <- mg_fun(6, "V1", c(0, 2.1)) + theme(legend.position = "none") # GWCR2
plot8 <- mg_fun(6, "V2", c(0, 2.2)) + theme(legend.position = "none")
plot9 <- mg_fun(6, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")
plot10 <- mg_fun(8, "V1", c(0, 2.1)) + theme(legend.position = "none") # MGWCR2
plot11 <- mg_fun(8, "V2", c(0, 2.2)) + theme(legend.position = "none")
plot12 <- mg_fun(8, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")

legend1 <- get_legend(plot1 + theme(legend.position = "right", legend.direction = "vertical"))
legend2 <- get_legend(plot2 + theme(legend.position = "right", legend.direction = "vertical"))
legend3 <- get_legend(plot3 + theme(legend.position = "right", legend.direction = "vertical"))

legends_combined <- plot_grid(legend1, legend2, legend3, ncol = 1)
combined_plot <- plot_grid(
  plot1, plot2, plot3, plot4,
  plot5, plot6, plot7, plot8,
  plot9, plot10, plot11, plot12,
  ncol = 4, byrow = FALSE,
  align="hv",axis="tblr",rel_widths=c(1,1,1,1),rel_heights=c(1,1,1)
)
final_plot <- plot_grid(combined_plot,legends_combined,ncol = 2,rel_widths = c(4, 0.4))
print(final_plot)



##############################
###### Figure 5 ##############
##############################

# Empirical Standard Error (sd) list
se_empirical_list <- list(
  sdbetas_mgwcr2, # 1
  sdbetas_mgwcr1, # 2
  sdbetas_gwcr2,  # 3
  sdbetas_gwcr1,  # 4
  sdbetas_glb2,   # 5
  sdbetas_glb1    # 6
)

# Average Standard Error (se) list
se_average_list <- list(
  seavg_mgwcr2,   # 1
  seavg_mgwcr1,   # 2
  seavg_gwcr2,    # 3
  seavg_gwcr1,    # 4
  seavg_glb2,     # 5
  seavg_glb1      # 6
)

mg_se <- function(idx, var, lim){
    griddf <- as.data.frame(cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                                      lngcoords = seq(1, 25, 1))), loc = 1:625))
    values1 <- se_empirical_list[[idx]]
  values2 <- se_average_list[[idx]]
  
  griddf1 <- merge(griddf, values1, by = "loc")
  griddf2 <- merge(griddf, values2, by = "loc")
  
  # Empirical SE plot (p1)
  p1 <- ggplot(griddf1, aes(x = latcoords, y = lngcoords, fill = get(var))) +
    geom_tile() +
    scale_fill_gradient(low = "snow", high = "darkred", limits = lim, labels = number_format(accuracy = 0.01)) +
    theme_void() +
    theme(legend.title = element_blank(), legend.position = "none")
  
  # Average SE plot (p2)
  p2 <- ggplot(griddf2, aes(x = latcoords, y = lngcoords, fill = get(var))) +
    geom_tile() +
    scale_fill_gradient(low = "snow", high = "darkred", limits = lim, labels = number_format(accuracy = 0.01)) +
    theme_void() +
    theme(legend.title = element_blank(), legend.position = "none")
  
  print(c(max(griddf1[,var]), max(griddf2[,var]), min(griddf1[,var]), min(griddf2[,var])))
  return(grid.arrange(p1, p2, ncol = 2))
}

## Simulation 1 
p1_se <- mg_se(2, "V1", c(0, 0.1))
p2_se <- mg_se(2, "V2", c(0, 0.1))
grid.arrange(p1_se, p2_se, ncol = 1)


## Simulation 2 
p1_se2 <- mg_se(1, "V1", c(0, 0.2))
p2_se2 <- mg_se(1, "V2", c(0, 0.2))
p3_se2 <- mg_se(1, "V3", c(0, 0.2))
grid.arrange(p1_se2, p2_se2, p3_se2, ncol = 1)


