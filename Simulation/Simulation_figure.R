############################################
###### MGWCR
###### Bandwidth figure - box plot
###### MAB, MSE figure - box plot
###### Standard error figure
###### Parameters, estimates surface
############################################

library(ggplot2)
library(openxlsx)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(dplyr)
select <- dplyr::select
library(tidyr)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(scales)


####################################
######## Figure 2  #################
####################################
# Bandwidth


## Simulation 1
bdw2 <- GWCR_bws_fix %>% left_join(MGWCR_bws_fix, by='r')
bdw2$r <- NULL

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

png(filename="~/bws_s1.png",width=700,height=300,unit="px",bg="transparent")
p
dev.off()

## Simulation 2
bdw2 <- GWCR_bws_fix %>% left_join(MGWCR_bws_fix, by='r')
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


png(filename="C:/Users/User/Dropbox/MGWCR/Results/Simulation/Multi_normal/figure/bws_s2.png",width=700,height=300,unit="px",bg="transparent")
p
dev.off()


####################################
######## Figure3  ##################
####################################
# MAB, MSE 
df_all <- cbind(mse_gwcr1,mse_gwcr2,mmse_gwcr1,mse_mgwcr2)
df.mab <- df_all %>% select(GLB_MAB1, GWCR_MAB1, MGWCR_MAB1, GLB_MAB2, GWCR_MAB2, MGWCR_MAB2)
df.mse <- df_all %>% select(GLB_MSE1, GWCR_MSE1, MGWCR_MSE1, GLB_MSE2, GWCR_MSE2, MGWCR_MSE2)

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
df_all <- cbind(mse_gwcr1,mse_gwcr2,mmse_gwcr1,mse_mgwcr2)
df.mab <- df_all %>% select(GLB_MAB1, GWCR_MAB1, MGWCR_MAB1, GLB_MAB2, GWCR_MAB2, MGWCR_MAB2)
df.mse <- df_all %>% select(GLB_MSE1, GWCR_MSE1, MGWCR_MSE1, GLB_MSE2, GWCR_MSE2, MGWCR_MSE2)

df.mab <- df_all %>% select(GLB_MAB1, GWCR_MAB1, MGWCR_MAB1, GLB_MAB2, GWCR_MAB2, MGWCR_MAB2, 
                            GLB_MAB3, GWCR_MAB3, MGWCR_MAB3)
df.mse <- df_all %>% select(GLB_MSE1, GWCR_MSE1, MGWCR_MSE1, GLB_MSE2, GWCR_MSE2, MGWCR_MSE2, 
                            GLB_MSE3, GWCR_MSE3, MGWCR_MSE3)

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
#  estimates

mg_fun <- function(sht, var, lim){
  griddf <- as.data.frame(cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                                      lngcoords = seq(1, 25, 1))), loc = 1:625))
  values <- read.xlsx("~/simulation_estimates.xlsx", sheet = sht)
  griddf0 <- merge(griddf, values, by = "loc")
  plot <- ggplot(griddf0, aes(x = latcoords, y = lngcoords, fill = get(var))) +
    geom_tile() +
    scale_fill_gradientn(colors = color,limits = lim) + 
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_blank(),
      axis.line =  element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) 
  return(plot)
}

#### Simulation 1
griddf <- as.data.frame(cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                                    lngcoords = seq(1, 25, 1))), loc = 1:625))
values <- read.xlsx("~/simulation_estimates.xlsx", sheet = 1)
griddf0 <- merge(griddf, values, by = "loc")

color <- rev(brewer.pal(10, "RdBu"))
plot1 <- mg_fun(7, "V1", c(-0.7, 3.6)) + theme(legend.position = "none")
plot2 <- mg_fun(7, "V2", c(-0.7, 3.6)) + theme(legend.position = "none")
plot3 <- mg_fun(7, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")
plot4 <- mg_fun(5, "V1", c(-0.7, 3.6)) + theme(legend.position = "none")
plot5 <- mg_fun(5, "V2", c(-0.7, 3.6)) + theme(legend.position = "none")
plot6 <- mg_fun(5, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")
plot7 <- mg_fun(3, "V1", c(-0.7, 3.6)) + theme(legend.position = "none")
plot8 <- mg_fun(3, "V2", c(-0.7, 3.6)) + theme(legend.position = "none")
plot9 <- mg_fun(3, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")
plot10 <- mg_fun(1, "V1", c(-0.7, 3.6)) + theme(legend.position = "none")
plot11 <- mg_fun(1, "V2", c(-0.7, 3.6)) + theme(legend.position = "none")
plot12 <- mg_fun(1, "V3", c(-0.7, 3.6)) + theme(legend.position = "none")
legend <- get_legend(
  plot1 + theme(legend.position = "right", legend.direction = "vertical")
)
combined_plot <- plot_grid(
  plot1, plot2, plot3, plot4,
  plot5, plot6, plot7, plot8,
  plot9, plot10, plot11, plot12,
  ncol = 4, byrow = FALSE,
  align="hv",axis="tblr",rel_widths=c(1,1,1,1),rel_heights=c(1,1,1)
)
final_plot <- plot_grid(combined_plot,legend,ncol = 2,rel_widths = c(4, 0.4))
print(final_plot)



#### Simulation 2
plot1 <- mg_fun(8, "V1", c(0, 2.1)) + theme(legend.position = "none")
plot2 <- mg_fun(8, "V2", c(-1, 1)) + theme(legend.position = "none")
plot3 <- mg_fun(6, "V1", c(0, 2.1)) + theme(legend.position = "none")
plot4 <- mg_fun(6, "V2", c(-1, 1)) + theme(legend.position = "none")
plot5 <- mg_fun(4, "V1", c(0, 2.1)) + theme(legend.position = "none")
plot6 <- mg_fun(4, "V2", c(-1, 1)) + theme(legend.position = "none")
plot7 <- mg_fun(2, "V1", c(0, 2.1)) + theme(legend.position = "none")
plot8 <- mg_fun(2, "V2", c(-1, 1)) + theme(legend.position = "none")
legend1 <- get_legend(
  mg_fun(8, "V1", c(0, 2.1)) + 
    theme(legend.position = "right", legend.direction = "vertical")
)
legend2 <- get_legend(
  mg_fun(8, "V2", c(-1, 1)) + 
    theme(legend.position = "right", legend.direction = "vertical")
)
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


##############################
###### Figure 5 ##############
##############################
# Standard errors

griddf <- as.data.frame(cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                                    lngcoords = seq(1, 25, 1))), loc = 1:625))

mg_se <- function(sht, var, lim){
  
  griddf <- as.data.frame(cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                                      lngcoords = seq(1, 25, 1))), loc = 1:625))
  values1 <- read.xlsx("~/simulation_empirical_standard_error.xlsx", sheet = sht)
  values2 <- read.xlsx("~/simulation_average_standard_error.xlsx", sheet = sht)
  griddf1 <- merge(griddf, values1, by = "loc")
  griddf2 <- merge(griddf, values2, by = "loc")
  
  # scale_fill_gradientn(colors = rainbow(5), limits = lim) + 
  p1 <- ggplot(griddf1, aes(x = latcoords, y = lngcoords, fill = get(var))) +
    geom_tile() +
    scale_fill_gradient(low = "snow", high = "darkred", limits = lim,labels = number_format(accuracy = 0.01)) +
    # scale_fill_gradientn(colors = color,limits = lim) + 
    # scale_fill_viridis_c(option = color, limits = lim) +
    theme_void()+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")
  
  
  p2 <- ggplot(griddf2, aes(x = latcoords, y = lngcoords, fill = get(var))) +
    geom_tile() +
    scale_fill_gradient(low = "snow", high = "darkred", limits = lim,labels = number_format(accuracy = 0.01)) +
    # scale_fill_gradientn(colors = color,limits = lim) + 
    # scale_fill_viridis_c(option = color, limits = lim) +
    theme_void() +
    theme(legend.title = element_blank())+
    theme(legend.position = "none")
  
  print(c(max(griddf1[,var]),max(griddf2[,var]),min(griddf1[,var]),min(griddf2[,var])))
  return(grid.arrange(p1, p2, ncol = 2))
  
}
+
  theme(legend.position = "none")


color = rainbow(7)
# MGWCR
p1 <- mg_se(1,"V1", c(0,0.2))
p2 <- mg_se(1,"V2", c(0,0.2))
p3 <- mg_se(1,"V3", c(0,0.2))
grid.arrange(p1, p2, p3, ncol = 1)
# GWCR
mg_se(3,"V1", c(0,0.2))
mg_se(3,"V2", c(0,0.3))
mg_se(3,"V3", c(0,0.2))

# Glb
mg_se(5,"V1", c(0,0.1))
mg_se(5,"V2", c(0,0.1))
mg_se(5,"V3", c(0,0.1))

## Simulation 2
# MGWCR
p1 <- mg_se(2,"V1", c(0,0.1))
p2 <- mg_se(2,"V2", c(0,0.1))
grid.arrange(p1, p2, ncol = 1)
# GWCR
mg_se(4,"V1", c(0,0.1))
mg_se(4,"V2", c(0,0.1))
# Glb
mg_se(6,"V1", c(0,0.1))
mg_se(6,"V2", c(0,0.1))


# GWCR_pre
griddf <- as.data.frame(cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                                    lngcoords = seq(1, 25, 1))), loc = 1:625))
values1 <- read.xlsx("C:/Users/User/Dropbox/MGWCR/Results/Simulation/Multi_normal/simulation_empirical_standard_error.xlsx", sheet = 1)
values2 <- read.xlsx("C:/Users/User/Dropbox/MGWCR/Results/Simulation/Multi_normal/simulation_average_standard_error.xlsx", sheet = 1)
griddf1 <- merge(griddf, values1, by = "loc")
griddf2 <- merge(griddf, values2, by = "loc")

# scale_fill_gradientn(colors = rainbow(5), limits = lim) + 
p1 <- ggplot(griddf1, aes(x = latcoords, y = lngcoords, fill = V1)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0,0.2)) +
  theme_minimal() 

p2 <- ggplot(griddf2, aes(x = latcoords, y = lngcoords, fill = V1)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0,0.2)) +
  theme_minimal() 

grid.arrange(p1, p2, ncol = 2)






