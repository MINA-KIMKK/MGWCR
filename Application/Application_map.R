###################################################################
## MGWCR
## Application 
###################################################################

library(ggplot2)
library(sf)
library(tigris)
library(dplyr)
library(openxlsx)


##################################
plot_state_counties <- function(variable, state,
                                title = "", 
                                subtitle = "", 
                                palette,
                                option,
                                filepath,
                                sheet,
                                limits=NULL) { 
  
  florida_counties <- counties(state = state, cb = TRUE, class = "sf")
  
  coef.app <- read.xlsx(filepath, sheet = sheet)
  coef.app <- coef.app %>% mutate(GEOID = as.factor(GEOID))
  
  florida_map <- florida_counties %>%
    left_join(coef.app, by = "GEOID")
  
  # automatically setting limitations
  if (is.null(limits)){
    variable_min <- floor(min(florida_map[[variable]], na.rm = TRUE) * 10) / 10
    variable_max <- ceiling(max(florida_map[[variable]], na.rm = TRUE) * 10) / 10
    limits <- c(variable_min, variable_max)
  } 
  breaks <- seq(limits[1], limits[2], length.out = 4)
  
  # Depending on the option, set the fill differently
  if (option == 3) {
    breaks <- seq(limits[1], limits[2], length.out = 3)
    p <- ggplot(data = florida_map) +
      geom_sf(aes_string(fill = paste0("ifelse(", variable, " >0.05, NA, ", variable, ")")), color = "white") +
      scale_fill_gradient(low = "#79AF97FF", high = "grey90", limits = limits, labels = scales::number_format(accuracy = 0.01), breaks = breaks, na.value = "grey90")
  } else {
    breaks <- seq(limits[1], limits[2], length.out = 4)
    p <- ggplot(data = florida_map) +
      geom_sf(aes_string(fill = variable), color = "white")
    
    if (option == 1) {
      p <- p + scale_fill_gradient2(low = "#003c67ff", mid = "snow", high = "#b24745ff", midpoint = 0, limits = limits, labels = scales::number_format(accuracy = 0.1), breaks = breaks)
    } else if (option == 2) {
      p <- p + scale_fill_viridis_c(option = palette, limits = limits, labels = scales::number_format(accuracy = 0.1), breaks = breaks)
    }
  }
  
  p <- p + labs(title = title, subtitle = subtitle) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.position = c(0.16, 0.15), 
          legend.direction = "horizontal",
          legend.text = element_text(size = 13))
  
  return(p)
}



########################
## Setting 
state="NY"
sheet=1

########################
## Coefficients 
filepath1="~/coefficients.xlsx"
plot_state_counties("agemo",state=state,palette="magma",filepath=filepath1,sheet=sheet,option=1)
plot_state_counties("marriage",state=state,palette="magma",filepath=filepath1,sheet=sheet,option=1)
plot_state_counties("income",state=state,palette="magma",filepath=filepath1,sheet=sheet,option=1)
plot_state_counties("pm25",state=state,palette="magma",filepath=filepath1,sheet=sheet,option=1)
plot_state_counties("NDVI270",state=state,palette="magma",filepath=filepath1,sheet=sheet,option=1)

## p-value
filepath2="~/pvalue.xlsx"
plot_state_counties("agemo",state=state,palette=color,filepath=filepath2,sheet=sheet,option=3,limits=c(0,0.1))
plot_state_counties("marriage",state=state,palette=color,filepath=filepath2,sheet=sheet,option=3,limits=c(0,0.1))
plot_state_counties("income",state=state,palette=color,filepath=filepath2,sheet=sheet,option=3,limits=c(0,0.1))
plot_state_counties("pm25",state=state,palette=color,filepath=filepath2,sheet=sheet,option=3,limits=c(0,0.1))
plot_state_counties("NDVI270",state=state,palette=color,filepath=filepath2,sheet=sheet,option=3,limits=c(0,0.1))




########################
## Distribution 
plot_ny_dist <- function(variable){
  
  filepath="~/distribution.xlsx"
  ny_counties <- counties(state = "NY", cb = TRUE, class = "sf")
  
  coef.app <- read.xlsx(filepath, sheet = 1)
  coef.app <- coef.app %>% mutate(GEOID = as.factor(GEOID))
  
  ny_map <- ny_counties %>%
    left_join(coef.app, by = "GEOID")
  
  variable_min <- floor(min(ny_map[[variable]], na.rm = TRUE) * 10) / 10
  variable_max <- ceiling(max(ny_map[[variable]], na.rm = TRUE) * 10) / 10
  limits <- c(variable_min, variable_max)
  breaks <- seq(limits[1], limits[2], length.out = 4)
  
  breaks <- seq(limits[1], limits[2], length.out = 3)
  p <- ggplot(data = ny_map) +
    geom_sf(aes_string(fill = variable)) +
    scale_fill_gradient(low = "white", high = "#00a1d5ff", limits = limits, labels = scales::number_format(accuracy = 0.01), breaks = breaks, na.value = "grey90")
  
  p <- p +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.position = c(0.16, 0.15), 
          legend.direction = "horizontal",
          legend.text = element_text(size = 13))
  
  print(p)
  
}
plot_ny_dist("agemo")
plot_ny_dist("marriage")
plot_ny_dist("income2")
plot_ny_dist("pm25")
plot_ny_dist("NDVI270")




