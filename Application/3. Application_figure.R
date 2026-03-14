###################################################################
## MGWCR
## Figure6
## coefficients&pvalue map
###################################################################

library(ggplot2)
library(sf)
library(tigris)
library(dplyr)


##################################
plot_state_counties <- function(data, 
                                variable, 
                                state,
                                title = "", 
                                subtitle = "", 
                                palette,
                                option,
                                limits=NULL) { 
  
  florida_counties <- counties(state = state, cb = TRUE, class = "sf")
  
  coef.app <- data %>% mutate(GEOID = as.factor(loc))
  
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



####################################
######## Figure 6  #################
####################################
color = "viridis" # Change freely
betas_app <- MGWCR_betas_app # GWCR_betas_app, MGWCR_betas_app
pvals_app <- MGWCR_pvals_app # GWCR_pvals_app, MGWCR_pvals_app

## Coefficients 
plot_state_counties(data = betas_app, variable = "agemo", state = "NY", palette = "magma", option = 1)
plot_state_counties(data = betas_app, variable = "marriage", state = "NY", palette = "magma", option = 1)
plot_state_counties(data = betas_app, variable = "income", state = "NY", palette = "magma", option = 1)
plot_state_counties(data = betas_app, variable = "pm25", state = "NY", palette = "magma", option = 1)
plot_state_counties(data = betas_app, variable = "NDVI270", state = "NY", palette = "magma", option = 1)

## p-value
plot_state_counties(data = GWCR_pvals_app, variable = "agemo", state = "NY", palette = color, option = 3, limits = c(0, 0.1))
plot_state_counties(data = GWCR_pvals_app, variable = "marriage", state = "NY", palette = color, option = 3, limits = c(0, 0.1))
plot_state_counties(data = GWCR_pvals_app, variable = "income", state = "NY", palette = color, option = 3, limits = c(0, 0.1))
plot_state_counties(data = GWCR_pvals_app, variable = "pm25", state = "NY", palette = color, option = 3, limits = c(0, 0.1))
plot_state_counties(data = GWCR_pvals_app, variable = "NDVI270", state = "NY", palette = color, option = 3, limits = c(0, 0.1))

