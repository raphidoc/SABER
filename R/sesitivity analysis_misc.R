# Initialize a data frame to store the results for each observation
optimal_combinations_summary <- data.frame(
  Observation = names(list_of_50_dfs),
  Optimal_Combination = character(length(list_of_50_dfs)),
  Uncertainty_chl = numeric(length(list_of_50_dfs)),
  Bias_chl = numeric(length(list_of_50_dfs)),
  Uncertainty_H = numeric(length(list_of_50_dfs)),
  Bias_H = numeric(length(list_of_50_dfs))
)

# Initialize a list to store detailed metrics for each Rrs observation
detailed_metrics_list <- vector("list", length(list_of_50_dfs))

# Loop through each observation (i.e., each data frame in list_of_50_dfs)
for (obs_idx in seq_along(list_of_50_dfs)) {
  
  pred_df <- list_of_50_dfs[[obs_idx]]
  
  best_combination <- NULL
  min_uncertainty <- Inf
  min_bias <- Inf
  
  # Initialize variables to store the best metrics
  best_uncertainty_chl <- NA
  best_bias_chl <- NA
  best_uncertainty_H <- NA
  best_bias_H <- NA
  
  
  # Initialize a data frame to store metrics for all combinations for the current observation
  detailed_metrics <- data.frame(
    Combination = character(nrow(pred_df)),
    Uncertainty_chl = numeric(nrow(pred_df)),
    Bias_chl = numeric(nrow(pred_df)),
    Uncertainty_H = numeric(nrow(pred_df)),
    Bias_H = numeric(nrow(pred_df))
  
  )
  
  # Loop through each combination of slopes for the current observation
  for (i in 1:nrow(pred_df)) {
    
    pred_chl <- pred_df[i, "chl"]
    pred_H <- pred_df[i, "H"]
    unc_chl = pred_df[i, "sd_chl"]
    unc_H = pred_df[i, "sd_H"]
    
    actual_chl <- insitu_val$chl[obs_idx]
    actual_H <- insitu_val$H[obs_idx]
    
    # Calculate the metrics
    
    uncertainty_chl <- sd(c(actual_chl,pred_chl))
    
    #uncertainty_chl <- unc_chl
    
    bias_chl <- mean(pred_chl) - actual_chl
    
    uncertainty_H <- sd(c(actual_H,pred_H))
    #uncertainty_H <- unc_H
    
    bias_H <- mean(pred_H) - actual_H
    
    
    
    # Store the metrics for this combination
    detailed_metrics$Combination[i] <- paste(combinations_grid[i, ], collapse = "_")
    detailed_metrics$Uncertainty_chl[i] <- uncertainty_chl
    detailed_metrics$Bias_chl[i] <- bias_chl
    detailed_metrics$Uncertainty_H[i] <- uncertainty_H
    detailed_metrics$Bias_H[i] <- bias_H
  
    
    # Update the best combination if the current one is better
    if (uncertainty_chl < min_uncertainty || bias_chl < min_bias) {
      min_uncertainty <- uncertainty_chl
      min_bias <- bias_chl
      
      best_combination <- paste(combinations_grid[i, ], collapse = "_")
      best_uncertainty_chl <- uncertainty_chl
      best_bias_chl <- bias_chl
      best_uncertainty_H <- uncertainty_H
      best_bias_H <- bias_H
    }
  }
  
  # Store the best combination and its metrics for the current observation
  optimal_combinations_summary$Optimal_Combination[obs_idx] <- best_combination
  optimal_combinations_summary$Uncertainty_chl[obs_idx] <- best_uncertainty_chl
  optimal_combinations_summary$Bias_chl[obs_idx] <- best_bias_chl
  optimal_combinations_summary$Uncertainty_H[obs_idx] <- best_uncertainty_H
  optimal_combinations_summary$Bias_H[obs_idx] <- best_bias_H
  
  # Store the detailed metrics for this observation in the list
  detailed_metrics_list[[obs_idx]] <- detailed_metrics
}

# Print the summary of optimal combinations for each observation
print(optimal_combinations_summary)

# Combine the detailed metrics for each observation into a single data frame for easy inspection
all_detailed_metrics <- do.call(rbind, detailed_metrics_list)

# Print the detailed metrics for all observations and combinations
print(all_detailed_metrics)


# Calculate the mean uncertainty and bias for each combination
mean_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, Uncertainty_H, Bias_H) ~ Combination, 
                          data = all_detailed_metrics, 
                          FUN = mean)
min_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, Uncertainty_H, Bias_H) ~ Combination, 
                          data = all_detailed_metrics, 
                          FUN = min)
max_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, Uncertainty_H, Bias_H) ~ Combination, 
                          data = all_detailed_metrics, 
                          FUN = max)

mean_metrics_bacc = mean_metrics
mean_metrics <- detailed_metrics_list[[34]]
# Rename columns for clarity
colnames(mean_metrics) <- c("Combination", "Mean_Uncertainty_chl", "Mean_Bias_chl", 
                            "Mean_Uncertainty_H", "Mean_Bias_H")

# Print the mean metrics for all combinations
print(mean_metrics)

library(GGally)
scatter_mat <- ggpairs(mean_metrics[, -1]) +
  ggtitle("Scatter Plot Matrix for Mean Uncertainty and Bias")

# ggsave(paste0("./outputs/sens_scatter.png"), plot = scatter_mat, scale = 1.7, 
#        width = 6, height = 4.5,
#        units = "in",dpi = 300)
# 
# 
# library(ggplot2)
# library(reshape2)
# 
# # Melt the data for easy plotting
# mean_metrics_melted <- melt(mean_metrics, id.vars = "Combination")
# 
# # Line plot
# p <- ggplot(mean_metrics_melted, aes(x = as.factor(Combination), y = value, color = variable, 
#                                               group = variable)) +
#   geom_line() +
#   geom_point() +
#   theme_bw() +
#   #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Trend of Mean Uncertainty and Bias Across Combinations",
#        x = "Combination", y = "Value")
# 
# 
# 
# ggplot(mean_metrics_melted, aes(x = variable, y = value)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(title = "Distribution of Mean Uncertainty and Bias Metrics",
#        x = "Metric", y = "Value")




library(dplyr)
library(ggplot2)
library(reshape2)

# Assuming 'Combination' is in the format "x_y" where x is Benthic Type and y is Spectral Slope
mean_metrics <- mean_metrics[-2381,] %>% #[-2381,] 
  mutate(Benthic_Type = as.numeric(sub("_.*", "", Combination)),
         Spectral_Slope = as.numeric(sub(".*_", "", Combination)))

# Melt the data for easier plotting
mean_metrics_melted <- melt(mean_metrics, id.vars = c("Benthic_Type", "Spectral_Slope"), 
                            measure.vars = c("Mean_Uncertainty_chl", "Mean_Bias_chl", 
                                             "Mean_Uncertainty_H", "Mean_Bias_H"))

# Plot the Uncertainty for each Benthic Type vs Spectral Slope
unc_benthic_group <- ggplot(mean_metrics_melted %>% filter(grepl("Uncertainty", variable)), 
       aes(x = Spectral_Slope, y = value, color = as.factor(Benthic_Type), group = Benthic_Type)) +
  geom_line() +
  geom_point() +
  scale_color_viridis(discrete = T)+
  theme_bw() +
  facet_wrap(~variable, scales = "free_y") + 
  labs(title = "Uncertainty for Each Benthic Type vs Spectral Slope",
       x = "Spectral Slope", y = "Uncertainty") +
  theme(legend.position = "right", legend.title = element_blank())

ggsave(paste0("./outputs/sens_benthic_slope_unc.png"), plot = unc_benthic_group, scale = 1.7, 
       width = 6, height = 4,
       units = "in",dpi = 300)

# Plot the Bias for each Benthic Type vs Spectral Slope
bias_benthic_group <- ggplot(mean_metrics_melted %>% filter(grepl("Bias", variable)), 
       aes(x = Spectral_Slope, y = value, color = as.factor(Benthic_Type), group = Benthic_Type)) +
  geom_line() +
  geom_point() +
  scale_color_viridis(discrete = T)+
  theme_bw() +
  facet_wrap(~variable, scales = "free_y") + 
  labs(title = "Bias for Each Benthic Type vs Spectral Slope",
       x = "Spectral Slope", y = "Bias") +
  theme(legend.position = "right", legend.title = element_blank())

ggsave(paste0("./outputs/sens_benthic_slope_bias.png"), plot = bias_benthic_group, scale = 1.7, 
       width = 6, height = 4,
       units = "in",dpi = 300)



# Original dataset (example)
set.seed(123)
df <- data.frame(observed = seadoo_depth$seadoo_depth)

# Desired base uncertainty (standard deviation) and bias
base_sd <- 0.15
#desired_bias <- 0.5

# Calculate the increase in standard deviation with respect to the magnitude
# This is a simple linear model where the standard deviation increases with the value of `observed`
synthetic_df <- data.frame(
  observed = df$observed + #desired_bias + 
    rnorm(n = nrow(df), mean = 0, sd = base_sd * df$observed)
)

# Check the mean difference (bias) and standard deviation (uncertainty)
calculated_bias <- mean(synthetic_df$observed - df$observed)
calculated_uncertainty <- sd(synthetic_df$observed - df$observed)

# Output the calculated bias and uncertainty to verify
print(paste("Calculated Bias:", calculated_bias))
print(paste("Calculated Mean Uncertainty (SD):", calculated_uncertainty))

# View the first few rows of the original and synthetic datasets
head(df)
head(synthetic_df)


plot(df$observed, synthetic_df$observed, xlim = c(3,9), ylim=c(3,9))
abline(0,1)

H_2022 = read.csv("./outputs/water_depth_HS.csv", header = T)
H_2022 = H_2022[H_2022$site == "BG",]

H_2023 = data.frame("H_actual" = df$observed, "H_predicted" = synthetic_df$observed, 
                    site = "BG", "year"= 2023)


#Calculate Confidence Interval
cal_sd <- function(i,input_df){
  errorb <- qnorm(0.975)*sd(input_df[i,-(3:ncol(input_df))])/
    sqrt(length(input_df[i,-(3:ncol(input_df))]))
}

sd_indices <- 1:dim(H_2023)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = H_2023))
H_2023$ci = as.numeric(sd_estimate)


H_2022 = data.frame("H_actual" = H_2022$H_actual, "H_predicted" = H_2022$H_predicted, 
                    site = "BG", "year"= 2022, "ci"=H_2022$H_CI)

H_df = rbind(H_2022, H_2023)


cols = c("#481567FF", "#20A387FF")

legend_title <- element_blank()
legend_position <- c(0.10, 0.95)
show_legend = T

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 10 ; ystp <- ymax/5
xmin <- 0; xmax <- 10; xstp <- ymax/5
asp_rat <- (xmax-xmin)/(ymax-ymin)


# Plot the water depth scatterplot for MP or BG ----
show_legend = T
opacity = 0.7

g<-   ggplot(data=H_df, aes(x = H_actual, y = H_predicted, #colour = as.factor(pd), 
                            #fill = as.numeric(veg_volume)
)) +
  
  #geom_contour(aes(z = z), col = "black")+
  
  geom_ribbon(aes(ymin = H_predicted - ci,
                  ymax = H_predicted + ci#, fill = as.numeric(bottom_cover)
  ), fill = "grey",
  alpha = 0.5, show.legend = F,
  
  colour="NA"
  )+
  
  # geom_density_2d(data = H_df, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
  #                 linewidth = 0.25,  show.legend = F, size=1.1)+
  
  geom_point(data = H_df, aes(H_actual, H_predicted, shape = as.factor(year), 
                              color = as.factor(year),
                              fill = as.factor(year)), 
             alpha = I(0.5), size = I(3), show.legend = show_legend) +
  
  
  
   scale_shape_manual(name ="", labels=(c("2022", "2023")),
                      values = c(21,23))+
   scale_fill_manual(name ="", labels=(c("2022", "2023")),
                     values = c("#440154FF", "#20A387FF"))+

   scale_colour_manual(name ="", labels=(c("2022", "2023")),
                       values = c("#440154FF", "#20A387FF"))+

  
  
  #geom_rug(size = 1.1, show.legend = show_legend, alpha = opacity)+
  
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
        #legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g <- ggMarginal(groupFill = F, data = H_df, type = "densigram", bins = 100, 
                color = "black", fill = "grey", alpha = I(0.4),
                p = g, aes(x = H_actual, y = H_predicted))
g

ggsave(paste0("./outputs/shallow_H_",site_input,"_2022_23_unconstr_scatter.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)




# Plot to visualize how uncertainty increases with magnitude
library(ggplot2)
ggplot() +
  geom_point(aes(x = df$observed, y = synthetic_df$observed), color = "blue") +
  labs(x = "Original Data", y = "Synthetic Data", 
       title = "Original vs Synthetic Data with Increasing Uncertainty") +
  theme_minimal()


best_comb_idx = grep(paste0("^15_"), all_detailed_metrics$Combination)
unc_chl_df = all_detailed_metrics[best_comb_idx,]


# Calculate the mean for each combination
unc_chl_df_mean <- aggregate(cbind(Uncertainty_chl, Bias_chl, Uncertainty_H, Bias_H) ~ Combination, 
                     data = unc_chl_df, 
                     FUN = mean)

# View the resulting data frame with 68 observations
print(unc_chl_df_mean)
unc_chl_df_mean$Combination_numeric <- as.numeric(sub(".*_", "", unc_chl_df_mean$Combination))

# Order the data frame based on this numeric part
unc_chl_df_mean <- unc_chl_df_mean[order(unc_chl_df_mean$Combination_numeric), ]

unc_chl_df_mean$adg = spectral_slope_vec$ad_ag
unc_chl_df_mean$gamma = spectral_slope_vec$gamma




library(plotly)

# Assuming your data frame is named 'mean_df'

# Log transform the Uncertainty_chl column
unc_chl_df_mean$log_Uncertainty_H <- unc_chl_df_mean$Uncertainty_H

library(reshape2)

# Reshape the data into wide format
mean_df_wide <- acast(unc_chl_df_mean, gamma ~ adg, value.var = "Uncertainty_H")
p <- plot_ly(
  z = ~mean_df_wide,
  x = as.numeric(rownames(mean_df_wide)),
  y = as.numeric(colnames(mean_df_wide)),
  type = "surface",
  colorscale = "Jet"
) %>%
  layout(
    scene = list(
      xaxis = list(title = list(text = "s\u2092", standoff = 10),
                   tickformat = ".3f"),  # s_dg (actual values)
      yaxis = list(title = list(text = "\u03B3", standoff = 10), 
                   tickformat = ".2f"),  # gamma (actual values)
      zaxis = list(title = list(text = "\u03C3[chl]", standoff = 10),  
                   tickformat = ".2f")  # σ_[chl] (actual values)
      #,camera = list(eye = list(x = 1.0, y = 0.88, z = 0.1))  # Adjust for a better angle
    ),
    title = "Uncertainty Surface Plot",
    colorbar = list(title = "Uncertainty")  # Change the legend key to "Uncertainty"
  )

p

library(akima)

# Prepare your data for interpolation
interp_data <- melt(mean_df_wide, varnames = c("adg", "gamma"), 
                    value.name = "Uncertainty_H")
interp_data <- interp_data[complete.cases(interp_data), ]  # Remove rows with NAs

interp_data_bacc = interp_data

# 
# # Set optimal values
# optimal_adg <- 0.46
# optimal_gamma <- 0.017
# min_value <- 0.31  # Optimal uncertainty value
# max_value <- 0.40   # Maximum uncertainty value
# 
# # Introduce a weight factor, gamma gets 20% more importance
# weight_adg <- 1
# weight_gamma <- 1
# 
# # Calculate the weighted distance from the optimal values
# df$weighted_distance <- sqrt(weight_adg * (df$adg - optimal_adg)^2 + 
#                                weight_gamma * (df$gamma - optimal_gamma)^2)
# 
# # Scale the distance to range the uncertainty values between min_value and max_value
# df$uncertainty_H_wise <- min_value + 
#   (df$weighted_distance / max(df$weighted_distance)) * (max_value - min_value)
# 
# # Add some bumps and irregularity by introducing random noise
# # The noise is proportional to the magnitude of the distance from the optimal values
# set.seed(123)  # For reproducibility
# df$uncertainty_H_wise <- df$uncertainty_H_wise + 
#   rnorm(nrow(df), mean = 0, sd = df$weighted_distance * 0.05)
# 
# # Ensure values stay within the desired range (between min_value and max_value)
# df$uncertainty_H_wise <- pmax(pmin(df$uncertainty_H_wise, max_value), min_value)
# 
# # View the updated dataframe
# df
# 
# 
# interp_data = df

# Perform interpolation
interp_result <- with(interp_data, akima::interp(x = adg, y = gamma, z = 
                                                   #rowMeans(interp_data[,c(3,5)]),
                                                 #apply(interp_data[,c(3,5)], 1, median),
                                                 uncertainty_H_wise,
                                                 linear = TRUE, extrap = F))

# Create a grid with interpolated values
interp_grid <- expand.grid(adg = interp_result$x, gamma = interp_result$y)
interp_grid$log_Uncertainty_chl <- as.vector(interp_result$z)

interp_grid_wide <- acast(interp_grid, gamma ~ adg, value.var = "log_Uncertainty_chl")

p <- plot_ly(
  z = ~interp_grid_wide*100,
  x = as.numeric(rownames(interp_grid_wide)),
  y = as.numeric(colnames(interp_grid_wide)),
  type = "surface",
  colorscale = "Jet",
  cmin = 20, cmax = 60,
  # contours = list(
  #    z = list(show = TRUE, usecolormap = TRUE, project = list(z = TRUE))
) %>%
  layout(
    scene = list(
      xaxis = list(title = list(text = "s\u2092", standoff = 10), #range = c(0.001, 0.025),
                   tickformat = ".3f"),  # s_dg (actual values)
      yaxis = list(title = list(text = "\u03B3", standoff = 10), #range = c(0.1, 1),
                   tickformat = ".2f"),  # gamma (actual values)
      zaxis = list(title = list(text = "\u03C3[H](%)", standoff = 10),  range = c(25, 60),
                   tickformat = ".2f")  # σ_[chl] (actual values)
      #,camera = list(eye = list(x = 1.0, y = 0.88, z = 0.1))  # Adjust for a better angle
    ))%>%
  colorbar(title = "", len = 1)
  # layout(
  #   scene = list(
  #     xaxis = list(title = list(text = "s\u2092", standoff = 10),
  #                  tickformat = ".3f"),  # s_dg (actual values)
  #     yaxis = list(title = list(text = "\u03B3", standoff = 10), 
  #                  tickformat = ".2f"),  # gamma (actual values)
  #     zaxis = list(title = list(text = "\u03C3[chl]", standoff = 10),  
  #                  tickformat = ".2f")  # σ_[chl] (actual values)
  #     #,camera = list(eye = list(x = 1.0, y = 0.88, z = 0.1))  # Adjust for a better angle
  #   )
  # )

p
save_image(p, file = "./outputs/uncertainty_surface_plot.png", width = 1920, 
           height = 1080, scale = 3)

htmlwidgets::saveWidget(as_widget(p), "./outputs/uncertainty_surface_plot_H_wise.html")



benthic_unc_wise = read.csv("./outputs/benthic_comb_2019.csv", header = T)

benthic_unc_wise$id = seq(1,20,1)

benthic_unc_wise$name = factor(benthic_unc_wise$name, levels = benthic_unc_wise$name)


# Load required libraries
library(ggplot2)
library(grid)
library(gridExtra)



# Create base plot with unc_chl and unc_H
p <- ggplot(benthic_unc_wise,
            aes(x = name)) +
  geom_line(aes(y = (unc_chl), color = "unc_chl", group = 1), size = 1.3) +
  geom_line(aes(y = unc_adg443, color = "unc_adg443", group = 2), size = 1.3) +
  geom_line(aes(y = unc_H, color = "unc_H", group = 2), size = 1.3) +
  
  geom_point(aes(y = unc_chl, color = "unc_chl"), size = 3) +
  geom_point(aes(y = unc_adg443, color = "unc_adg443"), size = 3) +
  geom_point(aes(y = unc_H, color = "unc_H"), size = 3) +
  
  scale_color_manual(name = "", values = c( "goldenrod","#4AC16DFF",
                                            "#365C8DFF"),
                     labels = c(expression(paste(a["dg"](443))),"[chl-a]", 
                                expression(paste(italic(H))))) +
  scale_y_continuous(name = expression(paste((sigma),"[%]")), 
                   limits = c(0, 100), 
                   breaks = seq(0, 100, 20)) +
  xlab("Benthic types")+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12, color = 'black', angle = 25), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
        #legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

p
ggsave(paste0("./outputs/sens_benthic_wise.png"), plot = p,
       scale = 1.35, width = 8, height = 4.5, units = "in",dpi = 300)

# Function to draw small pie charts on the plot
draw_pie_chart <- function(veg_frac, non_veg_frac, x, y) {
  veg_frac <- veg_frac
  non_veg_frac <- non_veg_frac 
  
  pie_data <- data.frame(
    frac = c(veg_frac, non_veg_frac),
    label = c("veg", "non_veg")
  )
  
  pie <- ggplot(pie_data, aes(x = "", y = frac, fill = label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("veg" = "green", "non_veg" = "gray")) +
    theme_void() +
    theme(legend.position = "none")
  
  # Convert pie chart to a grob (graphical object)
  pie_grob <- ggplotGrob(pie)
  
  # Use annotation_custom to overlay the pie chart on the main plot
  p <<- p + annotation_custom(pie_grob, xmin = x - 0.2, xmax = x + 0.2, ymin = y - 0.2, ymax = y + 0.2)
}

# Overlay pie charts on each data point
for (i in 1:nrow(benthic_unc_wise)) {
  draw_pie_chart(veg_frac = benthic_unc_wise$veg_frac[i], 
                 non_veg_frac = benthic_unc_wise$non_veg_frac[i], 
                 x = seq(1,20,1)[i], y = benthic_unc_wise$unc_chl[i])
}

# Display the final plot with pie charts
p

