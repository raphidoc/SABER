#===================================================================================
# Code to create and Plot the spectral reflectance of different benthic end-members
#===================================================================================


#======================================================
# Generate WISE-Man Rb end-members
#======================================================

# Specify the main directory path
main_dir <- "C:/R/SABER/data/WISE_Rb/"

# Get the subdirectories
sub_dirs <- list.dirs(main_dir, recursive = FALSE)

# Create an empty list to store the results
sub_dir_list <- vector("list", length(sub_dirs))

# Iterate over the subdirectories
for (i in seq_along(sub_dirs)) {
  sub_dir <- sub_dirs[i]
  file_paths <- list.files(sub_dir, full.names = TRUE)
  sub_dir_list[[i]] <- list(sub_dir = sub_dir, file_paths = file_paths)
}

# Create an empty data frame
result_df <- read.table(file_paths[1], header = TRUE)[1]

# Iterate over the subdirectories and file paths
for (i in seq_along(sub_dir_list)) {
  sub_dir <- sub_dir_list[[i]]$sub_dir
  file_paths <- sub_dir_list[[i]]$file_paths
  
  print(paste0(sub_dir, " is read"))
  eelgrass <- read.table(file_paths[1], header = TRUE)[1]
  
  for(j in 1:length(file_paths))
  {
    station <- unlist(strsplit(file_paths[j], "/"))
    fname <- station[length(station)]
    
    eeltemp=read.table(file_paths[j], header = TRUE)
    eelgrass[[j+1]] <- eeltemp$rho
  }
  
  asdmean = colMeans(t(as.matrix(eelgrass[,2:ncol(eelgrass)]))) #calculate mean
  
  # Add the mean values as a column to the data frame
  col_name <- unlist(strsplit(sub_dir, "/"))  # Generate a column name
  col_name = col_name[length(col_name)]
  result_df[[col_name]] <- asdmean
}
#result_df$Zostera = result_df$Zostera/10
result_df_filter = data.frame("wavelength"=result_df$wavelength, "Mud"=result_df$Mud,
                              "sand" = colMeans(t(as.matrix(result_df$`New-Sand`, result_df$Sand))),
                              "Eelgrass" = result_df$Zostera)

# Remove all variables except `result_df`
rm(list = ls()[!ls() %in% c('result_df', 'result_df_filter')])

Rb_long = reshape2::melt(result_df_filter, id.vars="wavelength")
names(Rb_long) = c("wavelength","class","rb")

xmin = 400; xmax= 800; xstp=100
ymin= 0; ymax=0.5;ystp= signif(ymax/5, digits = 2)
asp_rat <- (xmax-xmin)/(ymax-ymin)

g1 <- ggplot(data = Rb_long, aes(x = wavelength, y = rb, color = class)) + 
  geom_line(show.legend = T, lwd=1.5)+
  scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")")) , limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0.10, 0.9),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        #legend.spacing.y = unit(2.0, 'cm'),
        plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g1

rb = result_df_filter
names(rb) = c("wavelength", "class1", "class2", "class3")
write.csv(rb, file = "./data/Rb_WISE-Man_endmembers.csv", col.names = T, row.names = F, 
          quote = FALSE)
save(rb, file = "./data/WISE-Man.RData")

#======================================================
# Generate Algae_WISE Rb-endmembers
#======================================================
#1.7 Load Algae-WISE bottom reflectance
year = 2023

if(year == 2022){
  
  algae_rb = read.csv("./data/Rb_AlgaeWISE.csv", header = T) #This data has been accidentally replaced with wrong one
  algae_rb$type = paste0(algae_rb$Target, algae_rb$Background)
  
  bottom_type = unique(algae_rb$Target)
  bottom_type = bottom_type[-(c(2,4))]
  
} else {
  
  algae_rb = read.csv("./data/bottom_reflectance_2023-08-25.csv", header = T)
  
  
  bottom_type = unique(algae_rb$target)
  bottom_type = bottom_type[-(c(2,4))]
  
}


xmin = 400; xmax= 800; xstp=100
ymin= 0; ymax=100;ystp= signif(ymax/5, digits = 2)
asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(data = algae_rb, aes(x = algae_rb$lambda, y = algae_rb$rho, color = algae_rb$target)) + 
  geom_line(show.legend = T, lwd=1.1)+
  scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("a"),{}["t-w"],"(",lambda,",", 0^"-",")[", m^-1,"]")) , limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0.10, 0.9),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        #legend.spacing.y = unit(2.0, 'cm'),
        plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

#Split Rb data
bottom_type_list = split(algae_rb, algae_rb$target)

# Create a list where each entry contains a data frame with wavelength and Rb
Rb_algaeWISE <- lapply(bottom_type_list, function(df_observation) {
  wavelength <- unique(df_observation$lambda)
  Rb_df <- matrix(df_observation$rho, nrow=512, byrow = F)
  Rb = rowMeans(Rb_df, na.rm = T)
  
  # Create a list for wavelength and Rb
  df_wavelength_Rrs <- data.frame(wavelength, Rb)
})

Rb_df = data.frame("wave" = Rb_algaeWISE[[1]]$wavelength)

for (i in seq_along(Rb_algaeWISE)) {
  Rb_df <- cbind(Rb_df, Rb_algaeWISE[[i]]$Rb)
}

# Set column names of the result dataframe as the names of each dataframe in the list
names(Rb_df)[-1] <- names(Rb_algaeWISE)
Rb_df[Rb_df < 0 ] <- NA
Rb_df[-1][Rb_df[-1] > 100 & Rb_df[-1] <= 150] <- 100
Rb_df[-1][Rb_df[-1] > 100 ] <- NA
Rb_df <- Rb_df[Rb_df$wave >400 & Rb_df$wave <=800,]
Rb_df[-1] <- Rb_df[-1]/100

Rb_df[-1] = apply(Rb_df[-1], 2,FUN = zoo::na.approx)

# Apply spline interpolation on each column except the first one
Rb_df_smooth <- as.data.frame(apply(Rb_df[-1], 2, 
                                    function(x) smooth.spline(Rb_df$wave, x)$y))

Rb_df_interp <- as.data.frame(apply(Rb_df_smooth, 2, 
                                    function(x) spline(Rb_df$wave, x, xout = wavelength)$y))
Rb_df_interp[Rb_df_interp > 1] = 1
Rb_df_interp$wavelength = wavelength
Rb_df_interp <- Rb_df_interp %>% dplyr::select(wavelength, everything())

Rb_set = Rb_df_interp
rm(algae_rb, bottom_type_list, bottom_type, Rb_df_interp)

Rb_long = reshape2::melt(Rb_set, id.vars="wavelength")
names(Rb_long) = c("wavelength","class","rb")

xmin = 400; xmax= 750; xstp=50
ymin= 0; ymax=1;ystp= signif(ymax/5, digits = 2)
asp_rat <- (xmax-xmin)/(ymax-ymin)

legend_position <- c(0.35, 0.98)

g1 <- ggplot(data = Rb_long[Rb_long$class != "Touffue",], 
             aes(x = wavelength, y = rb, color = class)) + 
  geom_line(show.legend = T, lwd=1.5)+
  scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
                     limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")")) , 
                     limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  scale_color_viridis(name = "", discrete = T
                      # , 
                      # labels = c(expression(paste(italic("Lithothamnium sp."))), 
                      #                                     "Rock", 
                      #            expression(paste(italic("Saccharina latissima"))) 
                      #)
)+
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 25, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
        legend.direction = "vertical",
        legend.title = element_text(colour = "black", size = 15, face = "plain"),
        legend.text = element_text(colour = "black", size = 25, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        legend.text.align = 0,
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.0,0.9,0.0,0.0), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g1
ggsave(paste0("./outputs/rb_algaewise_",year,".png"), plot = g1, scale = 1.7, 
       width = 4.5, height = 4.5, 
       units = "in",dpi = 300)

if(year = 2022) {
  
  rb = Rb_set[-(5)]
  names(rb) = c("wavelength", "class1", "class2", "class3")
  write.csv(rb, file = "./data/Rb_algaeWISE_endmembers.csv", col.names = T, row.names = F, 
            quote = FALSE)
  save(rb, file = "./data/algae-WISE.RData")
  
} else {
  
  
  benthic_classes <- colnames(Rb_set)[-1]
  
  combinations <- combn(benthic_classes, 3, simplify = FALSE)
  
  # Function to get initials of benthic class names
  get_initials <- function(name) {
    initials <- toupper(substr(name, 1, 1))
    paste(initials, collapse = "")
  }
  
  combination_Rbs <- lapply(combinations, function(cols) {
    #Rb_set %>% dplyr::select(wavelength, all_of(cols))
    rb <- Rb_set %>% dplyr::select(c("wavelength", cols))
    initials <- paste(cols[1], cols[2], cols[3], collpase = "_")#get_initials(cols)
    filename <- paste0("./data/Rb_algaewise_2023_", initials, ".csv")
    
    write.csv(rb, file = filename, row.names = FALSE, quote = F)
    save(rb, file = paste0("./data/Rb_algaewise_2023_", initials, ".RData"))
  })
  
  
  
  
  rb = Rb_set[-(5)]
  names(rb) = c("wavelength", "class1", "class2", "class3")
  write.csv(rb, file = "./data/Rb_algaeWISE_endmembers.csv", col.names = T, row.names = F, 
            quote = FALSE)
  save(rb, file = "./data/algae-WISE_2023.RData")
  
}

