


library(tidyverse)
library(vegan)
library(nlme)


cover_ppt <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/cover_ppt_2024-12-19.csv")
info_df <- cover_ppt%>%
          dplyr::select(site_code, map, habitat.type)%>%
          unique()

#number of replicates, start by requiring at least 5 replicates
num_replicates <- cover_ppt%>%
                  dplyr::select(site_code, trt, block, plot, subplot)%>%
                  unique()%>%
                  group_by(site_code, trt)%>%
                  dplyr::summarise(replicates = n())%>%
                  pivot_wider(names_from = "trt", values_from = "replicates")%>%
                  subset(Control >= 5 & Drought >=5)

cover_ppt <- subset(cover_ppt, site_code %in% num_replicates$site_code)%>%
            subset(n_treat_years ==1) #starting with just the third treatment year. To add more treatment years you need to tinker with the loop below so that it contrains calculations for each treatment year
            




site_vector <- unique(cover_ppt$site_code)

distances_master <- {}

for(i in 1:length(site_vector)) {
  temp.df <- subset(cover_ppt, site_code == site_vector[i])
  
  temp.wide <- temp.df%>%
    pivot_wider(names_from = Taxon, values_from = max_cover, values_fill = 0)
  temp.distances <- vegdist(temp.wide[30:ncol(temp.wide)], method = "bray")#CHECK THE NUMBER OF THE STARTING COLUMN, SUBJECT TO CHANGE WITH DATA UPDATES
  temp.mod <- betadisper(temp.distances, group = temp.wide$trt, type = "centroid")
  distances_temp <- data.frame(site_code = site_vector[i], trt = temp.wide$trt, dist = temp.mod$dist)
  #distances_temp <- subset(distances_temp, dist > 0.00000000001) #not necessary when cO2 treatment excluded
  #distances_temp$dist <- ifelse(distances_temp$dist > 0.00000000001, distances_temp$dist, 0.001) #changes value for single serc experiment where distance equals essentially 0 which doesn't work with response ratios
  distances_master <- rbind(distances_master, distances_temp )
  rm(temp.df, temp.wide, temp.distances, temp.mod, distances_temp)
}

mean.dist.df <- distances_master%>%
                group_by(site_code, trt)%>%
                dplyr::summarize(mean_dist = mean(dist))



##quick look
dist.df <- left_join(mean.dist.df, info_df, by = "site_code")

mod <- lme(mean_dist~trt*map, random = ~1|site_code, data = dist.df)
summary(mod)

