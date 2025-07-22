

library(ggeffects)
library(tidyverse)
library(vegan)
library(nlme)
library(fixest)


#soil variables for soil variation
soil <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/IDE_soil_2024-12-16.csv")
summary(soil)#determine which soil variables have most complete data: ph, p, k, c,n

soil_variation <- soil%>%
  mutate(p = as.numeric(p))%>%
  group_by(site_code)%>%
  dplyr::summarise(ph_var = quantile(ph, 0.95, na.rm = TRUE)-quantile(ph, 0.05, na.rm = TRUE),
                  p_var = quantile(p, 0.95, na.rm = TRUE)-quantile(p, 0.05, na.rm = TRUE),
                  k_var = quantile(k, 0.95, na.rm = TRUE)-quantile(k, 0.05, na.rm = TRUE),
                  c_var = quantile(c, 0.95, na.rm = TRUE)-quantile(c, 0.05, na.rm = TRUE),
                  n_var = quantile(n, 0.95, na.rm = TRUE)-quantile(n, 0.05, na.rm = TRUE))
  


#cover data
cover_ppt.1 <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/cover_ppt_2024-12-19.csv")%>%
          subset(habitat.type == "Grassland")
info_df <- cover_ppt.1%>%
          dplyr::select(site_code, map, habitat.type)%>%
          unique()

#number of replicates, start by requiring at least 5 replicates
num_replicates <- cover_ppt.1%>%
                  dplyr::select(site_code, trt, block, plot, subplot)%>%
                  unique()%>%
                  group_by(site_code, trt)%>%
                  dplyr::summarise(replicates = n())%>%
                  pivot_wider(names_from = "trt", values_from = "replicates")%>%
                  subset(Control >= 5 & Drought >=5)

cover_ppt <- subset(cover_ppt.1, site_code %in% num_replicates$site_code)%>%
             subset(n_treat_years >= 1 & n_treat_years <= 4) #starting with just the first treatment year. To add more treatment years you need to tinker with the loop below so that it constrains calculations for each treatment year
            

relprecip_by_site_n_trt_year <- cover_ppt.1%>%
                                  dplyr::select(site_code, year, n_treat_years, trt, ppt.1, ppt.2, ppt.3, ppt.4, map)%>%
                                  unique()%>%
                                  mutate(relprecip.1 = (ppt.1-map)/map,
                                         relprecip.2 = (ppt.2-map)/map,
                                         relprecip.3 = (ppt.3-map)/map,
                                         relprecip.4 = (ppt.4-map)/map)


#calculate alpha diversity
alpha_richness_by_site <- cover_ppt.1%>%
  group_by(site_code, trt, n_treat_years, block, plot, subplot)%>%
  dplyr::summarize(richness = n())%>%
  group_by(site_code, trt, n_treat_years)%>%
  dplyr::summarise(mean_richness = mean(richness))



#calculate gamma diversity
gamma_by_site <- cover_ppt.1%>%
  subset(n_treat_years == 0)%>%
  dplyr::select(site_code, block, plot, subplot)%>%
  unique()%>%
  group_by(site_code)%>%
  dplyr::summarise(n_plots = n())%>%
  left_join(cover_ppt.1, by = "site_code")%>%
  subset(n_treat_years == 0)%>%
  dplyr::select(site_code, Taxon, n_plots)%>%
  unique()%>%
  group_by(site_code, n_plots)%>%
  dplyr::summarise(gamma_rich = n())%>%
  mutate(gamma_rich_relative = gamma_rich/n_plots) #you have to relatvize by number of plots for different sampling effort
  



#loop to calculate beta diversity by site, treatment, year

site_vector <- unique(cover_ppt$site_code)

distances_master <- {}

for(i in 1:length(site_vector)) {
  temp.df <- subset(cover_ppt, site_code == site_vector[i])
  
  temp.wide <- temp.df%>%
    dplyr::select(site_code, year, n_treat_years, trt, block, plot, subplot, Taxon, max_cover)%>%
    pivot_wider(names_from = Taxon, values_from = max_cover, values_fill = 0)
  temp.distances.bray <- vegdist(temp.wide[8:ncol(temp.wide)], method = "bray")#CHECK THE NUMBER OF THE STARTING COLUMN, SUBJECT TO CHANGE WITH DATA UPDATES
  temp.mod.bray <- betadisper(temp.distances.bray, group = temp.wide$trt, type = "centroid")
  
  temp.distances.jaccard <- vegdist(temp.wide[8:ncol(temp.wide)], method = "jaccard")#CHECK THE NUMBER OF THE STARTING COLUMN, SUBJECT TO CHANGE WITH DATA UPDATES
  temp.mod.jaccard <- betadisper(temp.distances.jaccard, group = temp.wide$trt, type = "centroid")
  
  distances_temp <- data.frame(site_code = site_vector[i], trt = temp.wide$trt, dist.bray = temp.mod.bray$dist, dist.jaccard = temp.mod.jaccard$dist, year = temp.wide$year, n_treat_years = temp.wide$n_treat_years )
  #distances_temp <- subset(distances_temp, dist > 0.00000000001) #not necessary when cO2 treatment excluded
  
  

  
  
  
  distances_master <- rbind(distances_master, distances_temp )
  rm(temp.df, temp.wide, temp.distances, temp.mod, distances_temp)
}

mean.dist.df <- distances_master%>%
                group_by(site_code, trt, year, n_treat_years)%>%
                dplyr::summarize(mean_dist.bray = mean(dist.bray), mean_dist.jaccard = mean(dist.jaccard))

dist.df <- left_join(mean.dist.df, info_df, by = "site_code")%>%
            left_join(soil_variation, by = "site_code")%>%
            left_join(relprecip_by_site_n_trt_year, by = c("site_code", "year", "n_treat_years", "trt", "map"))%>%
            mutate(multyear.relprecip = (relprecip.1+relprecip.2+relprecip.3+relprecip.4)/4)



##bray vs jaccard
#simple treatment
mod <- lme(mean_dist.bray ~ trt, random = list(n_treat_years=~1,site_code=~1), data = dist.df)
summary(mod)
mod <- feols(mean_dist.bray ~ trt|n_treat_years+site_code, cluster = ~site_code, data = dist.df)
summary(mod)
mod <- feols(mean_dist.bray ~ trt|n_treat_years, panel.id = ~site_code+n_treat_years, data = dist.df)
summary(mod)
mod <- feols(mean_dist.bray ~ trt, panel.id = c("site_code","n_treat_years"), data = dist.df)
summary(mod)

x <- ggpredict(mod, "trt")

ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))


mod <- feols(mean_dist.jaccard ~ trt|n_treat_years+site_code, cluster = ~site_code, data = dist.df)
summary(mod)

##relprecip
mod <- feols(mean_dist.bray ~ relprecip.1|n_treat_years+site_code, cluster = ~site_code, data = dist.df)
summary(mod)

mod <- feols(mean_dist.bray ~ relprecip.1, panel.id = c("site_code","n_treat_years"), data = dist.df)
summary(mod)
mod<- lme(mean_dist.bray ~ relprecip.1, random = list(site_code= ~1, n_treat_years=~1), data = dist.df)
summary(mod)

mod <- feols(mean_dist.jaccard ~ relprecip.1|n_treat_years+site_code, cluster = ~site_code, data = dist.df)
summary(mod)


##Over time
mod <- feols(mean_dist.bray ~ trt * n_treat_years|site_code, cluster = ~site_code, data = dist.df)
summary(mod, vcov = "cluster")
mod <- feols(mean_dist.bray ~ trt * n_treat_years, panel.id = ~site_code, data = dist.df)
summary(mod)
mod <- lme(mean_dist.bray ~ trt * n_treat_years, random =  list(site_code=~1), data = dist.df)
summary(mod)

mod <- feols(mean_dist.jaccard ~ trt * n_treat_years|site_code, cluster = ~site_code, data = dist.df)
summary(mod)

mod <- feols(mean_dist.bray ~ multyear.relprecip|site_code, cluster = ~site_code, data = dist.df)
summary(mod)

mod <- feols(mean_dist.jaccard ~ multyear.relprecip|site_code, cluster = ~site_code, data = dist.df)
summary(mod)



























##quick look

mod <- lme(mean_dist.bray~trt*map, random = ~1|site_code, data = dist.df)
summary(mod)
mod <- lme(mean_dist.jaccard~trt*map, random = ~1|site_code, data = dist.df)
summary(mod)


##alpha and gamma
a.g.df <- dist.df%>%
          left_join(alpha_richness_by_site, by = c("site_code", "trt", "n_treat_years"))%>%
          left_join(gamma_by_site, by = "site_code")

mod <- feols(mean_dist ~ mean_richness + relprecip.1 * gamma_rich_relative #+ ph_var + p_var + k_var + c_var + n_var 
             | n_treat_years, cluster = ~site_code, data = a.g.df) #need to add more information to run this. include interaction of site code and n_treat_years
summary(mod)



#minimalist
mod <- lm(mean_dist ~ relprecip.1 + map + ph_var + p_var + k_var + c_var + n_var+ site_code, data = dist.df) #need just a little more data for this
summary(mod)

mod <- feols(mean_dist ~ relprecip.1 * map 
             #+ ph_var + p_var + k_var + c_var + n_var + 
               ,cluster = ~site_code
                , data = dist.df)
summary(mod)




library(PerformanceAnalytics)
mat <- dist.df[,c("map", "ph_var", "p_var", "k_var", "c_var", "n_var", "relprecip.1")]
chart.Correlation(mat, histogram = TRUE)

mod <- lm(mean_dist ~ multyear.relprecip * map + ph_var + p_var + k_var + c_var + n_var + site_code, data = dist.df) #need just a little more data for this
summary(mod)

mod <- feols(mean_dist ~ multyear.relprecip * map 
             #+ ph_var + p_var + k_var + c_var + n_var 
             ,cluster = ~site_code, data = dist.df)
summary(mod)

mod <- feols(mean_dist ~ multyear.relprecip |  n_treat_years, cluster = ~site_code, data = dist.df)
summary(mod)

mod <- feols(mean_dist ~ multyear.relprecip * n_treat_years | site_code , data = dist.df)
summary(mod)

mod <- feols(mean_dist ~ relprecip.1 * n_treat_years ,cluster = ~site_code , data = dist.df)
summary(mod)

#nestedness vs turnover
#mod <- feols(`Beta diversity` ~ Nestedness + Turnover + `Soil variation vars`, data = ) #this one would require a different calculation of beta diversity that incorporates nestedness and turnover