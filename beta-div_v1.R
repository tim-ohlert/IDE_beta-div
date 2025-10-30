

library(ggeffects)
library(tidyverse)
library(vegan)
library(nlme)
library(fixest)
library(sandwich)
library(RcmdrMisc)
library(lmtest)
library(ggthemes)


#site info
Site_Elev.Disturb <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/Site_Elev-Disturb.csv")

#climate info
climate <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/climate/climate_mean_annual_by_site_v3.csv")


prop <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/community_comp/Prc_LifeHistory_Controls_Oct2023.csv")

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
cover_ppt.1 <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/cover_ppt_2025-10-20.csv")%>%
          #subset(habitat.type == "Grassland")%>%
          subset(habitat.type != "Forest")%>%
          mutate(block = ifelse(site_code == "allmendo.ch", "2", block))%>% #grouping allmend sites
          mutate(site_code = ifelse(site_code == "allmendo.ch", "allmendb.ch", site_code))

info_df <- cover_ppt.1%>%
          dplyr::select(site_code, map, habitat.type)%>%
          unique()

nreps <- cover_ppt.1%>%
  subset(n_treat_years == 0)%>%
  dplyr::select(site_code, block,plot, subplot)%>%
  unique()%>%
  group_by(site_code)%>%
  dplyr::summarize(reps = n())
  
gamma <- cover_ppt.1%>%
          subset(n_treat_years == 0)%>%
          dplyr::select(site_code, Taxon)%>%
          unique()%>%
          group_by(site_code)%>%
          dplyr::summarize(richness = n())%>%
          left_join(nreps, by = "site_code")%>%
          mutate(gamma_rich = richness/reps)%>%
          dplyr::select(site_code, gamma_rich)

dominance <- cover_ppt.1%>%
          subset(n_treat_years == 0)%>%
          dplyr::select(site_code, block, plot, subplot, Taxon, max_cover)%>%
          unique()%>%
          group_by(site_code, block, plot, subplot)%>%
          dplyr::summarize(dominance = max(max_cover)/sum(max_cover) )%>%
          group_by(site_code)%>%
          dplyr::summarize(bp_dominance = mean(dominance))
          
  
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
            

relprecip_by_site_n_trt_year <- cover_ppt%>%
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
            mutate(multyear.relprecip = (relprecip.1+relprecip.2+relprecip.3+relprecip.4)/4, site.year.dummy = paste(site_code, year, sep = "-"))#, n_treat_years = as.character(n_treat_years))



##bray vs jaccard
#simple treatment
mod <- feols(mean_dist.bray ~ trt|as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.df)
summary(mod, vcov = vcovHAC, cluster = ~site_code + n_treat_years)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray <- data.frame(metric = "Bray", mean = coeftest(mod, vcov = NeweyWest(mod))[1], se = coeftest(mod, vcov = NeweyWest(mod))[2])



#x <- ggpredict(mod, "trt")

#ggplot(x, aes(x, predicted))+
#  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))


mod <- feols(mean_dist.jaccard ~ trt|as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.df)
summary(mod)
summary(mod, vcov = vcovHAC, cluster = ~site_code + n_treat_years)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac <- data.frame(metric = "Jaccard", mean = coeftest(mod, vcov = NeweyWest(mod))[1], se = coeftest(mod, vcov = NeweyWest(mod))[2])


bray%>%
  rbind(jac)%>%
ggplot(aes(metric, mean))+
  geom_pointrange(aes(ymin = mean-se, ymax = mean+se))+
  geom_hline(yintercept = 0)+
  theme_base()






##relprecip##relprecipcluster = 
mod <- feols(mean_dist.bray ~ relprecip.1|as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray <- data.frame(metric = "Bray", intercept = mean(mod$sumFE), slope = coeftest(mod, vcov = NeweyWest(mod))[1], se = sd(mod$sumFE)/sqrt(40))


mod <- feols(mean_dist.jaccard ~ relprecip.1|as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.df)
#summary(mod)coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac <- data.frame(metric = "Jaccard", intercept = mean(mod$sumFE), slope = coeftest(mod, vcov = NeweyWest(mod))[1], se = sd(mod$sumFE)/sqrt(40))

#bray
bray%>%
  ggplot(aes( ))+
  ylim(0,0.7)+
  xlim(-1,1)+
  geom_abline(aes(slope = slope, intercept = intercept), color = "blue")+
  geom_abline(aes(slope = slope, intercept = intercept - se), linetype = "dotted", color = "blue")+
  geom_abline(aes(slope = slope, intercept = intercept + se), linetype = "dotted", color = "blue")+
  geom_point(data = dist.df, aes(x=relprecip.1, y=mean_dist.bray, color = habitat.type))+
  geom_vline(xintercept = 0)+
  xlab("Relative precipitation")+
  ylab("Beta diversity (Bray-Curtis)")+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/relprecip_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 5.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

#jaccard
jac%>%
  ggplot(aes( ))+
  #ylim(0,0.8)+
  xlim(-1,1)+
  geom_abline(aes(slope = slope, intercept = intercept), color = "blue")+
  geom_abline(aes(slope = slope, intercept = intercept - se), linetype = "dotted", color = "blue")+
  geom_abline(aes(slope = slope, intercept = intercept + se), linetype = "dotted", color = "blue")+
  geom_point(data = dist.df, aes(x=relprecip.1, y=mean_dist.jaccard, color = habitat.type))+
  geom_vline(xintercept = 0)+
  xlab("Relative precipitation")+
  ylab("Beta diversity (Jaccard)")+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/relprecip_jac.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 5.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

##Over time
mod <- feols(mean_dist.bray ~ trt * n_treat_years|site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.df, aes(n_treat_years, mean_dist.bray, color = trt))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Treatment year")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/time_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

mod <- feols(mean_dist.jaccard ~ trt * n_treat_years|site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.df, aes(n_treat_years, mean_dist.jaccard, color = trt))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Treatment year")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/time_jac.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

mod <- feols(mean_dist.bray ~ relprecip.1 * n_treat_years|site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.jaccard ~ relprecip.1 * n_treat_years|site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.bray ~ relprecip.1 * relprecip.2 |site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.jaccard ~ relprecip.1 * relprecip.2 |site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.jaccard ~ relprecip.1 * relprecip.2 * relprecip.3 |site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.bray ~ relprecip.1 * relprecip.2 * relprecip.3 |site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.bray ~ multyear.relprecip, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.jaccard ~ multyear.relprecip, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))


##################
##site attribute moderators
#mean annual precipitation
dist.climate <- dist.df%>%
  left_join(climate, by = "site_code")

mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*MAP |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.climate)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray_map_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*MAP |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.climate)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac_map_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])


mod <- feols(mean_dist.bray ~ trt + trt*MAP |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.climate)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.climate, aes(MAP,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("MAP")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/map_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

mod <- feols(mean_dist.jaccard ~ trt + trt*MAP |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.climate)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.climate, aes(MAP,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("MAP")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/map_jac.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

#gamma diversity
dist.gamma <- dist.climate%>%
  left_join(gamma, by = "site_code")

mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*gamma_rich |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.gamma)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray_gam_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*gamma_rich |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.gamma)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac_gam_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.bray ~ trt + trt*gamma_rich |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.gamma)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.gamma, aes(gamma_rich,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Gamma richness")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/gamma_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

mod <- feols(mean_dist.jaccard ~ trt + trt*gamma_rich |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.gamma)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.gamma, aes(gamma_rich,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Gamma richness")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/gamma_jac.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

#percent of annual species cover
dist.prop <- dist.gamma%>%
  left_join(prop, by = "site_code")

mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*PctAnnual |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.prop)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray_ann_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*PctAnnual |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.prop)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac_ann_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.bray ~ trt + trt*PctAnnual |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.prop)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.prop, aes(PctAnnual,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Percent Annual")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/annual_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

mod <- feols(mean_dist.jaccard ~ trt + trt*PctAnnual |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.prop)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.prop, aes(PctAnnual,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Percent Annual")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/annual_jac.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

#nestedness of the community pretreatment
#calculate that nestedness component
library(betapart)

pretreatment <- cover_ppt.1%>%
                subset(n_treat_years == 0)

site_vector <- unique(pretreatment$site_code)

nestedness_master <- {}

for(i in 1:length(site_vector)) {
  
  temp.df <- subset(pretreatment, site_code == site_vector[i])
  
  temp.wide <- temp.df%>%
    dplyr::select(site_code, block, plot, subplot, Taxon, max_cover)%>%
    mutate(max_cover = ifelse(max_cover>0, 1, 0))%>%
    pivot_wider(names_from = Taxon, values_from = max_cover, values_fill = 0)
  
  temp.beta <- beta.multi(temp.wide[5:ncol(temp.wide)]) 
  #SNE is nestedness w/ sorensen
  #SOR is total sorensen
  
  nestedness_temp <- data.frame(site_code = site_vector[i], nestedness = temp.beta$beta.SNE, total = temp.beta$beta.SOR, proportion.nestedness = temp.beta$beta.SNE/temp.beta$beta.SOR)
                               
                              
  nestedness_master <- rbind(nestedness_master, nestedness_temp )
  rm(temp.df, temp.wide, temp.beta, nestedness_temp)
  
}

head(nestedness_master)

dist.nest <- dist.prop%>%
  left_join(nestedness_master, by = "site_code")

mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*proportion.nestedness |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.nest)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray_nest_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*proportion.nestedness |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.nest)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac_nest_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.bray ~ trt + trt*proportion.nestedness |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.nest)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.nest, aes(proportion.nestedness,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Nestedness")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/nest_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)

mod <- feols(mean_dist.jaccard ~ trt + trt*proportion.nestedness |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.nest)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.nest, aes(proportion.nestedness,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Nestedness")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/nest_jac.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)


#dominance of community pre-treatment
dist.dom <- dist.nest%>%
  left_join(dominance, by = "site_code")


mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*bp_dominance |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray_dom_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*bp_dominance |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac_dom_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.bray ~ trt + trt*bp_dominance |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(bp_dominance,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Dominance")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/dom_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)


mod <- feols(mean_dist.jaccard ~ trt + trt*bp_dominance |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(bp_dominance,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("Dominance")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/dom_jac.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 4.5,
  height = 3,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)



#soil heterogeneity
mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*ph_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*ph_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.bray ~ trt + trt*ph_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(ph_var,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("ph variation")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/ph_bray.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)


mod <- feols(mean_dist.jaccard ~ trt + trt*ph_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(ph_var,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("ph variation")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/ph_jac.pdf",
        plot = last_plot(),
        device = "pdf",
        path = NULL,
        scale = 1,
        width = 4.5,
        height = 3,
        units = c("in"),
        dpi = 600,
        limitsize = TRUE
)


###Figures of moderator effects
bray_map_stats%>%
  rbind(bray_gam_stats)%>%
  rbind(bray_ann_stats)%>%
  rbind(bray_nest_stats)%>%
  rbind(bray_dom_stats)%>%
ggplot(aes(moderator, estimate))+
  geom_pointrange(aes(ymin = estimate-std.error, ymax = estimate+std.error))+
  geom_hline(yintercept = 0)+
  theme_base()


jac_map_stats%>%
  rbind(jac_gam_stats)%>%
  rbind(jac_ann_stats)%>%
  rbind(jac_nest_stats)%>%
  rbind(jac_dom_stats)%>%
  ggplot(aes(moderator, estimate))+
  geom_pointrange(aes(ymin = estimate-std.error, ymax = estimate+std.error))+
  geom_hline(yintercept = 0)+
  theme_base()

