

library(ggeffects)
library(tidyverse)
library(vegan)
library(nlme)
library(fixest)
library(sandwich)
library(RcmdrMisc)
library(lmtest)
library(ggthemes)
library(kernelshap)   #  General SHAP
library(shapviz)      #  SHAP plots
library(grf)


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
          subset(site_code != "wildflower.us")%>%#can be added in later once precip data is incorporated for wildflower.us
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
  geom_point(data = dist.df, aes(x=relprecip.1, y=mean_dist.bray, color = habitat.type),shape = 1)+
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
  geom_point(data = dist.df, aes(x=relprecip.1, y=mean_dist.jaccard, color = habitat.type),shape = 1)+
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
  geom_point(shape = 1)+
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
  geom_point(shape = 1)+
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



###Precip lags with causal forest
drt.sev <- dist.df%>%
  ungroup()%>%
            subset(trt == "Drought")%>%
            dplyr::select(site_code, n_treat_years, relprecip.1, relprecip.2, relprecip.3, relprecip.4)%>%
            group_by(site_code, n_treat_years)%>%
            dplyr::summarize(drtsev.1 = mean(relprecip.1),
                             drtsev.2 = mean(relprecip.2),
                             drtsev.3 = mean(relprecip.3),
                             drtsev.4 = mean(relprecip.4))


# As outcomes we'll look at the number of correct answers.
Y <- dist.df$mean_dist.bray
W <- dist.df%>%
  ungroup()%>%
  mutate(trt_num = ifelse(trt=="Control", 0, 1))%>%
  pull(trt_num)
#W <- dist.df$relprecip.1
#X <- dist.df%>%
#  ungroup()%>%
#  dplyr::select(relprecip.1,relprecip.2,relprecip.3,relprecip.4)
X <- dist.df%>%
  ungroup()%>%
  left_join(drt.sev, by = c("site_code", "n_treat_years"))%>%
  dplyr::select(drtsev.1,drtsev.2,drtsev.3,drtsev.4)

rf <- regression_forest(X, W, num.trees = 5000)
p.hat <- predict(rf)$predictions

hist(p.hat)

Y.forest <- regression_forest(X, Y) #this function doesn't know the treatment but that's the whole point
Y.hat <- predict(Y.forest)$predictions

varimp.Y <- variable_importance(Y.forest)

# Keep the top 10 variables for CATE estimation
keep <- colnames(X)[order(varimp.Y, decreasing = TRUE)[1:4]]
keep
#"relprecip.1" "relprecip.2" "relprecip.3" "relprecip.4"     

X.cf <- X[, keep]
W.hat <- 0.5

# Set aside the first half of the data for training and the second for evaluation.
# (Note that the results may change depending on which samples we hold out for training/evaluation)
train <- sample(1:nrow(X.cf), size = floor(0.5 * nrow(X.cf)))#random sample of data to train instead of jut the first half of the dataset

train.forest <- causal_forest(X.cf[train, ], Y[train], W[train], Y.hat = Y.hat[train], W.hat = W.hat)
tau.hat.eval <- predict(train.forest, X.cf[-train, ])$predictions

eval.forest <- causal_forest(X.cf[-train, ], Y[-train], W[-train], Y.hat = Y.hat[-train], W.hat = W.hat)

average_treatment_effect(eval.forest)
#estimate     std.err 
#0.003091148 0.016592967 

varimp <- variable_importance(eval.forest)
ranked.vars <- order(varimp, decreasing = TRUE)
colnames(X.cf)[ranked.vars[1:4]]
#relprecip.4" "relprecip.3" "relprecip.2" "relprecip.1"

rate.cate <- rank_average_treatment_effect(eval.forest, list(cate = -1 *tau.hat.eval))
#rate.age <- rank_average_treatment_effect(eval.forest, list(map = X[-train, "map"]))

plot(rate.cate, ylab = "Number of correct answers", main = "TOC: By most negative CATEs")
#plot(rate.age, ylab = "Number of correct answers", main = "TOC: By decreasing map")

#xvars <- c("ppt.1", "ppt.2", "ppt.3", "ppt.4", "n_treat_days", "n_treat_years", "map", "arid", "PctAnnual", "PctGrass", "sand_mean", "AI", "cv_ppt_inter", "richness", "seasonality_index", "r_monthly_t_p")
imp <- sort(setNames(variable_importance(eval.forest), keep))
#par(mai = c(0.7, 2, 0.2, 0.2))
barplot(imp, horiz = TRUE, las = 1, col = "orange")

pred_fun <- function(object, newdata, ...) {
  predict(object, newdata, ...)$predictions
}
library(hstats)
pdps <- lapply(colnames(X.cf[-train, ]), function(v) plot(partial_dep(eval.forest, v=v, X = X.cf[-train, ], pred_fun = pred_fun
)))
library(patchwork)
wrap_plots(pdps, guides = "collect", ncol = 5) &
  #  ylim(c(-0.11, -0.06)) &
  ylab("Treatment effect")

H <- hstats(eval.forest, X = X, pred_fun = pred_fun, verbose = FALSE)
plot(H)
#partial_dep(eval.forest, v = "map", X = X, pred_fun = pred_fun) |> 
 # plot()

partial_dep(eval.forest, v = colnames(X.cf[-train, ]), X = X.cf[-train, ])

# Explaining one CATE
kernelshap(eval.forest, X = X.cf[-train, ], bg_X = X, 
           pred_fun = pred_fun) |> 
  shapviz() |> 
  sv_waterfall() +
  xlab("Prediction")

# Explaining all CATEs globally
system.time(  # 13 min
  ks <- kernelshap(eval.forest, X = X.cf[-train, ], pred_fun = pred_fun)  
)
shap_values <- shapviz(ks)
sv_importance(shap_values)
sv_importance(shap_values, kind = "bee")
#sv_dependence(shap_values, v = xvars) +
#  plot_layout(ncol = 3) &
#  ylim(c(-0.04, 0.03))


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


#####Site ANPP
anpp <- read.csv("C:/Users/ohler/Dropbox/IDE/data_processed/anpp_ppt_2025-10-20.csv")%>%
  subset(n_treat_years == 0)%>%
  dplyr::select(site_code, mass)%>%
  group_by(site_code)%>%
  dplyr::summarise(anpp = mean(mass))

dist.anpp <- dist.dom%>%
  left_join(anpp, by = "site_code")


mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*anpp |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.anpp)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray_anpp_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*anpp |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.anpp)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
jac_anpp_stats <- data.frame(moderator = rownames(coeftest(mod, vcov = NeweyWest(mod)))[2], estimate = coeftest(mod, vcov = NeweyWest(mod))[2,1], std.error = coeftest(mod, vcov = NeweyWest(mod))[2,2], t.value = coeftest(mod, vcov = NeweyWest(mod))[2,3], p.value = coeftest(mod, vcov = NeweyWest(mod))[2,4])

mod <- feols(mean_dist.bray ~ trt + trt*anpp |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.anpp)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.anpp, aes(anpp,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("ANPP")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/anpp_bray.pdf",
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


mod <- feols(mean_dist.jaccard ~ trt + trt*anpp |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.anpp)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.anpp, aes(anpp,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("ANPP")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/anpp_jac.pdf",
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

###soil heterogeneity
#first ph
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


#second N variability
mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*n_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*n_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.bray ~ trt + trt*n_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(n_var,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("N variation")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/n_bray.pdf",
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


mod <- feols(mean_dist.jaccard ~ trt + trt*n_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(n_var,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("N variation")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/n_jac.pdf",
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


#Third P variability
mod <- feols(mean_dist.bray ~ relprecip.1 + relprecip.1*p_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.jaccard ~ relprecip.1 + relprecip.1*p_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

mod <- feols(mean_dist.bray ~ trt + trt*p_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(p_var,mean_dist.bray,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("P variation")+
  ylab("Bray distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/p_bray.pdf",
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


mod <- feols(mean_dist.jaccard ~ trt + trt*p_var |as.factor(n_treat_years)+site_code, cluster = ~site_code, data = dist.dom)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))

ggplot(dist.dom, aes(p_var,mean_dist.jaccard,color = trt))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm")+
  xlab("N variation")+
  ylab("Jaccard distance")+
  scale_color_manual(values = c("black", "#D35721"))+
  theme_base()

ggsave( "C:/Users/ohler/Dropbox/Tim+Laura/Beta diversity/figures/p_jac.pdf",
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



###Causal forest for moderators
# As outcomes we'll look at the number of correct answers.
Y <- dist.anpp$mean_dist.bray
W <- dist.anpp%>%
  ungroup()%>%
  mutate(trt_num = ifelse(trt=="Control", 0, 1))%>%
  pull(trt_num)
#W <- dist.anpp$relprecip.1
#X <- dist.anpp%>%
#  ungroup()%>%
#  dplyr::select(relprecip.1,relprecip.2,relprecip.3,relprecip.4)
X <- dist.anpp%>%
  ungroup()%>%
  dplyr::select(map,cv_ppt_intra, cv_ppt_inter, seasonality_index, aridity_index,MAT, gamma_rich, PctAnnual, PctGrass, proportion.nestedness, bp_dominance, anpp)

rf <- regression_forest(X, W, num.trees = 5000)
p.hat <- predict(rf)$predictions

hist(p.hat)

Y.forest <- regression_forest(X, Y) #this function doesn't know the treatment but that's the whole point
Y.hat <- predict(Y.forest)$predictions

varimp.Y <- variable_importance(Y.forest)

# Keep the top 10 variables for CATE estimation
keep <- colnames(X)[order(varimp.Y, decreasing = TRUE)[1:10]]
keep
#[1] "PctAnnual"             "MAT"                   "PctGrass"             
#[4] "cv_ppt_inter"          "bp_dominance"          "proportion.nestedness"
#[7] "gamma_rich"            "aridity_index"         "anpp"                 
#[10] "map"         

X.cf <- X[, keep]
W.hat <- 0.5

# Set aside the first half of the data for training and the second for evaluation.
# (Note that the results may change depending on which samples we hold out for training/evaluation)
train <- sample(1:nrow(X.cf), size = floor(0.5 * nrow(X.cf)))#random sample of data to train instead of jut the first half of the dataset

train.forest <- causal_forest(X.cf[train, ], Y[train], W[train], Y.hat = Y.hat[train], W.hat = W.hat)
tau.hat.eval <- predict(train.forest, X.cf[-train, ])$predictions

eval.forest <- causal_forest(X.cf[-train, ], Y[-train], W[-train], Y.hat = Y.hat[-train], W.hat = W.hat)

average_treatment_effect(eval.forest)
#estimate    std.err 
#0.01858178 0.01068278 

varimp <- variable_importance(eval.forest)
ranked.vars <- order(varimp, decreasing = TRUE)
colnames(X.cf)[ranked.vars[1:10]]
#[1] "PctGrass"              "gamma_rich"            "bp_dominance"         
#[4] "map"                   "aridity_index"         "cv_ppt_inter"         
#[7] "PctAnnual"             "anpp"                  "MAT"                  
#[10] "proportion.nestedness"

rate.cate <- rank_average_treatment_effect(eval.forest, list(cate = -1 *tau.hat.eval))
#rate.age <- rank_average_treatment_effect(eval.forest, list(map = X[-train, "map"]))

plot(rate.cate, ylab = "Number of correct answers", main = "TOC: By most negative CATEs")
#plot(rate.age, ylab = "Number of correct answers", main = "TOC: By decreasing map")

#xvars <- c("ppt.1", "ppt.2", "ppt.3", "ppt.4", "n_treat_days", "n_treat_years", "map", "arid", "PctAnnual", "PctGrass", "sand_mean", "AI", "cv_ppt_inter", "richness", "seasonality_index", "r_monthly_t_p")
imp <- sort(setNames(variable_importance(eval.forest), keep))
#par(mai = c(0.7, 2, 0.2, 0.2))
barplot(imp, horiz = TRUE, las = 1, col = "orange")

pred_fun <- function(object, newdata, ...) {
  predict(object, newdata, ...)$predictions
}
library(hstats)
pdps <- lapply(colnames(X.cf[-train, ]), function(v) plot(partial_dep(eval.forest, v=v, X = X.cf[-train, ], pred_fun = pred_fun
)))
library(patchwork)
wrap_plots(pdps, guides = "collect", ncol = 5) &
  #  ylim(c(-0.11, -0.06)) &
  ylab("Treatment effect")

#H <- hstats(eval.forest, X = X, pred_fun = pred_fun, verbose = FALSE)
#plot(H)
#partial_dep(eval.forest, v = "map", X = X, pred_fun = pred_fun) |> 
# plot()

partial_dep(eval.forest, v = colnames(X.cf[-train, ]), X = X.cf[-train, ])

# Explaining one CATE
kernelshap(eval.forest, X = X.cf[-train, ], bg_X = X, 
           pred_fun = pred_fun) |> 
  shapviz() |> 
  sv_waterfall() +
  xlab("Prediction")

# Explaining all CATEs globally
system.time(  # 13 min
  ks <- kernelshap(eval.forest, X = X.cf[-train, ], pred_fun = pred_fun)  
)
shap_values <- shapviz(ks)
sv_importance(shap_values)
sv_importance(shap_values, kind = "bee")
#sv_dependence(shap_values, v = xvars) +
#  plot_layout(ncol = 3) &
#  ylim(c(-0.04, 0.03))



################################
########Turnover and nestedness as response variables
#calculate these for each site/treatment/year
#loop
site_year_trt_vector <- cover_ppt%>%
                        unite(col = site_year_trt, site_code, n_treat_years, trt, sep = "::")%>%
  distinct(site_year_trt) %>%
  pull(site_year_trt)


#%>%
 #                       dplyr::select(site_year_trt)%>%
  #                      unique()%>%
   #                     as.vector()
  
  

turn_nest_master <- {}

for(i in 1:length(site_year_trt_vector)) {
  
  temp.df <- cover_ppt%>%
    unite(col = "site_year_trt", site_code, n_treat_years, trt, sep = "::", remove = TRUE)%>%
  subset( site_year_trt == site_year_trt_vector[i])
  
  temp.wide <- temp.df%>%
    dplyr::select(site_year_trt, block, plot, subplot, Taxon, max_cover)%>%
    mutate(max_cover = ifelse(max_cover>0, 1, 0))%>%
    pivot_wider(names_from = Taxon, values_from = max_cover, values_fill = 0)
  
  temp.beta <- beta.multi(temp.wide[5:ncol(temp.wide)]) 
  #SNE is nestedness w/ sorensen
  #SOR is total sorensen
  
  turn_nest_temp <- data.frame(site_year_trt = site_year_trt_vector[i], nestedness = temp.beta$beta.SNE, total = temp.beta$beta.SOR, proportion.nestedness = temp.beta$beta.SNE/temp.beta$beta.SOR)
  
  
  turn_nest_master <- rbind(turn_nest_master, turn_nest_temp )
  rm(temp.df, temp.wide, temp.beta, turn_nest_temp)
  
}


turn_nest <- turn_nest_master%>%
  separate(site_year_trt, c("site_code", "n_treat_years", "trt"), sep = "::")%>%
  mutate(n_treat_years = as.numeric(n_treat_years))%>%
  left_join(dist.df, by = c("site_code", "n_treat_years", "trt"))

mod <- feols(proportion.nestedness ~ trt|as.factor(n_treat_years)+site_code, cluster = ~site_code, data = turn_nest)
summary(mod, vcov = vcovHAC, cluster = ~site_code + n_treat_years)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray <- data.frame(metric = "Bray", mean = coeftest(mod, vcov = NeweyWest(mod))[1], se = coeftest(mod, vcov = NeweyWest(mod))[2])

##relprecip##relprecipcluster = 
mod <- feols(proportion.nestedness ~ relprecip.1|as.factor(n_treat_years)+site_code, cluster = ~site_code, data = turn_nest)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Quadratic Spectral"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Truncated"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Parzen"))
coeftest(mod, vcov = kernHAC(mod, kernel = "Tukey-Hanning"))
coeftest(mod, vcov = NeweyWest(mod))
bray <- data.frame(metric = "Bray", intercept = mean(mod$sumFE), slope = coeftest(mod, vcov = NeweyWest(mod))[1], se = sd(mod$sumFE)/sqrt(40))













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

