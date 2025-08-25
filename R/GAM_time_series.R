vgetwd()
setwd("C:/Users/USER/Documents/pancreas")

#load packages
library(gtools)
library(pROC)
library(ape)
library(ggdendro)
library(WGCNA)
library(stats)
library(flashClust)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(tidyverse)
library(gridExtra)
library(gplots)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(fitdistrplus)
library(MASS)
library(ggplot2)
library(mgcv)
library(DHARMa)
list.files(paste0(getwd(), "/input"))


#GLM model for time series eigengene 
eigen_data= read.csv("output/module_eigengenes.csv")[,1:2] %>% column_to_rownames("X") %>% 
            mutate(time_points= c(rep("Before", 3), rep("After.2d", 3), rep("After.7d", 3)),
            time_points_numeric= c(rep(1, 3), rep(2, 3), rep(3, 3)))

eigen_data$time_points= factor(eigen_data$time_points, level= c("Before","After.2d","After.7d" ))

modules=read.csv("output/module_genes.csv")
data= read.csv("input/pancreas_input333.csv" )
data_blue= data[data$X %in% modules$blue,]
#density plot
plot(density(eigengene$MEblue))
shapiro.test(eigengene$MEblue)

# QQ plot for Binomial distribution
# x = rbinom(1000, size=9, prob= .5)
# y= eigengene$MEblue
# qqplot(x,y, main = "QQ Plot for Binomial Distribution")


# Fit a Gaussian GLM (default)
#"NA" values in the coefficients of your glm model arises due to perfect multicollinearity among the predictors
glm_gaussian <- glm(MEblue ~ time_points,
                    data = eigen_data, family = gaussian())
AIC(glm_gaussian)
summary(glm_gaussian)

# Plot diagnostic plots for the GAM
par(mfrow = c(2, 2)) 
plot(glm_gaussian, pages = 1)

# Residuals vs Fitted
plot(glm_gaussian$fitted.values, residuals(glm_gaussian), 
     main = "Residuals vs Fitted", 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# Normal Q-Q plot
qqnorm(residuals(glm_gaussian), main = "Normal Q-Q")
qqline(residuals(glm_gaussian), col = "red")

#asses the model
#Let's check the model fit and residuals using DHARMa
simulationOutput<- simulateResiduals(fittedModel = glm_gaussian)
plot(simulationOutput)


testDispersion(glm_gaussian) 
testZeroInflation(glm_gaussian)
testOutliers(glm_gaussian) 
testUniformity(glm_gaussian)

dev.off()
#----------------------------------------------------------------------
#given the data not normally distributed, we revert to generalized models
# Fit a GAM
gam_model <- gam(MEblue ~ s(as.numeric(time_points_numeric), bs = "cs", k=3), data = eigen_data)
# Summarize the model
summary(gam_model)
plot(gam_model)
# Compare models using AIC
aic_values <- AIC(glm_gaussian, gam_model)
print(aic_values)
#_________________________________
#asses assumption of the model
# Plot the smooth terms
plot(gam_model, residuals = TRUE, pch = 19, cex = 0.5)
# Residuals vs Fitted values plot
plot(gam_model$fitted.values, residuals(gam_model), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

## Residuals vs Time
plot(as.numeric(eigengene$time_points_numeric), residuals(gam_model), 
     xlab = "Time Points", ylab = "Residuals",
     main = "Residuals vs Time Points")
abline(h = 0, col = "red")

#asses using dharma okay
# Alternatively, use anova for model comparison (if nested)
anova(glm_gaussian, gam_model, test = "Chisq")
#____________________________
## visualize
predictions <- predict(gam_model, newdata = eigen_data, se.fit = TRUE)
eigen_data= cbind(eigen_data, predictions )

# Plot the data and the GAM fit
p= ggplot() +
  geom_point(data = eigen_data, aes(x = time_points_numeric, y = MEblue)) +
  geom_line(data = eigen_data,
            aes(x = time_points_numeric, y = fit), color = "#557192", linewidth = 2) +
  geom_ribbon(data = eigen_data, aes(x = time_points_numeric,
                                     ymin = fit - 1.96 * se.fit, 
                                     ymax = fit + 1.96 * se.fit), alpha = 0.1) +
  scale_x_continuous(breaks = c(1,1.5,2,2.5,3),
                     labels = c("Before","","After.2d","","After.7d")
  )+
  labs(title = "Generalized Additive Model (GAM) Fit for ME", 
       x = "Time Points", y = "Module eigengene of blue module (ME)") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),  
    panel.grid.major = element_line(color = "grey90"),            
    panel.grid.minor = element_line(color = "grey90")           
  )

print(p)
ggsave("plots/GAM_plot.png", p,  height = 4, width = 5,   dpi = 600)
#-------------------------------------------------------------------------------
#generalized additive mixed model that allow for controlling  for multiple measurements/ autocorrelation with the previous time point
#check auto correlation
par(mfrow=c(1,2))
acf(resid(gam_model), lag.max = 36, main = "ACF")
pacf(resid(gam_model), lag.max = 36, main = "pACF")

gamm_model= gamm(MEblue ~ s(time_points_numeric, k = 3), data = eigen_data,
     correlation = corARMA(form = ~ subject, p = 1))


summary(gamm_model)

# data_blue= data[data$X %in% modules$blue,]
# row.names(data_blue)= NULL
# data_blue= data_blue %>% column_to_rownames("X")
# expdata= t(data_blue)
# expdata= as.data.frame(expdata) %>% 
#   mutate(time_points=  c(rep(1, 3), rep(2, 3), rep(3, 3)))
# 
# df_long=  expdata %>%
#   pivot_longer(
#     cols= !time_points,
#     names_to = "gene",
#     values_to = "Expression"
#   ) %>%
#   filter(gene != "" & !is.na(gene)) 


