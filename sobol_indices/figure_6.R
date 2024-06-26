###############################################################################
### Analysis of emulator built on sequential design draws from dfnWorks
### Using Global SA: A Primer by A. Saltelli & Saltelli 2010
## Author: Alexander C. Murph
## Date: September 2023
library(hetGP)
library(ggplot2)
library(reshape2)
library(gridExtra)
setwd("~/GitHub/sobol_dfnworks/sobol_indices")
source("MC_sobol_indices.R")

# Grab the final hetGP model:
models_loc = "~/GitHub/sobol_dfnworks/GP_analysis"
model      = NULL
load(file=paste(models_loc, "/hetGP_model_final.Rdata", sep = ""))

# Saltelli suggests at least 500 samples:
num_samples     = 1000
n_bootstraps    = 5000
input_dimension = 5


#########################
## Before we do the relative analysis, we need to examine the total effect
## of the seed variable:
setwd("~/GitHub/sobol_dfnworks/sobol_indices")
models_loc = "~/GitHub/sobol_dfnworks/GP_analysis"
model = NULL
load(file=paste(models_loc, "/hetGP_model_final.Rdata", sep = "")) # model
load(file=paste(models_loc, "/sparse_lhs_normalized_data_woOutliers.Rdata", sep = "")) # norm_data
norm_data = orig_norm_data

# Get the data in a form to perform a prediction using the dispersion GP
full_X = NULL
for(val_i in 1:nrow(norm_data$X0)){
  temp_rows = matrix(rep(unlist(as.vector(norm_data$X0[val_i,])), times = norm_data$mult[val_i]),
                     ncol = ncol(norm_data$X0), nrow = norm_data$mult[val_i], byrow = TRUE)
  full_X    = rbind(full_X, temp_rows)
}
Z_train          = norm_data$Z
colnames(full_X) = names(norm_data$X0)

# Get the predictions from the dispersion GP
predict_values     = predict(object=model, x = full_X)
dispersion_GP_vals = predict_values$nugs
avg_dis_GP_vals    = mean(dispersion_GP_vals)

#######
# There are two ways to do this.  Way 1: we use the sample variance directly:
var_y = var(norm_data$full_outputs)

# Predict the total effect (option 1):
avg_dis_GP_vals / var_y
# 0.05656231

# There are two ways to do this.  Way 1: use the estimates from the joint GP:
predict_values     = predict(object=model, x = full_X)
mean_GP_vals       = predict_values$mean
var_mean_GP_vals   = var(mean_GP_vals)

# Predict the total effect (option 2):
avg_dis_GP_vals / (var_mean_GP_vals + avg_dis_GP_vals)


########################
## This bootstrap this boi!
# Load in the final GP model:
load(file=paste(models_loc, "/hetGP_model_final.Rdata", sep = "")) # model
# Load in the full training data:
load(file=paste(models_loc, "/sparse_lhs_normalized_data_woOutliers.Rdata", sep = "")) # norm_data
# Load in the testing data:
load(file = paste(models_loc, "/compiled_testing_data_woOutliers.Rdata", sep = "")) # testing_data

orig_model   = model

full_data_Z         = norm_data$full_outputs
full_data_X         = norm_data$full_inputs
full_data_X$sim_num = NULL
full_data_X$ind     = NULL
bootstrapped_FOmean = NULL
bootstrapped_TEmean = NULL
bootstrapped_FOsd2  = NULL
bootstrapped_TEsd2  = NULL
bootstrapped_TEseed = c()

n_obs = nrow(full_data_X)
for(boot_idx in 1:n_bootstraps){
  print(paste("Starting bootstrap iteration:", boot_idx))
  resample_idxs = sample(1:n_obs, size = n_obs, replace = TRUE)
  resampled_X = full_data_X[resample_idxs,]
  resampled_Z = full_data_Z[resample_idxs]
  reps_list = find_reps(X = as.matrix(resampled_X), Z = as.vector(resampled_Z))
  model = mleHetGP( X = list(X0 = reps_list$X0, Z0 = reps_list$Z0, mult = reps_list$mult),
                    noiseControl = list(g_min = 3),
                    Z = reps_list$Z, covtype = "Matern3_2", 
                    settings = list(checkHom = FALSE ) )
  temp_SAs = estimate_sobol_indices(model, num_samples, input_dimension, 
                              var_names = c("alpha_semi", "beta_semi", 
                                            "sigma_semi", "alpha_TPL", "p32"))
  
  temp_row = data.frame(alpha_semi = temp_SAs[['first_order_mean']]$alpha_semi,
                        beta_semi  = temp_SAs[['first_order_mean']]$beta_semi,
                        sigma_semi = temp_SAs[['first_order_mean']]$sigma_semi,
                        alpha_TPL  = temp_SAs[['first_order_mean']]$alpha_TPL,
                        p32        = temp_SAs[['first_order_mean']]$p32)
  bootstrapped_FOmean = rbind(bootstrapped_FOmean, temp_row)
  
  temp_row = data.frame(alpha_semi = temp_SAs[['total_effect_mean']]$alpha_semi,
                        beta_semi  = temp_SAs[['total_effect_mean']]$beta_semi,
                        sigma_semi = temp_SAs[['total_effect_mean']]$sigma_semi,
                        alpha_TPL  = temp_SAs[['total_effect_mean']]$alpha_TPL,
                        p32        = temp_SAs[['total_effect_mean']]$p32)
  bootstrapped_TEmean = rbind(bootstrapped_TEmean, temp_row)
  
  temp_row = data.frame(alpha_semi = temp_SAs[['first_order_sd2']]$alpha_semi,
                        beta_semi  = temp_SAs[['first_order_sd2']]$beta_semi,
                        sigma_semi = temp_SAs[['first_order_sd2']]$sigma_semi,
                        alpha_TPL  = temp_SAs[['first_order_sd2']]$alpha_TPL,
                        p32        = temp_SAs[['first_order_sd2']]$p32)
  bootstrapped_FOsd2 = rbind(bootstrapped_FOsd2, temp_row)
  
  temp_row = data.frame(alpha_semi = temp_SAs[['total_effect_sd2']]$alpha_semi,
                        beta_semi  = temp_SAs[['total_effect_sd2']]$beta_semi,
                        sigma_semi = temp_SAs[['total_effect_sd2']]$sigma_semi,
                        alpha_TPL  = temp_SAs[['total_effect_sd2']]$alpha_TPL,
                        p32        = temp_SAs[['total_effect_sd2']]$p32)
  bootstrapped_TEsd2 = rbind(bootstrapped_TEsd2, temp_row)
  
  ############################
  # Get the predictions from the dispersion GP
  predict_values     = predict(object=model, x = full_X)
  dispersion_GP_vals = predict_values$nugs
  mean_GP_vals       = predict_values$mean
  avg_dis_GP_vals    = mean(dispersion_GP_vals)
  var_mean_GP_vals   = var(mean_GP_vals)
  
  # Predict the total effect:
  bootstrapped_TEseed = c(bootstrapped_TEseed, avg_dis_GP_vals / (var_mean_GP_vals + avg_dis_GP_vals))
  
}

write.csv(bootstrapped_FOmean, file = "bootstrapped_FOmean.csv")
write.csv(bootstrapped_TEmean, file = "bootstrapped_TEmean.csv")
write.csv(bootstrapped_FOsd2,  file = "bootstrapped_FOsd2.csv")
write.csv(bootstrapped_TEsd2,  file = "bootstrapped_TEsd2.csv")

bootstrapped_TEseed_df = data.frame(TEseed = bootstrapped_TEseed)
write.csv(bootstrapped_TEseed_df,  file = "bootstrapped_TEseed_df.csv")


#########################
## Let's make some boxplots :)
setwd("~/GitLab/sobol_dfnworks/sobol_indices")

# First Order Indices on the Mean Response GP:
bootstrapped_FOmean = read.csv(file = "bootstrapped_FOmean.csv")
bootstrapped_FOmean$X = NULL
graph_data = melt(bootstrapped_FOmean)
levels(graph_data$variable) = c(r"($\alpha$)", r"($\beta$)", r"($\sigma$)", r"($\gamma$)", r"(p32)")
g1 = ggplot(graph_data, aes(x = variable, y = value)) + geom_boxplot() + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  scale_x_discrete(labels = TeX(c(r"($\alpha$)", r"($\beta$)", r"($\sigma$)", r"($\gamma$)", r"(p32)")) ) + 
  xlab("Simulation Input Parameters") + 
  ylab("Monte Carlo Samples of \nFirst-Order Sobol' Index") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15), plot.title = element_text(size = 15)) +
  ylim(0,0.5)

bootstrapped_TEmean    = read.csv(file = "bootstrapped_TEmean.csv")
bootstrapped_TEseed_df = read.csv(file = "bootstrapped_TEseed_df.csv")
bootstrapped_TEmean$X = NULL
bootstrapped_TEmean$seed = bootstrapped_TEseed_df$TEseed
graph_data = melt(bootstrapped_TEmean)
levels(graph_data$variable) = c(r"($\alpha$)", r"($\beta$)", r"($\sigma$)", r"($\gamma$)", r"(p32)", r"($\epsilon$)")
g2 = ggplot(graph_data, aes(x = variable, y = value)) + geom_boxplot() + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  scale_x_discrete(labels = TeX(c(r"($\alpha$)", r"($\beta$)", r"($\sigma$)", r"($\gamma$)", r"(p32)", r"($\epsilon$)")) ) + 
  xlab("Simulation Input Parameters and Seed Variable") + 
  ylab("Monte Carlo Samples of \nTotal-Effect Sobol' Index") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15), plot.title = element_text(size = 15)) +
  ylim(0,0.5)

lst_p = list(g1, g2)
grid.arrange(lst_p[[1]], lst_p[[2]],
             layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))

#########################
## Now see if we can get these boxplots all on the same plot:
setwd("~/GitLab/sobol_dfnworks/sobol_indices")
bootstrapped_FOmean           = NULL
bootstrapped_FOmean           = read.csv(file = "bootstrapped_FOmean.csv")
bootstrapped_FOmean$X         = NULL
graph_data_1                  = melt(bootstrapped_FOmean)
levels(graph_data_1$variable) = c(r"($\alpha$)", r"($\beta$)", r"($\sigma$)", r"($\gamma$)", r"(p32)")
graph_data_1$Sobol_Index      = rep("First-Order", times = nrow(graph_data_1))

bootstrapped_TEmean           = read.csv(file = "bootstrapped_TEmean.csv")
bootstrapped_TEseed_df        = read.csv(file = "bootstrapped_TEseed_df.csv")
bootstrapped_TEmean$X         = NULL
bootstrapped_TEmean$seed      = bootstrapped_TEseed_df$TEseed
graph_data_2                  = melt(bootstrapped_TEmean)
levels(graph_data_2$variable) = c(r"($\alpha$)", r"($\beta$)", r"($\sigma$)", r"($\gamma$)", r"(p32)", r"($\epsilon$)")
graph_data_2$Sobol_Index      = rep("Total-Effect", times = nrow(graph_data_2))

graph_data = rbind(graph_data_1, graph_data_2)
ggplot(graph_data, aes(x = variable, y = value, color = Sobol_Index), fill = 'white') + geom_boxplot(position = position_dodge2(preserve = "single")) + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  scale_x_discrete(labels = TeX(c(r"($\alpha$)", r"($\beta$)", r"($\sigma$)", r"($\gamma$)", r"(p32)", r"($\epsilon$)")) ) + 
  xlab("Simulation Input Parameters and Seed Variable") + 
  ylab("Monte Carlo Samples \nof Sobol' Indices") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15), 
        plot.title = element_text(size = 15),legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  # ylim(0,0.5) +
  scale_color_manual(values = c("First-Order" = "black",
  "Total-Effect"="grey")) #+
  # scale_fill_manual(values = c("First-Order" = "grey",
                                # "Total-Effect"="white"))


#########################
## Doing boxplots on the dispersion GP would be misleading.
## Let's instead do a table.
setwd("~/GitLab/sobol_dfnworks/sobol_indices")
f = function(x){quantile(x,prob=c(0.025,0.975))}

bootstrapped_FOsd2   = read.csv(file = "bootstrapped_FOsd2.csv")
bootstrapped_FOsd2$X = NULL
apply(bootstrapped_FOsd2, 2, f)

bootstrapped_TEsd2  = read.csv(file = "bootstrapped_TEsd2.csv")
bootstrapped_TEsd2$X = NULL
apply(bootstrapped_TEsd2, 2, f)


bootstrapped_TEseed_df = read.csv(file = "bootstrapped_TEseed_df.csv")
bootstrapped_TEseed_df$X = NULL
apply(bootstrapped_TEseed_df, 2, f)
summary(bootstrapped_TEseed_df)

f = function(x){sd(x)}
apply(bootstrapped_TEseed_df, 2, f)



