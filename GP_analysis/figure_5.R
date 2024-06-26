###############################################################################
### After fitting the final GP, this script compares the nominal to the empirical
### vs. nominal coverage.
## Author: Alexander C. Murph
## Date: October 2023
library(Matrix)
library(hetGP)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggExtra) 
setwd("~/GitHub/sobol_dfnworks/GP_analysis")

# Load in the full training data:
load(file="sparse_lhs_normalized_data_woOutliers.Rdata") # norm_data
norm_data = orig_norm_data
# Load in the testing data:
load(file = "compiled_testing_data_woOutliers.Rdata") # testing_data

# Create the GP
curr_min            = norm_data$new_min
curr_max            = norm_data$new_max
X0                  = as.matrix(norm_data$X0)
Z0                  = norm_data$Z0
mult                = norm_data$mult
Z                   = norm_data$Z
X                   = list(X0 = X0, Z0 = Z0, mult = mult)
model_temp          = mleHetGP(X = X, Z = Z, lower = 0.0001, upper = 1)
model               = mleHetGP( X = list(X0 = model_temp$X0, Z0 = model_temp$Z0, mult = model_temp$mult),
                                noiseControl = list(g_min = 3),
                                Z = model_temp$Z, 
                                settings = list(checkHom = FALSE ) )

####################################
## The predict function requires a transformation on the data:
## Transformation on the testing data:
full_X_testing = NULL
for(val_i in 1:nrow(testing_data$X0)){
  temp_rows = matrix(rep(unlist(as.vector(testing_data$X0[val_i,])), times = testing_data$mult[val_i]),
                     ncol = ncol(testing_data$X0), nrow = testing_data$mult[val_i], byrow = TRUE)
  full_X_testing    = rbind(full_X_testing, temp_rows)
}
Z_test = testing_data$Z
colnames(full_X_testing) = names(norm_data$X0)

## Transformation on the training data:
curr_min            = norm_data$new_min
curr_max            = norm_data$new_max
X0                  = as.matrix(norm_data$X0)
Z0                  = norm_data$Z0
mult                = norm_data$mult
Z                   = norm_data$Z
X                   = list(X0 = X0, Z0 = Z0, mult = mult)

# Use these, with the test data, to calculate the metrics I am interested in.
full_X = NULL
for(val_i in 1:nrow(norm_data$X0)){
  temp_rows = matrix(rep(unlist(as.vector(norm_data$X0[val_i,])), times = norm_data$mult[val_i]),
                     ncol = ncol(norm_data$X0), nrow = norm_data$mult[val_i], byrow = TRUE)
  full_X    = rbind(full_X, temp_rows)
}
Z_train = norm_data$Z
colnames(full_X) = names(norm_data$X0)

####################################
# There is one somewhat obvious outlier in the testing data.  For the time being, I think that I'll just
# ignore it.  I since I'm not changing the training data, I shouldn't have to refit the GP.
full_X_testing = NULL
for(val_i in 1:nrow(testing_data$X0)){
  if(testing_data$Z[val_i] < 3){
    temp_rows = matrix(rep(unlist(as.vector(testing_data$X0[val_i,])), times = testing_data$mult[val_i]),
                       ncol = ncol(testing_data$X0), nrow = testing_data$mult[val_i], byrow = TRUE)
    full_X_testing    = rbind(full_X_testing, temp_rows)
  }
}
Z_test = testing_data$Z[which(testing_data$Z < 3)]
colnames(full_X_testing) = names(norm_data$X0)

# Now that I've refit the GP, I should check the final RMSE again:
scores(model=model, Xtest = full_X_testing, Ztest = Z_test, return.rmse = TRUE)$rmse

####################################
## Get the predictions using the model:
n_train     = length(Z_train)
train_preds = predict(object=model, x = full_X, xprime = full_X)
train_true  = as.matrix(Z_train)

n_test      = length(Z_test)
test_preds  = predict(object=model, x = full_X_testing, xprime = full_X_testing)
test_true   = Z_test
var_terms   = test_preds$sd2 + test_preds$nugs

coverage_probs = c()
for(idx in 1:n_test){
  # At each point in the input space, I want to know the probability of the smalleest
  # possible ball around the predicted value that still contains the true value.
  #### Note that this is a BALL -- so we're not working with upper and lower probabilities.
  dist_from_mean = abs(test_preds$mean[idx] - test_true[idx])
  min_cov_prob   = 1 - 2*(1 - pnorm(dist_from_mean, mean = 0, sd = sqrt(var_terms[idx])))
  coverage_probs = c(coverage_probs, min_cov_prob)
}

# Now draw the empirical CDF
graph_data = data.frame(value = coverage_probs)
ggplot(graph_data, aes(x=value)) + stat_ecdf() + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab("Nominal Coverage") + 
  ylab("Empirical Coverage") + 
  theme(axis.title = element_text(size = 21), axis.text = element_text(size = 15))+
  geom_abline(linetype='dashed')







