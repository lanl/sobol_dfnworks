###############################################################################
### After fitting the final GP, this script graphs the pivoted cholesky 
### transformation on the prediction errors as a means to assess the model fit.
### I'll likely do this on both the training and testing data.
## Author: Alexander C. Murph
## Date: September 2023
library(Matrix)
library(hetGP)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggExtra) 
setwd("~/GitHub/sobol_dfnworks/GP_analysis")
source("GP_fit_analysis_helpers.R")

num_of_experiments  = 200
num_of_replications = 10
orig_data_loc       = "~/GitHub/sobol_dfnworks/dfnworks_drivers"
test_data_loc       = "~/GitHub/sobol_dfnworks/test_data"
percentile          = 0.1

## Baseline data
# orig_norm_data           = get_normalized_full_data_woOutliers(num_of_replications,
#                                                                data_loc = orig_data_loc,
#                                                                num_of_experiments, percentile)
# save(orig_norm_data, file = "sparse_lhs_normalized_data_woOutliers.Rdata")
load("sparse_lhs_normalized_data_woOutliers.Rdata")
norm_data = orig_norm_data

## Testing Data
# testing_data = get_updated_normalized_test_data_woOutliers(norm_data, test_data_loc,
#                                                            percentile = percentile)
# save(testing_data, file = "compiled_testing_data_woOutliers.Rdata")
load("compiled_testing_data_woOutliers.Rdata")

## Build the GP
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

# Now I should have the model and the data from this point in the sequential design.
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
## Take a look at the outputs for a second...
graph_data = data.frame(bt = c(norm_data$Z, testing_data$Z), data_set = c(rep("train", times = length(norm_data$Z)),
                                                              rep("test", times = length(testing_data$Z))))
ggplot(graph_data, aes(x = bt, fill = data_set)) + geom_density(alpha = 0.5)

summary(Z_train)
summary(Z_test)

# Now that I've fit the GP, I should check the final RMSE:
scores(model=model, Xtest = full_X_testing, Ztest = Z_test, return.rmse = TRUE)$rmse

####################################
## Get the predictions using the model:
n_train     = length(Z_train)
train_preds = predict(object=model, x = full_X, xprime = full_X)
train_true  = as.matrix(Z_train)

n_test      = length(Z_test)
test_preds  = predict(object=model, x = full_X_testing, xprime = full_X_testing)
test_true   = Z_test

####################################
## Get the pivoted cholesky decomposition matrices and rescale prediction errors
cov_train = train_preds$cov + diag(train_preds$nugs)

# In the following, the decomposition on the cov matrix is done in the order
# or greatest to least greatest variance.
chol_train = chol(cov_train, pivot = TRUE)
pivot_train = attr(chol_train, "pivot")

# We get the prediction errors:
pivoted_error_terms_train = (train_true - as.matrix(train_preds$mean))
hist(pivoted_error_terms_train)

# Then order them according to the pivott in our cholesky decomposition:
pivoted_error_terms_train = pivoted_error_terms_train[pivot_train]
pivoted_cholesky_transformed_errors_train = solve(chol_train)%*%pivoted_error_terms_train
hist(pivoted_cholesky_transformed_errors_train)

# Finally, we do the same thing for the testing data:
cov_test = test_preds$cov + diag(test_preds$nugs)
chol_test = chol(cov_test, pivot = TRUE)
pivot_test = attr(chol_test, "pivot")
pivoted_error_terms_test = (test_true - as.matrix(test_preds$mean))
hist(pivoted_error_terms_test)

pivoted_error_terms_test = pivoted_error_terms_test[pivot_test]
pivoted_cholesky_transformed_errors_test = solve(chol_test)%*%pivoted_error_terms_test
hist(pivoted_cholesky_transformed_errors_test)

####################################
## Get the student-t bounds:
p = 5
df_train = n_train - p
upper_t_train = qt(0.975, df_train)
lower_t_train = qt(0.025, df_train)

df_test = n_test - p
upper_t_test = qt(0.975, df_test)
lower_t_test = qt(0.025, df_test)


####################################
## Create data and produce graphs.
lst_p = list()

pivot_data_for_train = data.frame(X = full_X[,'p32'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_train) ))
lst_p[[1]] = ggplot(pivot_data_for_train, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_train, upper_t_train), color='red')+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled p32")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

pivot_data_for_train = data.frame(X = full_X[,'alpha_semi'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_train) ))
lst_p[[2]] = ggplot(pivot_data_for_train, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_train, upper_t_train), color='red')+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\alpha$")) + 
  ylab(TeX("")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

pivot_data_for_train = data.frame(X = full_X[,'beta_semi'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_train) ))
lst_p[[3]] = ggplot(pivot_data_for_train, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_train, upper_t_train), color='red')+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\beta$")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

pivot_data_for_train = data.frame(X = full_X[,'sigma_semi'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_train) ))
lst_p[[4]] = ggplot(pivot_data_for_train, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_train, upper_t_train), color='red')+ 
  ggtitle(TeX("$")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\sigma$")) + 
  ylab(TeX("")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

pivot_data_for_train = data.frame(X = full_X[,'alpha_TPL'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_train) ))
lst_p[[5]] = ggplot(pivot_data_for_train, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_train, upper_t_train), color='red')+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\gamma$")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

# Convert to grobs
lst_p <- lapply(lst_p, ggplotGrob)

# Plot using gridExtra and grid
gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], lst_p[[3]], lst_p[[4]], grid::nullGrob(), lst_p[[5]], grid::nullGrob(),
                        layout_matrix = matrix(c(1,1,2,2,3,3,4,4,5,6,6,7), byrow = TRUE, ncol = 4))



### Now for the testing validation set:
lst_p = list()

pivot_data_for_test = data.frame(X = full_X_testing[,'p32'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_test) ))
lst_p[[1]] = NULL

p1 = ggplot(pivot_data_for_test, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_test, upper_t_test, 0), linetype=c('dashed', 'dashed', 'solid'))+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled p32")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15)) + 
  ylim(c(-3,3))
lst_p[[1]] = ggMarginal(p1, margins = "y", type = "density", size = 7, fill = 'gray', alpha = 0.2)

pivot_data_for_test = data.frame(X = full_X_testing[,'alpha_semi'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_test) ))
p2 = ggplot(pivot_data_for_test, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_test, upper_t_test, 0), linetype=c('dashed', 'dashed', 'solid'))+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\alpha$")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15)) + 
  ylim(c(-3,3))
lst_p[[2]] = ggMarginal(p2, margins = "y", type = "density", size = 7, fill = 'gray', alpha = 0.2)

pivot_data_for_test = data.frame(X = full_X_testing[,'beta_semi'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_test) ))
p3 = ggplot(pivot_data_for_test, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_test, upper_t_test, 0), linetype=c('dashed', 'dashed', 'solid'))+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\beta$")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15)) + 
  ylim(c(-3,3))
lst_p[[3]] = ggMarginal(p3, margins = "y", type = "density", size = 7, fill = 'gray', alpha = 0.2)

pivot_data_for_test = data.frame(X = full_X_testing[,'sigma_semi'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_test) ))
p4 = ggplot(pivot_data_for_test, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_test, upper_t_test, 0), linetype=c('dashed', 'dashed', 'solid'))+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\sigma$")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15)) + 
  ylim(c(-3,3))
lst_p[[4]] = ggMarginal(p4, margins = "y", type = "density", size = 7, fill = 'gray', alpha = 0.2)

pivot_data_for_test = data.frame(X = full_X_testing[,'alpha_TPL'], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_test) ))
p5 = ggplot(pivot_data_for_test, aes(x = X, y = Y)) + geom_point() + geom_hline(yintercept= c(lower_t_test, upper_t_test, 0), linetype=c('dashed', 'dashed', 'solid'))+ 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\gamma$")) + 
  ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15)) + 
  ylim(c(-3,3))
lst_p[[5]] = ggMarginal(p5, margins = "y", type = "density", size = 7, fill = 'gray', alpha = 0.2)

# Plot using gridExtra and grid
grid.arrange(lst_p[[1]], lst_p[[2]], lst_p[[3]], lst_p[[4]], grid::nullGrob(), lst_p[[5]], grid::nullGrob(),
                        layout_matrix = matrix(c(1,1,2,2,3,3,4,4,5,6,6,7), byrow = TRUE, ncol = 4))

save(model, file="hetGP_model_final.Rdata")

