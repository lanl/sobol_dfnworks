# Compare the performance of hetGP against joint GAMs and joint GLMs
## Author: Alexander C. Murph
## Date: June 2024

## Start by loading the data and performing the coverage experiment for the hetGP:
setwd("~/GitHub/sobol_dfnworks/GP_analysis")
source("figure_5.R")

## Now all the necessary data should be in the environment, so let's
setwd("~/GitHub/sobol_dfnworks/model_comparisons")
library(mgcv)
library(ggplot2)

# Get the JointModeling package from the CRAN archive:
url = "https://cran.r-project.org/src/contrib/Archive/JointModeling/JointModeling_1.0-2.tar.gz"
pkgFile = "JointModeling_1.0-2.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type = 'source', repos = NULL)
unlink(pkgFile)
library(JointModeling)


# GAMs example
full_data = data.frame(cbind(Z_train, full_X))
test_data = data.frame(cbind(Z_test, full_X_testing))

ajust <- fitjoint("gam", 'Z_train~alpha_semi+beta_semi+sigma_semi+alpha_TPL+p32', 
                  'd~alpha_semi+beta_semi+sigma_semi+alpha_TPL+p32', data = full_data, reml = FALSE)
# I need to manually grab the two GLMs:
mod.disp = ajust$mod.disp
mod.mean = ajust$mod.mean

n_test      = length(Z_test)
test_true   = Z_test

var_terms   = predict(object=mod.disp, newdata = data.frame(full_X_testing), 
                      type = 'response')
test_preds  = predict(object=mod.mean, data = data.frame(full_X_testing),
                      type = 'response')
coverage_probs = c()
for(idx in 1:n_test){
  # At each point in the input space, I want to know the probability of the smalleest
  # possible ball around the predicted value that still contains the true value.
  #### Note that this is a BALL -- so we're not working with upper and lower probabilities.
  dist_from_mean = abs(test_preds[idx] - test_true[idx])
  min_cov_prob   = 1 - 2*(1 - pnorm(dist_from_mean, mean = 0, sd = sqrt(var_terms)))
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



