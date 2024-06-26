# Verify heteroskedasticity with a linear model
## Author: Alexander C. Murph
## Date: June 2024
library(gridExtra)


## Start by loading the data and performing the coverage experiment for the hetGP:
setwd("~/GitHub/sobol_dfnworks/GP_analysis")
source("figure_5.R")

# GAMs example
full_data = data.frame(cbind(Z_train, full_X))
test_data = data.frame(cbind(Z_test, full_X_testing))

ajust <- lm('Z_train~alpha_semi+beta_semi+sigma_semi+alpha_TPL+p32', data = full_data)
resids = ajust$residuals

residual_data = data.frame(residuals = resids)
residual_data = cbind(residual_data, full_X)

####################################
## Create data and produce graphs.
lst_p = list()

lst_p[[1]] = ggplot(residual_data, aes(x = p32, y = residuals)) + geom_point() +
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled p32")) + 
  ylab(TeX("Residuals")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

lst_p[[2]] = ggplot(residual_data, aes(x = alpha_semi, y = residuals)) + geom_point() + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\alpha$")) + 
  ylab(TeX("")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

lst_p[[3]] = ggplot(residual_data, aes(x = beta_semi, y = residuals)) + geom_point() + 
  ggtitle(TeX("")) + 
  theme_bw() +  
  xlab(TeX("Scaled $\\beta$")) + 
  ylab(TeX("Residuals")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

lst_p[[4]] = ggplot(residual_data, aes(x = sigma_semi, y = residuals)) + geom_point() + 
  ggtitle(TeX("")) + 
  theme_bw() +  
  xlab(TeX("Scaled $\\sigma$")) + 
  ylab(TeX("")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

lst_p[[5]] = ggplot(residual_data, aes(x = alpha_TPL, y = residuals)) + geom_point() + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab(TeX("Scaled $\\gamma$")) + 
  ylab(TeX("Residuals")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

# Convert to grobs
lst_p <- lapply(lst_p, ggplotGrob)

# Plot using gridExtra and grid
gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], lst_p[[3]], lst_p[[4]], grid::nullGrob(), lst_p[[5]], grid::nullGrob(),
                        layout_matrix = matrix(c(1,1,2,2,3,3,4,4,5,6,6,7), byrow = TRUE, ncol = 4))