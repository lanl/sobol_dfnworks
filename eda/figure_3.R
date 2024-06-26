###############################################################################
### Script to look at the training data to establish heteroskedasticity.
## Author: Alexander C. Murph
## Date: September 2023
library(hetGP)
library(ggplot2)
library(latex2exp)
library(gridExtra)
setwd("~/GitLab/sobol_dfnworks/eda")
source("../GP_analysis/GP_fit_analysis_helpers.R")

# Make sure that they following is the same as in the sequential_draws file:
num_of_replications               = 10
num_of_experiments                = 200
total_num_from_initial_experiment = num_of_replications*num_of_experiments
percentile                        = 0.1

## Data locations
test_data_loc = "~/GitHub/sobol_dfnworks/test_data"
models_loc    = "~/GitHub/sobol_dfnworks/sequential_design_10th_percentile"
orig_data_loc = "~/GitHub/sobol_dfnworks/dfnworks_drivers"

## Baseline data
# setwd("~/GitHub/sobol_dfnworks/GP_analysis")
# orig_norm_data           = get_normalized_full_data_woOutliers(num_of_replications,
#                                                                data_loc = orig_data_loc,
#                                                                num_of_experiments, percentile)
# save(orig_norm_data, file = "../GP_analysis/sparse_lhs_normalized_data_woOutliers.Rdata")
load("../GP_analysis/sparse_lhs_normalized_data_woOutliers.Rdata")


########
# Let's simply look at the marginals of the input data against the response to see
# whether or not we have any evidence of heterskedasticity.
graph_data    = orig_norm_data$full_inputs
graph_data$BT = orig_norm_data$full_outputs

# Let's look at the sample (log) variances over subsequent intervals.
increment      = 0.1
lower_interval = -increment
upper_interval = 0
sample_log_var = c()
lower_bounds   = c()
while(upper_interval < 1){
  lower_interval = lower_interval + increment
  upper_interval = upper_interval + increment
  temp_data      = graph_data[which( (graph_data$alpha_TPL>=lower_interval)&(graph_data$alpha_TPL<upper_interval) ), ]
  sample_log_var = c(sample_log_var, log(var(temp_data$BT))) # this is log(var()) in the paper.
  lower_bounds   = c(lower_bounds, lower_interval)
}

graph_list      = list()
graph_data_1    = data.frame(TPL = lower_bounds, varBT = sample_log_var)
graph_data_2    = orig_norm_data$full_inputs
graph_data_2$BT = orig_norm_data$full_outputs

# lower_bounds    = lower_bounds[1:(length(lower_bounds)-1)]
g1 = ggplot(graph_data_1, aes(x = TPL + 0.05, y = varBT)) + geom_point() + 
                      ggtitle(TeX("Log Variance of 10th Percentile of Breakthrough Time over Intervals of $\\gamma$")) + 
                      theme_bw() + 
                      xlab(TeX("Scaled $\\gamma$")) + 
                      ylab(TeX("Sample Log Variance")) + 
                      theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15)) +
                      geom_vline(xintercept = lower_bounds, linetype = 'dashed', alpha = 0.8) +
                      xlim(0,1)

g2 = ggplot(graph_data_2, aes(x = alpha_TPL, y = BT)) + geom_point() + # this is just regular BT in the paper.
                      ggtitle(TeX("Radius TPL exponent ($\\gamma$) against 10th Percentile of Breakthrough Time")) + 
                      theme_bw() + 
                      xlab("") + 
                      theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))+
                      geom_vline(xintercept = lower_bounds, linetype = 'dashed', alpha = 0.8)+
                      xlim(0,1)

gA <- ggplotGrob(g2)
gB <- ggplotGrob(g1)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))



