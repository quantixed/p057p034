# Using the output from CreateSummary.R we will load and crunch the numbers

CSVFolder <- "Output/Data/" 

Summary <- read.csv(paste0(CSVFolder,"output-summary.csv"))
Summary$Average_surface_area <- Summary$Total_length_um * 2 * pi * (Summary$Average_width_um / 2)
Summary$Node_surface_area <- Summary$Total_Nodes * (2 * pi * (Summary$Average_width_um / 2)^2)
Summary$Node_ratio_surface_area <- Summary$Node_surface_area / (Summary$Average_surface_area + Summary$Node_surface_area)
Summary$Free_end_surface_area <- (Summary$Total_Nodes * Summary$Free_Ends) * (2 * pi * (Summary$Average_width_um / 2)^2)
Summary$Free_end_ratio_surface_area <- Summary$Free_end_surface_area / (Summary$Average_surface_area + Summary$Free_end_surface_area)
median(Summary$Free_end_ratio_surface_area)
# frequency of free ends to edges (from counts)
412 / (412 + 155)
# so corrected for s.a. gives
(412 / median(Summary$Free_end_ratio_surface_area)) / ( (412 / median(Summary$Free_end_ratio_surface_area)) + (155 / (1 - median(Summary$Free_end_ratio_surface_area))) )

data <- matrix(c(412,155,283.5,283.5), nrow = 2)
fisher.test(data)

# Fisher's Exact Test for Count Data
# 
# data:  data
# p-value = 4.432e-15
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.060314 3.431244
# sample estimates:
# odds ratio 
#     2.6557 

chisq.test(data,correct=F)

# Pearson's Chi-squared test
# 
# data:  data
# X-squared = 61.398, df = 1, p-value = 4.663e-15

# after correction
# ratio is 0.053 ends to edges surface area
# so our expected value should be 0.053 vs 0.947
# as a proportion of 567 this is 30 vs 537
data <- matrix(c(412,155,30,537), nrow = 2)
fisher.test(data)

# Fisher's Exact Test for Count Data
# 
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  31.14031 74.11907
# sample estimates:
# odds ratio 
#   47.34494 

chisq.test(data, correct = F)
# Pearson's Chi-squared test
# 
# data:  data
# X-squared = 541.02, df = 1, p-value < 2.2e-16