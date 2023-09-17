## Problem Set 2 ##

rm(list=ls())
library(ellipse)
library(MASS)

#=================================================================================================================================================================================================================================
#Implemented functions:

correlations_grid = function(X){
  
  R = cor(X)
  
  lowers = (dim(R)[1]**2 - dim(R)[1])/2
  
  R[!lower.tri(R)] = NA
  sort(R, decreasing = T)
  order.cor = order(R, decreasing = T)[1:lowers]
  stack.m = matrix(1:dim(R)[1]^2, ncol = dim(R)[1])
  
  Order.cor = matrix(rep(0, lowers*2), ncol = 2)
  for(element in seq(1, lowers)){
    Order.cor[element, ] = which(stack.m == order.cor[element], arr.ind = 1)
  }
  variables = matrix(names(X)[Order.cor], ncol = 2)
  
  values = c()
  for (index in seq(1, dim(variables)[1])){
    values[index] = paste(variables[index, 1], variables[index, 2], sep = " : ")
  }
  
  output = cbind(values,  R[Order.cor])
  colnames(output) = c("Variable's names", "Correlations")
  output = data.frame(output)
  output[, 2] = as.numeric(output[, 2])
  
  return(output)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

order_corrs = function(X, p){
  R = cor(X)
  
  m = matrix(sort(R[p, -p], decreasing = T), ncol = 1)
  rownames(m) = names(sort(R[p, -p], decreasing = T))
  colnames(m) = paste("Ordered Correlations: ", colnames(X)[p])
  
  return(m)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

corrs_report = function(X){
  
  corrs = list()
  for(p in seq(1, dim(X)[2])){
    corrs[[p]] = order_corrs(X, p)
  }
  
  names(corrs) = colnames(X)
  return(corrs)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

trace = function(matrix){
  summation = 0
  for(row in seq(1, dim(matrix)[1])){
    for(column in seq(1, dim(matrix)[2])){
      if(row == column){
        summation = summation + matrix[row, column]
      }
    }
  }
  
  return(summation)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

prop_var = function(L, R){
  variances = c()
  props = c()
  names = c()
  
  for(factor in seq(1, dim(L)[2])){
    variances[factor] = sum(L[, factor]**2)
    props[factor] = sum(L[, factor]**2)/trace(R)
    names[factor] = paste("Factor ", as.character(factor))
  }
  
  cum_var = cumsum(variances)
  cum_props = cumsum(props)
  
  output = cbind(props, cum_props, variances, cum_var)
  colnames(output) = c("Proportion of Variance", "Cumulative Proportion of Variance", "Variance", "Cumulative Variance")
  rownames(output) = names
  
  return(output)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

communalities = function(L, R){
  communalities = c()
  spec_var = c()
  
  for(p in seq(1, dim(L)[1])){
    communalities[p] = sum(L[p, ]**2)
    spec_var[p] = R[p, p] - communalities[p]
  }
  
  output = rbind(communalities, spec_var)
  colnames(output) = rownames(L)
  rownames(output) = c("Estimated Communality", "Specific Variance")
  
  return(output)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

frobenious = function(X){
  return(sqrt(trace(X%*%t(X))))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mean_squared_error_covariances = function(m){
  m[(!lower.tri(m)) & (!upper.tri(m))] = NA
  return(sqrt(apply((m^2), 2, mean, na.rm = T)))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

interpretation_report = function(L, comms){
  loads = c()
  factors = c()
  
  for(p in seq(1, dim(L)[1])){
    if((L[p, 1] < 0) & (L[p, 2] < 0)){
      loads[p] = L[p, which(L[p, ] == min(L[p, ]))]
    }
    else{
      loads[p] = L[p, which(L[p, ] == max(L[p, ]))]
    }
    
    if(abs(loads[p]) > 0.7){
      factors[p] = paste("Factor ", as.character(which(L[p, ] == loads[p])))
    }
    else{
      factors[p] = "Both"
    }
  }
  
  output = data.frame(loads, comms[1, ], factors)
  colnames(output) = c("Highest Loading", "Communality", "Factor")
  rownames(output) = colnames(comms)
  
  return(output)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

regression = function(xi, L, S, x_bar = NULL, scaled = T){
  if(scaled){
    return(t(L)%*%solve(S)%*%xi)
  }
  else{
    return(t(L)%*%solve(S)%*%t(xi - x_bar))
  }
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

self_scores = function(X, L, scaled = T){
  elements = rep(0, dim(X)[1] * dim(L)[2])
  scores = matrix(elements, nrow = dim(X)[1], ncol = dim(L)[2])
  
  if(scaled){
    X_scaled = scale(X)
    R = cor(X)
    for(i in seq(1, dim(X)[1])){
      scores[i, ] = regression(xi = X_scaled[i, ], L = L, S = R)
    }
  }
  else{
    S = cov(X)
    x_bar = colMeans(X)
    for(i in seq(1, dim(X)[1])){
      scores[i, ] = regression(xi = X[i, ], L = L, S = S, x_bar = x_bar, scaled = F)
    }
  }
  
  return(scores)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

self_mahalanobis = function(X, plot = F){
  S = cov(X)
  x.bar = apply(X, 2, mean)
  data = as.matrix(X)
  rownames(data) = seq(1, dim(X)[1])
  
  mahalanobis_squared = c()
  for(i in seq(1, dim(data)[1])){
    mahalanobis_squared[i] = t(data[i, ] - x.bar)%*%solve(S)%*%(data[i, ] - x.bar)
  }
  
  if(plot == T){
    probabilities = ppoints(length(mahalanobis_squared))
    theoretical_quantiles = qchisq(probabilities, df = dim(X)[2])
    sample_quantiles = sort(mahalanobis_squared)
    
    plot(y = sample_quantiles, x = theoretical_quantiles, pch = 16, main = "Q-Q plot for Squared Mahalanobis Distribution", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
    abline(a = 0, b = 1)
  }
  return(mahalanobis_squared)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

normality_criterion = function(components, chosen = c(1, 2), alphas = c(0.5, 0.9, 0.95)){
  data = components[, chosen]
  distances = self_mahalanobis(components[, chosen])
  
  elements = c()
  position = 1
  for(alpha in alphas){
    summation = 0
    bound = dim(data)[1]*(1 - alpha)
    for(distance in distances){
      if(distance > qchisq(alpha, df = length(chosen))){
        summation = summation + 1
      }
    }
    elements[position] = summation
    position = position + 1
  }
  
  percentages = c()
  position = 1
  for(element in elements){
    percentages[position] = paste(toString(round((element/dim(data)[1])*100, 2)), "%", sep = "")
    position = position + 1
  }
  
  expected = c()
  position = 1
  for(alpha in alphas){
    expected[position] = paste(toString(100*(1-alpha)), "%", sep = "")
    position = position + 1
  }
  
  names_ = c()
  position = 1
  for(alpha in alphas){
    names_[position] = paste("Quantile ", toString(alpha), ":", sep = "")
    position = position + 1
  }
  output = data.frame(elements, percentages, expected)
  colnames(output) = c("Number of elements", "Observed Percentage", "Expected Percentage")
  rownames(output) = names_
  
  return(output)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

self_qqplot = function(X){
  for(p in seq(1, dim(X)[2])){
    sample_quantiles = sort(X[, p])
    sample_mean = mean(X[, p])
    sample_sd = sd(X[, p])
    
    probabilities = ppoints(n)
    theoretical_quantiles = qnorm(probabilities)
    
    plot(y = sample_quantiles, x = theoretical_quantiles, pch = 16, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = paste("Q-Q plot for variable: ", colnames(X)[p], sep = ""))
    abline(a = sample_mean, b = sample_sd, col = "red")
    qqline(X[, p])
  }
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pooled = function(X, K, classes, target = "type"){
  p = dim(X)[2] - 1
  sigma_p = matrix(0, nrow = p, ncol = p)
  sigmas = list()
  
  for(index in seq(1, K)){
    sigmas[[index]] = cov(X[which(X[, target] == classes[index]), ][, -(p+1)])*(dim(X[which(X[, target] == classes[index]), ])[1] - 1)
  }
  
  for(k in seq(1, K)){
    sigma_p = sigma_p + sigmas[[k]]
  }
  sigma_p = (1/(dim(X)[1] - K))*sigma_p
  
  return(sigma_p)
}

#==================================================================================================================================================================================================================================

# EXERCISE 1: Pulp paper data

#Load and Display Data:

pulp_paper = read.table("C:/Users/39392/Documents/Luca/Università/Laurea Magistrale Stochastics and Data Science/I° Year/II° Semester/Multivariate Statistical Analysis/Problem Sets/Problem Set 2/data/pulp_paper.txt") 
names(pulp_paper) = c("BL", "EM", "SF", "BS", "AFL", "LFF", "FFF", "ZST")
head(pulp_paper)
dim(pulp_paper)

X = pulp_paper
n = dim(X)[1]
p = dim(X)[2]
S = cov(X)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 1.1) Obtain the maximum likelihood solution for m = 2 and m = 3 common factors on the standardize observations and compute the proportion of total sample variance due to each factor. List the estimated communalities, 
#specific variances, and the residual matrix Compare the results. Which choice of m do you prefer? Why?

#We look at the Correlation Structure:
R = cor(X); R
correlations_grid(X)
report = corrs_report(X); report
eigen(R)$values

#We can also make use of a pairplot to enlight the relationships among variables:
pairs(X, lower.panel = NULL, pch = 16)

#m = 2:
#We perform Factor Analysis with Maximum Likelyhood Estimation on the standardized observations:
X_s = scale(X)
fa_ml = factanal(covmat = R, factors = 2, rotation = "none"); fa_ml
L_z = fa_ml$loadings[,]; L_z

#We compute the proportion of total Sample Variance explained by each of the factor and the cumulative proportion:
props = prop_var(L_z, R); props

#We compute the Estimated Communalities and the Specific Variances for each of the variables:
comms = communalities(L_z, R); comms

#We compute the Residual Matrix and the approximation in gives of "R" (Frobenious Norm):
psi = diag(comms[2, ], ncol = p); psi
L_z%*%t(L_z) + psi; R
residual_matrix = R - (L_z%*%t(L_z) + psi); residual_matrix
frobenious(residual_matrix)
sum(eigen(R)$values[3:dim(R)[1]])

#m = 3:
#We perform Factor Analysis with Maximum Likelyhood Estimation on the standardized observations:
fa_ml_ = factanal(covmat = R, factors = 3, rotation = "none"); fa_ml_
L_z_ = fa_ml_$loadings[,]; L_z_

#We compute the proportion of total Sample Variance explained by each of the factor and the cumulative proportion:
props_ = prop_var(L_z_, R); props_

#We compute the Estimated Communalities and the Specific Variances for each of the variables:
comms_ = communalities(L_z_, R); comms_

#We compute the Residual Matrix and the approximation in gives of "R" (Frobenious Norm):
psi_ = diag(comms_[2, ], ncol = p); psi_
L_z_%*%t(L_z_) + psi_; R
residual_matrix_ = R - (L_z_%*%t(L_z_) + psi_); residual_matrix_
frobenious(residual_matrix_)
sum(eigen(R)$values[3:dim(R)[1]])

#We compare the two methods of extraction by considering the proportion of variance explained and the approximation of the Residual Matrix achieved:
prop = c(props[2, 2], props_[3, 2])
frobs = c(frobenious(residual_matrix), frobenious(residual_matrix_))
comparison = rbind(prop, frobs)
rownames(comparison) = c("Proportions of Variance Explained", "Frobenious Norms")
colnames(comparison) = c("M = 2", "M = 3")
comparison

#Reply steps with "Varimax":

#m = 2:
#We perform Factor Analysis with Maximum Likelyhood Estimation on the standardized observations:
fa_ml_v = factanal(covmat = R, factors = 2); fa_ml_v
L_z_v = fa_ml_v$loadings[,]; L_z_v

#We compute the proportion of total Sample Variance explained by each of the factor and the cumulative proportion:
props_v = prop_var(L_z_v, R); props_v

#We compute the Estimated Communalities and the Specific Variances for each of the variables:
comms_v = communalities(L_z_v, R); comms_v

#We compute the Residual Matrix and the approximation in gives of "R" (Frobenious Norm):
psi_v = diag(comms_v[2, ], ncol = p); psi_v
L_z_v%*%t(L_z_v) + psi_v; R
residual_matrix_v = R - (L_z_v%*%t(L_z_v) + psi_v); residual_matrix_v
frobenious(residual_matrix_v)
sum(eigen(R)$values[3:dim(R)[1]])

#m = 3:
#We perform Factor Analysis with Maximum Likelyhood Estimation on the standardized observations:
fa_ml__v = factanal(covmat = R, factors = 3); fa_ml__v
L_z__v = fa_ml__v$loadings[,]; L_z__v

#We compute the proportion of total Sample Variance explained by each of the factor and the cumulative proportion:
props__v = prop_var(L_z__v, R); props__v

#We compute the Estimated Communalities and the Specific Variances for each of the variables:
comms__v = communalities(L_z__v, R); comms__v

#We compute the Residual Matrix and the approximation in gives of "R" (Frobenious Norm):
psi__v = diag(comms__v[2, ], ncol = p); psi__v
L_z__v%*%t(L_z__v) + psi__v; R
residual_matrix__v = R - (L_z__v%*%t(L_z__v) + psi__v); residual_matrix__v
frobenious(residual_matrix__v)
sum(eigen(R)$values[3:dim(R)[1]])

#We compare the two methods of extraction by considering the proportion of variance explained and the approximation of the Residual Matrix achieved:
prop_v = c(props_v[2, 2], props__v[3, 2])
frobs_v = c(frobenious(residual_matrix_v), frobenious(residual_matrix__v))
comparison_v = rbind(prop_v, frobs_v)
rownames(comparison_v) = c("Proportions of Variance Explained", "Frobenious Norms")
colnames(comparison_v) = c("M = 2", "M = 3")
comparison_v

#We analyze the mean squared error in reproducing covariances with different number of factors:

mean_squared_error_covariances(residual_matrix)
mean_squared_error_covariances(residual_matrix_)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.2) Give an interpretation to the common factors in the m = 2 solution:

#We look at the loading's structure and communalities for m = 2 factors under a "Varimax" rotation:
L_z_v
interpretation_report(L_z_v, comms_v)

#We plot the loadings in the Factor Space:
plot(Factor1 ~ Factor2, data = L_z_v, pch = 16, main = "Projection of the Loadings in the Factor Space", 
     ylab = "Factor 1", xlab = "Factor 2")
abline(h = 0, v = 0, lty = 2)
points(Factor1 ~ Factor2, data = L_z_v, pch = 16, 
       col = c("red", "red", "red", "red", "blue", "blue", "orange", "green"))
text(Factor1 ~ Factor2, data = L_z_v, labels = colnames(X), pos = c(4, 4, 1, 4, 1, 4, 4, 4), 
     col = c("red", "red", "red", "red", "blue", "blue", "orange", "green"))

#We interpret taking care of the correlation structure:
report

#Attempt to change the position of "FFF":
modified = X
modified[, "FFF"] = - modified[, "FFF"]; head(modified)

fff_fa = factanal(covmat = cov(modified), factors = 2); fff_fa
L_fff = fff_fa$loadings[,]; L_fff

plot(Factor1 ~ Factor2, data = L_fff, pch = 16, main = "Projection of the Loadings in the Factor Space", 
     ylab = "Factor 1", xlab = "Factor 2", xlim = c(-0.5, 0.9), ylim = c(-0.2, 0.9))
abline(h = 0, v = 0, lty = 2)
points(Factor1 ~ Factor2, data = L_fff, pch = 16, 
       col = c("red", "red", "red", "red", "blue", "blue", "orange", "green"))
text(Factor1 ~ Factor2, data = L_fff, labels = colnames(X), pos = c(4, 4, 1, 4, 1, 4, 4, 4), 
     col = c("red", "red", "red", "red", "blue", "blue", "orange", "green"))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.3) Make a scatterplot of the factor scores for m = 2 obtained by the regression method. Is their correlation equal to zero? Should we expect so? Comment:

#We compute the Factor Scores with the built-in R method:
scores = factanal(x = X, factors = 2, scores = "regression")$scores; head(scores)

#We try to obtain the same result making personally the computations:
#Recall that we're working on scaled data, so mean is zero and data are divided for their variances:
head(self_scores(X, L_z_v))

#We plot the Factor Scores in order to check uncorrelation between factors:
plot(x = scores[, 1], y = scores[, 2], main = "Scatterplot of Factor 1 vs Factor 2", 
     xlab = "Factor 1", ylab = "Factor 2", pch = 16)
abline(h = 0, v = 0)
abline(a = 0, b = 1, col = "red")
legend("bottomright", fill = "red", legend = "Bisector")

cor(scores[, 1], scores[, 2]) #Compute Correlation

#We Check if points forms an elliptical shape for Bivariate Normality:
plot(x = scores[, 1], y = scores[, 2], main = "Scatterplot of Factor 1 vs Factor 2", 
     xlab = "Factor 1", ylab = "Factor 2", pch = 16, xlim = c(-2.4, 2.4))

f_cor = cor(scores); f_cor
f_cov = cov(scores); f_cov
f_bar = colMeans(scores); f_bar
lines(ellipse(x = f_cov, centre = f_bar, level = 0.5), col = "red", lwd = 3)
lines(ellipse(x = f_cov, centre = f_bar, level = 0.9), col = "red", lwd = 3)
lines(ellipse(x = f_cov, centre = f_bar, level = 0.95), col = "red", lwd = 3)

normality_criterion(scores, alphas = c(0.25, 0.5, 0.75, 0.9, 0.95)) #We also use the numerical criterion for Multivariate Normality

#we can also have a look to Univariate normality:
par(mfrow = c(2,1))

#Factor 1:
plot(density(scores[, 1]), lwd = 3, col = "red", main =  "Histogram for Factor 1", xlab = "Factor 1", ylim = c(0, 0.4))
color = rgb(red = 0, green = 0, blue = 0.5, alpha = 0.05)
hist(scores[, 1], probability = T, add = T, col = color)

sample_quantiles = sort(scores[, 1])
sample_mean = mean(scores[, 1])
sample_sd = sd(scores[, 1])

probabilities = ppoints(dim(scores)[1])
theoretical_quantiles = qnorm(probabilities)

plot(y = sample_quantiles, x = theoretical_quantiles, pch = 16, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = "Q-Q plot Factor 1")
abline(a = sample_mean, b = sample_sd, col = "red")
qqline(scores[, 1])

#Factor 2:
par(mfrow = c(2,1))

plot(density(scores[, 2]), lwd = 3, col = "red", main =  "Histogram for Factor 2", xlab = "Factor 2")
color = rgb(red = 0, green = 0, blue = 0.5, alpha = 0.05)
hist(scores[, 2], probability = T, add = T, col = color)

sample_quantiles = sort(scores[, 2])
sample_mean = mean(scores[, 2])
sample_sd = sd(scores[, 2])

probabilities = ppoints(dim(scores)[1])
theoretical_quantiles = qnorm(probabilities)

plot(y = sample_quantiles, x = theoretical_quantiles, pch = 16, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = "Q-Q plot Factor 2")
abline(a = sample_mean, b = sample_sd, col = "red")
qqline(scores[, 2])

#We have a look at the densities of the original variables:
par(mfrow = c(1,1))

color = rgb(red = 0, green = 0, blue = 0.5, alpha = 0.05)
for(p in seq(1, dim(X)[2])){
  plot(density(X[, p]), lwd = 3, col = "red", main = paste("Histogram for the Variable ", colnames(X)[p], sep = ""), xlab = paste("Variable ", colnames(X)[p], sep = ""))
  hist(X[, p], probability = T, add = T, col = color)
}

#Finally we perform the Shapiro-Wilks test:
shapiro.test(scores[, 1])
shapiro.test(scores[, 2])

#We check the other assumptions of the model about the Factors:
f_bar
f_cov

#We try to obtain the Error Terms variables and we check if assumptions are satisfied: THERE'RE PROBLEMS SINCE THE VARIANCES OF THE ERROR TERM VARIABLES DOESN'T MATCH WITH THE ESTIMATED SPECIFIC VARIANCES
epsilon = X_s - t(L_z_v%*%t(scores)); head(epsilon)
sum((X_s - (t(L_z_v%*%t(scores)) + epsilon))) #We check that we exactly reproduced the original Data Matrix
epsilon_bar = colMeans(epsilon); epsilon_bar
e_cov = cov(epsilon); e_cov
diag(e_cov)

#We summarize the assumptions for factors:
summary = cbind(f_bar, diag(f_cov), rep(f_cor[1, 2], 2), 
                c(shapiro.test(scores[, 1])$statistic, shapiro.test(scores[, 2])$statistic),
                c(shapiro.test(scores[, 1])$p.value, shapiro.test(scores[, 2])$p.value))
rownames(summary) = c("Factor 1", "Factor 2")
colnames(summary) = c("Expectations", "Variances", "Correlation", "Shapiro-Wilks statistics", "Shapiro-Wilks p-value")
summary

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.4) Suppose we have a new observation (15.5, 5.5, 2, -0.55, 0.6, 65, -5, 1.2). Calculate the corresponding m = 2 factor scores and add this bivariate point to the plot in 1.3). 
#How is it placed compared to the rest of the n = 62 points? Could you tell without computing the factor scores? Comment.

#We add the new observation to the original dataframe and we rescale data:
new = c(15.5, 5.5, 2, -0.55, 0.6, 65, -5, 1.2); new
X_n = X
X_n[63, ] = new; tail(X_n)
X_n_s = scale(X_n); tail(X_n_s)

#We estimate the new Factor Model:
f_n = factanal(x = X_n_s, factors = 2); f_n
L_n = f_n$loadings[,]; L_n

#We compute the Factor Scores and plot the results in the Factor space:
scores_n = factanal(x = X_n_s, factors = 2, scores = "regression")$scores; head(scores_n)

plot(x = scores_n[, 1], y = scores_n[, 2], main = "Scatterplot of New Factor 1 vs New Factor 2", 
     xlab = "New Factor 1", ylab = "New Factor 2", pch = 16)
abline(h = 0, v = 0, lty = 2)
abline(a = 0, b = 1, lty = 2, col = "red")
points(x = scores_n[63, 1], y = scores_n[63, 2], pch = 16, col = "red")
text(x = scores_n[63, 1], y = scores_n[63, 2], labels = "\n      New\nObservation", col = "red", pos = 4)

#We see how the new observation produces changes:

#How the point appears in the ordered distributions of the variables:
#Variables of the First Factor:
par(mfrow = c(2,4))

boxplot(X[, 1], main = colnames(X)[1], ylim = c(new[1], 26))
points(x = 1, y = X_n[63, 1], pch = 16, col = "red")
text(x = 1, y = X_n[63, 1], labels = "New Observation", col = "red", pos = 1)

boxplot(X[, 2], main = colnames(X)[2], ylim = c(new[2], 8.5))
points(x = 1, y = X_n[63, 2], pch = 16, col = "red")
text(x = 1, y = X_n[63, 2], labels = "New Observation", col = "red", pos = 1)

boxplot(X[, 3], main = colnames(X)[3], ylim = c(new[3], 8))
points(x = 1, y = X_n[63, 3], pch = 16, col = "red")
text(x = 1, y = X_n[63, 3], labels = "New Observation", col = "red", pos = 1)

boxplot(X_n[, 4], main = colnames(X)[4])
points(x = 1, y = X_n[63, 4], pch = 16, col = "red")
text(x = 1, y = X_n[63, 4], labels = "New Observation", col = "red", pos = 1)

#Variables AFL:
boxplot(X[, 5], main = colnames(X)[5], ylim = c(-0.7, new[5]))
points(x = 1, y = X_n[63, 5], pch = 16, col = "red")
text(x = 1, y = X_n[63, 5], labels = "New Observation", col = "red", pos = 4)

#Variables LFF and FFF:
boxplot(X[, 6], main = colnames(X)[6], ylim = c(0, new[6]))
points(x = 1, y = X_n[63, 6], pch = 16, col = "red")
text(x = 1, y = X_n[63, 7], labels = "New Observation", col = "red", pos = c(3, 1))

boxplot(X[, 7], main = colnames(X)[7], ylim = c(new[7], 85))
points(x = 1, y = X_n[63, 7], pch = 16, col = "red")
text(x = 1, y = X_n[63, 7], labels = "New Observation", col = "red", pos = c(3, 1))

#Variables ZST:
boxplot(X[, 8], main = colnames(X)[8], ylim = c(1, new[8]))
points(x = 1, y = X_n[63, 8], pch = 16, col = "red")
text(x = 1, y = X_n[63, 8], labels = "New Observation", col = "red", pos = 4)

par(mfrow = c(1,1))

which(X_n[,1] == min(X_n[,1]))
which(X_n[,2] == min(X_n[,2]))
which(X_n[,3] == min(X_n[,3]))
which(X_n[,4] == min(X_n[,4]))
which(X_n[,5] == max(X_n[,5]))
which(X_n[,6] == max(X_n[,6]))
which(X_n[,7] == min(X_n[,7]))
which(X_n[,8] == max(X_n[,8]))

#Mean Variations:
variations_m = round(((colMeans(X_n) - colMeans(X))/colMeans(X))*100, 1)
variations_m[5] = -variations_m[5]
for(position in seq(1, length(variations_m))){
  variations_m[position] = paste(as.character(variations_m[position]), " %")
}
mean_variations = data.frame(colMeans(X), colMeans(X_n), variations_m)
rownames(mean_variations) = colnames(X)
colnames(mean_variations) = c("Original Means", "New Means", "Percentage of Variation (absolute values)")
mean_variations

#Correlation Variations for variable "AFL":
variations_r = round(((corrs_report(X_n)$AFL[,] - report$AFL[,])/report$AFL[,])*100, 3)
for(position in seq(1, length(variations_r))){
  variations_r[position] = paste(as.character(variations_r[position]), " %")
}
corrs_variations = data.frame(report$AFL[,], corrs_report(X_n)$AFL[,], variations_r)
colnames(corrs_variations) = c("Original Correlations", "New Correlations", "Percentage of Variation (absolute values)")
corrs_variations

#Correlation Variations for variable "ZST":
variations_r_zst = round(((corrs_report(X_n)$ZST[,] - report$ZST[,])/report$ZST[,])*100, 3)
for(position in seq(1, length(variations_r_zst))){
  variations_r_zst[position] = paste(as.character(variations_r_zst[position]), " %")
}
corrs_variations_zst = data.frame(report$ZST[,], corrs_report(X_n)$ZST[,][c(3, 5, 4, 6, 2, 1, 7)], variations_r_zst[c(3, 5, 4, 6, 2, 1, 7)])
colnames(corrs_variations_zst) = c("Original Correlations", "New Correlations", "Percentage of Variation (absolute values)")
corrs_variations_zst

#Correlation Variations for variable "FFF":
variations_r_fff = round(((corrs_report(X_n)$FFF[,] - report$FFF[,])/report$FFF[,])*100, 3)
for(position in seq(1, length(variations_r_fff))){
  variations_r_fff[position] = paste(as.character(variations_r_fff[position]), " %")
}
corrs_variations_fff = data.frame(report$FFF[,], corrs_report(X_n)$FFF[,], variations_r_fff)
colnames(corrs_variations_fff) = c("Original Correlations", "New Correlations", "Percentage of Variation (absolute values)")
corrs_variations_fff

#Correlation Variations for variable "LFF":
variations_r_lff = round(((corrs_report(X_n)$LFF[,] - report$LFF[,])/report$LFF[,])*100, 3)
for(position in seq(1, length(variations_r_lff))){
  variations_r_lff[position] = paste(as.character(variations_r_lff[position]), " %")
}
corrs_variations_lff = data.frame(report$LFF[,], corrs_report(X_n)$LFF[,], variations_r_lff)
colnames(corrs_variations_lff) = c("Original Correlations", "New Correlations", "Percentage of Variation (absolute values)")
corrs_variations_lff[2, 2] = 0.6872388; corrs_variations_lff[2, 3] = "-13.307  %"
corrs_variations_lff[3, 2] = 0.7790708; corrs_variations_lff[3, 3] = "-2.158  %"
corrs_variations_lff

summary_corr = list(corrs_variations, corrs_variations_zst, corrs_variations_fff, corrs_variations_lff); names(summary_corr) = c("AFL", "ZST", "FFF", "LFF")
summary_corr

#Loadings changes:
variations_l1 = round(((L_n[, 1] - L_z_v[, 1])/L_z_v[, 1])*100, 3)
for(position in seq(1, length(variations_l1))){
  variations_l1[position] = paste(as.character(variations_l1[position]), " %")
}
loads1_variations = data.frame(L_z_v[, 1], L_n[, 1], variations_l1)
colnames(loads1_variations) = c("Original Loadings", "New Loadings", "Percentage of Variation")
loads1_variations

variations_l2 = round(((L_n[, 2] - L_z_v[, 2])/L_z_v[, 2])*100, 3)
for(position in seq(1, length(variations_l2))){
  variations_l2[position] = paste(as.character(variations_l2[position]), " %")
}
loads2_variations = data.frame(L_z_v[, 2], L_n[, 2], variations_l2)
colnames(loads2_variations) = c("Original Loadings", "New Loadings", "Percentage of Variation")
loads2_variations

plot(Factor1 ~ Factor2, data = L_n, pch = 16, main = "Projection of the New Loadings in the Factor Space", 
     ylab = "New Factor 1", xlab = "New Factor 2")
abline(h = 0, v = 0, lty = 2)
points(Factor1 ~ Factor2, data = L_n, pch = 16, 
       col = c("red", "red", "red", "red", "blue", "blue", "orange", "blue"))
text(Factor1 ~ Factor2, data = L_n, labels = colnames(X), pos = c(4, 4, 1, 4, 1, 4, 4, 4), 
     col = c("red", "red", "red", "red", "blue", "blue", "orange", "blue"))

#Add the new observation to the plot:
n_s = t(L_z_v)%*%solve(cor(X))%*%t((X_n[63, ] - colMeans(X))/(diag(cov(X))^(1/2))); n_s
plot(x = scores[, 1], y = scores[, 2], main = "Scatterplot of Factor 1 vs Factor 2", 
     xlab = "Factor 1", ylab = "Factor 2", pch = 16, xlim = c(-4.5, 2), ylim = c(-2.5, 4.5))
points(x = n_s[1,], y = n_s[2,], col = "red", pch = 16)
abline(h = 0, v = 0)
text(x = n_s[1,], y = n_s[2,], col = "red", labels = "New Observation", pos = 4)

#==================================================================================================================================================================================================================================

# EXERCISE 2: Glass Data

#Load and Display Data:

glass = read.table("C:/Users/39392/Documents/Luca/Università/Laurea Magistrale Stochastics and Data Science/I° Year/II° Semester/Multivariate Statistical Analysis/Problem Sets/Problem Set 2/data/glass.txt", header = T)
glass$type = factor(glass$type)
levels(glass$type) = c("WinF","WinNF","Veh","Con","Tabl","Head")
head(glass)
dim(glass)

K = length(unique(glass[, "type"])); K
classes = c("WinF","WinNF","Veh","Con","Tabl","Head"); classes

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2.1) Use linear discriminant analysis to predict the glass type. Look at the first two discriminant directions: what are the most important variables in separating the classes? Comment:

#Check if Hypothesis of LDA are confirmed:

#We look at differences in means, variances and correlations:
colMeans(glass[, -dim(glass)[2]])
diag(cov(glass[, -dim(glass)[2]]))
correlations_grid(glass[, -dim(glass)[2]])


#We compute Covariance Matrices, Centroids and Prior Probabilities for different classes of the response variable:
subsets = list()
centroids = matrix(rep(0, (dim(glass)[2] - 1)*length(classes)), nrow = 6)
sigmas = list()
probs = c()

for(index in seq(1, length(classes))){
  subsets[[index]] = glass[which(glass[, "type"] == classes[index]), ]
  centroids[index, ] = colMeans(subsets[[index]][, -dim(glass)[2]])
  sigmas[[index]] = cov(subsets[[index]][, -dim(glass)[2]])
  probs[index] = dim(subsets[[index]])[1]/dim(glass)[1]
}

names(subsets) = classes
names(sigmas) = classes
rownames(centroids) = classes
colnames(centroids) = colnames(glass)[-dim(glass)[2]]
names(probs) = classes

subsets
sigmas
probs
centroids

#We evaluate how much the Covariance Matrices differs from the Pooled Sample Covariance Matrix:
sigma_p = pooled(glass, K, classes); sigma_p

differences = c()
for(k in seq(1, K)){
  differences[k] = frobenious(sigma_p - sigmas[[k]])^2
}
names(differences) = names(sigmas)
differences

#We evaluate Multivariate Normality for the classes:
mahalanobis_squared = list()
for(k in seq(1, K)){
  if(k != 5){
    mahalanobis_squared[[k]] = self_mahalanobis(subsets[[k]][, -dim(glass)[2]])
  }
  else{
    mahalanobis_squared[[k]] = self_mahalanobis(subsets[[5]][,-c(6, 8, 9, dim(glass)[2])])
  }
}

par(mfrow = c(3, 2))

alphas = c(0.25, 0.5, 0.75, 0.95)

plot(x = seq(1, dim(subsets[[1]])[1]), y = mahalanobis_squared[[1]], pch = 16, main = classes[1], xlab = "Observations", ylab = "Squared Mahalanobis distances")
abline(h = qchisq(0.25, df = dim(glass)[2] - 1), col = "green", lty = "dashed")
abline(h = qchisq(0.5, df = dim(glass)[2] - 1), col = "red", lty = "dashed")
abline(h = qchisq(0.75, df = dim(glass)[2] - 1), col = "blue", lty = "dashed")
abline(h = qchisq(0.95, df = dim(glass)[2] - 1), col = "black", lty = "dashed")
legend("topleft", fill = c("green", "red", "blue", "black"), legend = c("Quantile: 0.25", "Quantile: 0.50", "Quantile: 0.75", "Quantile: 0.95"))
normality_criterion(subsets[[1]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)

plot(x = seq(1, dim(subsets[[2]])[1]), y = mahalanobis_squared[[2]], pch = 16, main = classes[2], xlab = "Observations", ylab = "Squared Mahalanobis distances")
abline(h = qchisq(0.25, df = dim(glass)[2] - 1), col = "green", lty = "dashed")
abline(h = qchisq(0.5, df = dim(glass)[2] - 1), col = "red", lty = "dashed")
abline(h = qchisq(0.75, df = dim(glass)[2] - 1), col = "blue", lty = "dashed")
abline(h = qchisq(0.95, df = dim(glass)[2] - 1), col = "black", lty = "dashed")
normality_criterion(subsets[[2]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)

plot(x = seq(1, dim(subsets[[3]])[1]), y = mahalanobis_squared[[3]], pch = 16, main = classes[3], xlab = "Observations", ylab = "Squared Mahalanobis distances")
abline(h = qchisq(0.25, df = dim(glass)[2] - 1), col = "green", lty = "dashed")
abline(h = qchisq(0.5, df = dim(glass)[2] - 1), col = "red", lty = "dashed")
abline(h = qchisq(0.75, df = dim(glass)[2] - 1), col = "blue", lty = "dashed")
abline(h = qchisq(0.95, df = dim(glass)[2] - 1), col = "black", lty = "dashed")
normality_criterion(subsets[[3]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)

plot(x = seq(1, dim(subsets[[4]])[1]), y = mahalanobis_squared[[4]], pch = 16, main = classes[4], xlab = "Observations", ylab = "Squared Mahalanobis distances")
abline(h = qchisq(0.25, df = dim(glass)[2] - 1), col = "green", lty = "dashed")
abline(h = qchisq(0.5, df = dim(glass)[2] - 1), col = "red", lty = "dashed")
abline(h = qchisq(0.75, df = dim(glass)[2] - 1), col = "blue", lty = "dashed")
abline(h = qchisq(0.95, df = dim(glass)[2] - 1), col = "black", lty = "dashed")
normality_criterion(subsets[[4]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)

plot(x = seq(1, dim(subsets[[5]])[1]), y = mahalanobis_squared[[5]], pch = 16, main = classes[5], xlab = "Observations", ylab = "Squared Mahalanobis distances")
abline(h = qchisq(0.25, df = dim(subsets[[5]])[2] - 4), col = "green", lty = "dashed")
abline(h = qchisq(0.5, df = dim(subsets[[5]])[2] - 4), col = "red", lty = "dashed")
abline(h = qchisq(0.75, df = dim(subsets[[5]])[2] - 4), col = "blue", lty = "dashed")
abline(h = qchisq(0.95, df = dim(subsets[[5]])[2] - 4), col = "black", lty = "dashed")
normality_criterion(subsets[[5]][, -c(6, 8, 9)], chosen = seq(1, dim(subsets[[5]][, -c(6, 8, 9)])[2] -1), alphas = alphas)

plot(x = seq(1, dim(subsets[[6]])[1]), y = mahalanobis_squared[[6]], pch = 16, main = classes[6], xlab = "Observations", ylab = "Squared Mahalanobis distances")
abline(h = qchisq(0.25, df = dim(glass)[2] - 1), col = "green", lty = "dashed")
abline(h = qchisq(0.5, df = dim(glass)[2] - 1), col = "red", lty = "dashed")
abline(h = qchisq(0.75, df = dim(glass)[2] - 1), col = "blue", lty = "dashed")
abline(h = qchisq(0.95, df = dim(glass)[2] - 1), col = "black", lty = "dashed")
normality_criterion(subsets[[6]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)

summary_criterions = data.frame(paste(as.character((rep(1, length(alphas)) - alphas)*100), "%"), normality_criterion(subsets[[1]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)[, 2],
                           normality_criterion(subsets[[2]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)[, 2], normality_criterion(subsets[[3]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)[, 2],
                           normality_criterion(subsets[[4]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)[, 2], normality_criterion(subsets[[5]][, -c(6, 8, 9)], chosen = seq(1, dim(subsets[[5]][, -c(6, 8, 9)])[2] -1), alphas = alphas)[, 2],
                           normality_criterion(subsets[[6]], chosen = seq(1, dim(glass)[2] -1), alphas = alphas)[, 2])
names = c()
for(position in seq(1, length(alphas))){
  names[position] = paste("Quantile: ", as.character(alphas[position]))
}
rownames(summary_criterions) = names
colnames(summary_criterions) = c("Expected", classes[1], classes[2], classes[3], classes[4], classes[5], classes[6])
summary_criterions

par(mfrow = c(1, 1))

#We produce pair scatterplots coloring observations according to their classes and marking centroids:
lookup = c("green",  "black", "yellow", "blue", "red", "cyan")
names(lookup) = classes; lookup
glass.col = lookup[glass$type]; head(glass.col)

pairs(glass[, -dim(glass)[2]], pch = 16, main = "Scatterplots", col = glass.col, lower.panel = NULL)
legend(x = "left", legend = classes, fill = lookup)

#We fit LDA Model in order to comment the first 2 discriminant directions:
lda.fit = lda(type ~ ., data = glass); lda.fit
discriminant_directions = lda.fit$scaling[,]; discriminant_directions

#We scale the first two discriminant coordinates:
scaled_ld = cbind(diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 1], diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 2])
rownames(scaled_ld) = colnames(glass)[-10]; colnames(scaled_ld) = c("Scaled LD1", "Scaled LD2")
scaled_ld

sort(scaled_ld[, 1], decreasing = T)
sort(scaled_ld[, 2], decreasing = T)

#We try to sphere data and estimates once again the discriminant directions:
lambda = diag(eigen(sigma_p)$values); lambda
P = eigen(sigma_p)$vectors; P
sphered = data.frame(t(diag(diag(lambda)^(-1/2))%*%t(P)%*%t(glass[, -dim(glass)[2]])))
sphered$type = glass[, "type"]
colnames(sphered) = colnames(glass); head(sphered)

lda.sphere = lda(type ~ ., data = sphered); lda.sphere
discriminant_directions_sp = lda.sphere$scaling[,]; discriminant_directions_sp

#We plot Data in the 2-D Reduced space determined by the first two discriminant directions in both the three cases of coefficients:

#Normal:
Z1 = t(discriminant_directions[, 1]%*%t(glass[, -dim(glass)[2]])); head(Z1)
Z2 = t(discriminant_directions[, 2]%*%t(glass[, -dim(glass)[2]])); head(Z2)

projected = matrix(rep(0, K*2), nrow = K)
for(row in seq(1, dim(centroids)[1])){
  projected[row, 1] = discriminant_directions[, 1]%*%centroids[row, ]
  projected[row, 2] = discriminant_directions[, 2]%*%centroids[row, ]
}
colnames(projected) = c("Z1", "Z2")
rownames(projected) = rownames(centroids)

plot(x = Z1, y = Z2, pch = 1, col = glass.col, main = "Data in the Reduced Subaspace", xlab = "LD1", ylab = "LD2")
for(row in seq(1, dim(centroids)[1])){
  points(x = projected[row, 1], y = projected[row, 2], 
         pch = 21, col = "black", bg = lookup[row], 
         cex = 2, lwd = 2)
}

#Standardized:
Z1_ = t(diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 1])%*%t(glass[, -dim(glass)[2]])
Z2_ = t(diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 2])%*%t(glass[, -dim(glass)[2]])

projected_ = matrix(rep(0, K*2), nrow = K)
for(row in seq(1, dim(centroids)[1])){
  projected_[row, 1] = t(diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 1])%*%centroids[row, ]
  projected_[row, 2] = t(diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 2])%*%centroids[row, ]
}
colnames(projected_) = c("Z1", "Z2")
rownames(projected_) = rownames(centroids)

plot(x = Z1_, y = Z2_, pch = 1, col = glass.col, main = "Data in the Reduced Subaspace", xlab = "LD1 scaled", ylab = "LD2 scaled")
for(row in seq(1, dim(centroids)[1])){
  points(x = projected_[row, 1], y = projected_[row, 2], 
         pch = 21, col = "black", bg = lookup[row], 
         cex = 2, lwd = 2)
}

#Sphered:
Z1_s = t(discriminant_directions_sp[, 1]%*%t(glass[, -dim(glass)[2]]))
Z2_s = t(discriminant_directions_sp[, 2]%*%t(glass[, -dim(glass)[2]]))

projected_s = matrix(rep(0, K*2), nrow = K)
for(row in seq(1, dim(centroids)[1])){
  projected_s[row, 1] = discriminant_directions_sp[, 1]%*%centroids[row, ]
  projected_s[row, 2] = discriminant_directions_sp[, 2]%*%centroids[row, ]
}
colnames(projected_s) = c("Z1", "Z2")
rownames(projected_s) = rownames(centroids)

plot(x = Z1_s, y = Z2_s, pch = 1, col = glass.col, main = "Data in the Reduced Subaspace", xlab = "LD1 sphered", ylab = "LD2 sphered")
for(row in seq(1, dim(centroids)[1])){
  points(x = projected_s[row, 1], y = projected_s[row, 2], 
         pch = 21, col = "black", bg = lookup[row], 
         cex = 2, lwd = 2)
}

#We reproduce the same plot as before in the scaled space and we draw the decision boundaries determined for the 2 Discriminant Variables:
lda.predict = predict(lda.fit)
means.hat = aggregate(lda.predict$x, by = list(glass$type), FUN = mean)
means.hat = means.hat[, -1]
rownames(means.hat) = rownames(centroids)
means.hat

L1 = as.matrix(lda.predict$x[, 1])
L2 = as.matrix(lda.predict$x[, 2])

plot(x = L1, y = L2, pch = 1, col = glass.col, main = "Data in the Reduced Subaspace (Scaled)", xlab = "LD1", ylab = "LD2")
for(row in seq(1, dim(means.hat)[1])){
  points(x = means.hat[row, 1], y = means.hat[row, 2], 
         pch = 21, col = "black", bg = lookup[row], 
         cex = 2, lwd = 2)
}

#Decision Boundaries:
len1<-80; len2<-100
delta<-0.2
grid.X1<-seq(from=min(lda.predict$x[,1])-delta,to=max(lda.predict$x[,1])+delta,length=len1)
grid.X2<-seq(from=min(lda.predict$x[,2])-delta,to=max(lda.predict$x[,2])+delta,length=len2)
dataT<-expand.grid(x.1=grid.X1,x.2=grid.X2)

lda.class<-rep(NA,length=len1*len2)
means.hat<-aggregate(lda.predict$x,by=list(glass$type),FUN=mean)
means.hat<-means.hat[,-1]

m<-as.matrix(means.hat[,1:2])
ones<-rep(1,11)
I.mat<-matrix(rep(c(1,0,0,1),6),byrow=T,ncol=2)

system.time(
  for(i in 1:(len1*len2) ){
    x<-matrix(I.mat%*%t(dataT[i,]),byrow=T,ncol=2)-m
    #x<-outer(ones,as.numeric(dataT[i,]))-m
    x<-diag(crossprod(t(x)))
    lda.class[i]<-order(x)[1]
    if ((i%%1000)==0) cat("iteration ",i," of ",len1*len2,"\n")
  }
)
lda.col<-lookup[lda.class]

Z<-class.ind(lda.class)
np<-len1
for(i in 1:length(lookup)){
  zp <- Z[, i] - apply(Z[,-i],1,max)
  contour(grid.X1, grid.X2, matrix(zp, np),
          add = T, levels = 0, labcex=0.1,lwd=3)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2.2) Compute the training error. Are there any groups less homogeneous than the others? Comment:

#Compute the predictions performed by the model and the Confusion Matrix and we divide the error rates with respect to the different classes of the response variable:
lda.predict = predict(lda.fit)
posterior = lda.predict$posterior; head(posterior)
pred_classes = lda.predict$class; head(pred_classes)

round(posterior[c(49, 64, 32, 194, 92), ], 3)
pred_classes[c(49, 64, 32, 194, 92)]

confusion_matrix = table(predicted = pred_classes, true = glass[,10]); confusion_matrix
training_rate = (sum(confusion_matrix) - trace(confusion_matrix))/214; training_rate

#Error rates labelled by classes:
rates = c()
for(class in seq(1, K)){
  rates[class] = (sum(confusion_matrix[, class]) - diag(confusion_matrix)[class])/sum(confusion_matrix[, class])
}
training_rates = matrix(rates); colnames(training_rates) = "Training Error Rates"; rownames(training_rates) = classes; training_rates

#Variances and generalized sample variances of the original predictors labelled by the classes:
variances = matrix(rep(0, K*(dim(glass)[2] - 1)), nrow = K)
colnames(variances) = colnames(glass)[-10]; rownames(variances) = classes
for(class in seq(1, K)){
  variances[class, ] = diag(sigmas[[class]])
}
variances

generalized = c()
for(class in seq(1, K)){
  generalized[class] = det(sigmas[[class]])
}
names(generalized) = classes; generalized

#We do the same but considering means:
centroids

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2.3) Implement a 10-fold cross validation using the partition of the observations provided by the variable groupCV to estimate the error rate. Comment:

#Import the subset's partition for the Cross Validation:
groupCV<-scan(file = "C:/Users/39392/Documents/Luca/Università/Laurea Magistrale Stochastics and Data Science/I° Year/II° Semester/Multivariate Statistical Analysis/Problem Sets/Problem Set 2/data/groupCV.txt")
glass$group = groupCV; head(glass)

#Implement Cross Validation and extract the Error Rates:
cv_error_rates = c()
train_error_rates = c()

for(fold in 1:10){
  test_set = glass[which(glass[, 11] == fold), -11]
  train_set = glass[which(glass[, 11] != fold), -11]
  
  lda = lda(type ~ ., data = train_set)
  
  cv_predictions = predict(lda, newdata = test_set)$class
  cv_conf_matrix = table(predicted = cv_predictions, true = test_set[, 10])
  cv_error_rates[fold] = (sum(cv_conf_matrix) - sum(diag(cv_conf_matrix)))/sum(cv_conf_matrix)
  
  train_predictions = predict(lda, newdata = train_set)$class
  train_conf_matrix = table(predicted = train_predictions, true = train_set[, 10])
  train_error_rates[fold] = (sum(train_conf_matrix) - sum(diag(train_conf_matrix)))/sum(train_conf_matrix)
}

names(cv_error_rates) = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5", "Fold 6", "Fold 7", "Fold 8", "Fold 9", "Fold 10")
names(train_error_rates) = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5", "Fold 6", "Fold 7", "Fold 8", "Fold 9", "Fold 10")
overall_rates = c(mean(cv_error_rates), mean(train_error_rates))
summary_cv = rbind(cv_error_rates, train_error_rates); rownames(summary_cv) = c("Test Error Rates", "Training Error Rates")
summary_cv = cbind(summary_cv, overall_rates, c(max(cv_error_rates), max(train_error_rates)), c(min(cv_error_rates), min(train_error_rates))); colnames(summary_cv)[11:13] = c("Average", "Worst", "Best")
t(summary_cv)

#Additional details on classes's error rates and classes prior probabilities (to understand why in some folds we find higher values for the general error rates):
cv_class_error = matrix(0, nrow = K, ncol = 10)
train_class_error = matrix(0, nrow = K, ncol = 10)
priors_tr = matrix(0, nrow = K, ncol = 10)
priors_te = matrix(0, nrow = K, ncol = 10)

for(fold in 1:10){
  test_set = glass[which(glass[, 11] == fold), -11]
  train_set = glass[which(glass[, 11] != fold), -11]
  
  lda = lda(type ~ ., data = train_set)
  
  cv_predictions = predict(lda, newdata = test_set)$class
  cv_conf_matrix = table(predicted = cv_predictions, true = test_set[, 10])
  train_predictions = predict(lda, newdata = train_set)$class
  train_conf_matrix = table(predicted = train_predictions, true = train_set[, 10])

  for(class in seq(1, K)){
    cv_class_error[class, fold] = (sum(cv_conf_matrix[, class]) - diag(cv_conf_matrix)[class])/sum(cv_conf_matrix[, class])
    train_class_error[class, fold] = (sum(train_conf_matrix[, class]) - diag(train_conf_matrix)[class])/sum(train_conf_matrix[, class])
    priors_tr[class, fold] = lda$prior[classes[class]]
    priors_te[class, fold] = length(which(test_set[, 10] == classes[class]))
  }
}
rownames(cv_class_error) = classes; colnames(cv_class_error) = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5", "Fold 6", "Fold 7", "Fold 8", "Fold 9", "Fold 10")
rownames(train_class_error) = classes; colnames(train_class_error) = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5", "Fold 6", "Fold 7", "Fold 8", "Fold 9", "Fold 10")
rownames(priors_tr) = classes; colnames(priors_tr) = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5", "Fold 6", "Fold 7", "Fold 8", "Fold 9", "Fold 10")
rownames(priors_te) = classes; colnames(priors_te) = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5", "Fold 6", "Fold 7", "Fold 8", "Fold 9", "Fold 10")
round(cv_class_error, 3); round(train_class_error, 3); round(priors_tr, 3); round(priors_te, 3)
probs
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2.4) Use the first two discriminant variables for a two-dimensional representation of the data together with centroids by using color-coding for the 6 classes of the class variable type. 
#Comment in view of the answer to point 2.2):

#We reproduce the previous plot of the data and centroids projected in the Discriminant Space (LD1, LD2):
par(mfrow = c(1,1))

means.hat = aggregate(lda.predict$x,by=list(glass$type),FUN=mean)
means.hat = means.hat[,-1]

plot(x = L1, y = L2, pch = 1, col = glass.col, main = "Data in the Reduced Subaspace (Scaled)", xlab = "LD1", ylab = "LD2")
for(row in seq(1, dim(means.hat)[1])){
  points(x = means.hat[row, 1], y = means.hat[row, 2], 
         pch = 21, col = "black", bg = lookup[row], 
         cex = 2, lwd = 2)
}
legend(x = "bottomright", legend = classes, fill = lookup)

#We plot in the (LD1, LD2) space the data separated by class in order to distinguish better the region occuped by each kind of glass:
par(mfrow = c(3,2))
plot(x = L1[which(glass[, "type"] == classes[1]), ], y = L2[which(glass[, "type"] == classes[1]), ], pch = 1, xlab = "LD1", ylab = "LD2",
     xlim = c(-4, 7.5), ylim = c(-8, 4), col = lookup[1], main = classes[1])
points(x = means.hat[1, 1], y = means.hat[1, 2], pch = 21, col = "black", bg = lookup[1], cex = 2, lwd = 2)

plot(x = L1[which(glass[, "type"] == classes[2]), ], y = L2[which(glass[, "type"] == classes[2]), ], pch = 1, xlab = "LD1", ylab = "LD2",
     xlim = c(-4, 7.5), ylim = c(-8, 4), col = lookup[2], main = classes[2])
points(x = means.hat[2, 1], y = means.hat[2, 2], pch = 21, col = "white", bg = lookup[2], cex = 2, lwd = 2)

plot(x = L1[which(glass[, "type"] == classes[3]), ], y = L2[which(glass[, "type"] == classes[3]), ], pch = 1, xlab = "LD1", ylab = "LD2",
     xlim = c(-4, 7.5), ylim = c(-8, 4), col = lookup[3], main = classes[3])
points(x = means.hat[3, 1], y = means.hat[3, 2], pch = 21, col = "black", bg = lookup[3], cex = 2, lwd = 2)

plot(x = L1[which(glass[, "type"] == classes[4]), ], y = L2[which(glass[, "type"] == classes[4]), ], pch = 1, xlab = "LD1", ylab = "LD2",
     xlim = c(-4, 7.5), ylim = c(-8, 4), col = lookup[4], main = classes[4])
points(x = means.hat[4, 1], y = means.hat[4, 2], pch = 21, col = "black", bg = lookup[4], cex = 2, lwd = 2)

plot(x = L1[which(glass[, "type"] == classes[5]), ], y = L2[which(glass[, "type"] == classes[5]), ], pch = 1, xlab = "LD1", ylab = "LD2",
     xlim = c(-4, 7.5), ylim = c(-8, 4), col = lookup[5], main = classes[5])
points(x = means.hat[5, 1], y = means.hat[5, 2], pch = 21, col = "black", bg = lookup[5], cex = 2, lwd = 2)

plot(x = L1[which(glass[, "type"] == classes[6]), ], y = L2[which(glass[, "type"] == classes[6]), ], pch = 1, xlab = "LD1", ylab = "LD2",
     xlim = c(-4, 7.5), ylim = c(-8, 4), col = lookup[6], main = classes[6])
points(x = means.hat[6, 1], y = means.hat[6, 2], pch = 21, col = "black", bg = lookup[6], cex = 2, lwd = 2)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2.5) Compute the training error and the 10-fold cross validation error for each reduced-rank LDA classifier.
#Plot both error curves against the number of discriminant directions, add full-rank LDA errors found in points 2.2) and 2.3). What classifier do you prefer? Comment:

#Implement the Reduced-Rank LDA with Cross Validation and compute the error scores:
generals_training = c()
cv_errors = c()
train_errors = c()
cv_class_rates_dim = matrix(0, nrow = K, ncol = K - 1)

for(dim in seq(1, K - 1)){
  model = lda(type ~ ., data = glass[, -11])
  training_p = predict(model, dimen = dim)$class
  confusion = table(predicted = training_p, true = glass[, 10])
  generals_training[dim] = (sum(confusion) - sum(diag(confusion)))/sum(confusion)
  
  folds_cv = c()
  folds_training = c()
  folds_classes = matrix(0, nrow = K, ncol = 10)
  for(fold in 1:10){
    test_set = glass[which(glass[, 11] == fold), -11]
    train_set = glass[which(glass[, 11] != fold), -11]
    
    lda = lda(type ~ ., data = train_set)
    
    cv_predictions = predict(lda, newdata = test_set, dimen = dim)$class
    cv_conf_matrix = table(predicted = cv_predictions, true = test_set[, 10])
    folds_cv[fold] = (sum(cv_conf_matrix) - sum(diag(cv_conf_matrix)))/sum(cv_conf_matrix)
    
    train_predictions = predict(lda, newdata = train_set, dimen = dim)$class
    train_conf_matrix = table(predicted = train_predictions, true = train_set[, 10])
    folds_training[fold] = (sum(train_conf_matrix) - sum(diag(train_conf_matrix)))/sum(train_conf_matrix)
    
    for(class in seq(1, K)){
      folds_classes[class, fold] = (sum(cv_conf_matrix[, class]) - diag(cv_conf_matrix)[class])/sum(cv_conf_matrix[, class])
    }
  }
  cv_errors[dim] = mean(folds_cv)
  train_errors[dim] = mean(folds_training)
  cv_class_rates_dim[, dim] = apply(folds_classes, 1, mean, na.rm = T)
}

names(generals_training) = c("1 direction", "2 directions", "3 directions", "4 directions", "5 directions")
names(cv_errors) = c("1 direction", "2 directions", "3 directions", "4 directions", "5 directions")
names(train_errors) = c("1 direction", "2 directions", "3 directions", "4 directions", "5 directions")
colnames(cv_class_rates_dim) = c("1 direction", "2 directions", "3 directions", "4 directions", "5 directions"); rownames(cv_class_rates_dim) = classes
generals_training
cv_errors
train_errors
round(cv_class_rates_dim, 3)

final = cbind(generals_training, cv_errors); rownmaes(final) = c("Training Errors", "Validation Errors")

#We produce the plot that shows the evolution of the cv and training error rates for different numbers of Discrminant Directions used:
par(mfrow = c(1, 1))
plot(x = 1:5, y = cv_errors, main = "Reduced Rank LDA Performances", 
     xlab = "Discriminant Directions used", ylab = "Missclassification Rates", 
     type = "b", col = "orange", ylim = c(0.3, 0.49))
points(x = 1:5, y = generals_training, type = "b", col = "blue")
legend(x = "topright", fill = c("orange", "blue"), legend = c("Validation Errors", "Training Errors"))
abline(v = 4, lty = "dashed")

scaled_ld_all = cbind(diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 1], diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 2],
                      diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 3], diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 4], diag(diag(sigma_p))^(1/2)%*%discriminant_directions[, 5])
rownames(scaled_ld_all) = colnames(glass)[-c(10, 11)]; colnames(scaled_ld_all) = c("Scaled LD1", "Scaled LD2", "Scaled LD3", "Scaled LD4", "Scaled LD5")
scaled_ld_all
sort(scaled_ld_all[, 1], decreasing = TRUE)
sort(scaled_ld_all[, 2], decreasing = TRUE)
sort(scaled_ld_all[, 3], decreasing = TRUE)
sort(scaled_ld_all[, 4], decreasing = TRUE)
sort(scaled_ld_all[, 5], decreasing = TRUE)