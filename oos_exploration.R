library(fbi) # this is the FredMD package
library(imputeTS) # for imputing mean with na_mean

## IDEA:  compare to mean model for out-of-sample

## DATA UPLOAD / PREPROCESSING ##
setwd("C:/Users/Jordan/Desktop/MS&E349/spca_play")
fred_md <- fredmd("current.csv", date_start = as.Date("1960-01-01"), date_end = as.Date("2020-01-01"))


colnames(fred_md) = colnames(read.csv("current.csv")) # fill in column names

na_idx=apply(fred_md,2,function(x)!any(is.na(x)))  # find indices that are not NA

data = fred_md[,na_idx]  #obtain 123 columns (not including the date), as in the paper


data <- rm_outliers.fredmd(data)
data <- na_mean(data)  # impute NA with mean
data = data[,2:124]  #obtain 123 columns

fred_md <- rm_outliers.fredmd(fred_md)
fred_md <- na_mean(fred_md)

cut_idx <- 300
spca_opt = TRUE


N <- dim(data)[2]


target <- fred_md[2:721, "CPIAUCSL"] # obtain one-month ahead values
target_train <- target[1:cut_idx]
r2_values <- c()

for (k in 1:8){
num_factors = k
num_factors = 4
target_train_factor <- target_train
target_train_ar <- target_train


for (i in cut_idx: 719){
  
  armodel <- arima(target, order = c(1, 0, 0), include.mean = FALSE)  # 9 lags for inflation
  
  
  target_train_ar <- append(target_train_ar, as.numeric(predict(armodel, n.ahead = 1)[1]))
  
  scaled_data <- scale(data[1:i,])
  ext_scaled_data <- scale(data[1:(i+1),]) 
  
  if (spca_opt == TRUE){
    
    
    gamma = c()
    
    for (j in colnames(data)){
      macro_var <- as.numeric(data[1:i,j])
      lm_model <- lm(target[1:i] ~ scale(macro_var))
      #scaler <- lm_model$coefficients[2] / (sum((lm_model$residuals)^2) / 718) # normalize by variance
      gamma = append(gamma, as.numeric(lm_model$coefficients[2]))
      #gamma = append(gamma, as.numeric(scaler))
    }
    
    gamma_X <- diag(gamma) %*% t(scaled_data)
    lambda_hat <-  svd(gamma_X %*% t(gamma_X))$v[,1:num_factors]
    Fhat <- 1/N * t(gamma_X) %*% lambda_hat
    
    spca_fit <- lm(target[1:i] ~ as.matrix(Fhat)[,1:num_factors])
    gamma_X_ext <- diag(gamma)%*%t(ext_scaled_data)
    lambda_hat_ext <-  svd(gamma_X_ext %*% t(gamma_X_ext))$v[,1:num_factors]
    
    Fhat_ext <- 1/N * t(gamma_X_ext) %*% lambda_hat_ext
    
    spca_pred <-  spca_fit$coefficients[2:(num_factors + 1)]%*% as.matrix(Fhat_ext)[(i+1),] + spca_fit$coefficients[1]
    target_train_factor <- append(target_train_factor, spca_pred)
    
  }
  else{
    #normal pca
    #conduct pca on range [1, i] to obtain factors and regress on target to obtain loadings
    lambda_hat <-  svd(t(scaled_data) %*% scaled_data)$v[,1:num_factors]
    Fhat <- 1/N * scaled_data %*% lambda_hat
    
    pca_fit <- lm(target[1:i] ~ as.matrix(Fhat)[,1:num_factors])
    
    #now that we have the loadings, conduct PCA on [1, i+1] to obtain factors and form a prediction of the target
    
    lambda_hat_ext <-  svd(t(ext_scaled_data) %*% ext_scaled_data)$v[,1:num_factors]
    Fhat_ext <- 1/N * ext_scaled_data %*% lambda_hat_ext
    
    pca_pred <-  pca_fit$coefficients[2:(num_factors + 1)]%*% as.matrix(Fhat_ext)[(i+1),] + pca_fit$coefficients[1]
    target_train_factor <- append(target_train_factor, pca_pred)
    
  }
  
}

r2 = 1 - sum((target - target_train_factor)^2) / sum((target - target_train_ar)^2)
print(r2)
r2_values <- append(r2_values, r2)
}
plot(r2_values)

r2_spca_values <- r2_values
r2_pca_values <- r2_values
r2_scaler_values <- r2_values

plot(r2_pca_values, type = "l", ylim = c(0.4, 0.50))
lines(r2_spca_values, col = "red")
lines(r2_scaler_values, col = "blue")
