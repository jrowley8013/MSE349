library(fbi) # this is the FredMD package
library(imputeTS) # for imputing mean with na_mean



## DATA UPLOAD / PREPROCESSING ##
setwd("C:/Users/Jordan/Desktop/MS&E349/spca_play")
fred_md <- fredmd("current.csv", date_start = as.Date("1960-01-01"), date_end = as.Date("2020-01-01"))

colnames(fred_md) = colnames(read.csv("current.csv")) # fill in column names

na_idx=apply(fred_md,2,function(x)!any(is.na(x)))  # find indices that are not NA

data = fred_md[,na_idx]  #obtain 123 columns (not including the date), as in the paper


data <- rm_outliers.fredmd(data)
data <- na_mean(data)  # impute NA with mean
data = data[,2:124]  #obtain 123 columns




## OBTAIN OUT-OF-SAMPLE R2 ##

oos_R2 <- function(cut_time, target_name, num_factors, spca_opt = TRUE) {


  cut_idx <- which.max(as.Date(cut_time) <= fred_md[,'sasdate'])

  N <- dim(data)[2]


  target <- fred_md[2:721, target_name] # obtain one-month ahead values
  target_train <- target[1:cut_idx]


  target_train_ar <- target_train
  target_train_factor <- target_train




  for (i in cut_idx: 719){
    
    # set up AR model and use it to fill in a forecast
    armodel <- arima(target_train_ar, order = c(4, 0, 0), include.mean = FALSE)  # fix at 4 lags for now (will change later)
    target_train_ar <- append(target_train_ar, as.numeric(predict(armodel, n.ahead = 1)[1]))
    
    scaled_data <- scale(data[1:i,])
    ext_scaled_data <- scale(data[1:(i+1),]) #ppp change this?
    
    if (spca_opt == TRUE){
      # we conduct SPCA and append the forecasted value
      
      gamma = c()
    
      for (j in colnames(data)){
        macro_var <- as.numeric(data[1:i,j])
        lm_model <- lm(target_train_factor ~ scale(macro_var))
        gamma = append(gamma, as.numeric(lm_model$coefficients[2]))
      }
      
      
      gamma_X <- diag(gamma) %*% t(scaled_data)
      lambda_hat <-  svd(gamma_X %*% t(gamma_X))$v[,1:num_factors]  #lambdas agree
      Fhat <- 1/N * t(gamma_X) %*% lambda_hat
      
      spca_fit <- lm(target_train_factor ~ as.matrix(Fhat)[,1:num_factors])
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
      
      pca_fit <- lm(target_train_factor ~ as.matrix(Fhat)[,1:num_factors])
      
      #now that we have the loadings, conduct PCA on [1, i+1] to obtain factors and form a prediction of the target

      lambda_hat_ext <-  svd(t(ext_scaled_data) %*% ext_scaled_data)$v[,1:num_factors]
      Fhat_ext <- 1/N * ext_scaled_data %*% lambda_hat_ext

      pca_pred <-  pca_fit$coefficients[2:(num_factors + 1)]%*% as.matrix(Fhat_ext)[(i+1),] + pca_fit$coefficients[1]
      target_train_factor <- append(target_train_factor, pca_pred)
      
    }
    
  }
  

  r2 = 1 - sum((target - target_train_factor)^2) / sum((target - target_train_ar)^2)
  
  return(r2)
}

#oos R2 using vanilla PCA
oos_R2("1984-12-01", "UNRATE", 6, spca_opt = FALSE)  # we expect 0.08, but obtain 0.0625

#oos R2 using SPCA
oos_R2("1984-12-01", "UNRATE", 6, spca_opt = TRUE)  # expect 0.1, but obtain 0.08154
