thames <- function(lps = NULL,params,n_samples = NULL,d = NULL, radius = NULL,
                   p = 0.025, q = 1-p, lp_func = NULL,
                   bound = NULL, n_simuls = 1e5){
  
  # dimension of parameter space
  if(is.null(d)){
    d <- dim(params)[2]
    if(is.null(d)){
      d <- 1
    }
  }
  
  # radius of A
  if(is.null(radius)){
    radius <- sqrt(d+1)
  }
  
  # number of posterior samples
  if(is.null(n_samples)){
    if(d==1){
      n_samples <- length(params)
    }else{
      n_samples <- dim(params)[1]
    }
  }
  
  # calculate unnormalized log posterior values
  if(is.null(lps)){
    lps <- lp_func(params)
  }
  
  if(length(lps) != n_samples){
    return('Error: number of unnormalized log posterior values does not match posterior sample size.')
  }
  
  # split the sample
  n1 <- n_samples %/% 2
  n2 <- n_samples - n1
  if(d==1){
    params1 <- params[1:n1]
    params2 <- params[(n1+1):n_samples]
  }else{
    params1 <- params[1:n1,]
    params2 <- params[(n1+1):n_samples,]
  }
  lps1 <- lps[1:n1]
  lps2 <- lps[(n1+1):n_samples]
  
  # calculate posterior mode and covariance from first half of sample
  if(d==1){
    theta_hat <- mean(params1)
    sigma_hat <- var(params1)
    sigma_svd = sigma_hat # one-dimensional svd is just the variance
    
    log_det_sigma_hat <- log(sigma_hat)
    
    # which samples are in A?
    inA <- sapply(params2,function(theta){
      # calculate distance from theta_hat
      theta_tilde <- (theta-theta_hat)/sqrt(sigma_hat)
      
      # is distance of theta less than the radius?
      return(sum(theta_tilde^2) < radius^2)
    })
  }else{
    theta_hat <- colMeans(params1)
    sigma_hat <- cov(params1)
    
    # calculate SVD of sigma_hat
    sigma_svd = svd(sigma_hat)
    
    # calculate log(det(sigma_hat))
    log_det_sigma_hat = sum(log(sigma_svd$d))
    
    # which samples are in A?
    inA <- apply(params2,1,function(theta){
      # calculate distance from theta_hat
      theta_tilde <-  sigma_svd$d^(-1/2) * (t(sigma_svd$v) %*% (theta-theta_hat))
      
      # is distance of theta less than the radius?
      return(sum(theta_tilde^2) < radius^2)
    })
  }
  num_inA <- as.numeric(inA)
  
  # log volume of A
  logvolA = d*log(radius)+(d/2)*log(pi)+log_det_sigma_hat/2-lgamma(d/2+1)
  
  log_zhat_inv <- log(mean(exp(-(lps2 - max(lps2))) * num_inA )) - logvolA - max(lps2)
  
  
  sample_colors <- ifelse(inA, "navy", "grey")
  plot(params2, col = sample_colors, pch = 19, cex = 0.3,
       xlab = expression(theta[1]), ylab = expression(theta[2]),
       main = "THAMES")
  
  
  draw_ellipse(center = theta_hat, cov = sigma_hat, radius = radius, col = "red", lwd = 1.5)
  
  
  legend("topright", legend = c("Samples in Ellipsoids", "Samples outside Ellipsoids", "Ellipsoid Boundaries"),
         col = c("navy", "grey", "red"),
         pch = c(19, 19, NA), lty = c(NA, NA, 1),cex = 0.5)
  
  return(list(log_zhat_inv = -log_zhat_inv,
              zhat_inv = exp(-log_zhat_inv)  ))
}