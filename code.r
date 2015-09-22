# comparing r and rcpp for linear regression bootstrapping
require(Rcpp)
require(inline)
require(RcppArmadillo)

boot_code = "

  arma::colvec y = Rcpp::as<arma::colvec> (y_);
  arma::mat    X = Rcpp::as<arma::mat>(X_) ;

  int p = X.n_cols;
  int n = X.n_rows;
  int m = as<int>(m_);

  arma::mat    H   = arma::inv(arma::trans(X)*X) * arma::trans(X);
  arma::colvec b   = H*y;
  arma::colvec Xb  = X*b;
  arma::colvec res = y-Xb;

  double sig2 = arma::as_scalar( arma::trans(res)*res/(n-p) );
  double sig  = sqrt(sig2);

  NumericMatrix beta(m,p);
  NumericVector bb(p);

  for (int i = 0; i < m; i++){  
    res       = rnorm(n,0,sig); 
    bb        = H*(Xb + res);
    beta(i,_) = bb;
  }
  return Rcpp::wrap(beta);
"

boot <- cxxfunction(signature(y_="numeric", X_="matrix",m_="integer"),
                    body=boot_code, 
                    plugin="RcppArmadillo")
                    
n    <- 10^2
 reps <- 10^5

 time        <- rep(0,2)
 names(time) <- c("R","Rcpp")

 X <- cbind(1,rnorm(n))
 y <- 10*X[,2] + rnorm(n,0,5)
 
  tick      <- proc.time()[3]
 beta1     <- matrix(0,reps,2)
 H         <- solve(t(X)%*%X)%*%t(X) 
 b         <- H%*%y
 Xb        <- X%*%b
 res       <- y-Xb
 sig       <- sqrt(sum(res^2)/(n-2))
 for(rep in 1:reps){
   y_rep       <- rnorm(n,Xb,sig)
   beta1[rep,] <- H%*%y_rep
 }
 tock      <- proc.time()[3]
 time[1]   <- tock-tick
 
  tick      <- proc.time()[3]
 beta2     <- boot(y,X,reps)
 tock      <- proc.time()[3]
 time[2]   <- tock-tick
 
 summary(lm(y~X-1))$coef
##      Estimate Std. Error    t value     Pr(>|t|)
## X1 -0.1836637  0.4306924 -0.4264382 6.707238e-01
## X2 10.2556490  0.4244457 24.1624513 4.553260e-43
apply(beta1,2,mean)
## [1] -0.1820817 10.2546089
apply(beta1,2,sd)
## [1] 0.4313784 0.4240668
apply(beta2,2,mean)
## [1] -0.1835478 10.2553940
apply(beta2,2,sd)
## [1] 0.4309970 0.4253691
time
##    R Rcpp 
## 1.14 0.69
