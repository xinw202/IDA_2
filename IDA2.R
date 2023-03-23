library(maxLik)
############Q2
load('dataex2.Rdata')
# set log likelihood function max
log_like = function(param,data){
  x = data[,1];r = data[,2]
  mu = param
  sum(r*dnorm(x, mean = mu, sd = 1.5, log = TRUE) +
        (1-r)*pnorm(x, mean = mu, sd = 1.5, log.p = TRUE))
}
mle = maxLik(logLik = log_like, data = dataex2, start = c(mu = 0))
summary(mle)


############Q4
load("dataex4.Rdata")
#Q function
Q_fn <- function(beta,input=dataex4){
  #Replace NA in y_obs
  y_obs <- which(!is.na(input$Y))
  y_mis <- which(is.na(input$Y))
  beta_0 <- beta[1]
  beta_1 <- beta[2]
  S <- rep(0,length(input$Y))
  for(i in y_obs){
    S[i] <- input$Y[i]*(beta_0+beta_1*input$X[i])-log(1+exp(beta_0+beta_1*input$X[i]))
  }
  beta_new <- sum(S)
  return(beta_new)
}
beta_o <- c(0,0);
beta_new <- maxLik(Q_fn,start=beta_o);
beta_new

#############Q5
load("dataex5.Rdata")
y <- dataex5
em.mixture.two.pareto <- function(y, theta_0, eps){
  n <- length(y)
  theta <- theta_0
  p <- theta[1]
  lambda <- theta[2]
  mu <- theta[3]
    
  diff <- 1
  while(diff > eps){
    theta.old <- theta
      
    #E-step
    ptilde_1 <- p*(lambda*y^(-lambda-1))
    ptilde_2 <- (1-p)*(mu*y^(-mu-1))
    ptilde <- ptilde_1/(ptilde_1 + ptilde_2)
    
    #M-step
    p <- mean(ptilde)
    lambda <- sum(ptilde)/sum(ptilde*log(y))
    mu <- sum(1-ptilde)/sum((1-ptilde)*log(y))
    theta <- c(p, lambda, mu)
    diff <- sum(abs(theta - theta.old))
  }
  return(theta)  }
res <- em.mixture.two.pareto(y = y, theta_0 = c(0.3, 0.3, 0.4), eps = 0.0001)
pest <- res[1]; lambdaest <- res[2]; muest <- res[3]
pest
lambdaest
muest
# Draw the histogram of the data with the estimated density superimposed.
hist(y, breaks = "Freedman-Diaconis", 
     main = " Histogram of the data",
     xlab = "y",
     ylab = "Density",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
     freq = F, xlim = c(0, 100), ylim = c(0,0.02))
curve((pest*(lambdaest*x^(-lambdaest-1)) + (1-pest)*(muest*x^(-muest-1))),
      add = TRUE, lwd = 2, col = "blue2")

