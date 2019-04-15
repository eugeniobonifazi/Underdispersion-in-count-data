q <- function(y, mu, ni) (mu^y/factorial(y))^ni

dcomp <- function(x, lambda, nu) 
{compoisson::dcom(x=x, lambda = lambda, nu = nu)}

rcomp <- function(n, lambda, nu) 
{compoisson::rcom(n=n, lambda=lambda, nu=lambda)}

a_ratio <- function(y, y_star, mu, nu, mu_cand, nu_cand, prior.mu_par=list(mu, sd), prior.nu_par=list(mu, sd)){
  
  a = sum(log(q(y, mu_cand, nu_cand))+log(q(y_star, mu, nu)) 
          + log(dnorm(mu_cand, prior.mu_par$mu, prior.mu_par$sd)) 
          + log(dnorm(nu_cand, prior.nu_par$mu,prior.nu_par$sd))
          -log(q(y, mu, nu)) - log(q(y_star, mu_cand, nu_cand)) 
          - log(dnorm(mu, prior.mu_par$mu, prior.mu_par$sd)) 
          - log(dnorm(nu, prior.nu_par$mu,prior.nu_par$sd)))
}


exchange.alg <- function(nsim, y, likelihood_par=list(n=1, mu=3, nu=1.5),
                          prior.mu_par=list(mu=0, sd=2), prior.nu_par=list(mu=0, sd=2)){
  
  acc = 0
  mu = rep(NA,nsim)
  nu = rep(NA,nsim)
  yres = rep(NA,nsim)
 yres[1] = rcomp(likelihood_par$n, likelihood_par$mu, likelihood_par$nu)
  mu_cand = 0
  nu_cand = 0
  mu[1] = mu_cand
  nu[1] = nu_cand
  
  for(i in 2:nsim){
    mu_cand = rnorm(n=1, prior.mu_par$mu, prior.mu_par$sd)
    nu_cand = rnorm(n=1, prior.nu_par$mu, prior.nu_par$sd)
    while (mu_cand<=0){
      mu_cand = rnorm(n=1, prior.mu_par$mu, prior.mu_par$sd)
    }
    while (nu_cand<=0){
      nu_cand <- rnorm(n=1, prior.nu_par$mu, prior.nu_par$sd)
    }
    
    y_star = rcomp(likelihood_par$n, likelihood_par$mu, likelihood_par$nu)
    
    a = a_ratio(y, y_star, mu[i-1], nu[i-1], mu_cand, nu_cand, 
                prior.mu_par=prior.mu_par, prior.nu_par=prior.nu_par)
    
    
    r = min(exp(a),1)
    u = runif(1)
    
    if(u<r){
      mu[i] = mu_cand
      nu[i] = nu_cand
      yres[i] = y_star
      acc = acc + 1/nsim
    }
    else{
      mu[i] = mu[i-1]
      nu[i] = nu[i-1]
      yres[i] = yres[i-1]
    }
  }
  result = list(mu_acc = mu, nu_acc = nu, yres = yres, acc = acc*100)
  return(result)
}

exchange.alg2 <- function(nsim, y, likelihood_par=list(n=1, mu=3, nu=1.8),
                         prior.mu_par=list(mu=0, sd=2), prior.nu_par=list(mu=0, sd=2)){
  
  n=length(y)
  acc = 0
  mu = rep(NA,nsim)
  nu = rep(NA,nsim)
  yres = rep(NA,nsim)
  yres[1] = rcomp(likelihood_par$n, likelihood_par$mu, likelihood_par$nu)
  mu_cand = 0
  nu_cand = 0
  mu[1] = mu_cand
  nu[1] = nu_cand
  
  for(i in 2:nsim){
    mu_cand = rnorm(n=1, prior.mu_par$mu, prior.mu_par$sd)
    nu_cand = rnorm(n=1, prior.nu_par$mu, prior.nu_par$sd)
    while (mu_cand<=0){
      mu_cand = rnorm(n=1, prior.mu_par$mu, prior.mu_par$sd)
    }
    while (nu_cand<=0){
      nu_cand <- rnorm(n=1, prior.nu_par$mu, prior.nu_par$sd)
    }
    
    y_star = rcomp(n, likelihood_par$mu, likelihood_par$nu)
    
    a = a_ratio(y, y_star, mu[i-1], nu[i-1], mu_cand, nu_cand, 
                prior.mu_par=prior.mu_par, prior.nu_par=prior.nu_par)
    
    
    r = min(exp(a),1)
    u = runif(1)
    
    if(u<a){
      mu[i] = mu_cand
      nu[i] = nu_cand
      yres[i] = y_star
      acc = acc + 1/nsim
    }
    else{
      mu[i] = mu[i-1]
      nu[i] = nu[i-1]
      yres[i] = yres[i-1]
    }
  }
  result = list(mu_acc = mu, nu_acc = nu, acc = acc*100)
  return(result)
}


exchange_alg.plot <- function(cmp_exc){
  
  par(mfrow=c(1,1))
  mu_acc = cmp_exc$mu_acc
  nu_acc = cmp_exc$nu_acc
  
  expr <- vector("expression", 2)
  mu_mean = round(mean(unique(mu_acc)),2)
  nu_mean = round(mean(unique(nu_acc)),2)
  
  expr[[1]] <- bquote(bar(mu)==.(mu_mean))
  expr[[2]] <- bquote(bar(nu)==.(nu_mean))
  
 if(length(mu_acc)==1243){ 
   plot(mu_acc, type = "l", main = bquote(mu ~ "accepted"), 
       xlab = "simulations", ylab = bquote(mu))
  abline(h=mu_mean, lwd = 2, col = 2, lty = "longdash")
  legend("topright", legend = expr[1], lwd = 2,
         col = 2, lty = "longdash", cex = .5)
  
  plot(nu_acc, type = "l", main = bquote(nu ~ "accepted"), 
       xlab = "simulations", ylab = bquote(nu))
  abline(h=nu_mean, lwd = 2, col = 2, lty = "longdash")
  legend("topright", legend = expr[2], lwd = 2, 
         col = 2, lty = "longdash", cex = .5)
 }
  else{
    last_1000 = (length(mu_acc)-1000):length(mu_acc)
    plot(mu_acc[last_1000], type = "l", main = bquote(mu ~ "accepted"), 
         xlab = "simulations", ylab = bquote(mu))
    abline(h=mu_mean, lwd = 2, col = 2, lty = "longdash")
    legend("topright", legend = expr[1], lwd = 2,
           col = 2, lty = "longdash", cex = .5)
    
    plot(nu_acc[last_1000], type = "l", main = bquote(nu ~ "accepted"), 
         xlab = "simulations", ylab = bquote(nu))
    abline(h=nu_mean, lwd = 2, col = 2, lty = "longdash")
    legend("topright", legend = expr[2], lwd = 2, 
           col = 2, lty = "longdash", cex = .5)
  }
}
