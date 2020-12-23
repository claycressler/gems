pick_individuals_multivariate <-function(N, traitmeans, traitsds, corr){
  ## Generate Lognormal random variables with zeros. 
  ## traitmeans - Vector of mean trait values
  ## traitsds - Vector of trait standard deviations
  ## R - Correlation matrix - should have 1s along the diagonal and off-diagonal elements reflecting correlations between different variables
  
  ## Number of traits
  p <- length(traitmeans);
  ## mu and sigma are the converted mean and covariance matrix of the log-transformed normal variable.
  mu <- matrix(0,p,1);
  sigma <- matrix(0,p,p);
  for (i in 1:p) {
      mu[i] = log(traitmeans[i]^2 / sqrt(traitsds[i]^2+traitmeans[i]^2))
      sigma[i,i] = log(1 + traitsds[i]^2/traitmeans[i]^2);
  }
  ## Use correlation matrix R to create the covariance matrix Cov
  Cov <- sqrt(sigma)%*%corr%*%sqrt(sigma);
  ## Following the protocol laid out for drawing multivariate normal random variables in
  ## https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
  ## Draw normal random variables with mean 0 and sd 1
  X <- matrix(rnorm(p*N,0,1), nrow = N, ncol = p);
  ## Convert to lognormal random variables with the appropriate covariance structure
  ## using Cholesky decomposition of the covariance matrix
  for (i in 1:N){
    X[i,] <- exp(mu + chol(Cov)%*%X[i,])
  }      
  return(X);
}

pick_individuals <- function(N0, traitmean, traitsd) {
  mu <- log(traitmean^2 / sqrt(traitsd^2+traitmean^2))
  sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
  ## record this initial distribution in the output
  return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}


## Comparing the two methods
traitmeans <- c(1.5,0.8)
traitsds <- c(0.15, 0.08)
## if correlations are equal to 0, then pick_individuals and pick_individuals_multivariate should look very similar
corr <- matrix(c(1,0,0,1), nrow=2, byrow=T)

trait1 <- pick_individuals(1e5, traitmean=traitmeans[1], traitsd=traitsds[1])
trait2 <- pick_individuals(1e5, traitmean=traitmeans[2], traitsd=traitsds[2])
traits <- pick_individuals_multivariate(1e5, traitmeans, traitsds, corr)

## Compare - methods are identical
mean(trait1)
sd(trait1)
mean(traits[,1])
sd(traits[,1])

mean(trait2)
sd(trait2)
mean(traits[,2])
sd(traits[,2])

## if correlations are non-zero (positive), how does that affect things?
corr <- matrix(c(1,0.5,0.5,1), nrow=2, byrow=T)
traits_pos <- pick_individuals_multivariate(1e5, traitmeans, traitsds, corr)
## means remain the same but sd for each trait have changed
apply(traits_pos, 2, mean)
apply(traits_pos, 2, sd)

## Graphical comparison shows that the correlations are working as expected
plot(traits_pos)

plot(traits, col=2)




