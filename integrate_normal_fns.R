### normal integral functions
# series of functions bulding towards function to compute integrals of x ^ k * dnorm(x, mu, sigma) from a to b
# used in VI and W update steps

# antiderivative of x ^ 8 * dnorm(x, mu, sigma)
E_x8_ind_a_b_normal_helper <- function(mu, sigma, x){
  
  retVal <- .5 * (mu ^ 8 + 28 * mu ^ 6 * sigma ^ 2 + 210 * mu ^ 4 * sigma ^ 4 + 420 * mu ^ 2 * sigma ^ 6 + 105 * sigma ^ 8)
  
  if (x == Inf){
    retVal
  }else{
    retVal * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      (sigma * (mu ^ 7 + mu ^ 6 * x + mu ^ 5 * (27 * sigma ^ 2 + x ^ 2) + mu ^ 4 * (25 * sigma ^ 2 * x + x ^ 3) +
                  mu ^ 3 * (185 * sigma ^ 4 + 22 * sigma ^ 2 * x ^ 2 + x ^ 4) +
                  mu ^ 2 * (141 * sigma ^ 4 * x + 18 * sigma ^ 2 * x^3 + x^5) +
                  mu * (279 * sigma ^ 6 + 87 * sigma ^ 4 * x ^ 2 + 13 * sigma ^ 2 * x ^ 4 + x ^ 6) +
                  105 * sigma ^ 6 * x + 35 * sigma ^ 4 * x ^ 3 + 7 * sigma ^ 2 * x ^ 5 + x ^ 7)) *
      exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}

# antiderivative of x ^ 7 * dnorm(x, mu, sigma)
E_x7_ind_a_b_normal_helper <- function(mu, sigma, x){
  
  retVal <- .5 * mu * (mu ^ 6 + 21 * mu ^ 4 * sigma ^ 2 + 105 * mu ^ 2 * sigma ^ 4 + 105 * sigma ^ 6)
  
  if (x == Inf){
    retVal
  }else{
    retVal * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      (sigma * (mu ^ 6 + mu ^ 5 * x + mu ^ 4 * (20 * sigma ^ 2 + x ^ 2) + mu ^ 3 * (18 * sigma ^ 2 * x + x ^ 3) +
                  mu ^ 2 * (87 * sigma ^ 4 + 15 * sigma ^ 2 * x ^ 2 + x ^ 4) +
                  mu * (57 * sigma ^ 4 * x + 11 * sigma ^ 2 * x ^ 3 + x^5) +
                  48 * sigma ^ 6 + 24 * sigma ^ 4 * x ^ 2 + 6 * sigma ^ 2 * x ^ 4 + x ^ 6)) *
      exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}


# antiderivative of x ^ 6 * dnorm(x, mu, sigma)
E_x6_ind_a_b_normal_helper <- function(mu, sigma, x){
  
  retVal <- .5 * (mu ^ 6 + 15 * mu ^ 4 * sigma ^ 2 + 45 * mu ^ 2 * sigma ^ 4 + 15 * sigma ^ 6)
  
  if (x == Inf){
    retVal
  }else{
    retVal * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      (sigma * (mu ^ 5 + mu ^ 4 * x + mu ^ 3 * (14 * sigma ^ 2 + x ^ 2) + mu ^ 2 * (12 * sigma ^ 2 * x + x ^ 3) +
                  mu * (33 * sigma ^ 4 + 9 * sigma ^ 2 * x ^ 2 + x ^ 4) + 15 * sigma ^ 4 * x + 5 * sigma ^ 2 * x ^ 3 + x ^ 5)) *
      exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}

# antiderivative of x ^ 5 * dnorm(x, mu, sigma)
E_x5_ind_a_b_normal_helper <- function(mu, sigma, x){
  
  retVal <- .5 * mu * (mu ^ 4 + 10 * sigma ^ 2 * mu ^ 2 + 15 * sigma ^ 4)
  
  if (x == Inf){
    retVal
  }else{
    retVal * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      (sigma * (mu ^ 4 + x * mu ^ 3 + mu ^ 2 * (9 * sigma ^ 2 + x ^ 2) +
                  mu * (7 * x * sigma ^ 2 + x ^ 3) + 8 * sigma ^ 4 + 4 * sigma ^ 2 * x ^ 2 + x ^ 4)) *
      exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}

# antiderivative of x ^ 4 * dnorm(x, mu, sigma)
E_x4_ind_a_b_normal_helper <- function(mu, sigma, x){
  
  retVal <- .5 * (3 * sigma ^ 4 + 6 * mu ^ 2 * sigma ^ 2 + mu ^ 4)
  
  if (x == Inf){
    retVal
  }else{
    retVal * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      sigma * (5 * mu * sigma ^ 2 + mu ^ 3 + (3 * sigma ^ 2 + mu ^ 2) * x + mu * x ^ 2 + x ^ 3) * 
      exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}

# antiderivative of x ^ 3 * dnorm(x, mu, sigma)
E_cube_ind_a_b_normal_helper <- function(mu, sigma, x){
  
  retVal <- .5 * mu * (3 * sigma ^ 2 + mu ^ 2)
  
  if (x == Inf){
    retVal
  }else{
    retVal * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      sigma * (2 * sigma ^ 2 + mu ^ 2 + mu * x + x ^ 2) *
      exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}

# antiderivative of x ^ 2 * dnorm(x, mu, sigma)
E_square_ind_a_b_normal_helper <- function(mu, sigma, x){
  if (x == Inf){
    .5 * (sigma ^ 2 + mu ^ 2)
  }else{
    .5 * (sigma ^ 2 + mu ^ 2) * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      sigma * (mu + x) *
      exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}

# antiderivative of x * dnorm(x, mu, sigma)
E_x_ind_a_b_normal_helper <- function(mu, sigma, x){
  if (x == Inf){
    .5 * mu
  }else{
    .5 * mu * pracma::erf((x - mu) / (sqrt(2) * sigma)) -
      sigma * exp(-(x - mu) ^ 2 / (2 * sigma ^ 2)) / sqrt(2 * pi)
  }
}

# antiderivative of dnorm(x, mu, sigma)
# cdf
E_ind_a_b_normal_helper <- function(mu, sigma, x){
  if (x == Inf){
    1
  }else{
    pnorm(x, mu, sigma)
  }
}

E_ind_a_b_normal_helper2 <- function(mu, sigma, x){
  if (x == Inf){
    .5
  }else{
    .5 * pracma::erf((x - mu) / (sqrt(2) * sigma))
  }
}

# list of helper functions / antiderivatives
# to access in function below
normalHelperFns <- list(
  E_ind_a_b_normal_helper,
  E_x_ind_a_b_normal_helper,
  E_square_ind_a_b_normal_helper,
  E_cube_ind_a_b_normal_helper,
  E_x4_ind_a_b_normal_helper,
  E_x5_ind_a_b_normal_helper,
  E_x6_ind_a_b_normal_helper,
  E_x7_ind_a_b_normal_helper,
  E_x8_ind_a_b_normal_helper
)

# integrate x ^ k * dnorm(x, mu, sigma) from a to b
# for k = 0,1,...,7
E_xk_ind_a_b_normal <- function(k, mu, sigma, a, b){
  
  normal_helper_fn <- normalHelperFns[[k + 1]]
  
  normal_helper_fn(mu, sigma, b) -
    normal_helper_fn(mu, sigma, a)
}

############ new functions for faster computation!!

get_R_1_mat <- function(mu, sigma){
  cbind(
    ## k = 0
    .5,
    ## k = 1
    .5 * mu,
    ## k = 2
    .5 * (sigma ^ 2 + mu ^ 2),
    ## k = 3
    .5 * mu * (3 * sigma ^ 2 + mu ^ 2),
    ## k = 4
    .5 * (3 * sigma ^ 4 + 6 * mu ^ 2 * sigma ^ 2 + mu ^ 4),
    ## k = 5
    .5 * mu * (mu ^ 4 + 10 * sigma ^ 2 * mu ^ 2 + 15 * sigma ^ 4),
    ## k = 6
    .5 * (mu ^ 6 + 15 * mu ^ 4 * sigma ^ 2 + 45 * mu ^ 2 * sigma ^ 4 + 15 * sigma ^ 6),
    ## k = 7
    .5 * mu * (mu ^ 6 + 21 * mu ^ 4 * sigma ^ 2 + 105 * mu ^ 2 * sigma ^ 4 + 105 * sigma ^ 6),
    ## k = 8
    .5 * (mu ^ 8 + 28 * mu ^ 6 * sigma ^ 2 + 210 * mu ^ 4 * sigma ^ 4 + 420 * mu ^ 2 * sigma ^ 6 + 105 * sigma ^ 8)
  )
}

get_R_2_mat <- function(x, sigma, mu){
  if(x == Inf){
    return(0)
  }else{
    return(
      cbind(
        ## k = 0
        0,
        ## k = 1
        sigma ^ 2,
        ## k = 2
        sigma ^ 2 * (mu + x),
        ## k = 3
        sigma ^ 2 * (2 * sigma ^ 2 + mu ^ 2 + mu * x + x ^ 2),
        ## k = 4
        sigma ^ 2 * (5 * mu * sigma ^ 2 + mu ^ 3 + (3 * sigma ^ 2 + mu ^ 2) * x + mu * x ^ 2 + x ^ 3),
        ## k = 5
        sigma ^ 2 * (mu ^ 4 + x * mu ^ 3 + mu ^ 2 * (9 * sigma ^ 2 + x ^ 2) +
                       mu * (7 * x * sigma ^ 2 + x ^ 3) + 8 * sigma ^ 4 + 4 * sigma ^ 2 * x ^ 2 + x ^ 4),
        ## k = 6
        sigma ^ 2 * (mu ^ 5 + mu ^ 4 * x + mu ^ 3 * (14 * sigma ^ 2 + x ^ 2) + mu ^ 2 * (12 * sigma ^ 2 * x + x ^ 3) +
                       mu * (33 * sigma ^ 4 + 9 * sigma ^ 2 * x ^ 2 + x ^ 4) + 15 * sigma ^ 4 * x + 5 * sigma ^ 2 * x ^ 3 + x ^ 5),
        ## k = 7
        sigma ^ 2 * (mu ^ 6 + mu ^ 5 * x + mu ^ 4 * (20 * sigma ^ 2 + x ^ 2) + mu ^ 3 * (18 * sigma ^ 2 * x + x ^ 3) +
                       mu ^ 2 * (87 * sigma ^ 4 + 15 * sigma ^ 2 * x ^ 2 + x ^ 4) +
                       mu * (57 * sigma ^ 4 * x + 11 * sigma ^ 2 * x ^ 3 + x^5) +
                       48 * sigma ^ 6 + 24 * sigma ^ 4 * x ^ 2 + 6 * sigma ^ 2 * x ^ 4 + x ^ 6),
        ## k = 8
        sigma ^ 2 * (mu ^ 7 + mu ^ 6 * x + mu ^ 5 * (27 * sigma ^ 2 + x ^ 2) + mu ^ 4 * (25 * sigma ^ 2 * x + x ^ 3) +
                       mu ^ 3 * (185 * sigma ^ 4 + 22 * sigma ^ 2 * x ^ 2 + x ^ 4) +
                       mu ^ 2 * (141 * sigma ^ 4 * x + 18 * sigma ^ 2 * x^3 + x^5) +
                       mu * (279 * sigma ^ 6 + 87 * sigma ^ 4 * x ^ 2 + 13 * sigma ^ 2 * x ^ 4 + x ^ 6) +
                       105 * sigma ^ 6 * x + 35 * sigma ^ 4 * x ^ 3 + 7 * sigma ^ 2 * x ^ 5 + x ^ 7)
      )
    )
  }
}

get_R_1_mat_6 <- function(mu, sigma){
  cbind(
    ## k = 0
    .5,
    ## k = 1
    .5 * mu,
    ## k = 2
    .5 * (sigma ^ 2 + mu ^ 2),
    ## k = 3
    .5 * mu * (3 * sigma ^ 2 + mu ^ 2),
    ## k = 4
    .5 * (3 * sigma ^ 4 + 6 * mu ^ 2 * sigma ^ 2 + mu ^ 4),
    ## k = 5
    .5 * mu * (mu ^ 4 + 10 * sigma ^ 2 * mu ^ 2 + 15 * sigma ^ 4),
    ## k = 6
    .5 * (mu ^ 6 + 15 * mu ^ 4 * sigma ^ 2 + 45 * mu ^ 2 * sigma ^ 4 + 15 * sigma ^ 6)
  )
}

get_R_2_mat_6 <- function(x, sigma, mu){
  if(x == Inf){
    return(0)
  }else{
    return(
      cbind(
        ## k = 0
        0,
        ## k = 1
        sigma ^ 2,
        ## k = 2
        sigma ^ 2 * (mu + x),
        ## k = 3
        sigma ^ 2 * (2 * sigma ^ 2 + mu ^ 2 + mu * x + x ^ 2),
        ## k = 4
        sigma ^ 2 * (5 * mu * sigma ^ 2 + mu ^ 3 + (3 * sigma ^ 2 + mu ^ 2) * x + mu * x ^ 2 + x ^ 3),
        ## k = 5
        sigma ^ 2 * (mu ^ 4 + x * mu ^ 3 + mu ^ 2 * (9 * sigma ^ 2 + x ^ 2) +
                       mu * (7 * x * sigma ^ 2 + x ^ 3) + 8 * sigma ^ 4 + 4 * sigma ^ 2 * x ^ 2 + x ^ 4),
        ## k = 6
        sigma ^ 2 * (mu ^ 5 + mu ^ 4 * x + mu ^ 3 * (14 * sigma ^ 2 + x ^ 2) + mu ^ 2 * (12 * sigma ^ 2 * x + x ^ 3) +
                       mu * (33 * sigma ^ 4 + 9 * sigma ^ 2 * x ^ 2 + x ^ 4) + 15 * sigma ^ 4 * x + 5 * sigma ^ 2 * x ^ 3 + x ^ 5)
      )
    )
  }
}

