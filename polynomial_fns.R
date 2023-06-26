# evaluate polynomial with given coefficients
f_poly <- function(x, poly_coefs){
  k <- length(poly_coefs) - 1
  out <- 0
  for (i in 0:k){
    out <-  out + poly_coefs[i + 1] * x ^ i
  }
  return(out)
}

# get coefs of derivative
poly_deriv <- function(poly_coefs){
  k <- length(poly_coefs) - 1
  if (k==0){
    return(0)
  }else{
    return(poly_coefs[-1] * (1:k))
  }
}

# multiply 2 polynomials with coefficients coefs1 and coefs2
# return coefficients of product polynomial
poly_mult <- function(coefs1, coefs2){
  K1 <- length(coefs1)
  K2 <- length(coefs2)
  
  coefs_prod <- numeric(K1 + K2 - 1)
  
  for (k1 in 1:length(coefs1)){
    for (k2 in 1:length(coefs2)){
      coefs_prod[k1 + k2 - 1] <- coefs_prod[k1 + k2 - 1] + coefs1[k1] * coefs2[k2]
    }
  }
  
  return(coefs_prod)
}

#add 2 polynomials with coefficients coefs1 and coefs2
# return coefficients of sum polynomial
poly_add <- function(coefs1, coefs2){
  K1 <- length(coefs1)
  K2 <- length(coefs2)
  K <- max(K1, K2)
  coefs_add <- ifelse(1:K <= K1, coefs1, 0) + ifelse(1:K <= K2, coefs2, 0)
  return(coefs_add)
}

# get expanded polynomial coefficients for the expression (x + c) ^ n
expand_coefs <- function(c, n){
  vapply(0:n, function(k) choose(n, k) * c ^ (n - k), 0)
}

# evaluate taylor approximation of given function f, at a, and with speficied degree
# given list of derivatives
taylor_approx <- function(x, f, f_deriv_list, a, degree){
  retValue <- f(a)
  
  for (d in 1:degree){
    retValue <- retValue + f_deriv_list[[d]](a) * (x - a) ^ d / factorial(d)
  }
  
  return(retValue)
}

# get polynomial coeffiecients for taylor approx
get_taylor_approx_coefs <- function(f, f_deriv_list, a, degree){
  taylor_approx_coefs <- f(a)
  
  for (d in 1:degree){
    taylor_approx_coefs <- poly_add(taylor_approx_coefs, f_deriv_list[[d]](a) * expand_coefs(c = -a, n = d) / factorial(d))
  }
  
  return(taylor_approx_coefs)
}

# evaluate at x cubic weight function between a and b
w_cubic <- function(x, a, b){
  (x < a) * 0 +
    (x >= a & x < b) * f_poly((x - a)/(b - a), c(0,0,3,-2)) +
    (x >= b) * 1
}

# mix 2 functions using specified weight function w (default to cubic)
mix_2_fns <- function(x, f_1, f_2, w=w_cubic, a, b){
  (1- w(x, a, b)) * f_1(x) + w(x, a, b) * f_2(x)
}

# get polynomial coeffiecients for cubic weight function between a and b
get_cubic_weight_coefs <- function(a, b){
  poly_add(3 * expand_coefs(c = -a, n = 2) / (b - a) ^ 2, -2 * expand_coefs(c = -a, n = 3) / (b - a) ^ 3)
}

# get polynomial coeffiecients for linear weight function between a and b
get_linear_weight_coefs <- function(a, b){
  c(-a / (b - a), 1 / (b - a))
}

# get polynomial coeffiecients for mix using specified weight function
mix_2_poly <- function(coefs1, coefs2, a, b, get_weight_fn){
  w <- get_weight_fn(a, b)
  poly_add(
    poly_mult(poly_add(1, -w), coefs1),
    poly_mult(w, coefs2)
  )
}

#curve(f_poly(x, mix_2_poly(c(0,1), c(1,-1), 0, 1, get_cubic_weight_coefs)))