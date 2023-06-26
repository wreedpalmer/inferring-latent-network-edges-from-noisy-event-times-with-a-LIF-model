
### functions for evaluating the f(x) = log(K * (x - 1)) and its derivatives
# inverse logit curve for soft threshold
st_inv_logit <- function(x, K=50){
  gtools::inv.logit(K * (x - 1))
}

# default to K = 50
f <- function(x, K=50){
  log(1 + exp(K * (x - 1)))
}

# first derivative
f_d1 <- function(x, K=50){
  K * exp(K * (x - 1)) / (1 + exp(K * (x - 1)))
}

# second derivative
f_d2 <- function(x, K=50){
  K ^ 2 * exp(K * (x - 1)) / (1 + exp(K * (x - 1))) ^ 2
}

# third derivative
f_d3 <- function(x, K=50){
  -K ^ 3 * exp(K * (x - 1)) * (exp(K * (x - 1)) - 1) / (1 + exp(K * (x - 1))) ^ 3
}

# fourth derivative
f_d4 <- function(x, K=50){
  K ^ 4 * exp(K * x - K) * (-4 * exp(K * (x - 1)) + exp(2 * (K * (x - 1))) + 1) / (1 + exp(K * (x - 1))) ^ 4
}

# list containing first 4 derivative functions of f
f_d_list <- list(f_d1, f_d2, f_d3, f_d4)

par(mfrow = c(2, 2))
curve(f(x, K=50), from=.8, to=1.2, main = "f", xlab = "", ylab = "")
curve(f_d1(x, K=50), from=.8, to=1.2, main = "f'", xlab = "", ylab = "")
curve(f_d2(x, K=50), from=.8, to=1.2, main = "f''", xlab = "", ylab = "")
curve(f_d3(x, K=50), from=.8, to=1.2, main = "f'''", xlab = "", ylab = "")
par(mfrow = c(1, 1))

# K1 = 50, K2 = 50
# c1, c2 set  at roots of 4th derivative
# where 'curvature' of f is highest -- values at which |f'''(x)| are maximized
c_1 <- uniroot(function(x) f_d4(x,  K=50), c(.9, 1))$root
c_2 <- uniroot(function(x) f_d4(x,  K=50), c(1, 1.1))$root

# c3, c4 set where second derivative of taylor approx around c_1, c_2, resp., is zero
# want to limit those approximations to the interval where they are convex.
c_3 <- (c_1 * f_d3(c_1) - f_d2(c_1))/f_d3(c_1)
c_4 <- (c_2 * f_d3(c_2) - f_d2(c_2))/f_d3(c_2)

# set the lower knots c_5 and c_7 sequentially at the points where the 2nd degree approximations stop being increasing functions
#c_5 <- (c_3 * f_d2(c_3) - f_d1(c_3))/f_d2(c_3)
#c_7 <- (c_5 * f_d2(c_5) - f_d1(c_5))/f_d2(c_5)

# set the lower knots c_5 and c_7 sequentially at the points where the 3nd degree approximations stop being convex
c_5 <- (c_3 * f_d3(c_3) - f_d2(c_3))/f_d3(c_3)
c_7 <- (c_5 * f_d3(c_5) - f_d2(c_5))/f_d3(c_5)

#set higher knots c_6 and c_8 at points where the 2nd degree approximation's 1st derivative creeps above 50
#c_6 <- (c_4 * f_d2(c_4) - f_d1(c_4) + 50)/f_d2(c_4)
#c_8 <- (c_6 * f_d2(c_6) - f_d1(c_6) + 50)/f_d2(c_6)

# set higher knots c_6 and c_8 at points where the 3nd degree approximations stop being convex
c_6 <- (c_4 * f_d3(c_4) - f_d2(c_4))/f_d3(c_4)
c_8 <- (c_6 * f_d3(c_6) - f_d2(c_6))/f_d3(c_6)

c_9 <- c_7 - f(c_7)/f_d1(c_7)
c_10 <- (50 - c_8 * f_d1(c_8) + f(c_8))/(50 - f_d1(c_8))

# here are the knots for our piecewise approximation
knots <- c(c_9, c_7, c_5, c_3, c_1, 1, c_2, c_4, c_6, c_8, c_10)

par(mfrow = c(2, 2))
curve(f(x), from=.8, to=1.2, main = "f", xlab = "", ylab = "")
points(knots, f(knots))
curve(f_d1(x), from=.8, to=1.2, main = "f'", xlab = "", ylab = "")
points(knots, f_d1(knots))
curve(f_d2(x), from=.8, to=1.2, main = "f''", xlab = "", ylab = "")
points(knots, f_d2(knots))
curve(f_d3(x), from=.8, to=1.2, main = "f'''", xlab = "", ylab = "")
points(knots, f_d3(knots))
par(mfrow = c(1, 1))

# a list giving the intervals where the approximation is non-zero
# along with the polynomial coefficients for the approximation on that interval
piecewise_approx_list <- list(
  lower_3 = list(interval=c(c_7, c_5)),
  lower_2 = list(interval=c(c_5, c_3)),
  lower_1 = list(interval=c(c_3, c_1)),
  middle = list(interval=c(c_1, c_2)),
  upper_1 = list(interval=c(c_2, c_4)),
  upper_2 = list(interval=c(c_4, c_6)),
  upper_3 = list(interval=c(c_6, c_8)),
  upper_4 = list(interval=c(c_8, Inf))
)

piecewise_approx_list <- list(
  lower_4 = list(interval=c(c_9, c_7)),
  lower_3 = list(interval=c(c_7, c_5)),
  lower_2 = list(interval=c(c_5, c_3)),
  lower_1 = list(interval=c(c_3, c_1)),
  middle1 = list(interval=c(c_1, 1)),
  middle2 = list(interval=c(1, c_2)),
  upper_1 = list(interval=c(c_2, c_4)),
  upper_2 = list(interval=c(c_4, c_6)),
  upper_3 = list(interval=c(c_6, c_8)),
  upper_4 = list(interval=c(c_8, c_10)),
  upper_5 = list(interval=c(c_10, Inf))
)

### getting the polynomial coefficients for the piecewise approximation

## set weight function for mixing polynomials
#get_weight_fn <- get_linear_weight_coefs
get_weight_fn <- get_cubic_weight_coefs

# for between c_1 and c_2, middle knots
# mix two third degree taylor polynomials, centered at c_1 and c_2, resp.
#piecewise_approx_list$middle$coefs <- mix_2_poly(
#  get_taylor_approx_coefs(f, f_d_list, c_1, 3),
#  get_taylor_approx_coefs(f, f_d_list, c_2, 3),
#  c_1, c_2,
#  get_weight_fn
#)

piecewise_approx_list$middle1$coefs <- mix_2_poly(
  get_taylor_approx_coefs(f, f_d_list, c_1, 3),
  get_taylor_approx_coefs(f, f_d_list, 1, 3),
  c_1, 1,
  get_weight_fn
)
piecewise_approx_list$middle2$coefs <- mix_2_poly(
  get_taylor_approx_coefs(f, f_d_list, 1, 3),
  get_taylor_approx_coefs(f, f_d_list, c_2, 3),
  1, c_2,
  get_weight_fn
)

# for between c_2 and c_4, first upper section
# mix third degree taylor polynomial centered at c_2 with second degree taylor polynomials centered at c_4
piecewise_approx_list$upper_1$coefs <- mix_2_poly(
  get_taylor_approx_coefs(f, f_d_list, c_2, 3),
  get_taylor_approx_coefs(f, f_d_list, c_4, 3), #THIRD
  #c(-50, 50),
  c_2, c_4,
  get_weight_fn
)

# for between c_2 and c_4, first upper section
# mix second degree taylor polynomial centered at c_4 with second degree taylor polynomials centered at c_6
piecewise_approx_list$upper_2$coefs <- mix_2_poly(
  get_taylor_approx_coefs(f, f_d_list, c_4, 3), ##3
  get_taylor_approx_coefs(f, f_d_list, c_6, 3), ##3
  #c(-50, 50),
  c_4, c_6,
  get_weight_fn
)

# for between c_6 and c_8, first upper section
# mix second degree taylor polynomials centered at c_6 with the limiting line 50x - 50
piecewise_approx_list$upper_3$coefs <- mix_2_poly(
  get_taylor_approx_coefs(f, f_d_list, c_6, 3), ##3
  #c(-50, 50),
  get_taylor_approx_coefs(f, f_d_list, c_8, 3), ##3
  c_6, c_8,
  get_weight_fn
)

piecewise_approx_list$upper_4$coefs <- get_taylor_approx_coefs(f, f_d_list, c_8, 1)

# for c_8 and higher set to the limiting line 50x - 50
piecewise_approx_list$upper_5$coefs <- c(-50, 50)

# for between c_3 and c_1, first lower section
# mix second degree taylor polynomials centered at c_3 with third degree taylor polynomials centered at c_1
piecewise_approx_list$lower_1$coefs <- mix_2_poly(
  get_taylor_approx_coefs(f, f_d_list, c_3, 3),
  get_taylor_approx_coefs(f, f_d_list, c_1, 3),
  c_3, c_1,
  get_weight_fn
)

# for between c_5 and c_3, second lower section
# mix second degree taylor polynomials centered at c_5 with second degree taylor polynomials centered at c_3
piecewise_approx_list$lower_2$coefs <- mix_2_poly(
  get_taylor_approx_coefs(f, f_d_list, c_5, 3),
  get_taylor_approx_coefs(f, f_d_list, c_3, 3),
  c_5, c_3,
  get_weight_fn
)

# for between c_7 and c_5, third lower section
# mix limiting constant 0 with second degree taylor polynomials centered at c_5
piecewise_approx_list$lower_3$coefs <- mix_2_poly(
  #0,
  get_taylor_approx_coefs(f, f_d_list, c_7, 3),
  get_taylor_approx_coefs(f, f_d_list, c_5, 3),
  c_7, c_5,
  get_weight_fn
)

piecewise_approx_list$lower_4$coefs <- get_taylor_approx_coefs(f, f_d_list, c_7, 1)

curve(f(x, K=50), from=.85, to=1.15, main = "f", xlab = "", ylab = "")
curve(0*x, from=0, to=c_9, col="red", add=T)
for (segment in piecewise_approx_list){
  if (segment$interval[2] == Inf){
    b=10
  }else{
    b=segment$interval[2]
  }
  curve(f_poly(x, segment$coefs), from=segment$interval[1], to=b, col="red", add=T)
}

points(knots, f(knots))

curve(f_d1(x, K=50), from=.8, to=1.2, main = "f", xlab = "", ylab = "")#, ylim = c(-400,800))
for (segment in piecewise_approx_list){
  if (segment$interval[2] == Inf){
    b=10
  }else{
    b=segment$interval[2]
  }
  curve(f_poly(x, segment$coefs %>% poly_deriv), from=segment$interval[1], to=b, col="green", add=T)
}

points(knots, f_d1(knots))

plot(knots, rep(0,length(knots)), ylim=c(-.01,.01), xlim=c(.8,1.2))
curve(-f(x), from=0, to=min(knots), col="green", add=T)
for (segment in piecewise_approx_list){
  if (segment$interval[2] == Inf){
    b=10
  }else{
    b=segment$interval[2]
  }
  curve(f_poly(x, segment$coefs) - f(x), from=segment$interval[1], to=b, col="green", add=T, n=1001)
}

