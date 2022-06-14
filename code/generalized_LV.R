generalized_LV <- function(N, r, s, dt = 0.0001, tmax = 75) {
  # simulation conditions
  t <- 0
  step_count <- 0
  write.step <- 1/dt/10
  line.out <- tmax/dt/write.step
  line.count <- 1
  
  # initialize storage vectors
  t.out <- rep(NA, times = line.out)
  N1.out <- rep(NA, times = line.out)
  N2.out <- rep(NA, times = line.out)
  
  # run the simulation
  N1 <- N[1]
  N2 <- N[2]
  
  while(t < tmax) {
    dN1 <- (r[1] + s[1,1]*N1 + s[1,2]*N2)*N1*dt
    dN2 <- (r[2] + s[2,1]*N1 + s[2,2]*N2)*N2*dt
    
    N1 <- N1 + dN1; if(N1 < 0) N1 <- 0
    N2 <- N2 + dN2; if(N2 < 0) N2 <- 0
    
    t <- t + dt; step_count <- step_count + 1
    if(step_count==write.step) {
      t.out[line.count] <- t
      N1.out[line.count] <- N1
      N2.out[line.count] <- N2
      line.count <- line.count + 1
      step_count = 0
    }
  }
  OUT <- cbind(t.out, N1.out, N2.out)
  colnames(OUT) <- c("t", "N1", "N2")
  return(OUT)
}

# This function approximates a time series of two-species' population growth
#   based on the generalized Lotka-Volterra equations. It is inspired by
# 
#   C. Lehman, S. Loberg, and A. Clark. 2019. Quantitative Ecology: A New
#   Unified Approach. University of Minnesota Libraries Publishing
#
#   who utilize Euler's method of very-small-time-step difference equations
#   to approximate the differential equations
# 
# The function works with a time step of dt (=0.0001 by default), but it records
#   data on a coarser scale (t = 0.1, 0.2, 0.3, etc. this is hard-coded in
#   the function)
# 
# The function requires the following arguments:
#   N: a length=2 vector of starting population sizes
#   r: a length=2 vector of per-capita growth rates
#   s: a 2x2 matrix of interaction coefficients of the form
#      [s11  s12]
#      [s21  s22]