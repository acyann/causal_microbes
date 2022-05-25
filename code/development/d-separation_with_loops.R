# HEADER ======
# 
# Purpose: this script is intended to explore the concept of d-separation
#   in systems containing feedback loops
# In particular, I am curious about whether it is possible to d-separate the
#   variables X and S2 in a system like this:
#     X -> S1 <--> S2
#     where the double-sided arrow represents an ecological interaction feedback
#       loop between S1 and S2
# While double-sided arrows are sometimes used in DAGs, they are generally used
#   to represent hidden common causes:
#     S1 <--> S2 is equivalent to S1 <- u -> S2
#   however, this is likely NOT the case for my purposes
#   in
#     X -> S1 <- u -> S2
#   S1 is a collider on the path between X and S2, meaning that X and S2 are
#      already d-separated, provided that we don't condition on S1
#   I suspect that for a true ecological interaction
#     X -> S1 <--> S2
#   S1 acts as both a collider and a chain/pipe, meaning that it may not be
#     possible to d-separate X and S2 at all

# Approach: because the d-separation rules were developed only for systems that
#   can be represented by acyclic graphs, I don't think we can 100% utilize the
#   d-separation rules for ecological interactions that result in loops
# I'm not interested in fully developing a set of d-separation rules for loops
# Instead, I'm going to simulate a competitive interaction between S1 and S2,
#   letting X be an additional factor in the population growth of S1
# Then I'll look to see how X and S2 are correlated with each other, both
#   marginally and conditional on S1

# Generalized Lotka-Volterra model ======
# I'm going to build on the generalized Lotka-Volterra modeling framework from
# 
#   C. Lehman, S. Loberg, and A. Clark. 2019. Quantitative Ecology: A New
#   Unified Approach. University of Minnesota Libraries Publishing
# 
# Here is there code from section 8.2, with only the following changes:
#   --I changed '=' to the standard R assignment operator '<-'
#   --I changed parameter 'step' to 'step_size' to avoid conflicts with the
#     the stats::step() function
#   --I added some comments

# parameters from section 8.2 (simulates competition)
# initial pop sizes
N1 <- 0.01
N2 <- 0.01
# per capita growth rates
r1 <- 0.5
r2 <- 0.8
# interaction coefficients
s11 <- -0.08
s12 <- -0.03
s21 <- -0.09
s22 <- -0.06

# simulation conditions
t <- 0
dt <- 0.0001
tmax <- 75
step_size <- 0
print( c(t, N1, N2) )

# run the simulation
while(t < tmax) {
  dN1 <- (r1 + s11*N1 + s12*N2)*N1*dt
  dN2 <- (r2 + s21*N1 + s22*N2)*N2*dt
  
  N1 <- N1 + dN1; if(N1 < 0) N1 <- 0
  N2 <- N2 + dN2; if(N2 < 0) N2 <- 0
  
  t <- t + dt; step_size <- step_size + 1
  if(step_size==1000) {
    print( c(t, N1, N2))
    step_size = 0
  }
}

## model notes ======
# 
# Instead of using the closed-form solution to the L-V equations, this approach
#   uses a difference equation with a very small time step (Euler's method)
# The parameters above have a time step of (dt = ) 0.0001, but it only prints
#   out data every (step_size = ) 1000 time steps; the end result is that you
#   wind up getting data in time increments of 0.1 (e.g. 1.1, 1.2, 1.3, etc.)
# The scaling in N1 and N2 indicates that we aren't tracking individuals
#   here, but some fraction there of ... i.e. we are modeling changes
#   for 'individuals per 100 individuals' or something like that; think through
#   how this affects rounding errors; I am getting pop. sizes reported to the
#   sixth decimal place
# The specific interaction parameters lead to stable coexistence
#   This should happen whenever -(r2/s22) > -(r1/s12) & -(r1/s11) > -(r2/s21)
# It is important to calculate the dN1 and dN2 prior to updating the population
#   sizes with the new changes (lines 66 - 70) ... that way we are using the
#   current pop. size of each species to calculate how the interaction works;
#   if we tried to do this all in one step (e.g. combining lines 66 and 69),
#   then N1 would update and this NEW value would be used for dN2
# In this form, the model just prints (every 1000 time steps) to the console,
#   so this will need to be modified to let me store the data

# Saving LV model runs ======
# initial pop sizes
N1 <- 0.01
N2 <- 0.01
# per capita growth rates
r1 <- 0.5
r2 <- 0.8
# interaction coefficients
s11 <- -0.08
s12 <- -0.03
s21 <- -0.09
s22 <- -0.06

# simulation conditions
t <- 0
dt <- 0.0001
tmax <- 75
step_size <- 0
write.step <- 1000
line.out <- tmax/dt/write.step
line.count <- 1

# initialize storage vectors
t.out <- rep(NA, times = line.out)
N1.out <- rep(NA, times = line.out)
N2.out <- rep(NA, times = line.out)

# run the simulation
while(t < tmax) {
  dN1 <- (r1 + s11*N1 + s12*N2)*N1*dt
  dN2 <- (r2 + s21*N1 + s22*N2)*N2*dt
  
  N1 <- N1 + dN1; if(N1 < 0) N1 <- 0
  N2 <- N2 + dN2; if(N2 < 0) N2 <- 0
  
  t <- t + dt; step_size <- step_size + 1
  if(step_size==write.step) {
    t.out[line.count] <- t
    N1.out[line.count] <- N1
    N2.out[line.count] <- N2
    line.count <- line.count + 1
    step_size = 0
  }
}

# data are now stored in individual vectors: t.out, N1.out, N2.out
plot(NULL, xlim = c(0, 75), xlab = "time", ylim = c(0, 12), ylab = "pop. size")
lines(x = t.out, y = N2.out, col = "orange", lwd = 2)
lines(x = t.out, y = N1.out, col = "blue", lwd = 2)
# this matches Lehman et al. 2019, Fig. 8.2 