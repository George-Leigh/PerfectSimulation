# -*- Text -*-
##############################################################################
# R code for demonstration of perfect simulation from unbiased simulation with
#   a standard normal distribution: This script uses the algorithm from the
#   paper, with maximal coupling using a uniform spherical jump distribution.
# This is random-walk MCMC.  We'll see how many dimensions we can do.  It may
#   fall over at a few dimensions.  [It actually works passably at d = 20, but
#   the probability that coalescence occurs within a particular number of
#   blocks shows signs of not decaying geometrically by d = 20.]
# Hamiltonian Monte Carlo, as far as I'm aware, has no limit to the number of
#   dimensions for which it will produce perfect simulation.  That is the
#   subject of another paper.

# Author: George Leigh, 202212, updated 202307
# Run under R 4.2.2, Linux version
# Copyright 2023 George Leigh
# Associated journal article: Leigh, G. M., Yang, W-H, Wickens, M. and
#   Northrop, A. R. (2023).  "Perfect simulation from unbiased simulation".
# Citation of the journal article in any publication that uses the algorithms
#   published in it would be greatly appreciated by the authors.

# The programs in this project are free software: you can redistribute them
#   and/or modify them under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 2 of the
#   License, or (at your option) any later version.
# The programs are distributed in the hope that they will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
# You should have received a copy of the GNU General Public License along with
#   this program.  If not, see <https://www.gnu.org/licenses/>.

# Apart from reading in the functions, this script is intended to be run
#   interactively, not by calling the "source" function.

##############################################################################
# We will simulate a standard normal distribution in d dimensions.  There are
#   two parameters in the simulation model: sigma, the standard deviation of
#   the jumps; and r, the radius of the solid spherical uniform distribution
#   used for maximal coupling.
# The code follows the algorithms in the paper as closely as feasible.  The
#   major difference is that the code carries along logical variables
#   indicating whether coalescence has been attained, so that repeated testing
#   of whether chains take equal values is not necessary.

######################################## Worker functions
Euclidean = function(x1, x2) sqrt(sum((x1 - x2)^2)) # Euclidean metric
Norm = function(x) sqrt(sum(x^2)) # Euclidean norm

#################### Function to find the element with minimum distance
MinInd = function(x, i1, i2) {
 # x: matrix with d columns (note: must be a matrix, even when d == 1)
 # i1: index of first row to consider
 # i2: index of row for which we are finding the closest row, of rows i1:(i2
 #   - 1); must have 1 <= i1 < i2.
 # First derive a version of x[i1:(i2 - 1),] which we can be sure is a matrix.
 d = ncol(x)
 if (d == 1) {
  xmat = cbind(x[i1:(i2 - 1),])
 } else {
  xmat = rbind(x[i1:(i2 - 1),])
 }
 Dist = apply(xmat, 1, Euclidean, x[i2,])
 which.min(Dist) + i1 - 1 # Index of row with minimum distance, starting at
 #   index 1 = row 1 of x (not row i1)
}

#################### Function to generate a jump for maximal coupling, using
#   the solid spherical uniform distribution
MaxCoupleJump = function(X, r, Rdir, Rmag) {
 # d: Number of dimensions (= length of X below)
 # X: Jump origin
 # r: Radius of solid spherical uniform distribution for jump
 # Rdir: Vector of length d of N(0, 1) variables for jump direction
 # Rmag: Uniform(0, 1) variable for jump magnitude
 d = length(X)
 uhat = Rdir / Norm(Rdir) # Unit vector
 JumpMag = r * Rmag^(1 / d) # Jump magnitude, ensures constant p.d.f. over the
 #   d-dimensional hypersphere.
 X + JumpMag * uhat
}

#################### Function to maximally couple a point Y with a previously
#   jumped point X which had jump destination Xstar
MaxCouple = function(X, Xstar, Y, r) {
 # X: Jump origin for X (the chain that jumps freely)
 # Xstar: Jump destination for X
 # Y: Jump origin for Y (the chain that we are maximally coupling with X)
 # r: Radius of solid spherical uniform distribution for jump
 if (Norm(Y - Xstar) <= r) {
  Ystar = Xstar # Successfully coupled
  lCoupled = TRUE # Mark as successfully coupled.
 } else {
  uX = Xstar - X # Jump vector for X
  L = Norm(Y - X) / 2 # Half distance from X to Y
  v = (Y - X) / (2 * L) # Unit vector in direction from X to Y
  c = L * v + uX - sum(uX * v) * v
  normc = Norm(c)
  if (normc < r) { # Need to translate Ystar to make sure it is not in the
   #   overlap zone.
   w = -L + sqrt(L^2 + r^2 - normc^2) # Distance from point C to D in figure
   uY = uX + 2 * w * v # Do the translation of Ystar (uY != uX).
   Ystar = Y + uY
  } else { # normc >= r: no translation of Ystar is needed.
   Ystar = Y + uX
  }
  lCoupled = FALSE # Not successfully coupled (applies to both cases above)
 }
 list(Ystar = Ystar, lCoupled = lCoupled) # Record whether successfully
 #   coupled.  Note that the coupling still won't be completed successfully
 #   unless both Metropolis-Hastings tests, to move from X to Xstar and Y to
 #   Ystar, are successful.
}

#################### Function to conduct a Metropolis-Hastings test
MHTest = function(X, Xstar, RMH, U) { # Metropolis-Hastings test
 # X: Jump origin
 # Xstar: Proposed jump destination
 # RMH: Uniform(0, 1) random number
 # U: Negative log-likelihood function
 if (RMH <= exp(U(X) - U(Xstar))) {
  Xdest = Xstar # Destination = proposal
  lJumped = TRUE # Successful jump
 } else {
  Xdest = X # Destination = origin
  lJumped = FALSE # Unsuccessful jump
 }
 list(Dest = Xdest, lJumped = lJumped)
}

#################### Standard normal negative log-likelihood
Unorm = function(X) # Negative log-likelihood for standard normal distribution
 0.5 * sum(X^2) # Scale factor for distribution is omitted.

#################### Function to generate starting points for MCMC
Start = function(K, d) { # Generate starting points (K x d matrix).  We will
 #   sample each coordinate independently from a uniform distribution on [-6,
 #   6].
 xlim = 6
 X = array(runif(K * d, -xlim, xlim), dim = c(K, d))
 X
}

#################### Function to generate random numbers for MCMC
Rand = function(d, M) # Generate random numbers for M MCMC iterations and
 #   one maximal coupling step.
 list(
  RMCMC = list(Rjump = array(rnorm(M * d), dim = c(M, d)), # MCMC jump vectors
   RMH = runif(M)), # MCMC Metropolis-Hastings test random variables
  Rdir = rnorm(d), # Maximal coupling direction
  Rmag = runif(1), # Maximal coupling magnitude
  RMH = runif(1)) # Maximal coupling Metropolis-Hastings test random variable

#################### Function to conduct MCMC: this has an inbuilt
#   Metropolis-Hastings algorithm, separate to the one used for maximal
#   coupling.
Mcmc = function(X, M, RMCMC, sigma, U) { # Conduct M MCMC iterations starting
 #   from X and using random numbers in RMCMC.
 # X: Starting point, also becomes finishing point
 # M: Number of iterations
 # RMCMC: Random numbers to use
 # sigma: standard deviation of jump
 # U: Negative log-likelihood function
 d = length(X) # Number of dimensions
 sigmascal = sigma / sqrt(d) # Scaled version of sigma, to account for d
 for (iit in 1:M) {
  Xnew = X + sigmascal * RMCMC$Rjump[iit,]
  if (RMCMC$RMH[iit] <= exp(U(X) - U(Xnew)))
   X = Xnew
 } # iit
 X
} # Mcmc

#################### Function to generate one sample set of simulations
# Note: This function assumes that coalescence will occur within K - 1 blocks.
#   It does not include the extra work to handle chains that don't coalesce.
SampleSet = function(K, B, M, d, sigma, r, Seed) { # Main function to run the
 #   algorithm, generates K slightly correlated perfect samples, if an error
 #   doesn't occur.
 # K: Number of blocks, also equal to number of chains
 # B: Number of MCMC iterations in a block
 # M: Number of MCMC iterations per maximal coupling step (should be a divisor
 #    of B)
 # d: Number of dimensions
 # sigma: Standard deviation of jump in MCMC
 # r: Radius of uniform spherical distribution for maximal coupling
 # Seed: Random number seed, for reproducibility
 if (Seed > 0) set.seed(Seed) # Zero implies don't reset seed.
 J = ceiling(B / M) # Number of maximal coupling steps per block
 B = J * M # Redefine B if not a multiple of M.
 Q = Start(K, d) # Generate starting points (K x d matrix) along diagonal of
 #   matrix in Fig. 2 of paper
 Qminus = Q # To store version of Q prior to jumping for maximal coupling
 Qstar = Q # To store jumped version of Q, prior to Metropolis-Hastings test
 Qjumped = rep(FALSE, K) # Records whether maximal-coupling M-H test is
 #    successful for Q.
 Q1 = array(0, dim = c(K, J, d)) # To store results after every maximal
 #    coupling step in the first row of the matrix in Fig. 2 of paper
 Q1minus = Q1 # Version of Q1 prior to jumping for maximal coupling
 Q1star = Q1 # Jumped version of Q1, prior to Metropolis-Hastings test
 Q1jumped = array(FALSE, dim = c(K, J)) # Records whether M-H test is
 #    successful for Q1.
 s = .Random.seed # Save current seed, to recover later.
 Error = FALSE
 BlocksCoal = rep(0, K) # Number of blocks taken to coalesce with the chain
 #   below it (zero mean not coalesced, which results in Error being true)
 a = rep(0, K) # Mark all chains as active (not coalesced with any other
 #   chain).
 # We'll store the whole Q matrix (whole of matrix from Fig. 2 of paper), for
 #   diagnostic purposes only.  It is not necessary for the operation of the
 #   algorithm, but without it a return value of true for Error won't provide
 #   any description of what has gone wrong.  If this storage is not desired,
 #   all lines of code containing "Qarray" can be removed or commented out.
 #Qarray = array(0, dim = c(K, K, d))

 for (j in 1:K) { # Loop over blocks (columns) in Fig. 2 upper triangle.
  for (j1 in 1:J) {
   R = Rand(d, M) # Generate random numbers for MCMC and maximal coupling.
   for (i in 1:j) { # Loop over chains (rows) in Fig. 2 upper triangle.
    if (a[i] > 0) { # Coalesced to a previous row, so we just copy that.
     Q[i,] = Q[a[i],]
    } else { # Not coalesced, so we run MCMC.
     Qminus[i,] = Mcmc(Q[i,], M, R$RMCMC, sigma, Unorm)
     if (i == 1) { # Row 1 of Fig. 2: Jump freely, M-H test and store.
      Qstar[1,] = MaxCoupleJump(Qminus[1,], r, R$Rdir, R$Rmag)
      MHResult = MHTest(Qminus[1,], Qstar[1,], R$RMH, Unorm)
      #print(MHResult)
      #print(Q)
      Q[1,] = MHResult$Dest
      Qjumped[1] = MHResult$lJumped
      Q1minus[j, j1,] = Qminus[1,]
      Q1star[j, j1,] = Qstar[1,]
      Q1[j, j1,] = Q[1,]
      Q1jumped[j, j1] = Qjumped[1]
     } else { # i > 1 in upper triangle of Fig. 2: Maximally couple with the
      #   closest previous row (m < i).
      m = MinInd(Qminus, 1, i)
      MaxCoupleResult = MaxCouple(Qminus[m,], Qstar[m,], Qminus[i,], r)
      Qstar[i,] = MaxCoupleResult$Ystar
      lCoupled = MaxCoupleResult$lCoupled
      MHResult = MHTest(Qminus[i,], Qstar[i,], R$RMH, Unorm)
      Q[i,] = MHResult$Dest
      Qjumped[i] = MHResult$lJumped
      # We don't need to test equality of Q[i,] with Q[m,].  We use the three
      #   logical variables instead.
      if (lCoupled & Qjumped[m] & Qjumped[i]) {
       if (a[m] > 0)
        m = a[m] # Reset m to the earliest chain with which chain i has
        #   coalesced.  This is the row with which row m coalesced, if any.
       a[i] = m # Mark row i as coalesced with the said row, so we can just
       #   copy the rest of row i from the previous row a[i].
       # Mark any lower rows that earlier coalesced with row i as having
       #   instead coalesced with row m (using the new value of m).
       a[a == i] = m
      } # Coalesced
     } # i > 1
    } # Not coalesced: may or may not be coalesced now!
   } # i
  } # j1
  # Extra loop to set BlocksCoal
  if (j >= 2)
   for (i in 1:(j - 1))
    if (BlocksCoal[i] == 0 & a[i + 1] > 0 & (a[i + 1] == i | a[i + 1] == a[i]))
     BlocksCoal[i] = j - i
  #Qarray[1:j, j,] = Q[1:j,] # Store upper triangle of matrix, diagnostic only.
 } # j

 .Random.seed <<- s # Restore seed to regenerate random numbers used for upper
 #   triangle.  We want to use the same random numbers in the part of each
 #   column j that is in the lower triangle of Fig. 2.
 for (j in 1:(K - 1)) { # Loop over blocks (columns) in lower triangle
  if (a[j + 1] != j) cat("Error", j, Q[j:(j + 1), 1], "\n")
  Error = Error | a[j + 1] != j # Error if any row j has not coalesced with
  #   the one below it after running for K blocks since it was initialised at
  #   the diagonal element of the matrix.  These are processes X and Y
  #   respectively in the theory of unbiased sampling.  It will actually be OK
  #   (no holes generated) if it coalesces in one more block, using new random
  #   numbers, but we're flagging it as an error here to keep the algorithm
  #   fairly simple.
  # Transfer all coalescence involving row j to coalescence with row j + 1,
  #   seeing as we have now finished with row j.  Row j has cycled all the way
  #   around to where it was initialised.  In transferring the coalescence, we
  #   are assuming that rows j and j + 1 have coalesced.  If they haven't, the
  #   variable "Error" will be set to true above.
  a[j + 1] = a[j]
  a[a == j] = j + 1
  for (j1 in 1:J) {
   # Reload saved copy of row 1 (all three versions); no need to recompute.
   Qminus[1,] = Q1minus[j, j1,]
   Qstar[1,] = Q1star[j, j1,]
   Q[1,] = Q1[j, j1,]
   Qjumped[1] = Q1jumped[j, j1]
   R = Rand(d, M) # Generate random numbers for MCMC and maximal coupling.
   for (i in (j + 1):K) {
    if (a[i] > 0) {
     Q[i,] = Q[a[i],] # Copy from chain with which row i has coalesced.
    } else { # No coalescence: run MCMC.
     Qminus[i,] = Mcmc(Q[i,], M, R$RMCMC, sigma, Unorm) # MCMC creates
     #   pre-jump version of Q.
     # Set the row with which we will maximally couple row i.  We are using
     #   rows 1 and j + 1, j + 2, ..., K.  The rows inbetween have cycled
     #   right around, back to their points of initialisation on the diagonal
     #   of Fig. 2, and we have finished with them.
     if (i == j + 1) {
      m = 1 # Only row j + 1 is maximally coupled with row 1.
     } else { # i > j + 1
      m = MinInd(Qminus, j + 1, i) # Choose the closest row from rows j + 1,
      #   ..., i - 1.
     }
     MaxCoupleResult = MaxCouple(Qminus[m,], Qstar[m,], Qminus[i,], r)
     Qstar[i,] = MaxCoupleResult$Ystar
     lCoupled = MaxCoupleResult$lCoupled
     MHResult = MHTest(Qminus[i,], Qstar[i,], R$RMH, Unorm)
     Q[i,] = MHResult$Dest
     Qjumped[i] = MHResult$lJumped
     # Again we don't need to test equality Q[i,] with Q[m,].  We use the three
     #   logical variables generated by the functions we've called.
     if (lCoupled & Qjumped[m] & Qjumped[i]) {
      if (a[m] > 1)
       m = a[m] # Reset m to earliest row with which it coalesced, if any,
       #   except that only row j + 1 is allowed to coalesce with row m.
      a[i] = m # Mark row i as coalesced with the said row.  Having marked row
      #   i as coalesced, we can just copy the rest of row i from the previous
      #   row a[i], until row a[i] finishes and cycles back to the diagonal of
      #   Fig. 2 where it was initialised.  After that, the coalescence will
      #   be transferred to row a[i] + 1 as above.
      # Mark any lower rows that earlier coalesced with row i as having
      #   instead coalesced with row a[i], if a[i] > 1.  As stated above, when
      #   i > j + 1 we prefer coalescence with row j + 1 over row 1.
      if (m > 1) a[a == i] = m
     } # Coalesced
    } # Not coalesced: may or may not be coalesced now!
   } # i
  } # j1
  # Extra loop to set BlocksCoal
  if (j + 1 < K) # Chain K is done separately below.
   for (i in (j + 1):(K - 1))
    if (BlocksCoal[i] == 0 & a[i + 1] > 0 & (a[i + 1] == i | a[i + 1] == a[i]))
     BlocksCoal[i] = j - i + K
  # The following line is executed for all values of j (1, ..., K - 1): i = K.
  # We are meant to record the block on which chain K coalesces with chain 1,
  #   but there's a trick!  We don't care so much about chain K, because it is
  #   guaranteed to coalesce at some stage with the one-off diagonal chain
  #   (i.e., one of the chains immediately below the diagonal), provided that
  #   the coalescence "Error" is not set for any of the chains 1, ..., K - 1.
  #   A more accurate measure of when chain K coalesces with chain 1 is
  #   actually when the one-off-diagonal chain coalesces with chain 1.  Given
  #   that chains j + 1, ..., K - 1 are in between chain K and chain 1, the
  #   block on which chain K coalesces with chain 1 would be an overestimate
  #   (overly pessimistic).
  if (BlocksCoal[K] == 0 & a[j + 1] == 1)
   BlocksCoal[K] = j
  #Qarray[(j + 1):K, j,] = Q[(j + 1):K,] # Store lower triangle of matrix, as
  #   diagnostic only.
 } # j
 if (a[K] != 1) cat("Error", K, Q[c(K, 1)], "\n")
 Error = Error | a[K] != 1 # Error if row K has not coalesced with row 1 after
 #   running for K blocks and ending up at position (K, K - 1) of the matrix
 #   in Fig. 2.  Now rows K and 1 are processes X and Y respectively in the
 #   theory of unbiased sampling.  Again it will be OK (no holes generated) if
 #   rows K and 1 coalesce in one more block, but we keep it simple and flag
 #   it as an error here.
 Q[1,] = Q1[K, J,] # Restore final value in first row of matrix, after we
 #    overwrote it.
 #list(Q = Q, Error = Error, Qarray = Qarray)
 list(Q = Q, Coal = BlocksCoal, Error = Error)
} # SampleSet

#################### Function to conduct perfect simulation using the unbiased
#   simulation algorithm, simulating a collection of independent sample sets
Unbiased = function(N, K, B, M, d, sigma, r, Seed = 1, nProg = 100) {
 # N: Total number of sample points, should be a multiple of K below
 # K: Number of blocks, also = number of chains, should be a divisor of N
 # B: Number of MCMC iterations in a block
 # M: Number of MCMC iterations per maximal coupling step, should be a divisor
 #    of B
 # d: Number of dimensions
 # sigma: Standard deviation of jump in MCMC
 # r: Radius of uniform spherical distribution for maximal coupling
 # Seed: Random number seed, for reproducibility
 # nProg: Number of progress records to display
 if (Seed > 0) set.seed(Seed) # Zero implies don't reset seed.
 NSet = ceiling(N / K) # Number of sample sets to simulate
 N = NSet * K # Redefine N if not a multiple of K.
 Qmat = array(0, dim = c(N, d)) # Matrix to hold final simulations
 CoalVec = rep(0, N) # Vector to hold number of blocks taken to reach
 #   coalescence for each sample point
 ErrorVec = rep(FALSE, NSet) # Vector to hold error indicator for each sample
 #   set
 iProg = 0 # Number of progress records written so far
 for (iSet in 1:NSet) {
  Set = SampleSet(K, B, M, d, sigma, r, Seed = 0) # The simulation
  Ind = (iSet - 1) * K + (1:K) # Indices for this sample set in whole result
  #   set
  Qmat[Ind,] = Set$Q
  CoalVec[Ind] = Set$Coal
  ErrorVec[iSet] = Set$Error
  if (max(Ind) * nProg / N >= iProg + 1) {
   iProg = floor(max(Ind) * nProg / N)
   cat(max(Ind), ":", table(CoalVec[1:max(Ind)]), "\n")
  }
 } # iSet
 list(Q = Qmat, Coal = CoalVec, Error = ErrorVec)
} # Unbiased

##############################################################################
# Code to run the above functions to generate the results in the paper

#*************************************** d = 1
# Preliminary runs to find a suitable value for B
Norm1 = Unbiased(N = 1000, K = 20, B = 20, M = 1, d = 1, sigma = 2, r = 3,
 Seed = 1)
# All coalesce in 1 block.
Norm1 = Unbiased(N = 1000, K = 20, B = 10, M = 1, d = 1, sigma = 2, r = 3,
 Seed = 1)
# 995 coalesce in 1 block, 5 in 2 blocks.
Norm1 = Unbiased(N = 1000, K = 20, B = 5, M = 1, d = 1, sigma = 2, r = 3,
 Seed = 1)
# 902 89 9

# Settled on B = 5, run a larger number of simulations.
Norm1 = Unbiased(N = 1e7, K = 20, B = 5, M = 1, d = 1, sigma = 2, r = 3,
 Seed = 1, nProg = 1e4)
save(Norm1, file = "Unbiased_1e7_5_1_3.RData")

# Tabulate the number of blocks needed to coalesce.  In the paper, we will
#   show the mean and maximum.
table(Norm1$Coal)
#       1       2       3       4       5       6       7       8       9 
# 9020142  862765  103848   11795    1282     144      18       4       2 
mean(Norm1$Coal) # 1.111185
# Check that the simulations follow the target distribution (takes time to
#   plot).
N = nrow(Norm1$Q)
K = 20
qqnorm(Norm1$Q)
c(mean(Norm1$Q), mean(Norm1$Q^2), 1 / sqrt(N)) # No problem
# Check serial correlation within sample sets.  Note that with d = 1, we have
#   B equal to only 5, so can expect some serial correlation.
Ind1 = 1:N
Ind2 = 2:(N + 1)
l = Ind1 %% K == 0 # Last simulation in sample set, needs to wrap around to
#   first one.
Ind2[l] = Ind2[l] - K
c(mean(Norm1$Q), mean(Norm1$Q * Norm1$Q[Ind2]), 1 / sqrt(N)) # Display
#   covariance and its standard deviation.
# 0.0094038867 0.0003162278
# It seems clear that there is a significant serial correlation, even though
#   the distribution of the covariance is not normal.  It's 30 standard
#   deviations.

#*************************************** d = 2
Norm1 = Unbiased(N = 1000, K = 20, B = 10, M = 1, d = 2, sigma = 2, r = 3,
 Seed = 1)
# 912 82 6

Norm1 = Unbiased(N = 1e6, K = 20, B = 10, M = 1, d = 2, sigma = 2, r = 3,
 Seed = 1, nProg = 1e3)
save(Norm1, file = "Unbiased_1e6_10_2_3.RData")

# Run similar statistics to d = 1 above.
table(Norm1$Coal)
#      1      2      3      4      5      6 
# 919249  77051   3514    179      5      2 
mean(Norm1$Coal) # 1.084646
# Check target distribution.
N = nrow(Norm1$Q)
K = 20
qqnorm(Norm1$Q[, 1])
qqnorm(Norm1$Q[, 2])
c(apply(Norm1$Q, 2, mean), 1 / sqrt(N)) # Mean of each component
apply(Norm1$Q^2, 2, mean) # Variance of each component
c(mean(Norm1$Q[, 1] * Norm1$Q[, 2]), 1 / sqrt(N)) # Covariance of components
# Check serial correlation within sample sets.  We have B = 10, so should
#   still have some detectable serial correlation, although not as much as
#   with d = 1 and B = 5.
Ind1 = 1:N
Ind2 = 2:(N + 1)
l = Ind1 %% K == 0 # Last simulation in sample set, needs to wrap around to
#   first one.
Ind2[l] = Ind2[l] - K
c(mean(Norm1$Q[, 1] * Norm1$Q[Ind2, 1]), 1 / sqrt(N)) # Display covariance and
#   its standard deviation.
# 0.002424157 0.001000000 # About 1/4 of what it was with d = 1, B = 5.  Maybe
#   statistically significant with N = 1e6, but it would be a close thing.
c(mean(Norm1$Q[, 2] * Norm1$Q[Ind2, 2]), 1 / sqrt(N))
# Try combining them.
c(mean(Norm1$Q * Norm1$Q[Ind2,]), 1 / sqrt(2 * N))
# 0.0026253881 0.0007071068, almost 4 standard deviations, bearing in mind
#   that it isn't normal.

#*************************************** d = 5
Norm1 = Unbiased(N = 1000, K = 20, B = 25, M = 1, d = 5, sigma = 2, r = 3,
 Seed = 1)
# 911 84 4 1

Norm1 = Unbiased(N = 1e6, K = 20, B = 25, M = 1, d = 5, sigma = 2, r = 3,
 Seed = 1, nProg = 1e3)
save(Norm1, file = "Unbiased_1e6_25_5_3.RData")

# Run similar statistics to d = 2 above.
table(Norm1$Coal)
#      1      2      3      4      5      6      7 
# 898587  95896   5082    397     36      1      1 
mean(Norm1$Coal) # 1.107406
# Check target distribution.
d = ncol(Norm1$Q)
N = nrow(Norm1$Q)
K = 20
qqnorm(Norm1$Q[, 1])
qqnorm(Norm1$Q[, 2])
qqnorm(Norm1$Q[, 5])
c(apply(Norm1$Q, 2, mean), 1 / sqrt(N)) # Mean of each component
apply(Norm1$Q^2, 2, mean) # Variance of each component
c(mean(Norm1$Q[, 1] * Norm1$Q[, 2]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 1] * Norm1$Q[, 5]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 2] * Norm1$Q[, 5]), 1 / sqrt(N)) # Covariance of components
# Check serial correlation within sample sets.  We have B = 25 now, so should
#   have much less serial correlation than above.
Ind1 = 1:N
Ind2 = 2:(N + 1)
l = Ind1 %% K == 0 # Last simulation in sample set, needs to wrap around to
#   first one.
Ind2[l] = Ind2[l] - K
# For d >= 5 we'll just combine all the components.
c(mean(Norm1$Q * Norm1$Q[Ind2,]), 1 / sqrt(d * N))
# 0.0020052451 0.0004472136, higher than I expected

#*************************************** d = 10
Norm1 = Unbiased(N = 1000, K = 20, B = 100, M = 1, d = 10, sigma = 2, r = 3,
 Seed = 1)
# 922 74 4
Norm1 = Unbiased(N = 1000, K = 20, B = 95, M = 1, d = 10, sigma = 2, r = 3,
 Seed = 1)
# 906 85 8 1
Norm1 = Unbiased(N = 1000, K = 20, B = 90, M = 1, d = 10, sigma = 2, r = 3,
 Seed = 1)
# 891 106 3

# Settled on B = 95, run the full set of simulations.
Norm1 = Unbiased(N = 1e6, K = 20, B = 95, M = 1, d = 10, sigma = 2, r = 3,
 Seed = 1, nProg = 1e4)
save(Norm1, file = "Unbiased_1e6_95_10_3.RData")
# We'll abbreviate the statistics for d >= 10.  We'll accept that the
#   methodology work and that the output follows the target distribution.
table(Norm1$Coal)
#      1      2      3      4      5      6      7 
# 903661  88523   6817    856    119     22      2 
mean(Norm1$Coal) # 1.105323
# Check the target distribution briefly.
d = ncol(Norm1$Q)
N = nrow(Norm1$Q)
K = 20
qqnorm(Norm1$Q[, 1])
qqnorm(Norm1$Q[, d])
c(apply(Norm1$Q, 2, mean), 1 / sqrt(N)) # Mean of each component
apply(Norm1$Q^2, 2, mean) # Variance of each component
c(mean(Norm1$Q[, 1] * Norm1$Q[, 2]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 1] * Norm1$Q[, d]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 2] * Norm1$Q[, d]), 1 / sqrt(N)) # Covariance of components
# Check serial correlation within sample sets.
Ind1 = 1:N
Ind2 = 2:(N + 1)
l = Ind1 %% K == 0 # Last simulation in sample set, needs to wrap around to
#   first one.
Ind2[l] = Ind2[l] - K
# For d >= 5 we just combine all the components.
c(mean(Norm1$Q * Norm1$Q[Ind2,]), 1 / sqrt(d * N))
# -0.0005468468 0.0003162278: this is looking undetectable at B = 95 with a
#   sample size 1e6, or 1e7 with all components combined.  It's actually
#   negative.

#*************************************** d = 15
Norm1 = Unbiased(N = 1000, K = 20, B = 100, M = 1, d = 15, sigma = 2, r = 3,
 Seed = 1)
# 185 361 200 100 59 36 22 18 9 6 3 1
Norm1 = Unbiased(N = 1000, K = 20, B = 200, M = 1, d = 15, sigma = 2, r = 3,
 Seed = 1)
# 600 282 71 30 13 2 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 400, M = 1, d = 15, sigma = 2, r = 3,
 Seed = 1)
# 890 92 12 5 1
Norm1 = Unbiased(N = 1000, K = 20, B = 425, M = 1, d = 15, sigma = 2, r = 3,
 Seed = 1)
# 911 81 7 1

# Settled on B = 425, run full set of simulations.
Norm1 = Unbiased(N = 1e6, K = 20, B = 425, M = 1, d = 15, sigma = 2, r = 3,
 Seed = 1, nProg = 1e4)
save(Norm1, file = "Unbiased_1e6_425_15_3.RData")
# Abbreviated statistics
table(Norm1$Coal) # Now looks like it's failing to decay geometrically, and
#   that the random-walk MCMC is starting to fall over, which is not
#   surprising.  I'm surprised it got this far.
#      1      2      3      4      5      6      7      8      9 
# 903480  82543  11341   2088    424     98     17      8      1 
mean(Norm1$Coal) # 1.113841
# Check the target distribution briefly.
d = ncol(Norm1$Q)
N = nrow(Norm1$Q)
K = 20
qqnorm(Norm1$Q[, 1])
qqnorm(Norm1$Q[, d])
c(apply(Norm1$Q, 2, mean), 1 / sqrt(N)) # Mean of each component
apply(Norm1$Q^2, 2, mean) # Variance of each component
c(mean(Norm1$Q[, 1] * Norm1$Q[, 2]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 1] * Norm1$Q[, d]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 2] * Norm1$Q[, d]), 1 / sqrt(N)) # Covariance of components
# Check serial correlation within sample sets.
Ind1 = 1:N
Ind2 = 2:(N + 1)
l = Ind1 %% K == 0 # Last simulation in sample set, needs to wrap around to
#   first one.
Ind2[l] = Ind2[l] - K
# For d >= 5 we just combine all the components.
c(mean(Norm1$Q * Norm1$Q[Ind2,]), 1 / sqrt(d * N))
# 0.0005548828 0.0002581989: about 2 standard deviations and not normal, so
#   highly doubtful that it's statistically significant

# Try a larger value for r.
Norm1 = Unbiased(N = 1000, K = 20, B = 200, M = 1, d = 15, sigma = 2, r = 4,
 Seed = 1)
# 857 108 23 8 2 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 250, M = 1, d = 15, sigma = 2, r = 4,
 Seed = 1)
# 922 67 9 2

# Higher value again
Norm1 = Unbiased(N = 1000, K = 20, B = 250, M = 1, d = 15, sigma = 2, r = 5,
 Seed = 1)
# 902 75 17 5 1
# Note that the problem of non-geometric decay becomes worse with larger
#   values for r.  It looks like this isn't the solution that might make
#   random-walk MCMC work in a higher number of dimensions.

#*************************************** d = 20
Norm1 = Unbiased(N = 1000, K = 20, B = 3000, M = 1, d = 20, sigma = 2, r = 3,
 Seed = 1)
# 860 113 20 6 1
Norm1 = Unbiased(N = 1000, K = 20, B = 3500, M = 1, d = 20, sigma = 2, r = 3,
 Seed = 1)
# 900 82 14 3 1

# Settled on B = 3500, run full set of simulations.
Norm1 = Unbiased(N = 1e5, K = 20, B = 3500, M = 1, d = 20, sigma = 2, r = 3,
 Seed = 1, nProg = 1e3)
save(Norm1, file = "Unbiased_1e5_3500_20_3.RData")

# Abbreviated statistics
table(Norm1$Coal) # Again not decaying geometrically
#     1     2     3     4     5     6     7     8     9 
# 88313  9772  1532   316    56     6     3     1     1 
mean(Norm1$Coal) # 1.14071
# Check the target distribution briefly.
d = ncol(Norm1$Q)
N = nrow(Norm1$Q)
K = 20
qqnorm(Norm1$Q[, 1])
qqnorm(Norm1$Q[, d])
c(apply(Norm1$Q, 2, mean), 1 / sqrt(N)) # Mean of each component
apply(Norm1$Q^2, 2, mean) # Variance of each component
c(mean(Norm1$Q[, 1] * Norm1$Q[, 2]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 1] * Norm1$Q[, d]), 1 / sqrt(N)) # Covariance of components
c(mean(Norm1$Q[, 2] * Norm1$Q[, d]), 1 / sqrt(N)) # Covariance of components
# Check serial correlation within sample sets.
Ind1 = 1:N
Ind2 = 2:(N + 1)
l = Ind1 %% K == 0 # Last simulation in sample set, needs to wrap around to
#   first one.
Ind2[l] = Ind2[l] - K
# For d >= 5 we just combine all the components.
c(mean(Norm1$Q * Norm1$Q[Ind2,]), 1 / sqrt(d * N))
# -0.0008390086  0.0007071068: negative again

# Try a larger value for r.
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 1, d = 20, sigma = 2, r = 5,
 Seed = 1)
# 898 77 17 7 1
Norm1 = Unbiased(N = 1e6, K = 20, B = 500, M = 1, d = 20, sigma = 2, r = 5,
 Seed = 1, nProg = 1e4)
save(Norm1, file = "Unbiased_1e6_500_20_5.RData")
