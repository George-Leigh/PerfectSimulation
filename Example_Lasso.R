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

# Author: George Leigh, 202212, updated 202307, 202312, 202402
# Run under R 4.2.2, Linux version
# Copyright 2023 George Leigh
# Associated journal article: Leigh, G. M., Yang, W-H, Wickens, M. E. and
#   Northrop, A. R. (2024).  "Algorithms for long run unbiased simulation and
#   maximal coupling with arbitrarily infrequent bias correction".
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
#   interactively, not by calling the "source" function of R.

# Note: This code is more general than the pseudocode in the paper.  It
#   includes the parameter M, which is the number of separate
#   (non-maximally-coupled) MCMC iterations to perform before each maximal
#   coupling step.  The setting M = 0 corresponds to Algorithms 3 and 4 of the
#   paper, whereby no separate MCMC is conducted.

##############################################################################
# We will estimate parameters in the Bayesian Lasso in 12 dimensions.
# The code follows the algorithms in the paper as closely as feasible.  The
#   major difference is that the code carries along logical variables
#   indicating whether coalescence has been attained, so that repeated testing
#   of whether chains take equal values is not necessary.

######################################## Settings for the Bayesian Lasso,
#   assigned to global variables
source("Prelim_Lasso.R")
dsetLasso = 12
# Values of lambda to use: 0, 0.237, 1, 2, 5, 10, 20
parsLasso = 0.237
names(parsLasso) = "Llambda"

######################################## Worker functions
Euclidean = function(x1, x2) sqrt(sum((x1 - x2)^2)) # Euclidean metric
Norm = function(x) sqrt(sum(x^2)) # Euclidean norm

#################### Function to find the element with minimum distance
MinInd = function(x) {
 # x: matrix with d columns (NOTE: must be a matrix, even when d == 1, and
 #   must have at least 2 rows)
 # We find the closest previous row to the final row of x.
 d = ncol(x)
 n = nrow(x)
 # First derive a version of x[1:(n - 1),] which we can be sure is a matrix.
 if (d == 1) {
  xmat = cbind(x[1:(n - 1),])
 } else {
  xmat = rbind(x[1:(n - 1),])
 }
 Dist = apply(xmat, 1, Euclidean, x[n,])
 which.min(Dist) # Index of row with minimum distance
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
  #if (is.na(normc) | is.na(r)) cat("Error: c r:", c, r, "\n")
  #if (is.na(normc) | is.na(r)) cat("Error: L v uX:", L, v, uX, "\n")
  #if (is.na(normc) | is.na(r)) cat("Error: X Y L:", X, Y, L, "\n")
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

#################### Bayesian Lasso with Efron et al. 2004 diabetes data
UfuncLasso = function(q) {
 # Additional global variables used that are not passed as parameters:
 # - LassoX: design matrix
 # - LassoY: y (dependent) variable
 # - betahat: OLS estimates of regression coefficients
 # - tauhat: estimate of parameter tau = 1 / sigma
 # - Ldiab: lower triangular de-scaling matrix from q[1:(d - 1)] to regression
 #   coefficients
 # - cvtauhat: coefficient of variation of tauhat; q[d] is the scaled log of
 #   tau (see calculation of Ltau below).
 # - LassoDf: number of degrees of freedom for the model, = number of data
 #   records + number of Lasso parameters.  We would subtract 1 if using tau
 #   as a model parameter but not if using log(tau) which is what we do.
 d = length(q)
 d1 = d - 1
 Llambda = parsLasso["Llambda"]
 Lbeta = betahat + Ldiab %*% cbind(q[1:d1]) # De-scale Lasso beta (regression
 #   coefficients)
 Ltau = tauhat * exp(cvtauhat * q[d])
 Ssum = sum((LassoY - LassoX %*% Lbeta)^2)
 Tsum = sum(abs(Lbeta[2:d1])) # Omit the intercept coefficient.
 # We could remove the log function below, but there is so much chance
 #   of introducing a bug that I haven't done it.
 Ltau * (0.5 * Ltau * Ssum + Llambda * Tsum) - LassoDf * log(Ltau)
}

UderivLasso = function(q) {
 d = length(q)
 d1 = d - 1
 Llambda = parsLasso["Llambda"]
 Lbeta = betahat + Ldiab %*% cbind(q[1:d1])
 signbeta = sign(Lbeta)
 signbeta[1] = 0 # Remove intercept parameter from Lasso.
 Ltau = tauhat * exp(cvtauhat * q[d])
 dtaudq = cvtauhat * Ltau # Derivative of Ltau w.r.t. q[d]
 Ssum = sum((LassoY - LassoX %*% Lbeta)^2)
 Tsum = sum(abs(Lbeta[2:d1])) # Omit the intercept coefficient.
 c(t((Llambda * Ltau) * signbeta - Ltau^2 * t(LassoX) %*%
  (LassoY - LassoX %*% Lbeta)) %*% Ldiab,
  dtaudq * (Ltau * Ssum + Llambda * Tsum - LassoDf / Ltau))
}

#################### Function to generate starting points for MCMC, from
#   a widely dispersed uniform distribution for each coordinate
# These starting points appear to be too widely dispersed.  The results of the
#   simulation depend too heavily on the dispersion of the starting points.
StartWide = function(K, d) { # Generate starting points (K x d matrix).  We
 #   will sample each coordinate independently from a uniform distribution on
 #   [-6, 6].
 xlim = 6
 X = array(runif(K * d, -xlim, xlim), dim = c(K, d))
 X
}

#################### Function to generate starting points for MCMC, from a
#   fairly widely dispersed uniform distributions on a sphere
Start = function(K, d) { # Generate starting points (K x d matrix).  We will
 #   sample from a d-dimensional hypersphere of radius 6.
 xlim = 6
 Dir = array(rnorm(K * d), dim = c(K, d)) # Normal variables to get direction
 Rad = xlim * runif(K)^(1 / d)
 Scal = function(x) x / Norm(x)
 # Need the rbind function here in the MS Windows version of R, but not in the
 #   Linux version.
 X = Rad * t(rbind(apply(Dir, 1, Scal))) # Rad is recycled after each column.
 X
}

#################### Function to generate starting points for MCMC, from the
#   target N(0, I) distribution, to find the best conditions that we can hope
#   for as d increases
StartTarg = function(K, d) { # Generate starting points (K x d matrix).  We
 #   will sample  each coordinate independently from a N(0, 1) distribution.
 X = array(rnorm(K * d), dim = c(K, d))
 X
}

#################### Function to generate random numbers for MCMC
# We allow the option M = 0, which means that we use only the maximal coupling
#   algorithm, without any separate MCMC process.  When the separate MCMC is
#   included, it uses a separate Metropolis-Hastings algorithm.
Rand = function(d, M) { # Generate random numbers for M MCMC iterations and
 #   one maximal coupling step.
 Rlist = list(
  Rdir = rnorm(d + 2), # 1:d: maximal coupling direction; 1:(d + 2): chi^2_{d
  #   + 2} distribution for r
  Rmag = runif(1), # Maximal coupling magnitude
  RMH = runif(1)) # Maximal coupling Metropolis-Hastings test random variable
 if (M > 0) {
  Rlist$RMCMC = list(
   Rjump = array(rnorm(M * d), dim = c(M, d)), # MCMC jump vectors
   RMH = runif(M))
 }
 Rlist
}

#################### Function to conduct MCMC: this has an inbuilt
#   Metropolis-Hastings algorithm, separate to the one used for maximal
#   coupling.
# This function is not used if M is set to zero.
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
# rho = scale parameter for r.
SampleSet = function(K, B, M, d, sigma = 1, rho, Seed, lStartTarg = FALSE) { #
 #   Main function to run the algorithm, generates K slightly correlated
 #   perfect samples, if an error doesn't occur.
 # K: Number of blocks, also equal to number of chains
 # B: Number of MCMC iterations in a block
 # M: Number of MCMC iterations per maximal coupling step (should be either
 #    zero or a divisor of B)
 # d: Number of dimensions
 # sigma: Standard deviation of jump in MCMC (not used if M = 0)
 # rho: Scale parameter for r, the radius of the uniform spherical
 #   distribution for maximal coupling.
 # Seed: Random number seed, for reproducibility
 # NOTE: The code includes an extra scale factor of 1 / sqrt(d) in sigma and
 #   rho, to prevent the Metropolis-Hastings acceptance probability from
 #   becoming very small as d increases.

 ##### Setup
 if (Seed > 0) set.seed(Seed) # Zero implies don't reset seed.
 if (M > 0) {
  J = ceiling(B / M) # Number of maximal coupling steps per block
  B = J * M # Redefine B if not a multiple of M.
 } else {
  J = B
 }
 # Generate starting points (K x d matrix) along diagonal of matrix in Fig. 2
 #   of paper.
 if (lStartTarg) { # Use N(0, I) distribution.
  Q = StartTarg(K, d)
 } else { # Use widely dispersed uniform distribution.
  Q = Start(K, d) # Generate starting points (K x d matrix) along diagonal of
 }
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
 BlocksCftp = rep(0, K) # Number of blocks into the past at which a chain
 #   coalesces with all the chains above it (wrapping around vertically as
 #   necessary): a value of K - 1 means it has not coalesced with the chain
 #   below it, which results in Error being true.  A value of zero means it
 #   has not shared any random numbers with any other chains.
 Coal1 = rep(0, K) # Column number at which a row coalesces with row 1 in the
 #   upper triangle, used in the calculation of BlocksCftp (first element is
 #   not used)
 a = rep(0, K) # Mark all chains as active (not coalesced with any other
 #   chain).
 # We can store the whole Q matrix (whole of matrix from Fig. 2 of paper), for
 #   diagnostic purposes only.  It is not necessary for the operation of the
 #   algorithm, but without it a return value of true for Error won't provide
 #   any description of what has gone wrong.  If this storage is not desired,
 #   all lines of code containing "Qarray" can be removed or commented out.
 #Qarray = array(0, dim = c(K, K, d))

 ##### Upper triangle
 for (j in 1:K) { # Loop over blocks (columns) in Fig. 2 upper triangle.
  for (j1 in 1:J) {
   R = Rand(d, M) # Generate random numbers for MCMC and maximal coupling.
   # Assign r immediately after calling Rand.  NOTE: R$Rdir has d + 2 elements.
   r = rho * Norm(R$Rdir) / sqrt(d)
   for (i in 1:j) { # Loop over chains (rows) in Fig. 2 upper triangle.
    if (a[i] > 0) { # Coalesced to a previous row, so we just copy that.
     Q[i,] = Q[a[i],]
    } else { # Not coalesced, so we run MCMC (if M > 0).
     if (M > 0) {
      Qminus[i,] = Mcmc(Q[i,], M, R$RMCMC, sigma, UfuncLasso) # MCMC creates
      #   pre-jump version of Q.
     } else {
      Qminus[i,] = Q[i,] # No MCMC separate from maximal coupling
     }
     if (i == 1) { # Row 1 of Fig. 2: Jump freely, M-H test and store.
      Qstar[1,] = MaxCoupleJump(Qminus[1,], r, R$Rdir[1:d], R$Rmag)
      MHResult = MHTest(Qminus[1,], Qstar[1,], R$RMH, UfuncLasso)
      #print(MHResult)
      #print(Q)
      Q[1,] = MHResult$Dest
      Qjumped[1] = MHResult$lJumped
      #if (length(Q1minus[j, j1,]) != length(Qminus[1,])) {
      # cat("Q1minus[j, j1,]", Q1minus[j, j1,], "\n")
      # cat("Qminus[1,]", Qminus[1,], "\n")
      #}
      Q1minus[j, j1,] = Qminus[1,]
      Q1star[j, j1,] = Qstar[1,]
      Q1[j, j1,] = Q[1,]
      Q1jumped[j, j1] = Qjumped[1]
     } else { # i > 1 in upper triangle of Fig. 2: Maximally couple with the
      #   closest previous row (m < i).
      # Different structure to pseudocode in paper: Subset the rows before
      #   calling MinInd, and include only rows with a[i] == 0 (i.e.,
      #   non-coalesced rows).  We then need only one argument to MinInd, and
      #   we are not passing values that won't be used.
      Ind = which(a[1:i] == 0)
      # We know that Ind has at least 2 elements, so we can use cbind with
      #   confidence.
      m = Ind[MinInd(cbind(Qminus[Ind,]))]
      #if (Norm(Qminus[m,] - Qminus[i,]) == 0)
      # cat("Error1: m i a[m] a[i]:", m, i, a[m], a[i], "\n")
      MaxCoupleResult = MaxCouple(Qminus[m,], Qstar[m,], Qminus[i,], r)
      Qstar[i,] = MaxCoupleResult$Ystar
      lCoupled = MaxCoupleResult$lCoupled
      MHResult = MHTest(Qminus[i,], Qstar[i,], R$RMH, UfuncLasso)
      Q[i,] = MHResult$Dest
      Qjumped[i] = MHResult$lJumped
      # We don't need to test equality of Q[i,] with Q[m,].  We use the three
      #   logical variables instead.
      if (lCoupled & Qjumped[m] & Qjumped[i]) {
       # It is not possible to have a[m] > 0, because we set Ind to be indices
       #   of nonzero a's above.  Therefore, comment out this code.
       #if (a[m] > 0)
       # m = a[m] # Reset m to the earliest chain with which chain i has
       # #   coalesced.  This is the row with which row m coalesced, if any.
       a[i] = m # Mark row i as coalesced with row m, so we can just copy the
       #   rest of row i from the previous row a[i].
       # Mark any lower rows that earlier coalesced with row i as having
       #   instead coalesced with row m (using the new value of m).
       a[a == i] = m
      } # Coalesced
     } # i > 1
    } # Not coalesced: may or may not be coalesced now!
    if (Coal1[i] == 0 & a[i] == 1)
     Coal1[i] = j
   } # i
  } # j1
  #Qarray[1:j, j,] = Q[1:j,] # Store upper triangle of matrix, diagnostic only.
 } # j

 ##### Lower triangle
 .Random.seed <<- s # Restore seed to regenerate random numbers used for upper
 #   triangle.  We want to use the same random numbers in the part of each
 #   column j that is in the lower triangle of Fig. 2.
 for (j in 1:(K - 1)) { # Loop over blocks (columns) in lower triangle
  # We will operate with the random numbers in column j, and will produce the
  #   final, CFTP-adjusted value of the chain in row j + 1.
  # Before doing that, however, we will check the coalescence of the final
  #   (CFTP-adjusted) element of the chain in row j, which we have just
  #   finished.

  ## Preliminary work on column j - 1 (wrapping around to K when j = 1), to
  #   determine coalescence status of row j (the CFTP adjusted version of
  #   diagonal element (j - 1, j - 1))
  # First check for an error in coalescence of row j with row j + 1, which
  #   could indicate the need for the BC (bias correction) term in row j,
  #   which we don't want.
  if (a[j + 1] != j) cat("Error", j, Q[j:(j + 1), 1], "\n")
  Error = Error | a[j + 1] != j # Error if any row j has not coalesced with
  #   the one below it after running for K blocks since it was initialised at
  #   the diagonal element of the matrix.  These are processes X and Y
  #   respectively in the theory of unbiased sampling.  It will actually be OK
  #   (no holes generated) if it coalesces in one more block, using new random
  #   numbers, but we're flagging it as an error here to keep the algorithm
  #   fairly simple.

  # Next, check all the rows to assign the value of BlocksCftp which will tell
  #   us how many blocks of random numbers from earlier chains have been
  #   shared in generating the CFTP-adjusted value in row j (and column j - 1,
  #   wrapping around where necessary).  We start from row j + 1 and stop when
  #   we find a value that is not equal to cell (j, j - 1).  Then we are
  #   taking the greatest number of blocks for which any cell in column j - 1
  #   is not equal to cell (j, j - 1).
  AllCoal = TRUE
  for (i in (j + 1):(j + K - 1)) {
   i1 = (i - 1) %% K + 1
   # Condition a[i1] == j below is for coalescence when i1 > j; a[j] == i1 is
   #   for i1 == 1; and Coal1[i] <= j - 1 is for 1 < i1 < j.
   if (!(a[i1] == j | a[j] == i1 | (i > K + 1 & Coal1[i1] <= j - 1))) {
    AllCoal = FALSE
    break
   }
  }
  if (!AllCoal)
   BlocksCftp[j] = (j - i) %% K

  ## Work on column j, from row j + 1 onwards
  # We have finished with row j.  It has cycled all the way around to where it
  #   was initialised.  Transfer all coalescence involving it to coalescence
  #   with row j + 1.  In transferring the coalescence, we are assuming that
  #   rows j and j + 1 have coalesced.  If they haven't, the variable "Error"
  #   will be set to true above, so we will detect an error.
  a[j + 1] = a[j]
  a[a == j] = j + 1
   
  for (j1 in 1:J) {
   # Restore saved copy of row 1 (all versions); no need to recompute.
   Qminus[1,] = Q1minus[j, j1,]
   Qstar[1,] = Q1star[j, j1,]
   Q[1,] = Q1[j, j1,]
   Qjumped[1] = Q1jumped[j, j1]
   # Copy Qminus, Qstar and Qjumped for row 1, if row j + 1 has coalesced with
   #   row 1.  This is needed because we later want to couple subsequent rows
   #   with row j + 1, not with row 1.  It wasn't needed for the upper
   #   triangle.  Note we copy Q lower down (and we do it for all values of
   #   i), so we don't do it here.
   if (a[j + 1] == 1) {
    Qminus[j + 1,] = Qminus[1,]
    Qstar[j + 1,] = Qstar[1,]
    Qjumped[j + 1] = Qjumped[1]
   }
   # We have to generate the random numbers even if we are not going to need
   #   them, because if we don't the next sample set will not start with fresh
   #   random numbers.
   R = Rand(d, M) # Generate random numbers for MCMC and maximal coupling.
   # Assign r immediately after calling Rand.  NOTE: R$Rdir has d + 2 elements.
   r = rho * Norm(R$Rdir) / sqrt(d)
   for (i in (j + 1):K) {
    if (a[i] > 0) {
     Q[i,] = Q[a[i],] # Copy from chain with which row i has coalesced.
    } else { # a[i] == 0: no coalescence: run MCMC.
     if (M > 0) {
      Qminus[i,] = Mcmc(Q[i,], M, R$RMCMC, sigma, UfuncLasso) # MCMC creates
      #   pre-jump version of Q.
     } else {
      Qminus[i,] = Q[i,] # No MCMC that is separate from maximal coupling
     }
     # Set the row with which we will maximally couple row i.  We are using
     #   rows 1 and j + 1, j + 2, ..., K.  The rows inbetween have cycled
     #   right around, back to their points of initialisation on the diagonal
     #   of Fig. 2, and we have finished with them.
     if (i == j + 1) {
      m = 1 # Only row j + 1 is maximally coupled with row 1.
     } else { # i > j + 1
      # Different structure to pseudocode in paper: Subset the rows before
      #   calling MinInd, and include only rows with a[i] == 0 (i.e.,
      #   non-coalesced rows).  We then need only one argument to MinInd, and
      #   we are not passing values that won't be used.
      # NOTE: We want j + 1 to be included in Ind, so condition on a <= 1
      #   here, not a == 0 as in the upper triangle.
      Ind = j + which(a[(j + 1):i] <= 1)
      # We know that Ind has at least 2 elements, so we can use cbind with
      #   confidence.
      m = Ind[MinInd(cbind(Qminus[Ind,]))] # Choose the closest row from rows
      #   j + 1, ..., i - 1.
     }
     #if (Norm(Qminus[m,] - Qminus[i,]) == 0)
     # cat("Error2: m i a[m] a[i]:", m, i, a[m], a[i], "\n")
     MaxCoupleResult = MaxCouple(Qminus[m,], Qstar[m,], Qminus[i,], r)
     Qstar[i,] = MaxCoupleResult$Ystar
     lCoupled = MaxCoupleResult$lCoupled
     MHResult = MHTest(Qminus[i,], Qstar[i,], R$RMH, UfuncLasso)
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
    } # a[i] == 0: may or may not be coalesced now!
   } # i
  } # j1
  #Qarray[(j + 1):K, j,] = Q[(j + 1):K,] # Store lower triangle of matrix, as
  #   diagnostic only.
 } # j

 ## Repeat the preliminary work, now for coalescence status of element (K, K -
 #   1), the CFTP adjusted version of diagonal element (K - 1, K - 1).
 # First check for an error in coalescence of row K with row 1, which could
 #   indicate the need for the BC (bias correction) term in row K, which we
 #   don't want.
 if (a[K] != 1) cat("Error", K, Q[c(K, 1)], "\n")
 Error = Error | a[K] != 1 # Error if row K has not coalesced with row 1 after
 #   running for K blocks and ending up at position (K, K - 1) of the matrix
 #   in Fig. 2.  Now rows K and 1 are processes X and Y respectively in the
 #   theory of unbiased sampling.  Again it will be OK (no holes generated) if
 #   rows K and 1 coalesce in one more block, but we keep it simple and flag
 #   it as an error here.

 # Next, check all the rows to assign the value of BlocksCftp which will tell
 #   us how many blocks of random numbers from earlier chains have been shared
 #   in generating the CFTP-adjusted element (K, K - 1).  We start from row 1
 #   and stop when we find a value that is not equal to cell (K, K - 1).  Then
 #   we are taking the greatest number of blocks for which any cell in column
 #   K - 1 is not equal to cell (K, K - 1).
 AllCoal = TRUE
 for (i in 1:(K - 1)) {
  # Condition a[K] == i is for i == 1; and Coal1[i] <= K - 1 is for 1 < i < K.
  if (!(a[K] == i | (i > 1 & Coal1[i] <= K - 1))) {
   AllCoal = FALSE
   break
  }
 }
 if (!AllCoal)
  BlocksCftp[K] = K - i

 ## Final collation of results
 Q[1,] = Q1[K, J,] # Restore final value in first row of matrix, after we
 #    overwrote it.
 #list(Q = Q, Error = Error, Qarray = Qarray)
 list(Q = Q, Cftp = BlocksCftp, Error = Error)
} # SampleSet

#################### Function to conduct long run unbiased simulation using
#   the unbiased simulation algorithm, simulating a collection of independent
#   sample sets
Unbiased = function(N, K, B, M, d, sigma, rho, Seed = 1, nProg = 100,
 lStartTarg = FALSE) {
 # N: Total number of sample points, should be a multiple of K below
 # K: Number of blocks, also = number of chains, should be a divisor of N
 # B: Number of MCMC iterations in a block
 # M: Number of MCMC iterations per maximal coupling step, should be a divisor
 #    of B
 # d: Number of dimensions
 # sigma: Standard deviation of jump in MCMC
 # rho: Scale factor for radius of uniform spherical distribution for maximal
 #   coupling
 # Seed: Random number seed, for reproducibility
 # nProg: Number of progress records to display
 # lStartTarg: Logical, whether to use the target distribution for starting
 #   points: obviously, this is cheating, but it provides the best conditions
 #   for the algorithm, for us to gauge what happens as d increases.
 if (Seed > 0) set.seed(Seed) # Zero implies don't reset seed.
 NSet = ceiling(N / K) # Number of sample sets to simulate
 N = NSet * K # Redefine N if not a multiple of K.
 Qmat = array(0, dim = c(N, d)) # Matrix to hold final simulations
 CftpVec = rep(0, N) # Vector to hold number of blocks taken to reach
 #   coalescence for each sample point
 ErrorVec = rep(FALSE, NSet) # Vector to hold error indicator for each sample
 #   set
 iProg = 0 # Number of progress records written so far
 for (iSet in 1:NSet) {
  Set = SampleSet(K, B, M, d, sigma, rho, Seed = 0, lStartTarg = lStartTarg)
  Ind = (iSet - 1) * K + (1:K) # Indices for this sample set in whole result
  #   set
  Qmat[Ind,] = Set$Q
  CftpVec[Ind] = Set$Cftp
  ErrorVec[iSet] = Set$Error
  if (max(Ind) * nProg / N >= iProg + 1) {
   iProg = floor(max(Ind) * nProg / N)
   cat(max(Ind), ":", table(CftpVec[1:max(Ind)]), "\n")
  }
 } # iSet
 list(Q = Qmat, Cftp = CftpVec, Error = ErrorVec)
} # Unbiased

##############################################################################
# Code to run the above functions to generate the results in the paper

#**************************************** Most of the runs use the target
#   distribution for the starting points: I think it is best to set rho this
#   way.  Then we do a final run with the fairly widely dispersed starting
#   points.

#******************** d = 12, lambda = 0 (rho = 3.25, B = 170, 420, 460)
parsLasso["Llambda"] = 0
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, lStartTarg = TRUE)
#  517 279 114 49 23 10 4 2 1 1 (10)
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#* 526 260 119 52 22 13 5 1 1 1 (10)
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  481 263 142 69 26 10 6 2 1 (9)

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 1000, K = 20, B = 450, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1)
#  884 106 9 1 (4)
Lass = Unbiased(N = 1000, K = 20, B = 460, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1)
#* 916 80 3 1 (4)
Lass0 = Unbiased(N = 1e6, K = 20, B = 460, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, nProg = 1000)
#  899915 91849 7553 627 50 6 (6)
save(Lass0, file = "Lass0.RData")

#******************** d = 12, lambda = 0.237 (rho = 3.25, B = 170, 420, 460)
parsLasso["Llambda"] = 0.237
setglobalLasso(parsLasso["Llambda"])

#***** Results (P ~ 0.1), dispersed starting points
# Confident that settings from lambda = 0 can be used here.
Lass0p237 = Unbiased(N = 1e6, K = 20, B = 460, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, nProg = 1000)
#  899530 92349 7473 593 52 3 (6)
save(Lass0p237, file = "Lass0p237.RData")

#******************** d = 12, lambda = 1 (rho = 3.25, B = 170, 420, 460)
parsLasso["Llambda"] = 1
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, lStartTarg = TRUE)
#  
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#* 
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1)
#  
Lass1 = Unbiased(N = 1e6, K = 20, B = 500, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, nProg = 1000)
#  914061 79909 5626 372 30 2 (6)
save(Lass1, file = "Lass1.RData")

#******************** d = 12, lambda = 2 (rho = 3.25, B = 170, 420, 460)
parsLasso["Llambda"] = 2
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, lStartTarg = TRUE)
#  415 306 159 68 31 13 4 2 1 1 (10)
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#* 477 295 134 66 18 7 3 (7)
Lass = Unbiased(N = 1000, K = 25, B = 170, M = 0, d = dsetLasso, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  463 282 144 59 25 12 6 4 3 2 (10)

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1)
#  888 106 5 1 (4)
Lass2 = Unbiased(N = 1e6, K = 20, B = 500, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, nProg = 1000)
#  895637 95545 8068 685 61 4 (6)
save(Lass2, file = "Lass2.RData")

#******************** d = 12, lambda = 5 (rho = 3.25, B = 170, 420, 460)
parsLasso["Llambda"] = 5
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 1000, K = 25, B = 250, M = 0, d = dsetLasso, rho = 2.75,
 Seed = 1, lStartTarg = TRUE)
#* 495 315 121 44 16 6 1 1 1 (9)
Lass = Unbiased(N = 1000, K = 25, B = 250, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, lStartTarg = TRUE)
#  488 313 128 42 18 8 3 (7)
Lass = Unbiased(N = 1000, K = 25, B = 250, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  468 308 134 54 24 9 2 1

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 1000, K = 20, B = 700, M = 0, d = dsetLasso, rho = 2.75,
 Seed = 1)
#  913 81 5 1 (4)
Lass = Unbiased(N = 1000, K = 20, B = 700, M = 0, d = dsetLasso, rho = 3,
 Seed = 1)
#* 917 76 7 (3) Nothing between rho = 2.75 and 3; I suggest that rho = 3 here
#   gives a smoother transition as lambda increases.
Lass = Unbiased(N = 1000, K = 20, B = 700, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1)
#  897 96 7 (3)
Lass5 = Unbiased(N = 1e6, K = 20, B = 700, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, nProg = 1000)
#  909457 83912 6165 438 27 1 (6)
save(Lass5, file = "Lass5.RData")

#******************** d = 12, lambda = 10 (rho = 3.25, B = 170, 420, 460)
parsLasso["Llambda"] = 10
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 1000, K = 25, B = 400, M = 0, d = dsetLasso, rho = 2.75,
 Seed = 1, lStartTarg = TRUE)
#  472 294 140 59 25 7 2 1
Lass = Unbiased(N = 1000, K = 25, B = 400, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, lStartTarg = TRUE)
#* 484 303 145 49 14 4 1 (7)
Lass = Unbiased(N = 1000, K = 25, B = 400, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  465 293 140 58 25 13 4 1 1 (9)
Lass = Unbiased(N = 1000, K = 25, B = 400, M = 0, d = dsetLasso, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  442 303 147 57 28 17 2 2 2 (9)  

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 1000, K = 20, B = 900, M = 0, d = dsetLasso, rho = 3,
 Seed = 1)
#* 847 132 19 2 (4)
Lass = Unbiased(N = 1000, K = 20, B = 900, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1)
#  888 106 5 1 (4) This result is anomalous.  I suggest that we want rho = 3.
Lass = Unbiased(N = 1000, K = 20, B = 1000, M = 0, d = dsetLasso, rho = 3,
 Seed = 1)
#  886 103 10 1 (4)
Lass = Unbiased(N = 1000, K = 20, B = 1100, M = 0, d = dsetLasso, rho = 3,
 Seed = 1)
#  910 86 4 (3)
#Lass10 = Unbiased(N = 1e6, K = 20, B = 900, M = 0, d = dsetLasso, rho = 3,
# Seed = 1, nProg = 1000)
##  850038 131822 15969 1918 219 29 5 (7)
#save(Lass10, file = "Lass10.RData")
# Above didn't work very well so think we'd better redo with greater B.
#   Grounds for setting that value of B were shaky anyway.
Lass10 = Unbiased(N = 1e6, K = 20, B = 1100, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, nProg = 1000)
#  908543 84704 6243 469 37 3 1 (7)
save(Lass10, file = "Lass10.RData")

#******************** d = 12, lambda = 20 (rho = 2.75, B = 900, 2500)
parsLasso["Llambda"] = 20
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 1000, K = 25, B = 900, M = 0, d = dsetLasso, rho = 2.5,
 Seed = 1, lStartTarg = TRUE)
#  537 274 111 44 21 9 3 1
Lass = Unbiased(N = 1000, K = 25, B = 900, M = 0, d = dsetLasso, rho = 2.75,
 Seed = 1, lStartTarg = TRUE)
#* 545 304 108 32 7 2 1 1
Lass = Unbiased(N = 1000, K = 25, B = 900, M = 0, d = dsetLasso, rho = 3,
 Seed = 1, lStartTarg = TRUE)
#  526 281 120 52 14 4 2 1 (8)
Lass = Unbiased(N = 1000, K = 25, B = 900, M = 0, d = dsetLasso, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  503 298 113 52 23 11 (6)
Lass = Unbiased(N = 1000, K = 25, B = 900, M = 0, d = dsetLasso, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 500, K = 20, B = 2500, M = 0, d = dsetLasso, rho = 2.75,
 Seed = 1)
#* 462 37 1 (3)
Lass = Unbiased(N = 500, K = 20, B = 2500, M = 0, d = dsetLasso, rho = 3,
 Seed = 1)
#  446 50 4
Lass20 = Unbiased(N = 1e6, K = 20, B = 2500, M = 0, d = dsetLasso, rho = 2.75,
 Seed = 1, nProg = 1000)
#  931097 64861 3798 228 14 2 (6)
save(Lass20, file = "Lass20.RData")

#******************** d = 12, lambda = 30 (rho = 2.75, B = 900, 2500)
parsLasso["Llambda"] = 30
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 1000, K = 25, B = 2000, M = 0, d = dsetLasso, rho = 2.5,
 Seed = 1, lStartTarg = TRUE)
#* 656 239 72 20 10 2 1
Lass = Unbiased(N = 1000, K = 25, B = 2000, M = 0, d = dsetLasso, rho = 2.75,
 Seed = 1, lStartTarg = TRUE)
#  608 250 93 34 10 4 1

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 500, K = 20, B = 4000, M = 0, d = dsetLasso, rho = 2.5,
 Seed = 1)
#  440 50 9 1 (4)
Lass = Unbiased(N = 500, K = 20, B = 4500, M = 0, d = dsetLasso, rho = 2.5,
 Seed = 1)
#* 464 32 3 1 (4) We'll choose the safer one.
Lass30 = Unbiased(N = 1e6, K = 20, B = 4500, M = 0, d = dsetLasso, rho = 2.5,
 Seed = 1, nProg = 1000)
#  925158 69942 4583 299 16 2 (6)
save(Lass30, file = "Lass30.RData")

#******************** d = 12, lambda = 40 (rho = 2.25, B = 7000)
parsLasso["Llambda"] = 40
setglobalLasso(parsLasso["Llambda"])

#***** Set rho (P ~ 0.5)
Lass = Unbiased(N = 500, K = 25, B = 2500, M = 0, d = dsetLasso, rho = 1.75,
 Seed = 1, lStartTarg = TRUE)
#  193 137 82 37 21 15 8 3 2 2
Lass = Unbiased(N = 500, K = 25, B = 2500, M = 0, d = dsetLasso, rho = 2,
 Seed = 1, lStartTarg = TRUE)
#  229 141 81 36 7 4 1 1 (8)
Lass = Unbiased(N = 500, K = 25, B = 2500, M = 0, d = dsetLasso, rho = 2.25,
 Seed = 1, lStartTarg = TRUE)
#* 228 147 69 33 16 6 1 (7) # Nothing to choose between 2 and 2.25, but 2.5 is
#   better than 1.75, so we'll choose 2.25.
Lass = Unbiased(N = 500, K = 25, B = 2500, M = 0, d = dsetLasso, rho = 2.5,
 Seed = 1, lStartTarg = TRUE)
#  210 149 76 38 15 8 3 1 (8)

#***** Results (P ~ 0.1), dispersed starting points
Lass = Unbiased(N = 1400, K = 20, B = 6500, M = 0, d = dsetLasso, rho = 2.25,
 Seed = 1, nProg = 7)
#  1235 149 16 (3) Looks too low, although we don't really have enough data to
#   say that.  We'll use a higher value of B for safety.
Lass40 = Unbiased(N = 1e6, K = 20, B = 7000, M = 0, d = dsetLasso, rho = 2.25,
 Seed = 1, nProg = 10000)
#  900224 91079 7901 723 67 6 (6)
save(Lass40, file = "Lass40.RData")

##############################################################################
# Reload to get results for table and plots.
LlambdaVec = c(0, 0.237, 1, 2, 5, 10, 20, 30, 40)
names(LlambdaVec) = c("0", "0p237", "1", "2", "5", "10", "20", "30", "40")

for (iLlambda in names(LlambdaVec)) {
 Llambda = LlambdaVec[iLlambda]
 cat("lambda =", Llambda, "\n")

 parsLasso["Llambda"] = Llambda
 setglobalLasso(parsLasso["Llambda"])
 load(paste0("Lass", iLlambda, ".RData"))
 Lass = get(paste0("Lass", iLlambda))
 remove(list = paste0("Lass", iLlambda))

 dSamp = 12
 d1Samp = 11
 nSamp = dim(Lass$Q)[1]
 LbetaSamp = rep(betahat, nSamp) +
  Ldiab %*% t(Lass$Q[, 1:d1Samp]) # De-scale Lasso beta (regression
  #   coefficients); produces a matrix d1Samp x nSamp
 LtauSamp = tauhat * exp(cvtauhat * Lass$Q[,dSamp])
 SsumSamp = apply((rep(LassoY, nSamp) - LassoX %*% LbetaSamp)^2, 2, sum)
 TsumSamp = apply(abs(LbetaSamp[2:d1Samp,]), 2, sum) # Omit the intercept
 #   coefficient.

 cat("Quantiles of S:", quantile(SsumSamp, probs = c(0, 0.05, 0.1, 0.25, 0.5,
  0.75, 0.9, 0.95, 1)), "\n")
 cat("Quantiles of T:", quantile(TsumSamp, probs = c(0, 0.05, 0.1, 0.25, 0.5,
  0.75, 0.9, 0.95, 1)), "\n")
 hist(TsumSamp, 100)
 #readline("Press enter to continue")
}

# lambda = 0 
# Quantiles of S: 1265733 1277034 1279898 1285678 1293667 1303501 1313983
#   1321033 1419508
# Quantiles of T: 78.92309 109.0791 116.8686 137.5688 168.7061 204.0939
#   236.8715 256.5723 434.0313

# lambda = 0.237
# Quantiles of S: 1265130 1277003 1279873 1285597 1293555 1303285 1313717
#   1320742 1408273
# Quantiles of T: 76.68341 106.6202 113.1078 130.5548 158.8255 192.8429
#   225.0672 244.779 427.0898

# lambda = 1 
# Quantiles of S: 1265366 1277280 1280165 1285938 1293869 1303553 1313960
#   1321006 1403081
# Quantiles of T: 72.50473 101.6540 106.2715 116.6983 136.2292 163.2307
#   191.9491 210.0534 371.1379

# lambda = 2
# Quantiles of S: 1264917 1278110 1281057 1286879 1294788 1304429 1314696
#   1321622 1412198
# Quantiles of T: 71.19344 97.86399 101.5957 109.0287 121.0395 139.9943
#   161.8083 176.8864 332.3045

# lambda = 5
# Quantiles of S: 1265607 1280434 1283390 1289122 1296897 1306410 1316595
#   1323485 1424960
# Quantiles of T: 68.67167 91.83866 94.86330 100.2426 107.0069 115.4711
#   125.8971 133.7971 256.3718

# lambda = 10
# Old results with B = 900 (too low; doesn't affect validity of results)
# Quantiles of S: 1268485 1282792 1285746 1291638 1299827 1309899 1320817
#   1328215 1435654
# Quantiles of T: 64.76576 85.62321 88.34821 93.02379 98.50021 104.4114
#   110.3865 114.4291 177.1026
# New results with B = 1100
# Quantiles of S: 1267452 1282765 1285737 1291612 1299788 1309854 1320770
#   1328136 1426756
# Quantiles of T: 65.90696 85.6462 88.35288 93.04503 98.51117 104.4215
#   110.3648 114.413 171.8842

# lambda = 20
# Quantiles of S: 1271660 1288507 1292409 1300227 1310910 1323916 1337542
#   1346533 1453375
# Quantiles of T: 58.66476 76.07183 78.52329 82.71053 87.54763 92.54476
#   97.21326 100.0908 133.9478

# lambda = 30
# Quantiles of S: 1273795 1299026 1304565 1315270 1329186 1345154 1361192
#   1371605 1493654
# Quantiles of T: 52.94148 68.53583 70.67433 74.44516 78.82135 83.38405
#   87.63918 90.24855 122.6319

# lambda = 40
# Quantiles of S: 1276272 1314423 1321543 1334719 1350919 1369062 1387152
#   1398949 1560792
# Quantiles of T: 48.01501 62.70826 64.60336 67.91123 71.78873 75.88623
#   79.76367 82.15969 104.271
