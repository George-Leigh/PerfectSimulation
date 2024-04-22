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
#   Northrop, A. R. (2024).  "Algorithms for maximal coupling and long run
#   unbiased simulation with arbitrarily infrequent bias correction".
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
      Qminus[i,] = Mcmc(Q[i,], M, R$RMCMC, sigma, Unorm) # MCMC creates
      #   pre-jump version of Q.
     } else {
      Qminus[i,] = Q[i,] # No MCMC separate from maximal coupling
     }
     if (i == 1) { # Row 1 of Fig. 2: Jump freely, M-H test and store.
      Qstar[1,] = MaxCoupleJump(Qminus[1,], r, R$Rdir[1:d], R$Rmag)
      MHResult = MHTest(Qminus[1,], Qstar[1,], R$RMH, Unorm)
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
      MHResult = MHTest(Qminus[i,], Qstar[i,], R$RMH, Unorm)
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
   # The following 5 lines are not in the pseudocode.  The pseudocode makes
   #   these assignments lower down when i = j + 1 and a[i] > 0, so they are
   #   not required here.  However, some of the assignments made in the
   #   pseudocode are not actually required.  This code makes only necessary
   #   assignments, so is structured slightly differently.  See also the
   #   comment below about "Different structure to pseudocode in paper", where
   #   we subset the rows, which means we don't need to follow the pseudocode
   #   in setting Qminus[i,] when a[i] > 0 below.  Instead, we set only Q.
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
      Qminus[i,] = Mcmc(Q[i,], M, R$RMCMC, sigma, Unorm) # MCMC creates
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

#******************** d = 1 (rho = 1.65, B = 2, 5, 10)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 1.25, Seed = 1,
 lStartTarg = TRUE)
#  60769 20249 8969 4550 2425 1359 764 415 229 115 68 40 22 12 6 3 2 2 1 (19)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 1.5, Seed = 1,
 lStartTarg = TRUE)
#  62471 21221 8642 3890 1854 932 498 238 122 65 30 16 11 6 3 1 (16)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 1.55, Seed = 1,
 lStartTarg = TRUE)
#  62560 21429 8677 3790 1765 883 467 215 106 58 25 10 7 4 3 1 (16)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 1.6, Seed = 1,
 lStartTarg = TRUE)
#  62747 21601 8579 3709 1715 831 418 199 96 49 26 15 8 4 3 (15)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 1.65, Seed = 1,
 lStartTarg = TRUE)
#* 62820 21711 8590 3675 1669 779 395 187 92 45 22 10 2 2 1 (15)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 1.7, Seed = 1,
 lStartTarg = TRUE)
#  62630 21921 8634 3670 1673 766 366 178 83 40 24 9 3 2 1 (15)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 1.75, Seed = 1,
 lStartTarg = TRUE)
#  62363 22105 8778 3652 1641 753 368 170 89 46 21 7 3 2 2 (15)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 2, Seed = 1,
 lStartTarg = TRUE)
#  61001 22918 9407 3846 1610 669 286 142 67 30 13 4 3 1 1 1 1 (17)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#  58804 23572 10047 4260 1857 833 344 151 74 37 12 6 2 1 (14)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  56756 23891 10688 4735 2137 995 445 190 78 39 23 12 6 3 1 1 (16)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#  54542 24105 11360 5274 2462 1226 545 248 122 56 31 15 7 4 2 1 (16)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  52497 24192 11904 5731 2830 1417 699 332 185 97 56 26 18 7 4 3 1 1 (18)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 3.25, Seed = 1,
 lStartTarg = TRUE)
#  50636 24088 12407 6231 3281 1651 857 413 218 110 55 24 13 8 5 2 1 (17)
Norm1 = Unbiased(N = 100000, K = 25, B = 2, M = 0, d = 1, rho = 3.5, Seed = 1,
 lStartTarg = TRUE)
#  48693 23953 12872 6711 3598 1934 1036 552 306 164 88 39 23 16 8 4 3 (17)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 100000, K = 20, B = 5, M = 0, d = 1, rho = 1.65, Seed = 1,
 lStartTarg = TRUE)
#* 91235 7538 1058 144 23 2 (6)
Norm1 = Unbiased(N = 100000, K = 20, B = 6, M = 0, d = 1, rho = 1.65, Seed = 1,
 lStartTarg = TRUE)
#  94218 5198 526 47 9 1 1 (7)

#***** Dispersed starting points
Norm1 = Unbiased(N = 100000, K = 20, B = 5, M = 0, d = 1, rho = 1.65, Seed = 1)
#  62479 26471 8427 2049 463 90 18 2 1
#  Note: Much too low!
Norm1 = Unbiased(N = 100000, K = 20, B = 10, M = 0, d = 1, rho = 1.65, Seed = 1)
#* 90205 9315 461 19 (4)
Norm1 = Unbiased(N = 1e8, K = 20, B = 10, M = 0, d = 1, rho = 1.65, Seed = 1,
 nProg = 1000)
#  90274783 9243138 464580 16916 560 21 2 (7)
save(Norm1, file = "Norm1.RData")

#******************** d = 2 (rho = 2.1, B = 3, 10, 18)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 1.6, Seed = 1,
 lStartTarg = TRUE)
#  44687 22104 13231 7946 4846 2959 1760 1011 595 337 195 125 82 46 35 23 10 3
#   1 1 2 1 (22)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 1.65, Seed = 1,
 lStartTarg = TRUE)
#  45313 22339 13210 7781 4606 2806 1640 947 562 316 186 112 74 42 25 19 12 6
#   2 1 1 (21)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 1.7, Seed = 1,
 lStartTarg = TRUE)
#  45842 22528 13099 7662 4454 2699 1548 884 523 310 173 104 71 41 25 16 11 5
#   2 2 1 (21)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 1.75, Seed = 1,
 lStartTarg = TRUE)
#  46355 22705 12999 7455 4353 2609 1479 828 489 277 162 105 61 45 29 17 11 8
#   6 3 1 1 1 1 (24)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 1.8, Seed = 1,
 lStartTarg = TRUE)
#  46681 22854 12956 7432 4294 2484 1423 796 463 258 142 83 57 39 18 10 5 3 1
#   1 (20)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 1.85, Seed = 1,
 lStartTarg = TRUE)
#  46850 23051 12902 7376 4237 2434 1363 745 455 266 139 78 49 24 13 7 3 3 3 1
#   1 (21)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 2, Seed = 1,
 lStartTarg = TRUE)
#  47806 23416 12877 7150 3936 2170 1199 631 363 190 118 70 31 17 10 7 6 3 (18)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 2.1, Seed = 1,
 lStartTarg = TRUE)
#* 47926 23620 12839 7046 3890 2154 1163 622 336 181 90 54 29 18 12 6 5 3 1 1
#   2 1 1 (23)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 2.15, Seed = 1,
 lStartTarg = TRUE)
#  47845 23637 12899 7038 3879 2170 1152 625 356 192 87 49 29 17 9 7 5 1 1 1 1
#   (21)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#  47524 23793 13139 7164 3906 2076 1094 594 323 177 98 58 25 16 9 3 1 (17)
Norm1 = Unbiased(N = 100000, K = 25, B = 3, M = 0, d = 2, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  46291 23753 13491 7449 4089 2216 1243 703 371 192 95 49 27 16 10 4 1 (17)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 100000, K = 25, B = 10, M = 0, d = 2, rho = 2.1, Seed = 1,
 lStartTarg = TRUE)
#* 89749 8953 1125 153 19 1 (6)
Norm1 = Unbiased(N = 100000, K = 25, B = 12, M = 0, d = 2, rho = 2.1, Seed = 1,
 lStartTarg = TRUE)
#  93408 6037 497 55 3 (5)

#***** Dispersed starting points
Norm1 = Unbiased(N = 100000, K = 25, B = 18, M = 0, d = 2, rho = 2.1, Seed = 1)
#* 90403 9307 278 12 (4)
Norm2 = Unbiased(N = 1e7, K = 25, B = 18, M = 0, d = 2, rho = 2.1, Seed = 1,
 nProg = 1000)
#  9045549 922010 31517 910 13 1 (6)
save(Norm2, file = "Norm2.RData")

#******************** d = 3 (rho = 2.25, B = 6, 17, 27)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 100000, K = 25, B = 5, M = 0, d = 3, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#  44534 24027 13931 7908 4386 2355 1294 679 394 220 121 69 37 20 13 7 3 2 (18)
Norm1 = Unbiased(N = 100000, K = 25, B = 6, M = 0, d = 3, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#* 51452 24453 12289 6084 2986 1413 668 322 169 83 41 16 10 6 4 2 2 (17)
Norm1 = Unbiased(N = 100000, K = 25, B = 7, M = 0, d = 3, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#  58244 23940 10351 4317 1797 748 330 157 70 31 9 3 2 1 (14)

Norm1 = Unbiased(N = 100000, K = 25, B = 6, M = 0, d = 3, rho = 2, Seed = 1,
 lStartTarg = TRUE)
#  50076 24369 12730 6427 3273 1653 797 359 173 69 36 17 8 6 5 2 (16)
Norm1 = Unbiased(N = 100000, K = 25, B = 6, M = 0, d = 3, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  51265 24608 12369 6115 2952 1425 673 312 147 68 33 13 9 5 4 2 (16)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 100000, K = 25, B = 17, M = 0, d = 3, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#* 89435 9207 1171 162 24 1
Norm1 = Unbiased(N = 100000, K = 25, B = 18, M = 0, d = 3, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#  90712 8208 935 127 17 1 (6)

#***** Dispersed starting points
Norm1 = Unbiased(N = 100000, K = 25, B = 27, M = 0, d = 3, rho = 2.25, Seed = 1)
#* 89354 10210 414 22 (4)
Norm3 = Unbiased(N = 1e7, K = 25, B = 27, M = 0, d = 3, rho = 2.25, Seed = 1,
 nProg = 1000)
#  8944150 1015407 38986 1401 55 1 (6)
save(Norm3, file = "Norm3.RData")

#******************** d = 4 (rho = 2.5, B = 9, 27, 39)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 50000, K = 25, B = 9, M = 0, d = 4, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#* 24013 12525 6669 3371 1688 876 444 227 96 53 25 6 4 1 1 1 (16)
Norm1 = Unbiased(N = 50000, K = 25, B = 10, M = 0, d = 4, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#  26144 12673 6051 2780 1262 592 266 124 57 25 14 5 2 2 2 1 (16)

Norm1 = Unbiased(N = 50000, K = 25, B = 9, M = 0, d = 4, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#*  24360 12671 6545 3292 1628 776 384 178 87 40 17 11 4 2 2 1 1 1 (18)
Norm1 = Unbiased(N = 50000, K = 25, B = 9, M = 0, d = 4, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  22949 12543 6941 3751 1905 940 461 247 136 68 30 17 6 4 2 (15)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 50000, K = 20, B = 25, M = 0, d = 4, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  43985 5297 630 74 14 (5)
Norm1 = Unbiased(N = 50000, K = 20, B = 27, M = 0, d = 4, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#* 45016 4465 462 52 5 (5)

#***** Dispersed starting points
Norm1 = Unbiased(N = 50000, K = 20, B = 39, M = 0, d = 4, rho = 2.5, Seed = 1)
#* 45003 4803 188 6 (4)
Norm4 = Unbiased(N = 1e7, K = 20, B = 39, M = 0, d = 4, rho = 2.5, Seed = 1,
 nProg = 1000)
#  9016044 943434 38897 1557 65 3 (6)
save(Norm4, file = "Norm4.RData")

#******************** d = 5 (rho = 2.5, B = 14, 40, 56)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 50000, K = 25, B = 14, M = 0, d = 5, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#* 24528 13125 6655 3104 1396 635 297 148 62 27 15 4 2 1 1 (15)
Norm1 = Unbiased(N = 50000, K = 25, B = 15, M = 0, d = 5, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  26237 13099 6082 2628 1165 488 187 74 30 8 2 (11)

Norm1 = Unbiased(N = 50000, K = 25, B = 14, M = 0, d = 5, rho = 2.25, Seed = 1,
 lStartTarg = TRUE)
#  23823 13360 6802 3195 1500 715 315 147 71 40 16 10 4 2 (14)
Norm1 = Unbiased(N = 50000, K = 25, B = 14, M = 0, d = 5, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#  24361 13093 6697 3179 1436 670 304 142 64 27 11 6 4 4 2

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 50000, K = 20, B = 40, M = 0, d = 5, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#* 44935 4539 473 48 5 (5)
Norm1 = Unbiased(N = 50000, K = 20, B = 41, M = 0, d = 5, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  45109 4435 418 33 5 (5)

#***** Dispersed starting points
Norm1 = Unbiased(N = 50000, K = 20, B = 52, M = 0, d = 5, rho = 2.5, Seed = 1)
#  39965 9524 489 21 1 (5)
Norm1 = Unbiased(N = 50000, K = 20, B = 55, M = 0, d = 5, rho = 2.5, Seed = 1)
#* 45304 4480 209 7 (4)
Norm1 = Unbiased(N = 50000, K = 20, B = 56, M = 0, d = 5, rho = 2.5, Seed = 1)
#  45570 4248 172 10 (4)
Norm5 = Unbiased(N = 1e7, K = 20, B = 55, M = 0, d = 5, rho = 2.5, Seed = 1,
 nProg = 1000)
#  9070272 889744 38263 1656 62 3 (6)
save(Norm5, file = "Norm5.RData")

#******************** d = 6 (rho = 2.75, B = 20, 60, 75)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 25000, K = 25, B = 20, M = 0, d = 6, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#* 11920 6687 3435 1616 750 330 148 68 31 9 3 1 1 1 (14)

Norm1 = Unbiased(N = 25000, K = 25, B = 20, M = 0, d = 6, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  11741 6884 3469 1612 730 308 149 63 25 9 4 3 1 1 1 (15)
Norm1 = Unbiased(N = 25000, K = 25, B = 20, M = 0, d = 6, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  11773 6676 3391 1645 781 376 186 92 38 22 11 6 2 1

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 25000, K = 20, B = 55, M = 0, d = 6, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#  22114 2554 299 30 2 1 (6)
Norm1 = Unbiased(N = 25000, K = 20, B = 60, M = 0, d = 6, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#* 22642 2134 206 17 1 (5)

#***** Dispersed starting points
Norm1 = Unbiased(N = 25000, K = 20, B = 75, M = 0, d = 6, rho = 2.75, Seed = 1)
#* 22683 2211 100 6 (4)
Norm6 = Unbiased(N = 1e7, K = 20, B = 75, M = 0, d = 6, rho = 2.75, Seed = 1,
 nProg = 1000)
#  9088313 867525 42083 1987 87 5 (6)
save(Norm6, file = "Norm6.RData")

#******************** d = 7 (rho = 2.75, B = 30, 80, 95)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 25000, K = 25, B = 30, M = 0, d = 7, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#* 12565 6815 3217 1419 595 228 89 39 21 8 2 1 1 (13)

Norm1 = Unbiased(N = 25000, K = 25, B = 30, M = 0, d = 7, rho = 2.5, Seed = 1,
 lStartTarg = TRUE)
#  12150 7009 3335 1469 616 242 105 45 18 9 1 1 (12)
Norm1 = Unbiased(N = 25000, K = 25, B = 30, M = 0, d = 7, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  12383 6826 3260 1462 618 245 107 49 28 12 7 2 1 (13)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 25000, K = 20, B = 80, M = 0, d = 7, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#* 22390 2357 228 21 4 (5)
Norm1 = Unbiased(N = 25000, K = 20, B = 85, M = 0, d = 7, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#  22711 2121 154 13 1 (5)

#***** Dispersed starting points
Norm1 = Unbiased(N = 25000, K = 20, B = 95, M = 0, d = 7, rho = 2.75, Seed = 1)
#* 22383 2457 150 10 (4)
Norm7 = Unbiased(N = 1e7, K = 20, B = 95, M = 0, d = 7, rho = 2.75, Seed = 1,
 nProg = 1000)
#  8920914 1011664 63204 3972 231 14 1 (7)
save(Norm7, file = "Norm7.RData")

#******************** d = 8 (rho = 3, B = 45, 115, 140)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 25000, K = 25, B = 45, M = 0, d = 8, rho = 2.75, Seed = 1,
 lStartTarg = TRUE)
#  12979 6986 3021 1252 470 173 73 26 10 5 3 1 1 (13)
Norm1 = Unbiased(N = 25000, K = 25, B = 45, M = 0, d = 8, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#* 13053 6880 3007 1256 499 184 71 34 10 3 2 1

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 25000, K = 20, B = 115, M = 0, d = 8, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#* 22517 2247 215 20 1 (5)
Norm1 = Unbiased(N = 25000, K = 20, B = 120, M = 0, d = 8, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  22722 2075 184 19 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 25000, K = 20, B = 140, M = 0, d = 8, rho = 3, Seed = 1)
#* 22746 2149 101 4 (4)
Norm8 = Unbiased(N = 1e7, K = 20, B = 140, M = 0, d = 8, rho = 3, Seed = 1,
 nProg = 1000)
#  9146662 806424 44332 2443 131 7 1 (7)
save(Norm8, file = "Norm8.RData")

#******************** d = 9 (rho = 3, B = 60, 160, 180)

#***** Set rho (P ~ 0.5)
gNorm1 = Unbiased(N = 25000, K = 25, B = 60, M = 0, d = 9, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#* 12542 6990 3169 1339 551 243 96 37 21 8 3 1 (12)
gNorm1 = Unbiased(N = 25000, K = 25, B = 65, M = 0, d = 9, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  13584 6916 2811 1072 373 147 62 26 5 2 1 1 (12)

gNorm1 = Unbiased(N = 25000, K = 25, B = 60, M = 0, d = 9, rho = 3.25, Seed = 1,
 lStartTarg = TRUE)
# 12403 7062 3187 1361 567 239 105 42 21 9 3 1

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 25000, K = 20, B = 155, M = 0, d = 9, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  22320 2410 244 24 2 (5)
Norm1 = Unbiased(N = 25000, K = 20, B = 160, M = 0, d = 9, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#* 22563 2207 213 17 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 25000, K = 20, B = 180, M = 0, d = 9, rho = 3, Seed = 1)
#* 22560 2281 145 11 2 1 (6)
Norm9 = Unbiased(N = 1e6, K = 20, B = 180, M = 0, d = 9, rho = 3, Seed = 1,
 nProg = 1000)
#* 902946 90542 6080 404 25 3 (6)
save(Norm9, file = "Norm9.RData")

#******************** d = 10 (rho = 3.25, B = 85, 230, 240)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 15000, K = 25, B = 85, M = 0, d = 10, rho = 3, Seed = 1,
 lStartTarg = TRUE)
#  7505 4269 1861 784 321 140 61 31 14 8 4 2 (12)
Norm1 = Unbiased(N = 15000, K = 25, B = 85, M = 0, d = 10, rho = 3.25, Seed = 1,
 lStartTarg = TRUE)
#* 7509 4181 1872 816 352 157 73 29 8 1 1 1 (12)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 15000, K = 20, B = 220, M = 0, d = 10, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  13408 1440 137 11 3 1 (6)
Norm1 = Unbiased(N = 15000, K = 20, B = 230, M = 0, d = 10, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#* 13569 1294 121 14 2 (5)

#***** Dispersed starting points
Norm1 = Unbiased(N = 15000, K = 20, B = 240, M = 0, d = 10, rho = 3.25,
 Seed = 1)
#* 13419 1463 112 6 (4)
Norm1 = Unbiased(N = 15000, K = 20, B = 250, M = 0, d = 10, rho = 3.25,
 Seed = 1)
#  13595 1305 88 11 1 (5)
Norm10 = Unbiased(N = 1e6, K = 20, B = 250, M = 0, d = 10, rho = 3.25,
 Seed = 1, nProg = 1000)
#* 905825 87490 6172 466 45 2 (6)
save(Norm10, file = "Norm10.RData")

#******************** d = 11 (rho = 3.25, B = 120, 310, 330)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 10000, K = 25, B = 120, M = 0, d = 11, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#* 5224 2778 1220 472 184 72 30 12 6 2 (10)
Norm1 = Unbiased(N = 10000, K = 25, B = 120, M = 0, d = 11, rho = 3.5, Seed = 1,
 lStartTarg = TRUE)
#  4998 2845 1235 520 232 96 43 19 7 2 1 1 1 (13)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 10000, K = 20, B = 300, M = 0, d = 11, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  8937 958 97 7 1 (5)
Norm1 = Unbiased(N = 10000, K = 20, B = 310, M = 0, d = 11, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#*  9032 883 77 8 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 10000, K = 20, B = 330, M = 0, d = 11, rho = 3.25,
 Seed = 1)
#* 8970 952 68 9 1 (5)
Norm11 = Unbiased(N = 1e6, K = 20, B = 330, M = 0, d = 11, rho = 3.25,
 Seed = 1, nProg = 1000)
#  897656 94385 7352 547 53 6 1 (7)
save(Norm11, file = "Norm11.RData")

#******************** d = 12 (rho = 3.25, B = 170, 420, 460)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 5000, K = 25, B = 170, M = 0, d = 12, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#* 2583 1394 600 259 105 38 16 4 1 (9)
Norm1 = Unbiased(N = 5000, K = 25, B = 170, M = 0, d = 12, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  2581 1391 612 249 94 49 16 5 2 1 (10)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 5000, K = 20, B = 400, M = 0, d = 12, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  4409 527 59 5 (4)
Norm1 = Unbiased(N = 5000, K = 20, B = 420, M = 0, d = 12, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#* 4470 475 48 7 (4)
Norm1 = Unbiased(N = 5000, K = 20, B = 480, M = 0, d = 12, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  4656 325 19 (3)
Norm1 = Unbiased(N = 5000, K = 20, B = 500, M = 0, d = 12, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  4666 310 23 1 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 5000, K = 20, B = 460, M = 0, d = 12, rho = 3.25,
 Seed = 1)
#* 4518 446 31 5 (4)
Norm12 = Unbiased(N = 1e6, K = 20, B = 460, M = 0, d = 12, rho = 3.25,
 Seed = 1, nProg = 1000)
#  901613 90569 7196 565 48 9 (6)
save(Norm12, file = "Norm12.RData")

#******************** d = 13 (rho = 3.5, B = 250, 650, 700)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 3000, K = 25, B = 250, M = 0, d = 13, rho = 3.25,
 Seed = 1, lStartTarg = TRUE)
#  1531 864 363 150 56 26 7 1 1 1
Norm1 = Unbiased(N = 3000, K = 25, B = 250, M = 0, d = 13, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#* 1557 863 360 135 50 23 9 3 (8)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 3000, K = 20, B = 600, M = 0, d = 13, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  2670 297 30 3 (4)
Norm1 = Unbiased(N = 3000, K = 20, B = 650, M = 0, d = 13, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#* 2718 253 28 1 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 3000, K = 20, B = 650, M = 0, d = 13, rho = 3.5,
 Seed = 1)
#* 2666 311 23 (3)
Norm1 = Unbiased(N = 3000, K = 20, B = 700, M = 0, d = 13, rho = 3.5,
 Seed = 1)
#  2752 226 19 3 (4)
Norm13 = Unbiased(N = 1e6, K = 20, B = 700, M = 0, d = 13, rho = 3.5,
 Seed = 1, nProg = 1000)
#  922451 72452 4782 294 21 (5)
save(Norm13, file = "Norm13.RData")

#******************** d = 14 (rho = 3.75, B = 330, 860, 920)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 3000, K = 25, B = 330, M = 0, d = 14, rho = 3.5,
 Seed = 1, lStartTarg = TRUE)
#  1491 799 381 183 78 35 20 9 2 1 1 (11)
Norm1 = Unbiased(N = 3000, K = 25, B = 330, M = 0, d = 14, rho = 3.75,
 Seed = 1, lStartTarg = TRUE)
#* 1532 868 361 150 60 19 7 2 1 (9)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 3000, K = 20, B = 860, M = 0, d = 14, rho = 3.75,
 Seed = 1, lStartTarg = TRUE)
#* 2706 271 19 3 1
Norm1 = Unbiased(N = 3000, K = 20, B = 900, M = 0, d = 14, rho = 3.75,
 Seed = 1, lStartTarg = TRUE)
#  2712 262 23 3 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 3000, K = 20, B = 920, M = 0, d = 14, rho = 3.75,
 Seed = 1)
#* 2714 265 20 1 (4)
Norm14 = Unbiased(N = 1e6, K = 20, B = 920, M = 0, d = 14, rho = 3.75,
 Seed = 1, nProg = 1000)
#  907008 85504 6882 558 42 6 (6)
save(Norm14, file = "Norm14.RData")

#******************** d = 15 (rho = 4, B = 450, 1300, 1400)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 2000, K = 25, B = 450, M = 0, d = 15, rho = 3.75,
 Seed = 1, lStartTarg = TRUE)
#  956 583 278 110 44 15 9 3 1 1 (10)
Norm1 = Unbiased(N = 2000, K = 25, B = 450, M = 0, d = 15, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#* 988 556 261 120 43 11 8 5 4 2 1 1 (12)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 2000, K = 20, B = 1200, M = 0, d = 15, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#  1789 191 18 2 (4)
Norm1 = Unbiased(N = 2000, K = 20, B = 1250, M = 0, d = 15, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#  1772 199 26 3 (4)
Norm1 = Unbiased(N = 2000, K = 20, B = 1300, M = 0, d = 15, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#* 1823 152 22 3 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 2000, K = 20, B = 1300, M = 0, d = 15, rho = 4,
 Seed = 1)
#* 1825 164 10 1 (4)
Norm1 = Unbiased(N = 2000, K = 20, B = 1400, M = 0, d = 15, rho = 4,
 Seed = 1)
#  1847 142 11 (3)
Norm15 = Unbiased(N = 1e6, K = 20, B = 1300, M = 0, d = 15, rho = 4,
 Seed = 1, nProg = 1000)
#  908482 84286 6679 507 42 3 1 (7)
save(Norm15, file = "Norm15.RData")

#******************** d = 16 (rho = 4, B = 650, 1700, 1700)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 2000, K = 25, B = 650, M = 0, d = 16, rho = 3.75,
 Seed = 1, lStartTarg = TRUE)
#  1002 550 264 113 43 16 6 3 1 1 1 (11)
Norm1 = Unbiased(N = 2000, K = 25, B = 650, M = 0, d = 16, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#* 1024 565 250 87 40 19 11 3 1 (9)
Norm1 = Unbiased(N = 2000, K = 25, B = 650, M = 0, d = 16, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  956 593 287 98 45 17 3 1 (8)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 2000, K = 20, B = 1600, M = 0, d = 16, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#  1758 212 25 5 (4)
Norm1 = Unbiased(N = 2000, K = 20, B = 1700, M = 0, d = 16, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#* 1776 201 20 3 (4)
Norm1 = Unbiased(N = 2000, K = 20, B = 1800, M = 0, d = 16, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#  1827 161 12 (3)

#***** Dispersed starting points
Norm1 = Unbiased(N = 2000, K = 20, B = 1700, M = 0, d = 16, rho = 4,
 Seed = 1)
#* 1785 192 22 1 (4)
Norm1 = Unbiased(N = 2000, K = 20, B = 1800, M = 0, d = 16, rho = 4,
 Seed = 1)
#  1830 160 10 (3)
Norm16 = Unbiased(N = 1e6, K = 20, B = 1700, M = 0, d = 16, rho = 4,
 Seed = 1, nProg = 1000)
#  892063 97738 9244 868 82 4 1 (7)
save(Norm16, file = "Norm16.RData")

#******************** d = 17 (rho = 4.5, B = 940, 2300, 2400)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 1500, K = 25, B = 940, M = 0, d = 17, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#  756 437 193 62 27 14 5 3 2 1 (10)
Norm1 = Unbiased(N = 1500, K = 25, B = 940, M = 0, d = 17, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  774 390 191 87 42 11 2 1 1 1 (10)
Norm1 = Unbiased(N = 1500, K = 25, B = 940, M = 0, d = 17, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#* 774 428 179 71 28 9 5 4 1 1 (10)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 1500, K = 20, B = 2300, M = 0, d = 17, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  1318 161 18 3 (4)
Norm1 = Unbiased(N = 1500, K = 20, B = 2400, M = 0, d = 17, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  1367 126 7 (3)
Norm1 = Unbiased(N = 1500, K = 20, B = 2300, M = 0, d = 17, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#* 1344 143 10 3

#***** Dispersed starting points
Norm1 = Unbiased(N = 1500, K = 20, B = 2400, M = 0, d = 17, rho = 4.5,
 Seed = 1)
#* 1350 142 6 1 1 (5)
Norm17 = Unbiased(N = 1e6, K = 20, B = 2400, M = 0, d = 17, rho = 4.5,
 Seed = 1, nProg = 1000)
#  892786 96942 9206 957 102 7 (6)
save(Norm17, file = "Norm17.RData")

#******************** d = 18 (rho = 4.5, B = 1350, 3400, 3600)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 1000, K = 25, B = 1350, M = 0, d = 18, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#  526 281 121 46 19 6 1 (7)
Norm1 = Unbiased(N = 1000, K = 25, B = 1350, M = 0, d = 18, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  518 287 123 46 16 8 2 (7)
#  Note: Results are anomalous for d = 18 and P ~ 0.5 (indicate too small a
#   value for rho, compared to other values of d).  Results are consistent for
#   P ~ 0.1, so we have gone with them.
Norm1 = Unbiased(N = 1000, K = 25, B = 1350, M = 0, d = 18, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#  509 291 125 45 20 8 1 1 (8)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 1500, K = 20, B = 3200, M = 0, d = 18, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  1307 173 20 (3)
Norm1 = Unbiased(N = 1500, K = 20, B = 3400, M = 0, d = 18, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  1347 140 12 1 (4)
Norm1 = Unbiased(N = 1500, K = 20, B = 3400, M = 0, d = 18, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#* 1373 116 9 1 1 (5)

#***** Dispersed starting points
Norm1 = Unbiased(N = 1500, K = 20, B = 3600, M = 0, d = 18, rho = 4.5,
 Seed = 1)
#* 1366 120 11 3 (4)
Norm18 = Unbiased(N = 1e6, K = 20, B = 3600, M = 0, d = 18, rho = 4.5,
 Seed = 1, nProg = 1000)
#  914191 79059 6242 474 32 2 (6)
save(Norm18, file = "Norm18.RData")

#******************** d = 19 (rho = 4.5, B = 1750, 4600, 5000)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 1000, K = 25, B = 1750, M = 0, d = 19, rho = 4,
 Seed = 1, lStartTarg = TRUE)
#  414 290 162 72 38 14 4 4 2 (9)
Norm1 = Unbiased(N = 1000, K = 25, B = 1750, M = 0, d = 19, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  489 282 130 57 26 12 4 (7)
Norm1 = Unbiased(N = 1000, K = 25, B = 1750, M = 0, d = 19, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#* 495 285 122 51 23 13 6 4 1 (9)
Norm1 = Unbiased(N = 1000, K = 25, B = 1750, M = 0, d = 19, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#  479 263 132 76 27 13 6 3 1 (9)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 1000, K = 20, B = 4000, M = 0, d = 19, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  842 131 22 4 1 (5)
Norm1 = Unbiased(N = 1000, K = 20, B = 4400, M = 0, d = 19, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#  882 103 14 1 (4)
Norm1 = Unbiased(N = 1000, K = 20, B = 4600, M = 0, d = 19, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#* 886 104 10 (3)

#***** Dispersed starting points
Norm1 = Unbiased(N = 1000, K = 20, B = 5000, M = 0, d = 19, rho = 4.5,
 Seed = 1)
#* 916 78 4 1 1 (5)
Norm19 = Unbiased(N = 100000, K = 20, B = 5000, M = 0, d = 19, rho = 4.5,
 Seed = 1, nProg = 1000)
#  91407 7948 604 37 3 1 (6)
save(Norm19, file = "Norm19.RData")

#******************** d = 20 (rho = 4.75, B = 2600, 6200, 6400)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 1000, K = 25, B = 2600, M = 0, d = 20, rho = 4.25,
 Seed = 1, lStartTarg = TRUE)
#  492 273 130 50 23 15 7 6 2 2 (10)
Norm1 = Unbiased(N = 1000, K = 25, B = 2600, M = 0, d = 20, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#  521 258 119 61 22 13 5 1 (8)
Norm1 = Unbiased(N = 1000, K = 25, B = 2600, M = 0, d = 20, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#* 527 263 123 58 25 3 1 (7)
Norm1 = Unbiased(N = 1000, K = 25, B = 2600, M = 0, d = 20, rho = 5,
 Seed = 1, lStartTarg = TRUE)
#  514 280 119 48 22 9 6 2 (8)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 1000, K = 20, B = 6200, M = 0, d = 20, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#  899 95 5 1 (4)
Norm1 = Unbiased(N = 1000, K = 20, B = 6400, M = 0, d = 20, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#  888 100 12 (3)
#  Note: Probability of coalescence has gone backwards from B = 6200 to 6400,
#   for rho = 4.5, which is nonsensical.  Both sample average probabilities of
#   coalescence are below the target value of 0.9 (although they are very
#   close).  Therefore would choose the higher value of B.  In the end,
#   though, we don't have to worry, because we want rho = 4.75!
Norm1 = Unbiased(N = 1000, K = 20, B = 6200, M = 0, d = 20, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#* 900 93 5 2 (4)

# Check that a smaller value of rho doesn't work.
Norm1 = Unbiased(N = 100, K = 20, B = 6200, M = 0, d = 20, rho = 2,
 Seed = 1, lStartTarg = TRUE)
#  56 25 12 4 2 1 (6)
#  All good, it doesn't work.

#***** Dispersed starting points
Norm1 = Unbiased(N = 1000, K = 20, B = 6400, M = 0, d = 20, rho = 4.75,
 Seed = 1)
#  883 107 10 (3)
Norm20 = Unbiased(N = 100000, K = 20, B = 6400, M = 0, d = 20, rho = 4.75,
 Seed = 1, nProg = 1000)
save(Norm20, file = "Norm20.RData")
#* 89487 9485 926 91 11 (5) 14.5 hours (run on MS Windows)

#******************** d = 21 (rho = 4.75, B = 3400, 8600, 9200)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 1000, K = 25, B = 3400, M = 0, d = 21, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#  470 285 147 64 21 8 4 1 (8)
Norm1 = Unbiased(N = 1000, K = 25, B = 3400, M = 0, d = 21, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#  488 289 144 52 21 3 2 1 (8)
Norm1 = Unbiased(N = 1000, K = 25, B = 3400, M = 0, d = 21, rho = 5,
 Seed = 1, lStartTarg = TRUE)
#* 506 271 119 60 28 9 3 1 2 1 (10)
Norm1 = Unbiased(N = 1000, K = 25, B = 3400, M = 0, d = 21, rho = 5.25,
 Seed = 1, lStartTarg = TRUE)
#  500 291 128 54 19 7 1 (7)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 1000, K = 20, B = 8600, M = 0, d = 21, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#* 902 84 12 2 (4)
Norm1 = Unbiased(N = 1000, K = 20, B = 9000, M = 0, d = 21, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#  901 93 6 (3)
Norm1 = Unbiased(N = 1000, K = 20, B = 8400, M = 0, d = 21, rho = 5,
 Seed = 1, lStartTarg = TRUE)
#  877 105 18 (3)

#***** Dispersed starting points
Norm1 = Unbiased(N = 1000, K = 20, B = 9200, M = 0, d = 21, rho = 4.75,
 Seed = 1)
#  894 92 10 3 1 (5)
Norm1 = Unbiased(N = 1000, K = 20, B = 9200, M = 0, d = 21, rho = 5,
 Seed = 1)
#  880 108 12 (3)
Norm21 = Unbiased(N = 100000, K = 20, B = 9200, M = 0, d = 21, rho = 4.75,
 Seed = 1, nProg = 1000)
save(Norm21, file = "Norm21.RData")
#* 90486 8686 761 59 7 1 (6) 20.5 hours (run on MS Windows)

#******************** d = 22 (rho = 5.25, B = 5500, 12000, 1300)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 1000, K = 25, B = 5500, M = 0, d = 22, rho = 4.5,
 Seed = 1, lStartTarg = TRUE)
#  552 269 107 42 17 6 4 2 1 (9)
Norm1 = Unbiased(N = 1000, K = 25, B = 5500, M = 0, d = 22, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#  585 253 97 36 20 5 2 2 (8)
Norm1 = Unbiased(N = 1000, K = 25, B = 5500, M = 0, d = 22, rho = 5,
 Seed = 1, lStartTarg = TRUE)
#  594 247 92 42 16 5 2 1 1
Norm1 = Unbiased(N = 1000, K = 25, B = 5500, M = 0, d = 22, rho = 5.25,
 Seed = 1, lStartTarg = TRUE)
#* 605 257 95 31 11 1
Norm1 = Unbiased(N = 1000, K = 25, B = 5500, M = 0, d = 22, rho = 5.5,
 Seed = 1, lStartTarg = TRUE)
#  559 269 105 44 19 2 2 (7)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 1000, K = 20, B = 13000, M = 0, d = 22, rho = 4.75,
 Seed = 1, lStartTarg = TRUE)
#  888 101 11 (3)
Norm1 = Unbiased(N = 1000, K = 20, B = 14000, M = 0, d = 22, rho = 5,
 Seed = 1, lStartTarg = TRUE)
#  934 60 6 (3)
Norm1 = Unbiased(N = 1000, K = 20, B = 12000, M = 0, d = 22, rho = 5.25,
 Seed = 1, lStartTarg = TRUE)
#* 910 84 6 (3)

#***** Dispersed starting points
Norm1 = Unbiased(N = 1000, K = 20, B = 13000, M = 0, d = 22, rho = 5.25,
 Seed = 1)
#* 919 76 5 (3)
Norm22 = Unbiased(N = 100000, K = 20, B = 13000, M = 0, d = 22, rho = 5.25,
 Seed = 1, nProg = 1000)
save(Norm22, file = "Norm22.RData")
#* 91176 8099 665 54 6 (5)

#******************** d = 23 (rho = 5.25, B = 8000, 18000, 18000)

#***** Set rho (P ~ 0.5)
Norm1 = Unbiased(N = 1000, K = 25, B = 8000, M = 0, d = 23, rho = 5.25,
 Seed = 1, lStartTarg = TRUE)
#* 642 247 79 23 8 1 (6)
Norm1 = Unbiased(N = 1000, K = 25, B = 8000, M = 0, d = 23, rho = 5.5,
 Seed = 1, lStartTarg = TRUE)
#  589 267 98 30 11 4 1 (7)
Norm1 = Unbiased(N = 1000, K = 25, B = 8000, M = 0, d = 23, rho = 5.75,
 Seed = 1, lStartTarg = TRUE)
#  590 272 87 34 8 6 2 1 (8)

#***** Results (P ~ 0.1)
Norm1 = Unbiased(N = 1000, K = 20, B = 18000, M = 0, d = 23, rho = 5.25,
 Seed = 1, lStartTarg = TRUE)
#* 919 78 3 (3)
Norm1 = Unbiased(N = 1000, K = 20, B = 18000, M = 0, d = 23, rho = 5.5,
 Seed = 1, lStartTarg = TRUE)
#  904 89 6 1 (4)

#***** Dispersed starting points
Norm1 = Unbiased(N = 1000, K = 20, B = 19000, M = 0, d = 23, rho = 5.25,
 Seed = 1)
#  927 69 4 (3)
Norm1 = Unbiased(N = 1000, K = 20, B = 18000, M = 0, d = 23, rho = 5.25,
 Seed = 1)
#* 914 80 5 1 (4)
Norm23 = Unbiased(N = 100000, K = 20, B = 18000, M = 0, d = 23, rho = 5.25,
 Seed = 1, nProg = 1000)
save(Norm23, file = "Norm23.RData")
#* 91378 7905 658 56 3 (5)

#**************************************** Results for the very widely
#   dispersed starting point, using the uniform (-6, 6) distribution for each
#   coordinate of the starting points: these are very extreme points in high
#   numbers of dimensions.  In many practical cases, likelihoods would be
#   undefined at these points.  I am not planning to put these into the paper.

#******************** d = 1 (rho = 3.25, B = 3, 8)
Norm1 = Unbiased(N = 100000, K = 20, B = 3, M = 0, d = 1, rho = 2, Seed = 1)
#  45413 28039 15177 6901 2819 1073 401 125 36 15 1 (11)
Norm1 = Unbiased(N = 100000, K = 20, B = 3, M = 0, d = 1, rho = 2.5, Seed = 1)
#  50635 27603 13125 5390 2077 758 265 93 37 14 2 1 (12)
Norm1 = Unbiased(N = 100000, K = 20, B = 3, M = 0, d = 1, rho = 3, Seed = 1)
#  52559 26843 12353 5116 1998 720 257 88 45 14 7 (11)
Norm1 = Unbiased(N = 100000, K = 20, B = 3, M = 0, d = 1, rho = 3.25, Seed = 1)
#* 52581 26648 12305 5173 2066 765 303 98 39 16 6 (11)
Norm1 = Unbiased(N = 100000, K = 20, B = 3, M = 0, d = 1, rho = 3.5, Seed = 1)
#  52344 26441 12296 5307 2173 863 345 142 57 25 4 2 1 (13)
Norm1 = Unbiased(N = 100000, K = 20, B = 3, M = 0, d = 1, rho = 4, Seed = 1)
#  51523 25803 12460 5653 2557 1115 500 217 94 42 19 11 3 1 1 1 (16)

Norm1 = Unbiased(N = 100000, K = 20, B = 8, M = 0, d = 1, rho = 3.25, Seed = 1)
#* 89619 9618 716 42 5 (5)

#******************** d = 2 (rho = 3.5, B = 8, 19)
Norm1 = Unbiased(N = 50000, K = 20, B = 5, M = 0, d = 2, rho = 3.25, Seed = 1)
#  16542 14396 9466 5139 2439 1140 503 223 98 37 9 5 2 1 (14)
Norm1 = Unbiased(N = 50000, K = 20, B = 7, M = 0, d = 2, rho = 3.25, Seed = 1)
#  23362 15958 7000 2458 840 269 78 24 8 2 1 (11)
Norm1 = Unbiased(N = 50000, K = 20, B = 8, M = 0, d = 2, rho = 3.25, Seed = 1)
#* 26515 15665 5586 1639 433 117 30 12 3 (9)

Norm1 = Unbiased(N = 50000, K = 20, B = 8, M = 0, d = 2, rho = 3, Seed = 1)
#  26466 15971 5483 1562 398 102 14 4 (8)
Norm1 = Unbiased(N = 50000, K = 20, B = 8, M = 0, d = 2, rho = 3.5, Seed = 1)
#* 26606 15329 5622 1736 501 150 39 9 4 3 1 (11)
Norm1 = Unbiased(N = 50000, K = 20, B = 8, M = 0, d = 2, rho = 3.75, Seed = 1)
#  26475 14896 5736 1938 612 223 77 29 8 4 1 1 (12)

Norm1 = Unbiased(N = 50000, K = 20, B = 18, M = 0, d = 2, rho = 3.5, Seed = 1)
#  44198 5411 360 29 2 (5)
Norm1 = Unbiased(N = 50000, K = 20, B = 19, M = 0, d = 2, rho = 3.5, Seed = 1)
#* 44951 4781 257 11 (4)
Norm1 = Unbiased(N = 50000, K = 20, B = 20, M = 0, d = 2, rho = 3.5, Seed = 1)
#  45740 4031 214 14 1 (5)

#******************** d = 3 (rho = 3.25, B = 14, 32)
Norm1 = Unbiased(N = 30000, K = 20, B = 12, M = 0, d = 3, rho = 3.5, Seed = 1)
#  12391 10616 4593 1591 558 171 51 21 8 (9)
Norm1 = Unbiased(N = 30000, K = 20, B = 14, M = 0, d = 3, rho = 3.5, Seed = 1)
#* 14524 10535 3538 1019 280 75 21 5 2 1 (10)
Norm1 = Unbiased(N = 30000, K = 20, B = 15, M = 0, d = 3, rho = 3.5, Seed = 1)
#  15749 10210 3045 761 179 44 9 2 1 (9)

Norm1 = Unbiased(N = 30000, K = 20, B = 14, M = 0, d = 3, rho = 3, Seed = 1)
#  14518 11215 3302 780 148 28 7 2 (8)
Norm1 = Unbiased(N = 30000, K = 20, B = 14, M = 0, d = 3, rho = 3.25, Seed = 1)
#* 14682 10711 3443 877 217 52 15 3 (8)

Norm1 = Unbiased(N = 30000, K = 20, B = 25, M = 0, d = 3, rho = 3.25, Seed = 1)
#  24392 5165 410 32 1 (5)
Norm1 = Unbiased(N = 30000, K = 20, B = 30, M = 0, d = 3, rho = 3.25, Seed = 1)
#  26522 3332 138 8 (4)
Norm1 = Unbiased(N = 30000, K = 20, B = 32, M = 0, d = 3, rho = 3.25, Seed = 1)
#* 27056 2835 103 6 (4)

#******************** d = 4 (rho = 3.25, B = 23, 47)
Norm1 = Unbiased(N = 30000, K = 20, B = 22, M = 0, d = 4, rho = 3.25, Seed = 1)
#  14402 11602 3192 660 114 27 3 (7)
Norm1 = Unbiased(N = 30000, K = 20, B = 23, M = 0, d = 4, rho = 3.25, Seed = 1)
#* 15251 11236 2799 565 109 33 6 1 (8)

Norm1 = Unbiased(N = 30000, K = 20, B = 23, M = 0, d = 4, rho = 3, Seed = 1)
#  15209 11568 2635 494 82 11 1 (7) 
Norm1 = Unbiased(N = 30000, K = 20, B = 23, M = 0, d = 4, rho = 3.5, Seed = 1)
#  15161 10912 3025 684 165 38 11 3 1 (9)

Norm1 = Unbiased(N = 30000, K = 20, B = 47, M = 0, d = 4, rho = 3.25, Seed = 1)
#* 26916 2974 105 5 (4)
Norm1 = Unbiased(N = 30000, K = 20, B = 48, M = 0, d = 4, rho = 3.25, Seed = 1)
#  27191 2713 93 3 (4)
Norm1 = Unbiased(N = 30000, K = 20, B = 55, M = 0, d = 4, rho = 3.25, Seed = 1)
#  28284 1681 34 1 (4)

#******************** d = 5 (rho = 3, B = 33, 65)
Norm1 = Unbiased(N = 20000, K = 20, B = 33, M = 0, d = 5, rho = 3.25, Seed = 1)
#* 9950 7843 1786 338 67 13 3 (7)

Norm1 = Unbiased(N = 20000, K = 20, B = 33, M = 0, d = 5, rho = 2.75, Seed = 1)
#  9968 8194 1535 246 47 10 (6)
Norm1 = Unbiased(N = 20000, K = 20, B = 33, M = 0, d = 5, rho = 3, Seed = 1)
#* 10161 7891 1622 265 52 8 1 (7)
Norm1 = Unbiased(N = 20000, K = 20, B = 33, M = 0, d = 5, rho = 3.5, Seed = 1)
#  9791 7658 1977 444 98 25 5 2 (8)

Norm1 = Unbiased(N = 20000, K = 20, B = 65, M = 0, d = 5, rho = 3, Seed = 1)
#* 18012 1930 53 5 (4)

#******************** d = 6 (rho = 3, B > 45, B = 87)
Norm1 = Unbiased(N = 20000, K = 20, B = 45, M = 0, d = 6, rho = 3, Seed = 1)
#* 9872 8144 1630 299 45 9 1 (7)

Norm1 = Unbiased(N = 20000, K = 20, B = 45, M = 0, d = 6, rho = 2.75, Seed = 1)
#  9768 8305 1596 280 47 3 1 (7)
Norm1 = Unbiased(N = 20000, K = 20, B = 45, M = 0, d = 6, rho = 3.25, Seed = 1)
#  9888 7965 1718 330 81 17 1 (7)
Norm1 = Unbiased(N = 20000, K = 20, B = 45, M = 0, d = 6, rho = 3.5, Seed = 1)
#  9743 7740 1949 443 97 21 5 2 (8)

Norm1 = Unbiased(N = 20000, K = 20, B = 85, M = 0, d = 6, rho = 3, Seed = 1)
#  17783 2123 90 4 (4)
Norm1 = Unbiased(N = 20000, K = 20, B = 87, M = 0, d = 6, rho = 3, Seed = 1)
#* 17956 1969 72 3 (4)
Norm1 = Unbiased(N = 20000, K = 20, B = 88, M = 0, d = 6, rho = 3, Seed = 1)
#  18063 1866 67 4 (4)

#******************** d = 7 (rho = 3.25, B = 60, B = 120)
Norm1 = Unbiased(N = 15000, K = 20, B = 60, M = 0, d = 7, rho = 3, Seed = 1)
#  7575 5914 1217 231 50 11 2 (7)
Norm1 = Unbiased(N = 15000, K = 20, B = 60, M = 0, d = 7, rho = 3.25, Seed = 1)
#* 7609 5780 1283 247 65 14 2 (7)
Norm1 = Unbiased(N = 15000, K = 20, B = 60, M = 0, d = 7, rho = 3.5, Seed = 1)
#  7285 5848 1455 323 74 13 2 (7)

Norm1 = Unbiased(N = 15000, K = 20, B = 115, M = 0, d = 7, rho = 3.25, Seed = 1)
#  13304 1627 65 4 (4)
Norm1 = Unbiased(N = 15000, K = 20, B = 120, M = 0, d = 7, rho = 3.25, Seed = 1)
#* 13514 1438 46 2 (4)

#******************** d = 8 (rho = 3, B = 78, B = 160)
Norm1 = Unbiased(N = 10000, K = 20, B = 78, M = 0, d = 8, rho = 2.75, Seed = 1)
#  4970 3937 856 190 42 5 (6)
Norm1 = Unbiased(N = 10000, K = 20, B = 78, M = 0, d = 8, rho = 3, Seed = 1)
#* 5048 3868 855 180 38 7 4 (7)
Norm1 = Unbiased(N = 10000, K = 20, B = 78, M = 0, d = 8, rho = 3.25, Seed = 1)
#  5012 3842 890 196 46 12 1 1 (8)

Norm1 = Unbiased(N = 10000, K = 20, B = 150, M = 0, d = 8, rho = 3, Seed = 1)
#  8921 1040 39 (3)
Norm1 = Unbiased(N = 10000, K = 20, B = 155, M = 0, d = 8, rho = 3, Seed = 1)
Norm1 = Unbiased(N = 10000, K = 20, B = 160, M = 0, d = 8, rho = 3, Seed = 1)





# Older results including non-CFTP and smaller sample sizes
#******************** d = 1
# Vary M, find B such that prob. of lack of coalescence = 0.5.
# Preliminary runs to find a suitable value for B
Norm1 = Unbiased(N = 1000, K = 20, B = 3, M = 0, d = 1, rho = 2, Seed = 1)
# 475 279 138 68 21 12 7
# 436 292 151 72 29 13 7
Norm1 = Unbiased(N = 1000, K = 20, B = 4, M = 0, d = 1, rho = 2, Seed = 1)
# 619 252 89 33 4 2 1
# 596 260 101 33 5 2 2 1

Norm1 = Unbiased(N = 1000, K = 20, B = 3, M = 0, d = 1, rho = 2.5, Seed = 1)
# 529 266 133 47 15 5 5
# 503 280 138 50 19 5 5
Norm1 = Unbiased(N = 1000, K = 20, B = 3, M = 0, d = 1, rho = 3, Seed = 1)
# 535 259 130 49 16 7 3 1
# 504 269 134 58 21 8 4 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 3, M = 0, d = 1, rho = 3.5, Seed = 1)
# *554 247 123 50 13 11 2
# *527 262 127 52 19 11 2
Norm1 = Unbiased(N = 1000, K = 20, B = 3, M = 0, d = 1, rho = 4, Seed = 1)
# 523 271 126 46 20 10 3 1
# 512 276 125 49 24 10 3 1

Norm1 = Unbiased(N = 1000, K = 20, B = 2, M = 0, d = 1, rho = 3.5, Seed = 1)
# 380 250 154 96 51 36 22 6 4 1
# 350 240 164 104 62 43 23 9 4 1

# Set B = ceiling(3 * log(10) / log(2)) = 10, for probability of coalescence
#   of 0.9.
Norm1 = Unbiased(N = 1000, K = 20, B = 10, M = 0, d = 1, rho = 3.5, Seed = 1)
# *945 51 4
#  943 53 4

# CFTP results only
Norm1 = Unbiased(N = 1000, K = 20, B = 9, M = 0, d = 1, rho = 3.5, Seed = 1)
# 928 69 3
Norm1 = Unbiased(N = 1000, K = 20, B = 8, M = 0, d = 1, rho = 3.5, Seed = 1)
# *891 105 3 1

#******************** d = 2
Norm1 = Unbiased(N = 1000, K = 20, B = 5, M = 0, d = 2, rho = 3.5, Seed = 1)
# 392 287 157 79 46 25 11 1 2
# 363 289 168 89 49 25 13 2 2
Norm1 = Unbiased(N = 1000, K = 20, B = 6, M = 0, d = 2, rho = 3.5, Seed = 1)
# 450 297 151 55 25 12 5 3 1 1
# 390 314 167 71 31 12 6 4 2 1 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 7, M = 0, d = 2, rho = 3.5, Seed = 1)
# *516 294 117 40 19 9 2 2 1
# *468 299 137 54 23 10 4 2 1 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 8, M = 0, d = 2, rho = 3.5, Seed = 1)
# 545 297 105 36 10 5 1 1 CFTP
Norm1 = Unbiased(N = 1000, K = 20, B = 10, M = 0, d = 2, rho = 3.5, Seed = 1)
# 663 256 61 17 1 2
# 635 272 71 18 2 2

Norm1 = Unbiased(N = 1000, K = 20, B = 7, M = 0, d = 2, rho = 2.5, Seed = 1)
# 480 312 144 42 16 5 1
# 444 329 153 50 17 6 1
Norm1 = Unbiased(N = 1000, K = 20, B = 7, M = 0, d = 2, rho = 3, Seed = 1)
# 510 297 120 50 15 4 2 2
# 463 314 138 55 18 7 3 2
Norm1 = Unbiased(N = 1000, K = 20, B = 7, M = 0, d = 2, rho = 3.25, Seed = 1)
# 462 308 135 58 25 5 2 2 1 1 1 CFTP
Norm1 = Unbiased(N = 1000, K = 20, B = 7, M = 0, d = 2, rho = 3.75, Seed = 1)
# *475 294 130 51 28 14 4 2 1 1 CFTP
Norm1 = Unbiased(N = 1000, K = 20, B = 7, M = 0, d = 2, rho = 4, Seed = 1)
# 516 283 121 44 18 9 4 3 1 1
# 468 288 138 55 25 13 7 3 2 1
Norm1 = Unbiased(N = 1000, K = 20, B = 7, M = 0, d = 2, rho = 4.5, Seed = 1)
# 495 268 127 56 20 16 9 3 2 3 1
# 450 268 138 67 34 21 11 3 3 4 1

7 * log(10) / log(2) # 23.25
Norm1 = Unbiased(N = 1000, K = 20, B = 23, M = 0, d = 2, rho = 3.75, Seed = 1)
# 940 57 3 CFTP
Norm1 = Unbiased(N = 1000, K = 20, B = 21, M = 0, d = 2, rho = 3.75, Seed = 1)
# 904 89 7 CFTP
Norm1 = Unbiased(N = 1000, K = 20, B = 20, M = 0, d = 2, rho = 3.75, Seed = 1)
# *899 91 9 1 CFTP

#******************** d = 3, CFTP only
Norm1 = Unbiased(N = 1000, K = 20, B = 12, M = 0, d = 3, rho = 3.75, Seed = 1)
# 459 334 137 49 13 4 2 2
Norm1 = Unbiased(N = 1000, K = 20, B = 13, M = 0, d = 3, rho = 3.75, Seed = 1)
# 451 343 139 46 15 4 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 14, M = 0, d = 3, rho = 3.75, Seed = 1)
# *472 355 124 40 6 3
Norm1 = Unbiased(N = 1000, K = 20, B = 15, M = 0, d = 3, rho = 3.75, Seed = 1)
# 539 325 105 23 6 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 20, M = 0, d = 3, rho = 3.75, Seed = 1)
# 687 257 45 9 2

Norm1 = Unbiased(N = 1000, K = 20, B = 14, M = 0, d = 3, rho = 3, Seed = 1)
# 481 381 101 30 6 1
Norm1 = Unbiased(N = 1000, K = 20, B = 14, M = 0, d = 3, rho = 3.25, Seed = 1)
# *495 371 100 26 7 1
Norm1 = Unbiased(N = 1000, K = 20, B = 14, M = 0, d = 3, rho = 3.5, Seed = 1)
# 493 358 110 29 9 1

Norm1 = Unbiased(N = 1000, K = 20, B = 28, M = 0, d = 3, rho = 3.25, Seed = 1)
# 852 139 8 1
Norm1 = Unbiased(N = 1000, K = 20, B = 29, M = 0, d = 3, rho = 3.25, Seed = 1)
# 873 120 7
Norm1 = Unbiased(N = 1000, K = 20, B = 30, M = 0, d = 3, rho = 3.25, Seed = 1)
# 892 106 2
Norm1 = Unbiased(N = 1000, K = 20, B = 31, M = 0, d = 3, rho = 3.25, Seed = 1)
# 886 108 6
Norm1 = Unbiased(N = 1000, K = 20, B = 32, M = 0, d = 3, rho = 3.25, Seed = 1)
# *899 100 1

#******************** d = 4, CFTP only
Norm1 = Unbiased(N = 1000, K = 20, B = 20, M = 0, d = 4, rho = 3.25, Seed = 1)
# 399 405 138 43 12 3
Norm1 = Unbiased(N = 1000, K = 20, B = 22, M = 0, d = 4, rho = 3.25, Seed = 1)
# 474 394 109 18 3 2
Norm1 = Unbiased(N = 1000, K = 20, B = 23, M = 0, d = 4, rho = 3.25, Seed = 1)
# *506 395 79 17 2 1
Norm1 = Unbiased(N = 1000, K = 20, B = 24, M = 0, d = 4, rho = 3.25, Seed = 1)
# 551 358 76 12 3
Norm1 = Unbiased(N = 1000, K = 20, B = 25, M = 0, d = 4, rho = 3.25, Seed = 1)
# 531 384 73 9 3
Norm1 = Unbiased(N = 1000, K = 20, B = 32, M = 0, d = 4, rho = 3.25, Seed = 1)
# 725 246 23 5 1

Norm1 = Unbiased(N = 1000, K = 20, B = 23, M = 0, d = 4, rho = 3, Seed = 1)
# 493 398 92 16 1
Norm1 = Unbiased(N = 1000, K = 20, B = 23, M = 0, d = 4, rho = 3.5, Seed = 1)
# 498 367 106 23 5 1
Norm1 = Unbiased(N = 1000, K = 20, B = 23, M = 0, d = 4, rho = 3.75, Seed = 1)
# 479 387 101 26 6 1

Norm1 = Unbiased(N = 1000, K = 20, B = 46, M = 0, d = 4, rho = 3.25, Seed = 1)
# 884 111 5
Norm1 = Unbiased(N = 1000, K = 20, B = 50, M = 0, d = 4, rho = 3.25, Seed = 1)
# *902 95 3
Norm1 = Unbiased(N = 1000, K = 20, B = 69, M = 0, d = 4, rho = 3.25, Seed = 1)
# 982 18

#******************** d = 5
Norm1 = Unbiased(N = 1000, K = 20, B = 30, M = 0, d = 5, rho = 3.25, Seed = 1)
# 458 405 111 24 2
Norm1 = Unbiased(N = 1000, K = 20, B = 33, M = 0, d = 5, rho = 3.25, Seed = 1)
# *501 380 90 24 5
Norm1 = Unbiased(N = 1000, K = 20, B = 40, M = 0, d = 5, rho = 3.25, Seed = 1)
# 635 302 54 8 1
Norm1 = Unbiased(N = 1000, K = 20, B = 50, M = 0, d = 5, rho = 3.25, Seed = 1)
# 771 210 18 1

Norm1 = Unbiased(N = 1000, K = 20, B = 33, M = 0, d = 5, rho = 2.5, Seed = 1)
# 476 433 75 12 3 1
Norm1 = Unbiased(N = 1000, K = 20, B = 33, M = 0, d = 5, rho = 2.75, Seed = 1)
# 513 391 80 14 2
Norm1 = Unbiased(N = 1000, K = 20, B = 33, M = 0, d = 5, rho = 3, Seed = 1)
# 507 406 72 13 2
Norm1 = Unbiased(N = 1000, K = 20, B = 33, M = 0, d = 5, rho = 3.5, Seed = 1)
# 479 385 102 26 7 1
Norm1 = Unbiased(N = 1000, K = 20, B = 33, M = 0, d = 5, rho = 3.75, Seed = 1)
# 476 372 114 26 8 3 1

Norm1 = Unbiased(N = 1000, K = 20, B = 50, M = 0, d = 5, rho = 3.25, Seed = 1)
# 771 210 18 1
Norm1 = Unbiased(N = 1000, K = 20, B = 65, M = 0, d = 5, rho = 3.25, Seed = 1)
# 875 120 4 1
Norm1 = Unbiased(N = 1000, K = 20, B = 67, M = 0, d = 5, rho = 3.25, Seed = 1)
# *898 100 2
Norm1 = Unbiased(N = 1000, K = 20, B = 70, M = 0, d = 5, rho = 3.25, Seed = 1)
# 927 70 3
Norm1 = Unbiased(N = 1000, K = 20, B = 100, M = 0, d = 5, rho = 3.25, Seed = 1)
# 984 16


# Non-CFTP results
Norm1 = Unbiased(N = 1000, K = 20, B = 23, M = 0, d = 5, rho = 3.5, Seed = 1)
# 372 366 170 63 16 6 2 3 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 30, M = 0, d = 5, rho = 3.5, Seed = 1)
# 475 395 98 22 9 1
Norm1 = Unbiased(N = 1000, K = 20, B = 31, M = 0, d = 5, rho = 3.5, Seed = 1)
# *492 371 97 29 8 1 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 32, M = 0, d = 5, rho = 3.5, Seed = 1)
# 535 351 88 17 6 3

Norm1 = Unbiased(N = 1000, K = 20, B = 31, M = 0, d = 5, rho = 3, Seed = 1)
# 473 410 90 23 4
Norm1 = Unbiased(N = 1000, K = 20, B = 31, M = 0, d = 5, rho = 4, Seed = 1)
# 443 367 130 40 11 3 2 1 2 1

31 * log(10) / log(2) # 102.98
Norm1 = Unbiased(N = 1000, K = 20, B = 103, M = 0, d = 5, rho = 3.5, Seed = 1)
# *980 20

#******************** d = 10
Norm1 = Unbiased(N = 1000, K = 20, B = 103, M = 0, d = 10, rho = 3.5, Seed = 1)
# 362 403 159 58 13 4 1
Norm1 = Unbiased(N = 1000, K = 20, B = 130, M = 0, d = 10, rho = 3.5, Seed = 1)
# *505 361 96 30 5 2 1

Norm1 = Unbiased(N = 1000, K = 20, B = 130, M = 0, d = 10, rho = 2, Seed = 1)
# 464 404 102 25 5
Norm1 = Unbiased(N = 1000, K = 20, B = 130, M = 0, d = 10, rho = 2.5, Seed = 1)
# 531 364 86 16 3
Norm1 = Unbiased(N = 1000, K = 20, B = 130, M = 0, d = 10, rho = 2.75, Seed = 1)
# *540 333 95 26 4 2
Norm1 = Unbiased(N = 1000, K = 20, B = 130, M = 0, d = 10, rho = 3, Seed = 1)
# 538 348 93 19 2
Norm1 = Unbiased(N = 1000, K = 20, B = 130, M = 0, d = 10, rho = 4, Seed = 1)
# 479 337 124 42 9 5 2 2

130 * log(10) / log(2) # 431.85
Norm1 = Unbiased(N = 1000, K = 20, B = 250, M = 0, d = 10, rho = 2.75, Seed = 1)
# 878 111 11
Norm1 = Unbiased(N = 1000, K = 20, B = 275, M = 0, d = 10, rho = 2.75, Seed = 1)
# *1000 : 906 90 4
Norm1 = Unbiased(N = 1000, K = 20, B = 300, M = 0, d = 10, rho = 2.75, Seed = 1)
# 930 65 5
Norm1 = Unbiased(N = 1000, K = 20, B = 432, M = 0, d = 10, rho = 2.75, Seed = 1)
# 989 11

#******************** d = 15
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 2.5,
 Seed = 1)
# 473 318 128 54 16 8 1 1 1
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 2.75,
 Seed = 1)
# 487 297 143 48 15 5 2 2 1
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 3, Seed = 1)
# 510 295 113 45 18 8 7 4
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 3.25, Seed = 1)
# 496 285 110 59 31 14 2 3
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 3.5, Seed = 1)
# 518 292 108 47 18 7 8 2
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 3.75, Seed = 1)
# *535 287 111 37 21 3 3 2 1
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 4, Seed = 1)
# 519 298 116 42 14 6 2 2 1
Norm1 = Unbiased(N = 1000, K = 20, B = 500, M = 0, d = 15, rho = 4.25, Seed = 1)
# 490 297 132 54 16 5 5 1

Norm1 = Unbiased(N = 1000, K = 20, B = 1200, M = 0, d = 15, rho = 3.75,
 Seed = 1)
# 871 118 10 1
Norm1 = Unbiased(N = 1000, K = 20, B = 1250, M = 0, d = 15, rho = 3.75,
 Seed = 1)
# *891 92 13 1 1 1 1
# 884 98 13 2 1 1 1
