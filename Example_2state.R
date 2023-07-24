# -*- Text -*-
##############################################################################
# R code for demonstration of perfect simulation from unbiased simulation with
#   a simple example of a discrete Markov chain with two states

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

# Apart from reading in the function "Unbiased", this script is intended to be
#   run interactively, not by calling the "source" function.

######################################## Function code
Unbiased = function(n, p, theta, k, Seed = 1) {
 # Parameters:
 # - n: Number of independent samples to generate
 # - p: Probability of moving from state 2 to state 1 in one iteration (0 < p
 #   < 1): small values are more challenging.
 # - theta: Multiplier for the probability of moving from state 1 to state 2,
 #   which determines the ergodic state (0 < theta <= 1); introduced to stop
 #   the ergodic distribution being 50-50, seeing as we want to start from
 #   50-50.  Small values make the ergodic distribution more asymmetric.
 # - k: Number of iterations before the unbiased simulation methodology takes
 #   effect
 # - Seed: Random number seed, provides reproducibility of results.
 # For a two-state process, the Markov transition matrix is cbind(c(1 - theta
 #   * p, theta * p), c(p, 1 - p)).
 # For this example, the distribution of starting points is (0.5, 0.5).  The
 #   ergodic distribution is (1 / (1 + theta), theta / (1 + theta))
 # The two chains X and Y are coupled by using the same uniform(0, 1) random
 #   number R to move between states.  The chain takes state 1 if R <= the
 #   probability of moving to or staying in state 1, given the current state.
 set.seed(Seed)
 nState = 2 # Make a variable, in case we ever want to change it.
 Pinit = 0.5 # Probability of starting from state 1
 Trans = rbind(c(1 - theta * p, theta * p), c(p, 1 - p)) # Transition matrix

 Count = array(0, dim = c(n, 2 * nState)) # Array to hold the answers: columns
 #   1:nState hold the unbiased sampling results, while columns (nState +
 #   1):(2 * nState) hold the results that we would get if we stopped after k
 #   iterations and didn't continue with the unbiased sampling algorithm.
 X = 1 + (runif(n) > Pinit) # Initial state of the vector of chains X, which
 #   have one additional iteration compared to Y
 Y = 1 + (runif(n) > Pinit) # Independent setting of the initial state of the
 #   vector of chains Y, which will have one less iteration than X
 Ind = 1:n # Vector of indices of independent samples

 # Apply the first iteration to X only.
 R = runif(n)
 X = 1 + (R > Trans[X, 1])

 # Apply iterations 2 to k to both X and Y, without testing for coalescence.
 if (k > 1) {
  for (iit in 2:k) {
   R = runif(n)
   # Update both X and Y, using the same random numbers R.
   X = 1 + (R > Trans[X, 1])
   Y = 1 + (R > Trans[Y, 1])
  } # iit
 } # k > 1

 # Put X into the results.  Use matrix subscripts.
 Count[cbind(Ind, X)] = 1 # Update this one below.
 Count[cbind(Ind, X + nState)] = 1 # Don't update.

 # The two chains have coalesced when X == Y.  Keep running the non-coalesced
 #   ones until they coalesce.  In the results, count +1 for X and -1 for Y.
 # For computational efficiency, we'll program it to use only the
 #   non-coalesced samples, on the assumption that k has been chosen large
 #   enough that most of the samples have already coalesced.
 lNonCoal = X != Y
 while (any(lNonCoal)) {

  # We need to generate new random numbers only for the non-coalesced samples,
  #   but doing that would affect the random numbers that we get, which could
  #   result in, for example, a *higher* number of iterations required to
  #   achieve coalescence for a *lower* value of k, which is undesirable when
  #   we wish to compare results for different values of k.
  # It is safer to regenerate all the random numbers, although it takes some
  #   extra computational effort.  I've found uniform random number generation
  #   very quick in the past, so I don't think that will be a big problem.
  R = runif(n)
  # Update X and Y for non-coalesced samples only.
  X[lNonCoal] = 1 + (R[lNonCoal] > Trans[X[lNonCoal], 1])
  Y[lNonCoal] = 1 + (R[lNonCoal] > Trans[Y[lNonCoal], 1])

  # Update lNonCoal.  We can do this before storing the results because, if X
  #   and Y have coalesced (i.e., X == Y), they will cancel each other out
  #   when we add a count of 1 for X and subtract a count of 1 for Y.
  lNonCoal = X != Y

  # Store the results.  As described above, we need to do this only for
  #   non-coalesced values.
  # We count +1 for X and -1 for Y.
  Count[cbind(Ind[lNonCoal], X[lNonCoal])] =
   Count[cbind(Ind[lNonCoal], X[lNonCoal])] + 1
  Count[cbind(Ind[lNonCoal], Y[lNonCoal])] =
   Count[cbind(Ind[lNonCoal], Y[lNonCoal])] - 1
 } # while (any(lNonCoal))
 Count
} # Unbiased

######################################## Code to run the function

# With theta = 1/9, the results should converge to ergodic probabilities of
#   (0.9, 0.1).  The standard deviation of independent, binomially distributed
#   results should be sqrt(0.9 * 0.1) = 0.3.

# n = 1e6 is a large enough sample for our purposes, and runs fairly quickly.
#   S.d. of sample mean in a binomial distribution is then 0.0003.

#################### 1 million sample strings for each value of k over a wide
#   range
krange = 5:120
CountList = list() # Save all results
Summ = array(0, dim = c(length(krange), 6)) # Summary matrix for different
#   values of k
dimnames(Summ) = list(krange, c("Mean3", "Mean1", "Prop", "Total", "Min",
 "Sd"))

for (k in krange) {
#for (k in c(5, 10, 20, 50, 100, 110)) {
 cat(k, "\n")
 Count = Unbiased(n = 1e6, p = 0.1, theta = 1 / 9, k = k)
 CountList[[as.character(k)]] = Count
 Summ[as.character(k), "Mean3"] = mean(Count[, 3]) # Mean of unbiased samples
 Summ[as.character(k), "Mean1"] = mean(Count[, 1]) # Mean before unbiased
 #   sampling methodology is executed
 Summ[as.character(k), "Prop"] = mean(Count[, 1] < 0 | Count[, 2] < 0) #
 #   Proportion of strings that contain holes
 Summ[as.character(k), "Total"] = -mean(pmin(Count[, 1], Count[, 2], 0)) #
 #   Total number of holes, as a proportion of n
 Summ[as.character(k), "Min"] = -min(Count) # Max number of holes
 Summ[as.character(k), "Sd"] = sqrt(var(Count[, 2])) # Sample s.d. (might be
 #   more accurate if we used column 2 instead of column 1, as values are
 #   smaller; answer is the same)
}

save(CountList, Summ, file = "Example_2state.RData") # Read back in to create
#   figure and table for paper.

#################### 10 million strings with k = 5
Count = Unbiased(n = 1e7, p = 0.1, theta = 1 / 9, k = 5)
apply(Count, 2, sum)
# 9046424  953576
apply(Count < 0, 2, sum)
# 1261960 1512221: 27.7% of counts are negative
max(abs(Count)[Count < 0]) # Up to 120
