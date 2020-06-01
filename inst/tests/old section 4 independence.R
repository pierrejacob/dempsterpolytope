### This scripts reproduces the examples of Section 7 of the paper on independence testing
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())


## For $K = 4$ and arranged as a 2 by 2 contingency table with proportions
## $\boldsymbol{\theta} = (\theta_{00}, \theta_{01}, \theta_{10}, \theta_{11})$,
## one is interested in testing the hypothesis of independence: $H_0:
## \theta_{00} \theta_{11} = \theta_{01} \theta_{10}$. Here are a few
## alternative methods to do so.

## 1. On the numerical example (10, 7, 22, 11). 

## We first illustrate with the observation ${\bf x} = (10, 7, 22, 11)$, where $n = \| {\bf x}\| = 50$.
counts <- c(10, 7, 22, 11)
K <- 4
counts_table = matrix(counts, byrow = T, nrow = 2)

### i) Pearson's Chi-squared test.

## Pearson's Chi-squared test statistic is X^2 = \sum_{i,j} \frac{(x_{ij} -
## e_{ij})^2}{e_{ij}}, where $e_{ij} = n\hat{p}_{i+}\hat{p}_{+j}$ is the expected
## number of counts in cell $ij$ under the null. The Pearson test statistic is
## asymptotically distributed as $\chi^2_1$ under the null. It is exactly
## equivalent to the score test. For this example, the test statistic is $X^2 =
## 0.3$ and its associated p-value is $0.584$. There is no clear evidence to
## reject the independence hypothesis.

e = matrix(rowSums(counts_table), nrow = 2) %*% matrix(colSums(counts_table), nrow = 1)/sum(counts_table)
(pearson = sum((counts_table - e)^2/e))
(pchisq(pearson, df = 1, lower.tail = F))

### ii) Likelihood ratio test. 

## The likelihood ratio test is asympotitically equivalent to the Pearson's
## Chi-squared test, but may differ in finite samples. The likelihood ratio test
## statistic, $$G^2 = -2\log (L_0/L_1) = 2\sum_{i,j} x_{ij} \log(x_{ij}/e_{ij}),
## $$ where $e_{ij}$ is as defined previously. Under the null, $G^2 \sim
## \chi_1^2$ asymptotically. For this example, the LR test statistic $G^2 =
## 0.297$, and its associated p-value is $0.586$. There is no clear evidence to
## reject the independence hypothesis.

(lr = 2*sum(counts_table*log(counts_table/e)))
(pchisq(lr, df = 1, lower.tail = F))

### iii) Bayesian analysis of association.

## The independence hypothesis is challenging for testing from the Bayesian
## perspective due to its reduced dimensions. One may resort to an analysis of
## Bayes Factors. An alternative way to think about it is, a strong evidence
## towards either positive or negative association is also strong evidence
## against independence. Below we compute the posterior probability of a
## positive vs. negative association in the two by two table. Suppose the prior
## distribution on $\boldsymbol{\theta} \sim Dir(\boldsymbol{\alpha})$, where
## $\boldsymbol{\alpha} = (1, 1, 1, 1)$ or something of the modeler's choice.
## The posterior distribution $$\boldsymbol{\theta} \mid {\bf x} \sim  Dir({\bf
## x} + \boldsymbol{\alpha}),$$ and the posterior probability of positive
## association is $P(H_+ \mid {\bf x})  = \int_{H_+} p(\boldsymbol{\theta} \mid
## {\bf x})\partial \boldsymbol{\theta},$ where $H_+ : \theta_{00} \theta_{11} >
## \theta_{01} \theta_{10}$. The posterior probability of negative association,
## $H_- : \theta_{00} \theta_{11} < \theta_{01} \theta_{10}$, can be computed
## analogously. For the observed dataset, the estimated posterior probabilities
## are respectively $$\hat{P}(H_+ \mid {\bf x}) = 0.286, \quad \hat{P}(H_- \mid
## {\bf x}) = 0.714.$$ This quantity can be easily estimated via Monte Carlo.
## There is no clear evidence to support either association hypotheses $H_+$ or
## $H_-$.

posterior_association <- function(data, prior = 1, nsim = 1e5){
  r_ = mapply(function(l){rgamma(nsim, l)}, data+prior)
  positive = mean(r_[,1]*r_[,4] > r_[,2]*r_[,3])
  return(data.frame(positive = positive, negative = 1-positive, 
             se = sqrt(positive*(1-positive)/nsim)))
}
posterior_association(data = counts_table)
  
## We compare this result with the DS Gibbs sampler output, which reports the
## $(p,q,r)$ for $H_+$ and $H_-$ respectively. The estimates from SMC output
## gives $$\hat{\texttt p}(H_+) = 0.148, \; \hat{\texttt q}(H_+) = 0.644, \;
## \hat{\texttt r}(H_+) = 0.209,$$ and on the reverse side (since $H_+ = \neg
## H_-$ with probability one), $$\hat{\texttt p}(H_-) = 0.644, \; \hat{\texttt
## q}(H_-) =  0.148, \; \hat{\texttt r}(H_-) = 0.209.$$ This result is congruent
## with the Bayesian analysis results, in the sense that for both hypotheses,
## the estimated *lower probabilities*, ${\texttt p}(H_+)$ and ${\texttt
## p}(H_-)$, are lower than the Bayes posterior probabilities, and the estimated
## *upper probabilities*, ${\texttt p}(H_+) + {\texttt r}(H_+)$ and ${\texttt
## p}(H_-)+{\texttt r}(H_+)$, are higher than their Bayesian counterparts. Note
## that this is not always guaranteed, due to Bayesian prior specifications and
## Monte Carlo errors. Nevertheless, the ${\texttt r}$ component adds an
## expression of ''don't know'' to the inference.


## test function
h <- function(etas) unlist(compare_with_independence(etas))
## that test function constructs the polytope of theta s.t. log(theta_j) - log(theta_d) <= log(etas[d,j])  for all d,j in [K]
## and returns four booleans:
## 1) whether that polytope intersects with the 'negative association polytope' of thetas 
## satisfying log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  <= 0
## 2) whether it is contained in it
## 3) whether that polytope intersects with the 'positive association polytope' of thetas 
## satisfying log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  >= 0
## 4) whether it is contained in it 
## first determine burn-in 
omega <- 0.9
lag <- 100
NREP <- 1e3
meetingtimes <- foreach(irep = 1:NREP, .combine = c) %dorng% {
  meeting_times(counts, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
}
## TV upper bounds
niterations <- 100
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes, lag, t))
g <- ggplot(data = data.frame(iteration = 1:niterations, ubounds = ubounds), aes(x = iteration, y = ubounds)) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
g 
## this suggests that a burn-in of 50 iterations should be enough, so we choose 100 to be safe
nchains <- 250
niterations <- 500
burnin <- 100
parallelestimators <- foreach(irep = 1:nchains, .combine = rbind) %dorng% {
  init <- rexp(K)
  init <- init/sum(init)
  samples_gibbs <- gibbs_sampler(niterations, counts, theta_0 = init)
  etas <- samples_gibbs$etas[(burnin+1):niterations,,]
  estimators <- t(apply(etas, 1, h))
  colMeans(estimators)[c(2,4)] ## extracts lower probabilities of negative and positive associations
}
pqr_estimators <- colMeans(parallelestimators)
pqr_estimators

p.negative <- pqr_estimators[1]
p.positive <- pqr_estimators[2]
r = 1 - p.negative - p.positive
DS_association = data.frame(positive = c(p.positive, p.negative, r),
                            negative = c(p.negative, p.positive, r))
rownames(DS_association) = c('p', 'q', 'r')
DS_association
## standard deviation of Monte Carlo estimates
apply(parallelestimators, 2, sd) / sqrt(nrow(parallelestimators))

## 2. An example of London underground incidents.

## This data comes from *Rosenbaum (2002, p. 191)* and was reanalyzied by *Ding
## and Miratrix (2019)*, who supplied the following description: ''In the London
## underground, some train stations have a drainage pit below the tracks. When an
## incident happens (i.e., a passenger falls, jumps, or is pushed from the
## station platform), such a pit is a place to escape contact with the wheels of
## the train. Researchers are interested in the mortality in stations with and
## without such a pit. In stations without a pit, only 5 lived out of 21 recorded
## incidents. For incidents in stations with a pit, 18 out of 32 lived. ''

## The observed data can be summarized by ${\bf x} = (16, 5, 14, 18)$, where the
## row variable is *no pit/0* versus *with pit/1*, and the column variable is
## *death/0* versus *survival/1*. The analysis supplied by *Ding and Miratrix
## (2019)* suggests that the existence of a pit significantly increases the
## chance of incident survival, reporting an improvement estimate $\hat{\tau} =
## 0.324$ with an associated confidence interval $(0.106, 0.543)$. They also
## reported a Neyman's confidence interval of $(0.072, 0.577)$, also excluding 0.

## We analyze the table to see whether the existence of pit is associated with
## the chance of survival. Pearson's Chi-squared test yields strong evidence
## against the null hypothesis of independence, with test statistic $X^2 = 5.43$
## and a p-value of $0.02$.

counts <- c(16, 5, 14, 18)
counts_table = matrix(counts, byrow = T, nrow = 2)

e = matrix(rowSums(counts_table), nrow = 2) %*% matrix(colSums(counts_table), nrow = 1)/sum(counts_table)
(pearson = sum((counts_table - e)^2/e))
(pchisq(pearson, df = 1, lower.tail = F))

## Similarly, the LR test statistic $G^2 = 5.63$ and p-value $0.017$.

(lr = 2*sum(counts_table*log(counts_table/e)))
(pchisq(lr, df = 1, lower.tail = F))

## The Bayesian analysis of association also suggests strong positive correlation of the lack of pit and death, with estimated posterior probabilities
## $$\hat{P}(H_+ \mid {\bf x}) = 0.99, \quad \hat{P}(H_- \mid {\bf x}) = 0.01.$$
  
posterior_association(data = counts_table)

## Finally, DS analyis
## first determine burn-in 
omega <- 0.9
lag <- 100
NREP <- 1e3
meetingtimes <- foreach(irep = 1:NREP, .combine = c) %dorng% {
  meeting_times(counts, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
}
## TV upper bounds
niterations <- 100
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes, lag, t))
g <- ggplot(data = data.frame(iteration = 1:niterations, ubounds = ubounds), aes(x = iteration, y = ubounds)) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
g 
## this suggests that a burn-in of 75 iterations should be enough, so we choose 150 to be safe
nchains <- 250
niterations <- 500
burnin <- 150
parallelestimators <- foreach(irep = 1:nchains, .combine = rbind) %dorng% {
  init <- rexp(K)
  init <- init/sum(init)
  samples_gibbs <- gibbs_sampler(niterations, counts, theta_0 = init)
  etas <- samples_gibbs$etas[(burnin+1):niterations,,]
  estimators <- t(apply(etas, 1, h))
  colMeans(estimators)[c(2,4)] ## extracts lower probabilities of negative and positive associations
}
pqr_estimators <- colMeans(parallelestimators)
pqr_estimators

p.negative <- pqr_estimators[1]
p.positive <- pqr_estimators[2]
r = 1 - p.negative - p.positive
DS_association = data.frame(positive = c(p.positive, p.negative, r),
                            negative = c(p.negative, p.positive, r))
rownames(DS_association) = c('p', 'q', 'r')
DS_association
## standard deviation of Monte Carlo estimates
apply(parallelestimators, 2, sd) / sqrt(nrow(parallelestimators))






