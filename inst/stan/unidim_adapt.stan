#include /pre/license.stan
functions {
#include /functions/pairwise.stan
}
data {
  real alphaScalePrior;
  // dimensions
  int<lower=1> NPA;             // worths or players or objects or things
  int<lower=1> NCMP;            // unique comparisons
  int<lower=1> N;               // observations
  int<lower=1> numRefresh;      // when change in item/pa1/pa2
  int<lower=1> NTHRESH;         // number of thresholds
  real varCorrection;
  // response data
  array[numRefresh] int<lower=1, upper=NPA> pa1;
  array[numRefresh] int<lower=1, upper=NPA> pa2;
  array[NCMP] int weight;
  array[NCMP] int pick;
  array[numRefresh] int refresh;
  array[numRefresh] int numOutcome;
}
transformed data {
  array[NCMP] int rcat;
  real alpha = 1.749;

  for (cmp in 1:NCMP) {
    rcat[cmp] = pick[cmp] + NTHRESH + 1;
  }
}
parameters {
  vector[NPA] rawTheta;
  vector<lower=0,upper=1>[NTHRESH] rawThreshold;
  real<lower=0> sigma;
}
transformed parameters {
  vector[NPA] theta = sqrt(sigma) * rawTheta;
  real scale = sd(theta) ^ varCorrection;
  vector[NTHRESH] threshold;
  vector[NTHRESH] rawCumTh;
  real maxSpan = max(theta) - min(theta);
  threshold = maxSpan * rawThreshold;
  rawCumTh = cumulative_sum(threshold);
}
model {

  sigma ~ inv_gamma(1, 1);
  rawTheta ~ std_normal();
  rawThreshold ~ beta(1.1, 2);

  {
    int cmpStart = 1;
    for (rx in 1:numRefresh) {
      target += pairwise_logprob(rcat, weight, cmpStart, refresh[rx],
                                 scale, alpha, theta[pa1[rx]], theta[pa2[rx]], rawCumTh);
      cmpStart += refresh[rx];
    }
  }
}
generated quantities {
  real thetaVar = variance(theta);
}
