#include /pre/license.stan
functions {
#include /functions/pairwise.stan
}
data {
  // dimensions
  int<lower=1> NPA;             // worths or players or objects or things
  int<lower=1> NCMP;            // unique comparisons
  int<lower=1> N;               // observations
  int<lower=1> numRefresh;      // when change in item/pa1/pa2
  real alphaShape;
  int<lower=1> NTHRESH;         // number of thresholds
  real varCorrection;
  // response data
  int<lower=1, upper=NPA> pa1[numRefresh];
  int<lower=1, upper=NPA> pa2[numRefresh];
  int weight[NCMP];
  int pick[NCMP];
  int refresh[numRefresh];
  int numOutcome[numRefresh];
}
transformed data {
  int rcat[NCMP];
  real alpha = 1.749;

  for (cmp in 1:NCMP) {
    rcat[cmp] = pick[cmp] + NTHRESH + 1;
  }
}
parameters {
  vector[NPA] theta;
  vector<lower=0,upper=1>[NTHRESH] rawThreshold;
  real<lower=0> sigma;
}
transformed parameters {
  real scale = sd(theta) ^ varCorrection;
  vector[NTHRESH] threshold;
  vector[NTHRESH] cumTh;
  real maxSpan = max(theta) - min(theta);
  threshold = maxSpan * rawThreshold;
  cumTh = cumulative_sum(threshold);
}
model {

  sigma ~ lognormal(1, 1);
  theta ~ normal(0, sigma);
  rawThreshold ~ beta(1.1, 1.1);

  {
    int cmpStart = 1;
    for (rx in 1:numRefresh) {
      target += pairwise_logprob(rcat, weight, cmpStart, refresh[rx],
                                 scale, alpha, theta[pa1[rx]], theta[pa2[rx]], cumTh);
      cmpStart += refresh[rx];
    }
  }
}
generated quantities {
  real thetaVar = variance(theta);
}
