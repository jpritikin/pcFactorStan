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
  int<lower=1> NTHRESH;         // number of thresholds
  real scale;
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

  for (cmp in 1:NCMP) {
    rcat[cmp] = pick[cmp] + NTHRESH + 1;
  }
}
parameters {
  vector[NPA] theta;
  vector[NTHRESH] threshold;
  real<lower=0> alpha;
}
transformed parameters {
  vector[NTHRESH] cumTh = cumulative_sum(threshold);
}
model {
  alpha ~ exponential(0.1);
  theta ~ normal(0, 1.0);
  threshold ~ normal(0, 2.0);

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
