#include /pre/license.stan
functions {
#include /functions/cmp_prob2.stan
}
data {
  // dimensions
  int<lower=1> NPA;             // number of players or objects or things
  int<lower=1> NCMP;            // number of unique comparisons
  int<lower=1> N;               // number of observations
  int<lower=1> NTHRESH;         // number of thresholds
  real scale;
  // response data
  int<lower=1, upper=NPA> pa1[NCMP];        // PA1 for observation N
  int<lower=1, upper=NPA> pa2[NCMP];        // PA2 for observation N
  int weight[NCMP];
  int pick[NCMP];
  int refresh[NCMP];
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
  vector[NTHRESH*2 + 1] prob;
  alpha ~ exponential(0.1);
  theta ~ normal(0, 1.0);
  threshold ~ normal(0, 2.0);
  for (cmp in 1:NCMP) {
    if (refresh[cmp]) {
      prob = cmp_probs(scale, alpha, theta[pa1[cmp]], theta[pa2[cmp]], cumTh);
    }
    if (weight[cmp] == 1) {
      target += categorical_lpmf(rcat[cmp] | prob);
    } else {
      target += weight[cmp] * categorical_lpmf(rcat[cmp] | prob);
    }
  }
}
generated quantities {
  real thetaVar = variance(theta);
  vector[NTHRESH*2 + 1] prob;
  vector[N] log_lik;
  int cur = 1;

  for (cmp in 1:NCMP) {
    real ll;
    if (refresh[cmp]) {
      prob = cmp_probs(scale, alpha, theta[pa1[cmp]], theta[pa2[cmp]], cumTh);
    }
    ll = categorical_lpmf(rcat[cmp] | prob);
    for (wx in 1:weight[cmp]) {
      log_lik[cur] = ll;
      cur = cur + 1;
    }
  }
}
