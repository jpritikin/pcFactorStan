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
  vector[NPA] rawTheta;
  vector<lower=0,upper=1>[NTHRESH] rawThreshold;
  real<lower=0> alpha;
}
transformed parameters {
  vector[NTHRESH] threshold;
  vector[NTHRESH] rawCumTh;
  real maxSpan = max(rawTheta) - min(rawTheta);
  threshold = maxSpan * rawThreshold;
  rawCumTh = cumulative_sum(threshold);
}
model {

  alpha ~ normal(1.749, alphaScalePrior) T[0,];
  rawTheta ~ std_normal();
  rawThreshold ~ beta(1.1, 2);

  {
    int cmpStart = 1;
    for (rx in 1:numRefresh) {
      target += pairwise_logprob(rcat, weight, cmpStart, refresh[rx],
                                 scale, alpha, rawTheta[pa1[rx]], rawTheta[pa2[rx]], rawCumTh);
      cmpStart += refresh[rx];
    }
  }
}
generated quantities {
  real log_lik[N];

  vector[NPA] theta = rawTheta;
  theta -= mean(theta);
  theta /= sd(theta);
{
    int cmpStart = 1;
    int cur = 1;

    for (rx in 1:numRefresh) {
      int nout = numOutcome[rx];
      int last = cur + nout - 1;
      log_lik[cur:last] =
        pairwise_loo(rcat, weight, nout, cmpStart, refresh[rx],
                     scale, alpha, theta[pa1[rx]], theta[pa2[rx]], rawCumTh);
      cmpStart += refresh[rx];
      cur += numOutcome[rx];
    }
  }
}
