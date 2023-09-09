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
  array[numRefresh] int<lower=1, upper=NPA> pa1;
  array[numRefresh] int<lower=1, upper=NPA> pa2;
  array[NCMP] int weight;
  array[NCMP] int pick;
  array[numRefresh] int refresh;
  array[numRefresh] int numOutcome;
}
transformed data {
  array[NCMP] int rcat;

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
  vector[NPA] theta = rawTheta;
  theta -= mean(theta);
  theta /= sd(theta);
}
