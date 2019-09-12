#include /pre/license.stan
functions {
#include /functions/cmp_prob2.stan
}
data {
  // dimensions
  int<lower=1> NPA;             // number of players or objects or things
  int<lower=1> NCMP;            // number of unique comparisons
  int<lower=1> N;               // number of observations
  int<lower=1> NITEMS;
  int<lower=1> NTHRESH[NITEMS];         // number of thresholds
  int<lower=1> TOFFSET[NITEMS];
  vector[NITEMS] scale;
  vector[NITEMS] alpha;
  // response data
  int<lower=1, upper=NPA> pa1[NCMP];        // PA1 for observation N
  int<lower=1, upper=NPA> pa2[NCMP];        // PA2 for observation N
  int weight[NCMP];
  int pick[NCMP];
  int refresh[NCMP];
  int item[NCMP];
}
transformed data {
  int totalThresholds = sum(NTHRESH);
  int rcat[NCMP];
  for (cmp in 1:NCMP) {
    rcat[cmp] = pick[cmp] + NTHRESH[item[cmp]] + 1;
  }
}
parameters {
  vector[totalThresholds] threshold;
  matrix[NPA,NITEMS] rawUniqueTheta; // do not interpret, see uniqueTheta
  vector[NPA] rawFactor;      // do not interpret, see factor
  vector[NITEMS] rawLoadings; // do not interpret, see factorLoadings
  vector[NITEMS] rawUnique;      // do not interpret, see unique
}
transformed parameters {
  vector[NITEMS] rawFactorProp;  // always positive
  vector[totalThresholds] cumTh;
  matrix[NPA,NITEMS] theta;
  for (pa in 1:NPA) {
    // theta must be normal(0,1) distributed
    theta[pa,] = ((rawFactor[pa] * rawLoadings)' +
                  rawUniqueTheta[pa,] .* rawUnique');
  }
  for (ix in 1:NITEMS) {
    int from = TOFFSET[ix];
    int to = TOFFSET[ix] + NTHRESH[ix] - 1;
    cumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
  for (fx in 1:NITEMS) {
    real resid = variance(rawUniqueTheta[,fx] * rawUnique[fx]);
    real pred = variance(rawFactor * rawLoadings[fx]);
    rawFactorProp[fx] = 0.5 + pred / (2 * (pred + resid));
  }
}
model {
  vector[NITEMS] rawLogitProp;
  vector[max(NTHRESH)*2 + 1] prob;
  int probSize;

  threshold ~ normal(0, 2.0);
  target += NITEMS * normal_lpdf(rawFactor | 0, 1.0);
  for (ix in 1:NITEMS) {
    rawLoadings[ix] ~ normal(0, 2.0);
    rawUnique[ix] ~ normal(0, 2.0);
  }
  for (pa in 1:NPA) {
    rawUniqueTheta[pa,] ~ normal(0, 1.0);
  }
  for (cmp in 1:NCMP) {
    if (refresh[cmp]) {
      int ix = item[cmp];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      probSize = (2*NTHRESH[ix]+1);
      prob[:probSize] = cmp_probs(scale[ix], alpha[ix],
               theta[pa1[cmp], ix],
               theta[pa2[cmp], ix], cumTh[from:to]);
    }
    if (weight[cmp] == 1) {
      target += categorical_lpmf(rcat[cmp] | prob[:probSize]);
    } else {
      target += weight[cmp] * categorical_lpmf(rcat[cmp] | prob[:probSize]);
    }
  }
  for (fx in 1:NITEMS) {
    rawLogitProp[fx] = logit(rawFactorProp[fx]);
  }
  target += normal_lpdf(rawLogitProp | 0, 0.5);
}
generated quantities {
  vector[NITEMS] factorProp = rawFactorProp;
  vector[NITEMS] sigma;
  vector[NITEMS] factorLoadings = rawLoadings;
  vector[NPA] factor = rawFactor;
  row_vector[NITEMS] unique = rawUnique';
  matrix[NPA,NITEMS] uniqueTheta = rawUniqueTheta;

  for (fx in 1:NITEMS) {
    sigma[fx] = sd(theta[,fx]);
    if (unique[fx] < 0) {
      unique[fx] = -unique[fx];
      uniqueTheta[,fx] = -uniqueTheta[,fx];
    }
  }
  if (factorLoadings[1] < 0) {
    factorLoadings = -factorLoadings;
    factor = -factor;
  }
  for (fx in 1:NITEMS) {
    if (factorLoadings[fx] < 0) factorProp[fx] = -factorProp[fx];
  }
}
