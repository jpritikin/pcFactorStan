#include /pre/license.stan
functions {
#include /functions/cmp_prob.stan
}
data {
  // dimensions
  int<lower=1> NPA;             // number of players or objects or things
  int<lower=1> NCMP;            // number of unique comparisons
  int<lower=1> N;               // number of observations
  int<lower=1> NITEMS;
  int<lower=1> NTHRESH[NITEMS];         // number of thresholds
  int<lower=1> TOFFSET[NITEMS];
  real scale;
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
  row_vector[NITEMS] rawUnique;      // do not interpret, see unique
  matrix[NPA,NITEMS] rawUniqueTheta; // do not interpret, see uniqueTheta
  vector[NPA] rawFactor;      // do not interpret, see factor
  vector[NITEMS] rawLoadings; // do not interpret, see factorLoadings
}
transformed parameters {
  vector[totalThresholds] cumTh;
  matrix[NPA,NITEMS] theta;
  for (pa in 1:NPA) {
    theta[pa,] = (rawFactor[pa] * rawLoadings)' +
      rawUniqueTheta[pa,] .* rawUnique;
  }
  for (ix in 1:NITEMS) {
    int from = TOFFSET[ix];
    int to = TOFFSET[ix] + NTHRESH[ix] - 1;
    cumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
}
model {
  vector[max(NTHRESH)*2 + 1] prob;
  int probSize;

  threshold ~ normal(0, 2.0);
  rawFactor ~ normal(0, 1.0);
  rawLoadings ~ normal(0, 1.0);
  for (pa in 1:NPA) {
    rawUniqueTheta[pa,] ~ normal(0, 1.0);
  }
  rawUnique ~ normal(1.0, 1.0);
  for (cmp in 1:NCMP) {
    if (refresh[cmp]) {
      int ix = item[cmp];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      probSize = (2*NTHRESH[ix]+1);
      prob[:probSize] = cmp_probs(scale,
               theta[pa1[cmp], ix],
               theta[pa2[cmp], ix], cumTh[from:to]);
    }
    if (weight[cmp] == 1) {
      target += categorical_lpmf(rcat[cmp] | prob[:probSize]);
    } else {
      target += weight[cmp] * categorical_lpmf(rcat[cmp] | prob[:probSize]);
    }
  }
}
generated quantities {
  vector[NITEMS] sigma;
  vector[NITEMS] factorLoadings = rawLoadings;
  vector[NITEMS] factorProp;
  vector[NPA] factor = rawFactor;
  row_vector[NITEMS] unique = rawUnique;
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
    // https://www.tandfonline.com/doi/full/10.1080/00031305.2018.1549100
    real resid = variance(rawUniqueTheta[,fx] * rawUnique[fx]);
    real pred = variance(rawFactor * rawLoadings[fx]);
    factorProp[fx] = pred / (pred + resid);
    if (factorLoadings[fx] < 0) factorProp[fx] = -factorProp[fx];
  }
}
