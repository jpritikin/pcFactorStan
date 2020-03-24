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
  int<lower=1> NITEMS;
  int<lower=1> NTHRESH[NITEMS];         // number of thresholds
  int<lower=1> TOFFSET[NITEMS];
  vector[NITEMS] scale;
  vector[NITEMS] alpha;
  int<lower=1> NFACTORS;
  real factorScalePrior[NFACTORS];
  int<lower=1> NPATHS;
  int factorItemPath[2,NPATHS];  // 1 is factor index, 2 is item index
  // response data
  int<lower=1, upper=NPA> pa1[numRefresh];
  int<lower=1, upper=NPA> pa2[numRefresh];
  int weight[NCMP];
  int pick[NCMP];
  int refresh[numRefresh];
  int numOutcome[numRefresh];
  int item[numRefresh];
}
transformed data {
  int totalThresholds = sum(NTHRESH);
  int rcat[NCMP];
  vector[NPATHS] pathScalePrior;
  {
    int cmpStart = 0;
    for (rx in 1:numRefresh) {
      int ix = item[rx];
      for (cmp in 1:refresh[rx]) {
        rcat[cmpStart + cmp] = pick[cmpStart + cmp] + NTHRESH[ix] + 1;
      }
      cmpStart += refresh[rx];
    }
  }
  for (px in 1:NPATHS) {
    int fx = factorItemPath[1,px];
    int ix = factorItemPath[2,px];
    if (fx < 1 || fx > NFACTORS) {
      reject("factorItemPath[1,","px","] names factor ", fx, " (NFACTORS=",NFACTORS,")");
    }
    if (ix < 1 || ix > NITEMS) {
      reject("factorItemPath[2,","px","] names item ", ix, " (NITEMS=",NITEMS,")");
    }
    pathScalePrior[px] = factorScalePrior[fx];
  }
}
parameters {
  vector[totalThresholds] threshold;
  matrix[NPA,NFACTORS] rawFactor;      // do not interpret, see factor
  vector[NPATHS] rawLoadings; // do not interpret, see factorLoadings
  matrix[NPA,NITEMS] rawUniqueTheta; // do not interpret, see uniqueTheta
  vector[NITEMS] rawUnique;      // do not interpret, see unique
}
transformed parameters {
  vector[totalThresholds] cumTh;
  matrix[NPA,NITEMS] theta;
  vector[NPATHS] rawPathProp;  // always positive
  real rawPerComponentVar[NITEMS,1+NFACTORS];
  for (ix in 1:NITEMS) {
    int from = TOFFSET[ix];
    int to = TOFFSET[ix] + NTHRESH[ix] - 1;
    cumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
  for (ix in 1:NITEMS) {
    theta[,ix] = rawUniqueTheta[,ix] * rawUnique[ix];
    rawPerComponentVar[ix, 1] = variance(theta[,ix]);
  }
  for (fx in 1:NFACTORS) {
    for (ix in 1:NITEMS) rawPerComponentVar[ix,1+fx] = 0;
  }
  for (px in 1:NPATHS) {
    int fx = factorItemPath[1,px];
    int ix = factorItemPath[2,px];
    vector[NPA] theta1 = rawLoadings[px] * rawFactor[,fx];
    rawPerComponentVar[ix,1+fx] = variance(theta1);
    theta[,ix] += theta1;
  }
  for (px in 1:NPATHS) {
    int fx = factorItemPath[1,px];
    int ix = factorItemPath[2,px];
    real resid = 0;
    real pred;
    for (cx in 1:(1+NFACTORS)) {
      if (cx == fx+1) {
        pred = rawPerComponentVar[ix,cx];
      } else {
        resid += rawPerComponentVar[ix,cx];
      }
    }
    rawPathProp[px] = pred / (pred + resid);
  }
}
model {
  threshold ~ normal(0, 2.0);
  rawFactor[,1] ~ std_normal();
  rawLoadings ~ normal(0, 5.0);
  rawUnique ~ normal(0, 5.0);
  for (ix in 1:NITEMS) {
    rawUniqueTheta[,ix] ~ std_normal();
  }
  {
    int cmpStart = 1;
    for (rx in 1:numRefresh) {
      int ix = item[rx];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      target += pairwise_logprob(rcat, weight, cmpStart, refresh[rx],
                                 scale[ix], alpha[ix], theta[pa1[rx], ix],
                                 theta[pa2[rx], ix], cumTh[from:to]);
      cmpStart += refresh[rx];
    }
  }
  target += normal_lpdf(logit(0.5 + rawPathProp/2.0) | 0, pathScalePrior);
}
generated quantities {
  vector[NPATHS] pathProp = rawPathProp;
  vector[NITEMS] sigma;
  vector[NPATHS] pathLoadings = rawLoadings;
  matrix[NPA,NFACTORS] factor = rawFactor;
  matrix[NPA,NITEMS] residual;
  matrix[NITEMS,NITEMS] residualItemCor;
  int rawSeenFactor[NFACTORS];
  int rawNegateFactor[NFACTORS];

  for (ix in 1:NITEMS) {
    residual[,ix] = rawUniqueTheta[,ix] * rawUnique[ix];
    residual[,ix] -= mean(residual[,ix]);
  }
  residualItemCor = crossprod(residual);
  residualItemCor = quad_form_diag(residualItemCor, 1.0 ./ sqrt(diagonal(residualItemCor)));

  for (fx in 1:NITEMS) {
    sigma[fx] = sd(theta[,fx]);
  }
  for (fx in 1:NFACTORS) rawSeenFactor[fx] = 0;
  for (px in 1:NPATHS) {
    int fx = factorItemPath[1,px];
    int ix = factorItemPath[2,px];
    if (rawSeenFactor[fx] == 0) {
      rawSeenFactor[fx] = 1;
      rawNegateFactor[fx] = rawLoadings[px] < 0;
    }
    if (rawNegateFactor[fx]) {
      pathLoadings[px] = -pathLoadings[px];
    }
  }
  for (fx in 1:NFACTORS) {
    if (!rawNegateFactor[fx]) continue;
    factor[,fx] = -factor[,fx];
  }
  for (fx in 1:NPATHS) {
    if (pathLoadings[fx] < 0) pathProp[fx] = -pathProp[fx];
  }
}
