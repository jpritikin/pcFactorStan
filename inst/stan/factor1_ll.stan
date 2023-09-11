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
  int<lower=1> NITEMS;
  array[NITEMS] int<lower=1> NTHRESH;         // number of thresholds
  array[NITEMS] int<lower=1> TOFFSET;
  vector[NITEMS] scale;
  real propShape;
  int<lower=1> NFACTORS;
  array[NFACTORS] real factorScalePrior;
  int<lower=1> NPATHS;
  array[2,NPATHS] int factorItemPath;  // 1 is factor index, 2 is item index
  // response data
  array[numRefresh] int<lower=1, upper=NPA> pa1;
  array[numRefresh] int<lower=1, upper=NPA> pa2;
  array[NCMP] int weight;
  array[NCMP] int pick;
  array[numRefresh] int refresh;
  array[numRefresh] int numOutcome;
  array[numRefresh] int item;
}
transformed data {
  int totalThresholds = sum(NTHRESH);
  array[NCMP] int rcat;
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
  array[NITEMS] real<lower=0> alpha;
  vector<lower=0,upper=1>[totalThresholds] rawThreshold;
  matrix[NPA,NFACTORS] rawFactor;      // do not interpret, see factor
  vector<lower=0,upper=1>[NPATHS] rawLoadings; // do not interpret, see factorLoadings
  matrix[NPA,NITEMS] rawUniqueTheta; // do not interpret, see uniqueTheta
  vector<lower=0,upper=1>[NITEMS] rawUnique;      // do not interpret, see unique
}
transformed parameters {
  vector[totalThresholds] threshold;
  vector[totalThresholds] rawCumTh;
  matrix[NPA,NITEMS] theta;
  vector[NPATHS] rawPathProp;  // always positive
  array[NITEMS,1+NFACTORS] real rawPerComponentVar;
  for (ix in 1:NITEMS) {
    theta[,ix] = rawUniqueTheta[,ix] * (2*rawUnique[ix]-1);
    rawPerComponentVar[ix, 1] = variance(theta[,ix]);
  }
  for (fx in 1:NFACTORS) {
    for (ix in 1:NITEMS) rawPerComponentVar[ix,1+fx] = 0;
  }
  for (px in 1:NPATHS) {
    int fx = factorItemPath[1,px];
    int ix = factorItemPath[2,px];
    vector[NPA] theta1 = (2*rawLoadings[px]-1) * rawFactor[,fx];
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
  for (ix in 1:NITEMS) {
    real maxSpan = max(theta[,ix]) - min(theta[,ix]);
    int from = TOFFSET[ix];
    int to = TOFFSET[ix] + NTHRESH[ix] - 1;
    threshold[from:to] = maxSpan * rawThreshold[from:to];
    rawCumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
}
model {
  for (ix in 1:NITEMS) alpha[ix] ~ normal(1.749, alphaScalePrior) T[0,];
  rawThreshold ~ beta(1.1, 2);
  rawFactor[,1] ~ std_normal();
  rawLoadings ~ beta(propShape, propShape);
  rawUnique ~ beta(propShape, propShape);
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
                                 theta[pa2[rx], ix], rawCumTh[from:to]);
      cmpStart += refresh[rx];
    }
  }
  // 1.0 excessive, 1.5 not enough
  target += normal_lpdf(logit(0.5 + rawPathProp/2.0) | 0, pathScalePrior);
}
generated quantities {
  array[N] real log_lik;

  vector[NPATHS] pathProp = rawPathProp;
  matrix[NPA,NFACTORS] factor = rawFactor;
  matrix[NITEMS,NITEMS] residualItemCor;

  {
    matrix[NPA,NITEMS] residual;
    for (ix in 1:NITEMS) {
      residual[,ix] = rawUniqueTheta[,ix] * (2*rawUnique[ix]-1);
      residual[,ix] -= mean(residual[,ix]);
      residual[,ix] /= sd(residual[,ix]);
    }
    residualItemCor = crossprod(residual);
    residualItemCor = quad_form_diag(residualItemCor, 1.0 ./ sqrt(diagonal(residualItemCor)));
  }

  {
    vector[NPATHS] pathLoadings = (2*rawLoadings-1);
    array[NFACTORS] int rawSeenFactor;
    array[NFACTORS] int rawNegateFactor;
    for (fx in 1:NFACTORS) rawSeenFactor[fx] = 0;
    for (px in 1:NPATHS) {
      int fx = factorItemPath[1,px];
      int ix = factorItemPath[2,px];
      if (rawSeenFactor[fx] == 0) {
        rawSeenFactor[fx] = 1;
        rawNegateFactor[fx] = rawLoadings[px] < 0.5;
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
  for (fx in 1:NFACTORS) {
    factor[,fx] -= mean(factor[,fx]);
    factor[,fx] /= sd(factor[,fx]);
  }

  {
    int cmpStart = 1;
    int cur = 1;
    for (rx in 1:numRefresh) {
      int ix = item[rx];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      int nout = numOutcome[rx];
      int last = cur + nout - 1;
      log_lik[cur:last] =
        pairwise_loo(rcat, weight, nout, cmpStart, refresh[rx],
                     scale[ix], alpha[ix],
                     theta[pa1[rx],ix], theta[pa2[rx],ix], rawCumTh[from:to]);
      cmpStart += refresh[rx];
      cur += numOutcome[rx];
    }
  }
}
