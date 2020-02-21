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
  int<lower=1> NFACTORS;
  real factorScalePrior[NFACTORS];
  int<lower=1> NPSI;  // = NFACTORS * (NFACTORS-1) / 2;
  real psiScalePrior[NPSI];
  int<lower=1> NPATHS;
  int factorItemPath[2,NPATHS];  // 1 is factor index, 2 is item index
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
  vector[NPATHS] pathScalePrior;
  for (cmp in 1:NCMP) {
    rcat[cmp] = pick[cmp] + NTHRESH[item[cmp]] + 1;
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
  corr_matrix[NFACTORS] Psi;
  matrix[NPA,NFACTORS] rawFactor;      // do not interpret, see factor
  vector[NPATHS] rawLoadings; // do not interpret, see factorLoadings
  matrix[NPA,NITEMS] rawUniqueTheta; // do not interpret, see uniqueTheta
  vector[NITEMS] rawUnique;      // do not interpret, see unique
}
transformed parameters {
  vector[totalThresholds] cumTh;
  cholesky_factor_corr[NFACTORS] CholPsi = cholesky_decompose(Psi);
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
  vector[max(NTHRESH)*2 + 1] prob;
  int probSize;
  int px=1;

  threshold ~ normal(0, 2.0);
  for (cx in 1:(NFACTORS-1)) {
    for (rx in (cx+1):NFACTORS) {
      target += normal_lpdf(logit(0.5 + Psi[rx,cx]/2.0) | 0, psiScalePrior[px]);
      px += 1;
    }
  }
  for (xx in 1:NPA) {
    rawFactor[xx,] ~ multi_normal_cholesky_lpdf(rep_vector(0, NFACTORS), CholPsi);
  }
  rawLoadings ~ normal(0, 5.0);
  rawUnique ~ normal(0, 5.0);
  for (ix in 1:NITEMS) {
    rawUniqueTheta[,ix] ~ std_normal();
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
