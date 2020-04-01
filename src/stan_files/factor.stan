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
  int<lower=1> NITEMS;
  int<lower=1> NTHRESH[NITEMS];         // number of thresholds
  int<lower=1> TOFFSET[NITEMS];
  vector[NITEMS] scale;
  real alpha[NITEMS];
  int<lower=1> NFACTORS;
  real factorScalePrior[NFACTORS];
  int<lower=1> NPSI;  // = NFACTORS * (NFACTORS-1) / 2;
  real psiScalePrior[NPSI];
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
  vector<lower=0,upper=1>[totalThresholds] rawThreshold;
  corr_matrix[NFACTORS] Psi;
  matrix[NPA,NFACTORS] rawFactor;      // do not interpret, see factor
  vector<lower=0,upper=1>[NPATHS] rawLoadings; // do not interpret, see factorLoadings
  matrix[NPA,NITEMS] rawUniqueTheta; // do not interpret, see uniqueTheta
  vector<lower=0,upper=1>[NITEMS] rawUnique;      // do not interpret, see unique
}
transformed parameters {
  vector[totalThresholds] threshold;
  vector[totalThresholds] cumTh;
  cholesky_factor_corr[NFACTORS] CholPsi = cholesky_decompose(Psi);
  matrix[NPA,NITEMS] theta;
  vector[NPATHS] rawPathProp;  // always positive
  real rawPerComponentVar[NITEMS,1+NFACTORS];
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
    cumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
}
model {
  rawThreshold ~ beta(1.1, 1.1);
  {
    int px=1;
    for (cx in 1:(NFACTORS-1)) {
      for (rx in (cx+1):NFACTORS) {
        target += normal_lpdf(logit(0.5 + Psi[rx,cx]/2.0) | 0, psiScalePrior[px]);
        px += 1;
      }
    }
  }
  for (xx in 1:NPA) {
    rawFactor[xx,] ~ multi_normal_cholesky_lpdf(rep_vector(0, NFACTORS), CholPsi);
  }
  rawLoadings ~ beta(3.0, 3.0);
  rawUnique ~ beta(3.0, 3.0);
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
  matrix[NPA,NFACTORS] factor = rawFactor;
  matrix[NITEMS,NITEMS] residualItemCor;

  {
    matrix[NPA,NITEMS] residual;
    for (ix in 1:NITEMS) {
      residual[,ix] = rawUniqueTheta[,ix] * (2*rawUnique[ix]-1);
      residual[,ix] -= mean(residual[,ix]);
    }
    residualItemCor = crossprod(residual);
    residualItemCor = quad_form_diag(residualItemCor, 1.0 ./ sqrt(diagonal(residualItemCor)));
  }

  {
    vector[NPATHS] pathLoadings = (2*rawLoadings-1);
    int rawSeenFactor[NFACTORS];
    int rawNegateFactor[NFACTORS];
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
    factor[,fx] /= sd(factor[,fx]);
  }
}
