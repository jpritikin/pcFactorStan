vector cmp_probs(real scale, real alpha, real pa1, real pa2, vector thr, array[] int want) {
  int nth = num_elements(thr);
  int nth2 = nth*2;
  array[1+nth2] real pr;
  vector[1 + nth2] out;
  real paDiff = alpha * scale * (pa1 - pa2);
  vector[nth] thrAlpha = thr * alpha;
  pr[1+nth2] = 1;
  for (tx in 1:nth2) {
    if (want[tx] || want[tx+1]) {
      real at;
      if (tx <= nth) {
        at = -thrAlpha[nth+1-tx];
      } else {
        at = thrAlpha[tx-nth];
      }
      pr[tx] = 1.0/(1.0+exp(-(paDiff + at)));
    } else {
      pr[tx] = 0; // not needed
    }
  }
  out[1] = pr[1];
  for (tx in 2:num_elements(out)) {
    out[tx] = pr[tx] - pr[tx-1];
  }
  for (tx in 1:num_elements(want)) {
    if (want[tx] && (out[tx] <= 0 || out[tx] >= 1)) {
      reject("out[",tx,"]=",out[tx],"; scale=",
             scale, " alpha=", alpha, " pa1=", pa1, " pa2=", pa2, " th=", thr);
    }
  }
  return out;
}

real pairwise_logprob(array[] int rcat, array[] int weight, int cmpStart, int len,
                      real scale, real alpha, real pa1, real pa2, vector cumTh)
{
  real lp = 0;
  int nth = num_elements(cumTh);
  vector[1+nth*2] prob;
  array[1+nth*2] int want;
  for (ox in 1:num_elements(want)) want[ox] = 0;
  for (ox in cmpStart:(cmpStart + len - 1)) {
    want[ rcat[ox] ] = 1;
  }
  prob = cmp_probs(scale, alpha, pa1, pa2, cumTh, want);
  for (ox in cmpStart:(cmpStart + len - 1)) {
    real lp1 = log(prob[rcat[ox]]);
    if (weight[ox] == 1) {
      lp += lp1;
    } else {
      lp += weight[ox] * lp1;
    }
  }
  return lp;
}

array[] real pairwise_loo(array[] int rcat, array[] int weight, int numOutcome, int cmpStart, int len,
                  real scale, real alpha, real pa1, real pa2, vector cumTh)
{
  array[numOutcome] real lp;
  int cur = 1;
  int nth = num_elements(cumTh);
  vector[1+nth*2] prob;
  array[1+nth*2] int want;
  for (ox in 1:num_elements(want)) want[ox] = 0;
  for (ox in cmpStart:(cmpStart + len - 1)) {
    want[ rcat[ox] ] = 1;
  }
  prob = cmp_probs(scale, alpha, pa1, pa2, cumTh, want);
  for (ox in cmpStart:(cmpStart + len - 1)) {
    real lp1 = log(prob[rcat[ox]]);
    for (wx in 1:weight[ox]) {
      lp[cur] = lp1;
      cur += 1;
    }
  }
  return lp;
}
