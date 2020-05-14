// [[Rcpp::interfaces(r, cpp)]]
//'

#include <boost/math/constants/constants.hpp> // for pi
#include <fftw3.h>
#include <Rcpp.h>
using namespace Rcpp;

// used to make fft computations more readable
#define REAL 0
#define IMAG 1

/******************************************************************/
/**   Helper Functions                                           **/
/******************************************************************/

// helper function for normalisation of PMFs (i.e. ensure that sum = 1)
void norm_dpb(NumericVector &pmf){
  // sums of PMF
  double new_sum = sum(pmf), old_sum = 0, older_sum = 0, oldest_sum = 0;
  //Rcout << ((new_sum < 1)?"l ":((new_sum == 1)?"e ":"g "));
  while(new_sum != 1){
    oldest_sum = older_sum;
    older_sum = old_sum;
    old_sum = new_sum;
    NumericVector old_pmf = pmf;
    pmf = pmf / new_sum;
    new_sum = sum(pmf);
    //Rcout << ((new_sum < 1)?"l ":((new_sum == 1)?"e ":"g "));
    if(new_sum >= 1 || new_sum == old_sum || new_sum == older_sum || new_sum == oldest_sum) break;
    if(new_sum < 1 && new_sum <= old_sum){
      pmf = old_pmf;
      break;
    }
  }
  //Rcout << "\n";
}

// "generic" function for computing some of the PMFs
NumericVector dpb_generic(IntegerVector obs, NumericVector probs, NumericVector cdf){
  // maximum observed value
  int max_q = max(obs);
  // results vector
  NumericVector results(max_q + 1);
  
  // compute masses
  results[0] = cdf[0];
  for(int i = 1; i <= max_q; i++)
    results[i] = cdf[i] - cdf[i - 1];
  
  // return final results
  return results[obs];
}

// "generic" function for computing some of the CDFs
NumericVector ppb_generic(IntegerVector obs, NumericVector probs, NumericVector pmf){
  // maximum observed value
  int max_q = max(obs);
  // results vector
  NumericVector results = NumericVector(max_q + 1);
  
  // compute cumulative probabilities
  results[0] = pmf[0];
  for(int i = 1; i <= max_q; i++)
    results[i] = pmf[i] + results[i - 1];
  
  // "correct" numerically too large results
  results[results > 1] = 1;
  
  // return final results
  return results[obs];
}

IntegerVector order(IntegerVector x){
  IntegerVector uni = unique(x).sort();
  IntegerVector order(x.length());
  int k = 0;
  for(int i = 0; i < uni.length(); i++){
    for(int j = 0; j < x.length(); j++){
      if(uni[i] == x[j]) order[k++] = j;
    }
  }
  return order;
}

/******************************************************************/
/**   Functions for "ordinary" Poisson binomial distribution     **/
/******************************************************************/

// PMFs
NumericVector dpb_conv(IntegerVector obs, NumericVector probs);
NumericVector dpb_dc(IntegerVector obs, NumericVector probs);
NumericVector dpb_dftcf(IntegerVector obs, NumericVector probs);
NumericVector dpb_rf(IntegerVector obs, NumericVector probs);
NumericVector dpb_mean(IntegerVector obs, NumericVector probs);
NumericVector dpb_gmba(IntegerVector obs, NumericVector probs, bool anti);
NumericVector dpb_pa(IntegerVector obs, NumericVector probs);
NumericVector dpb_na(IntegerVector obs, NumericVector probs, bool refined);

// CDFs
NumericVector ppb_conv(IntegerVector obs, NumericVector probs);
NumericVector ppb_dc(IntegerVector obs, NumericVector probs);
NumericVector ppb_dftcf(IntegerVector obs, NumericVector probs);
NumericVector ppb_rf(IntegerVector obs, NumericVector probs);
NumericVector ppb_mean(IntegerVector obs, NumericVector probs);
NumericVector ppb_gmba(IntegerVector obs, NumericVector probs, bool anti);
NumericVector ppb_pa(IntegerVector obs, NumericVector probs);
NumericVector ppb_na(IntegerVector obs, NumericVector probs, bool refined);

// Direct Convolution
// [[Rcpp::export]]
NumericVector dpb_conv(IntegerVector obs, NumericVector probs){
  // number of input probabilities
  double size = probs.length();
  
  NumericVector temp;
  NumericVector results(size + 1);
  results[0] = 1 - probs[0];
  results[1] = probs[0];
  
  for(int i = 1; i < size; i++){
    checkUserInterrupt();
    temp = results[Range(0, i)];
    results[Range(0, i)] = temp * (1 - probs[i]);
    results[Range(1, i + 1)] = results[Range(1, i + 1)] + temp * probs[i];
  }
  // make sure that probability masses sum up to 1
  //results = results / sum(results);
  norm_dpb(results);
  
  // return final results
  return results[obs];
}

// [[Rcpp::export]]
NumericVector ppb_conv(IntegerVector obs, NumericVector probs){
  // maximum observed value
  int max_q = max(obs);
  
  // probability masses
  NumericVector pmf = dpb_conv(IntegerVector(Range(0, max_q)), probs);
  
  // compute CDF
  NumericVector results = ppb_generic(obs, probs, pmf);
  
  // ensure that sum = 1, if last value = n (the highest oberservable value)
  if(max_q == probs.length())
    results[max_q] = 1;
  
  // return final results
  return results;
}

// Divide & Conquer FFT (DC-FFT)
NumericVector fft_probs(NumericVector probsA, NumericVector probsB){
  // sizes of input vectors and the result
  int sizeA = probsA.length();
  int sizeB = probsB.length();
  int sizeResult = sizeA + sizeB - 1;
  
  // results vector
  NumericVector result(sizeResult);
  
  // allocate memory for FFTs of the probs and the convolution result
  fftw_complex *probsA_fft, *probsB_fft, *result_fft;
  
  // 0-padding of probsA vector and perform FFT of it
  NumericVector padded_probsA(sizeResult);
  padded_probsA[Range(0, sizeA - 1)] = probsA;
  probsA_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)* sizeResult);
  fftw_plan planA = fftw_plan_dft_r2c_1d(sizeResult, padded_probsA.begin(), probsA_fft, FFTW_ESTIMATE);
  fftw_execute(planA);
  fftw_destroy_plan(planA);
  
  // 0-padding of probsB vector and perform FFT of it
  NumericVector padded_probsB(sizeResult);
  padded_probsB[Range(0, sizeB - 1)] = probsB;
  probsB_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)* sizeResult);
  fftw_plan planB = fftw_plan_dft_r2c_1d(sizeResult, padded_probsB.begin(), probsB_fft, FFTW_ESTIMATE);
  fftw_execute(planB);
  fftw_destroy_plan(planB);
  
  // convolute by complex multiplication of the transformed input probs
  result_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*sizeResult);
  for(int i = 0; i < sizeResult; i++){
    result_fft[i][REAL] = (probsA_fft[i][REAL]*probsB_fft[i][REAL] - probsA_fft[i][IMAG]*probsB_fft[i][IMAG])/sizeResult;
    result_fft[i][IMAG] = (probsA_fft[i][REAL]*probsB_fft[i][IMAG] + probsA_fft[i][IMAG]*probsB_fft[i][REAL])/sizeResult;
  }
  
  // inverse tranformation of the above multiplications
  fftw_plan planResult = fftw_plan_dft_c2r_1d(sizeResult, result_fft, result.begin(), FFTW_ESTIMATE);
  fftw_execute(planResult);
  fftw_destroy_plan(planResult);
  
  // garbage collection
  fftw_free(probsA_fft);
  fftw_free(probsB_fft);
  fftw_free(result_fft);
  
  // return final results
  return result;
}

// [[Rcpp::export]]
NumericVector dpb_dc(IntegerVector obs, NumericVector probs){
  // number of probabilities of success
  int size = probs.length();
  
  // direct convolution is sufficient, if size is below 600
  if(size <= 600) return dpb_conv(obs, probs);
  // automatically determine number of splits, if size is above 600
  int num_splits = (int)std::ceil(std::log(size / 600) / std::log(2.0));
  // number of groups
  int num_groups = (int)std::pow(2, num_splits);
  
  // loop variables
  int i, start, end;
  
  // compute group sizes with minimum size disparity
  IntegerVector group_sizes(num_groups, size / num_groups);
  int remainder = size % num_groups;
  for(i = 0; i < remainder; i++) group_sizes[i]++;
  
  // compute first and last indices of the groups
  IntegerVector starts(num_groups), ends(num_groups);
  starts[0] = 0;
  ends[0] = group_sizes[0] - 1;
  for(i = 1; i < num_groups; i++){
    starts[i] = starts[i - 1] + group_sizes[i - 1];
    ends[i] = ends[i - 1] + group_sizes[i];
  }
  
  // results vector; direct allocation will increase size of each group by 1
  NumericVector results(size + num_groups);
  
  // compute direct convolutions for each group
  for(int i = 0; i < num_groups; i++){
    checkUserInterrupt();
    // compute new starting and ending indices, because groups grow by 1
    start = starts[i] + i;
    end = ends[i] + i + 1;
    
    // target range
    Range target(start, end);
    
    // direct convolution
    results[target] = dpb_conv(IntegerVector(Range(0, end - start)), probs[Range(starts[i], ends[i])]);
    
    // update starting and ending indices
    starts[i] = start;
    ends[i] = end;
  }
  
  int num_groups_reduced = num_groups / 2;
  while(num_splits > 0){
    for(i = 0; i < num_groups_reduced; i++){
      checkUserInterrupt();
      // compute new starting and ending indices, because group sizes are
      // reduced by 1, due to FFT convolution
      start = starts[2*i] - i;
      end = ends[2*i + 1] - i - 1;
      
      //convolution
      results[Range(start, end)] = fft_probs(results[Range(starts[2*i], ends[2*i])], results[Range(starts[2*i + 1], ends[2*i + 1])]);
      
      // update starting and ending indices
      starts[i] = start;
      ends[i] = end;
    }
    num_groups_reduced /= 2;
    num_splits -= 1;
  }
  
  // select final results
  results = NumericVector(results[Range(0, size)]);
  
  // "correct" numerically false (and thus useless) results
  results[results < 5.55e-17] = 0;
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  //results = results / sum(results);
  norm_dpb(results);
  
  // return final results
  return results[obs];
}

// [[Rcpp::export]]
NumericVector ppb_dc(IntegerVector obs, NumericVector probs){
  // maximum observed value
  int max_q = max(obs);
  
  // probability masses
  NumericVector pmf = dpb_dc(IntegerVector(Range(0, max_q)), probs);
  
  // compute CDF
  NumericVector results = ppb_generic(obs, probs, pmf);
  
  // ensure that sum = 1, if last value = n (the highest oberservable value)
  if(max_q == probs.length())
    results[max_q] = 1;
  
  // return final results
  return results;
}

// Discrete Fourier Transformation of Characteristic Function (DFT-CF)
// [[Rcpp::export]]
NumericVector dpb_dftcf(IntegerVector obs, NumericVector probs){
  // number of probabilities of success
  int sizeIn = probs.length();
  // number of distribution
  int sizeOut = sizeIn + 1;
  
  // "initialize" DFT input vector
  fftw_complex *input_fft;
  input_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut);
  input_fft[0][REAL] = 1.0;
  input_fft[0][IMAG] = 0.0;
  
  // initialize complex numbers for "C" and "C to the power of i"
  std::complex<double> C = exp(std::complex<double>(0.0, 2.0) * boost::math::double_constants::pi / ((double)sizeOut));
  std::complex<double> C_power(1.0, 0.0);
  
  int mid = sizeIn / 2 + 1;
  /*double omega = 2 * boost::math::double_constants::pi / ((double)sizeOut);
  double d, x, c, s;
  std::complex<double> z;
  for(int l = 1; l <= mid; l++){
    d = 0.0;
    x = 0.0;
    for(int j = 0; j < sizeIn; j++){
      z.real(probs[j] * (std::cos(omega * l) - 1) + 1);
      z.imag(probs[j] * std::sin(omega * l));
      d += std::log(std::abs(z));
      x += std::arg(z);
    }
    d = std::exp(d);
    c = d * std::cos(x);
    s = d * std::sin(x);
    input_fft[l][REAL] = c;
    input_fft[l][IMAG] = s;
    input_fft[sizeOut - l][REAL] = c;
    input_fft[sizeOut - l][IMAG] = -s;
  }*/
  
  // compute closed-form expression of Hernandez and Williams
  for(int i = 1; i <= mid; i++){
    checkUserInterrupt();
    std::complex<double> product = 1.0;
    //double realCpow = C_power.real(), imagCpow = C_power.imag();
    //C_power.real(C.real() * realCpow - C.imag() * imagCpow);
    //C_power.imag(C.imag() * realCpow + C.real() * imagCpow);
    C_power *= C;
    
    for(int j = 0; j < sizeIn; j++) product *= 1.0 + (C_power - 1.0) * probs[j];
    
    input_fft[i][REAL] = product.real();
    input_fft[i][IMAG] = product.imag();
    input_fft[sizeOut - i][REAL] = product.real();
    input_fft[sizeOut - i][IMAG] = -product.imag();
  }
  
  // vector of DFT results
  fftw_complex *result_fft;
  result_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut);
  
  // perform DFT
  fftw_plan planDFT;
  planDFT = fftw_plan_dft_1d(sizeOut, input_fft, result_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(planDFT);
  
  // gather results
  NumericVector results(sizeOut);
  for(int i = 0; i < sizeOut; i++) results[i] = result_fft[i][REAL] / sizeOut;
  
  // garbage collection
  fftw_destroy_plan(planDFT);
  fftw_free(input_fft);
  fftw_free(result_fft);
  
  // "correct" numerically false (and thus useless) results
  results[results < 2.22e-16] = 0;
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  //results = results / sum(results);
  norm_dpb(results);
  
  // return final results
  return results[obs];
}

// [[Rcpp::export]]
NumericVector ppb_dftcf(IntegerVector obs, NumericVector probs){
  // maximum observed value
  int max_q = max(obs);
  
  // probability masses
  NumericVector pmf = dpb_dftcf(IntegerVector(Range(0, max_q)), probs);
  
  // compute CDF
  NumericVector results = ppb_generic(obs, probs, pmf);
  
  // ensure that sum = 1, if last value = n (the highest oberservable value)
  if(max_q == probs.length())
    results[max_q] = 1;
  
  // return final results
  return results;
}

// Recursive Formula
// [[Rcpp::export]]
NumericVector dpb_rf(IntegerVector obs, NumericVector probs){
  // number of input probabilities
  double size = probs.length();
  // maximum observed value
  int max_q = max(obs);
  
  NumericMatrix dist(size + 1, 2);
  NumericVector results(max_q + 1);
  int col_new = 0, col_old = 1;
  
  dist(0, col_new) = 1.0;
  dist(1, col_new) = 1 - probs[0];
  for(int j = 2; j <= size; j++) dist(j, col_new) = (1 - probs[j - 1]) * dist(j - 1, col_new);
  results[0] = dist(size, col_new);
  
  for(int i = 1; i <= max_q; i++){
    checkUserInterrupt();
    col_new -= std::pow(-1, i);
    col_old += std::pow(-1, i);
    
    for(int j = 0; j <= i - 1; j++)
      dist(j, col_new) = 0;
    
    for(int j = i; j <= size; j++){
      dist(j, col_new) = (1 - probs[j - 1]) * dist(j - 1, col_new) + probs[j - 1] * dist(j - 1, col_old);
    }
    
    results[i] = dist(size, col_new);
  }
  // make sure that probability masses sum up to 1
  //results = results / sum(results);
  norm_dpb(results);
  
  // return final results
  return results[obs];
}

// [[Rcpp::export]]
NumericVector ppb_rf(IntegerVector obs, NumericVector probs){
  // probability masses
  NumericVector pmf = dpb_rf(IntegerVector(Range(0, max(obs))), probs);
  
  // return final results
  return ppb_generic(obs, probs, pmf);
}

// Arithmetic Mean Binomial Approximation
// [[Rcpp::export]]
NumericVector dpb_mean(IntegerVector obs, NumericVector probs){
  // cumulative probabilities
  NumericVector cdf = ppb_mean(IntegerVector(Range(0, max(obs))), probs);
  
  // compute and return results
  return dpb_generic(obs, probs, cdf);
}

// [[Rcpp::export]]
NumericVector ppb_mean(IntegerVector obs, NumericVector probs){
  // mean of probabilities is the approximate binomial probability
  double bin_prob = mean(probs);
  
  // cumulative probabilities
  NumericVector results = pbinom(obs, (double)probs.length(), bin_prob);
  
  // "correct" numerically too large results
  results[results > 1] = 1;
  
  // compute and return results
  return results;
}

// Geometric Mean Binomial Approximations
// [[Rcpp::export]]
NumericVector dpb_gmba(IntegerVector obs, NumericVector probs, bool anti = false){
  // cumulative probabilities
  NumericVector cdf = ppb_gmba(IntegerVector(Range(0, max(obs))), probs, anti);
  
  // compute and return results
  return dpb_generic(obs, probs, cdf);
}

// [[Rcpp::export]]
NumericVector ppb_gmba(IntegerVector obs, NumericVector probs, bool anti = false){
  // number of probabilities of success
  int size = probs.length();
  
  // logarithms of 'probs' (sums of logarithms are numerically more stable than
  // products of probabilities, especially when the probabilities are small)
  NumericVector logs;
  double bin_prob;
  
  if(anti){
    logs = NumericVector(log(1 - probs));
    bin_prob = 1 - std::exp(mean(logs));
  }else{
    logs = NumericVector(log(probs));
    bin_prob = std::exp(mean(logs));
  }
  
  // cumulative probabilities
  NumericVector results = pbinom(obs, (double)size, bin_prob);
  
  // "correct" numerically too large results
  results[results > 1] = 1;
  
  // compute and return results
  return results;
}

// Poisson Approximation
// [[Rcpp::export]]
NumericVector dpb_pa(IntegerVector obs, NumericVector probs){
  // cumulative probabilities
  NumericVector cdf = ppb_pa(IntegerVector(Range(0, max(obs))), probs);
  
  // compute and return results
  return dpb_generic(obs, probs, cdf);
}

// [[Rcpp::export]]
NumericVector ppb_pa(IntegerVector obs, NumericVector probs){
  // sum of probability is the expectation of the Poisson approximation
  double lambda = sum(probs);
  
  // cumulative probabilities
  NumericVector results = ppois(obs, lambda);
  
  // "correct" numerically too large results
  results[results > 1] = 1;
  
  // make sure that largest observation has probability of 1
  results[obs == probs.length()] = 1;
  
  // compute and return results
  return results;
}

// Normal Approximations (NA, RNA)
// [[Rcpp::export]]
NumericVector dpb_na(IntegerVector obs, NumericVector probs, bool refined = true){
  // cumulative probabilities
  NumericVector cdf = ppb_na(IntegerVector(Range(0, max(obs))), probs, refined);
  
  // compute and return results
  return dpb_generic(obs, probs, cdf);
}

// [[Rcpp::export]]
NumericVector ppb_na(IntegerVector obs, NumericVector probs, bool refined = true){
  // maximum observed value
  int max_q = max(obs);
  // mu
  double mu = sum(probs);
  // p * q
  NumericVector pq = probs * (1 - probs);
  // sigma
  double sigma = std::sqrt(sum(pq));
  // standardised observations with continuity correction
  NumericVector obs_std = (NumericVector(obs) + 0.5 - mu)/sigma;
  // vector to store results
  NumericVector results = pnorm(obs_std);
  // cumulated probabilities
  if(refined){
    // gamma
    double gamma = sum(pq * (1 - 2 * probs));
    // probabilities
    results = results + gamma/(6 * std::pow(sigma, 3.0)) * (1 - pow(obs_std, 2.0)) * dnorm(obs_std);
  }
  // make sure that all probabilities do not exceed 1 and are at least 0
  results[results < 0] = 0;
  results[results > 1] = 1;
  
  // make sure largest possible value has cumulative probability of 1
  if(max_q == probs.length())
    results[obs == max_q] = 1;
  
  // return final results
  return results;
}


/******************************************************************/
/**   Functions for generalized Poisson binomial distribution    **/
/******************************************************************/

// PMFs
NumericVector dgpb_conv(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector dgpb_dc(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector dgpb_dftcf(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector dgpb_na(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool refined);

// CDFs
NumericVector pgpb_conv(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector pgpb_dc(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector pgpb_dftcf(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector pgpb_na(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool refined);


// Generalized Direct Convolution (G-DC)
// [[Rcpp::export]]
NumericVector dgpb_conv(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q){
  // number of probabilities of success
  int sizeIn = probs.length();
  // determine pairwise minimum and maximum
  NumericVector vp = pmin(val_p, val_q);
  NumericVector vq = pmax(val_p, val_q);
  // compute differences
  NumericVector d = vq - vp;
  // re-order 'probs' into 'pp' so that each probability is for the smaller value
  NumericVector pp(sizeIn);
  for(int i = 0; i < sizeIn; i++){
    if(val_p[i] > vp[i]) pp[i] = 1 - probs[i]; else pp[i] = probs[i];
  }
  
  // (maximum) output size
  int sizeOut = sum(d) + 1;
  
  // results vectors
  NumericVector results(sizeOut);
  
  // initialize results (first convolution step)
  results[0] = 1.0;
  int idx_end = 0;
  
  // perform convolution
  NumericVector temp;
  for(int i = 0; i < sizeIn; i++){
    checkUserInterrupt();
    temp = results[Range(0, idx_end)];
    for(int j = 0; j <= idx_end; j++){
      if(temp[j]) results[j] = temp[j] * pp[i];
    }
    for(int j = 0; j <= idx_end; j++){
      if(temp[j]) results[j + d[i]] += temp[j] * (1 - pp[i]);
    }
    idx_end += d[i];
  }
  
  // "correct" numerically false (and thus useless) results
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  norm_dpb(results);
  
  return results[obs - sum(vp)];
}

// [[Rcpp::export]]
NumericVector pgpb_conv(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q){
  // probability masses
  NumericVector pmf = dgpb_conv(IntegerVector(Range(sum(pmin(val_p, val_q)), max(obs))), probs, val_p, val_q);
  
  // return final results
  return ppb_generic(obs - sum(pmin(val_p, val_q)), probs, pmf);
}

// Generalized Divide & Conquer FFT Tree Convolution (G-DC-FFT)
// [[Rcpp::export]]
NumericVector dgpb_dc(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q){
  // number of probabilities of success
  int size = probs.length();
  // determine pairwise minimum and maximum
  NumericVector vp = pmin(val_p, val_q);
  NumericVector vq = pmax(val_p, val_q);
  // compute differences
  NumericVector d = vq - vp;
  // compute cumulative sums
  NumericVector cs = cumsum(d);
  // output size
  double sizeOut = cs[size - 1];
  
  // number of tree splits
  int num_splits = std::max<int>(0, (int)round(0.2 * std::log(sizeOut)/std::log(2.0) + 0.85 * std::log(size)/std::log(2.0) - 8.75));
  // direct convolution is sufficient, if no splits are necessary
  if(num_splits <= 0) return dgpb_conv(obs, probs, val_p, val_q);
  // number of groups
  int num_groups = (int)std::pow(2, num_splits);
  while(num_splits >= 0 && (num_groups > size || sizeOut/num_groups < max(d))){
    num_groups /= 2;
    num_splits -= 1;
  }
  // loop variables
  int i = 0, start, end;
  
  // compute group sizes with minimum size disparity
  IntegerVector group_sizes(num_groups);
  IntegerVector probs_start(num_groups);
  IntegerVector probs_end(num_groups);
  i = 0;
  for(int j = 1; j <= num_groups; j++){
    int sum = round(sizeOut * j/(double)num_groups);
    while(i < size && cs[i] < sum) i++;
    if(cs[i] - sum < sum - cs[i - 1]){
      //if(j > 1 && cs[i] == group_sizes[j - 2]) i++; 
      group_sizes[j - 1] = cs[i];
      probs_end[j - 1] = i;
      if(j < num_groups) probs_start[j] = i + 1;
    }else{
      //if(j > 1 && cs[i - 1] == group_sizes[j - 2]) i++; 
      group_sizes[j - 1] = cs[i - 1];
      probs_end[j - 1] = i - 1;
      if(j < num_groups) probs_start[j] = i;
    }
    i++;
  }
  
  for(int j = num_groups - 1; j > 0; j--) group_sizes[j] -= group_sizes[j - 1];
  
  // compute first and last indices of the groups
  IntegerVector group_starts(num_groups), group_ends(num_groups);
  group_starts[0] = 0;
  group_ends[0] = group_sizes[0] - 1;
  for(i = 1; i < num_groups; i++){
    group_starts[i] = group_starts[i - 1] + group_sizes[i - 1];
    group_ends[i] = group_ends[i - 1] + group_sizes[i];
  }
  
  // results vector; direct convolution will increase size of each group by 1
  NumericVector results(sizeOut + num_groups);
  
  // compute direct convolutions for each group
  for(int i = 0; i < num_groups; i++){
    checkUserInterrupt();
    // compute new starting and ending indices, because groups grow by 1
    start = group_starts[i] + i;//////////
    end = group_ends[i] + i + 1;
    
    // target range
    Range target(start, end);
    vp = val_p[Range(probs_start[i], probs_end[i])];
    vq = val_q[Range(probs_start[i], probs_end[i])];
    
    // direct convolution
    results[target] = dgpb_conv(IntegerVector(Range(sum(pmin(vp, vq)), sum(pmax(vp, vq)))), probs[Range(probs_start[i], probs_end[i])], vp, vq);
    
    // update starting and ending indices
    group_starts[i] = start;
    group_ends[i] = end;
  }
  
  int num_groups_reduced = num_groups / 2;
  while(num_splits > 0){
    for(i = 0; i < num_groups_reduced; i++){
      checkUserInterrupt();
      // compute new starting and ending indices, because group sizes are
      // reduced by 1, due to FFT convolution
      start = group_starts[2*i] - i;
      end = group_ends[2*i + 1] - i - 1;
      
      //convolution
      results[Range(start, end)] = fft_probs(results[Range(group_starts[2*i], group_ends[2*i])], results[Range(group_starts[2*i + 1], group_ends[2*i + 1])]);
      
      // update starting and ending indices
      group_starts[i] = start;
      group_ends[i] = end;
    }
    num_groups_reduced /= 2;
    num_splits -= 1;
  }
  
  // select final results
  results = results[Range(0, sizeOut)];
  
  // "correct" numerically false (and thus useless) results
  results[results < 5.55e-17] = 0;
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  //results = results / sum(results);
  norm_dpb(results);
  
  // return final results
  return results[obs - sum(pmin(val_p, val_q))];
}

// [[Rcpp::export]]
NumericVector pgpb_dc(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q){
  // maximum observed value
  int max_q = max(obs);
  
  // probability masses
  NumericVector pmf = dgpb_dc(IntegerVector(Range(sum(pmin(val_p, val_q)), max_q)), probs, val_p, val_q);
  
  // compute CDF
  NumericVector results = ppb_generic(obs - sum(pmin(val_p, val_q)), probs, pmf);
  
  // ensure that sum = 1, if last value equals the highest oberservable value
  if(max_q == sum(pmax(val_p, val_q)))
    results[obs == max_q] = 1;
  
  // return final results
  return results;
}

// Generalized Discrete Fourier Transformation of Characteristic Function (G-DFT-CF)
// [[Rcpp::export]]
NumericVector dgpb_dftcf(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q){
  // number of probabilities of success
  int sizeIn = probs.length();
  // determine minimum and maximum sum and their difference
  int sum_min = sum(pmin(val_p, val_q));
  int sum_max = sum(pmax(val_p, val_q));
  // output size
  int sizeOut = sum_max - sum_min + 1;
  
  // "initialize" DFT input vector
  fftw_complex *input_fft;
  input_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut);
  input_fft[0][REAL] = 1.0;
  input_fft[0][IMAG] = 0.0;
  
  // initialize complex numbers for "C" and "C to the power of i"
  std::complex<double> C = 1.0;
  if(sum_min) C = exp(std::complex<double>(0.0, -sum_min * 2.0) * boost::math::double_constants::pi / ((double)sizeOut));
  std::complex<double> C_power(1.0, 0.0);
  std::vector< std::complex<double> > Cp(sizeIn);
  std::vector< std::complex<double> > Cq(sizeIn);
  std::vector< std::complex<double> > Cp_power(sizeIn, 1.0);
  std::vector< std::complex<double> > Cq_power(sizeIn, 1.0);
  for(int i = 0; i < sizeIn; i++){
    if(val_p[i])
      Cp[i] = exp(std::complex<double>(0.0, val_p[i] * 2.0) * boost::math::double_constants::pi / ((double)sizeOut));
    else Cp[i] = 1.0;
    if(val_q[i])
      Cq[i] = exp(std::complex<double>(0.0, val_q[i] * 2.0) * boost::math::double_constants::pi / ((double)sizeOut));
    else Cq[i] = 1.0;
    
    //Cp_power[i] = 1.0;
    //Cq_power[i] = 1.0;
  }
  
  // compute closed-form expression of Hernandez and Williams
  for(int l = 1; l <= sizeOut / 2; l++){
    checkUserInterrupt();
    std::complex<double> product = 1.0;
    for(int k = 0; k < sizeIn; k++){
      if(val_p[k]) Cp_power[k] *= Cp[k];
      if(val_q[k]) Cq_power[k] *= Cq[k];
      product *= (probs[k] * Cp_power[k] + (1.0 - probs[k]) * Cq_power[k]);
    }
    
    if(sum_min){
      C_power *= C;
      product *= C_power;
    }
    
    input_fft[l][REAL] = product.real();
    input_fft[l][IMAG] = product.imag();
    input_fft[sizeOut - l][REAL] = product.real();
    input_fft[sizeOut - l][IMAG] = -product.imag();
  }
  
  // vector of DFT results
  fftw_complex *result_fft;
  result_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut);
  
  // perform DFT
  fftw_plan planDFT;
  planDFT = fftw_plan_dft_1d(sizeOut, input_fft, result_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(planDFT);
  
  // gather results
  NumericVector results(sizeOut);
  for(int i = 0; i < sizeOut; i++) results[i] = result_fft[i][REAL] / sizeOut;
  
  // garbage collection
  fftw_destroy_plan(planDFT);
  fftw_free(input_fft);
  fftw_free(result_fft);
  
  // "correct" numerically false (and thus useless) results
  results[results < 2.22e-16] = 0;
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  //results = results / sum(results);
  norm_dpb(results);
  
  // return final results
  return results[obs - sum_min];
}

// [[Rcpp::export]]
NumericVector pgpb_dftcf(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q){
  // maximum observed value
  int max_q = max(obs);
  
  // probability masses
  NumericVector pmf = dgpb_dftcf(IntegerVector(Range(sum(pmin(val_p, val_q)), max_q)), probs, val_p, val_q);
  
  // compute CDF
  NumericVector results = ppb_generic(obs - sum(pmin(val_p, val_q)), probs, pmf);
  
  // ensure that sum = 1, if last value equals the highest oberservable value
  if(max_q == sum(pmax(val_p, val_q)))
    results[obs == max_q] = 1;
  
  // return final results
  return results;
}

// Generalized Normal Approximations (G-NA, G-RNA)
// [[Rcpp::export]]
NumericVector dgpb_na(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool refined = true){
  // cumulative probabilities
  NumericVector cdf = pgpb_na(IntegerVector(Range(sum(pmin(val_p, val_q)), max(obs))), probs, val_p, val_q, refined);
  
  // compute and return results
  return dpb_generic(obs - sum(pmin(val_p, val_q)), probs, cdf);
}

// [[Rcpp::export]]
NumericVector pgpb_na(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool refined = true){
  // maximum observed value
  int max_q = max(obs);
  // mu
  double mu = sum(probs * val_p + (1 - probs) * val_q);
  // p * q
  NumericVector pq = probs * (1 - probs);
  // sigma
  double sigma = std::sqrt(sum(pq * pow(val_p - val_q, 2)));
  // standardized observations with continuity correction
  NumericVector obs_std = (NumericVector(obs) + 0.5 - mu)/sigma;
  // vector to store results
  NumericVector results = pnorm(obs_std);
  // cumulated probabilities
  double gamma = 0;
  if(refined && sigma){
    // gamma
    gamma = sum(pq * (1 - 2 * probs) * pow(val_p - val_q, 3))/std::pow(sigma, 3);
    // probabilities
    results = results + gamma * (1 - pow(obs_std, 2)) * dnorm(obs_std) / 6;
  }
  // make sure that all probabilities do not exceed 1 and are at least 0
  results[results < 0] = 0;
  results[results > 1] = 1;
  
  // make sure largest possible value has cumulative probability of 1
  if(max_q == sum(pmax(val_p, val_q)))
    results[obs == max_q] = 1;
  
  // return final results
  return results;
}
