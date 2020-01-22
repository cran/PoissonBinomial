// [[Rcpp::interfaces(r, cpp)]]
//'

#include <boost/math/constants/constants.hpp> // for pi
#include <fftw3.h>
#include <Rcpp.h>
using namespace Rcpp;

// used to make fft computations more readable
#define REAL 0
#define IMAG 1

// PMFs
NumericVector dpb_conv(IntegerVector obs, NumericVector probs);
NumericVector dpb_dc(IntegerVector obs, NumericVector probs);
NumericVector dpb_dftcf(IntegerVector obs, NumericVector probs);
NumericVector dpb_rf(IntegerVector obs, NumericVector probs);
NumericVector dpb_mean(IntegerVector obs, NumericVector probs);
NumericVector dpb_gmba(IntegerVector obs, NumericVector probs, bool anti = false);
NumericVector dpb_pa(IntegerVector obs, NumericVector probs);
NumericVector dpb_na(IntegerVector obs, NumericVector probs, bool refined = true);

// CDFs
NumericVector ppb_conv(IntegerVector obs, NumericVector probs);
NumericVector ppb_dc(IntegerVector obs, NumericVector probs);
NumericVector ppb_dftcf(IntegerVector obs, NumericVector probs);
NumericVector ppb_rf(IntegerVector obs, NumericVector probs);
NumericVector ppb_mean(IntegerVector obs, NumericVector probs);
NumericVector ppb_gmba(IntegerVector obs, NumericVector probs, bool anti = false);
NumericVector ppb_pa(IntegerVector obs, NumericVector probs);
NumericVector ppb_na(IntegerVector obs, NumericVector probs, bool refined = true);

// helper function for normalisation of PMFs (i.e. ensure that sum = 1)
void norm_dpb(NumericVector &pmf){
  // sums of PMF
  double new_sum = sum(pmf), old_sum = 0, older_sum = 0, oldest_sum = 0;
  //Rcout << ((new_sum < 1)?"l ":((new_sum == 1)?"e ":"g "));
  while(new_sum != 1){
    oldest_sum = older_sum;
    older_sum = old_sum;
    old_sum = new_sum;
    pmf = pmf / new_sum;
    new_sum = sum(pmf);
    //Rcout << ((new_sum < 1)?"l ":((new_sum == 1)?"e ":"g "));
    if(new_sum > 1 || new_sum == old_sum || new_sum == older_sum || new_sum == oldest_sum) break;
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
  
  // ensure that sum = 1, if last value = n (the highest oberservable value)
  if(max_q == probs.length())
    results[max_q] = 1;
  
  // "correct" numerically too large results
  results[results > 1] = 1;
  
  // return final results
  return results[obs];
}

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
  // probability masses
  NumericVector pmf = dpb_conv(IntegerVector(Range(0, max(obs))), probs);
  
  // return final results
  return ppb_generic(obs, probs, pmf);
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
  
  // direct convolution is sufficient, if size is small
  if(size <= 750) return dpb_conv(obs, probs);
  
  // number of tree splits
  int num_splits = (int)std::ceil(std::log(size / 750.0) / std::log(2.0));
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
  // probability masses
  NumericVector pmf = dpb_dc(IntegerVector(Range(0, max(obs))), probs);
  
  // return final results
  return ppb_generic(obs, probs, pmf);
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
  
  // compute closed-form expression of Hernandez and Williams
  int mid = sizeIn / 2 + 1;
  for(int i = 1; i <= mid; i++){
    std::complex<double> product = 1.0;
    double realCpow = C_power.real(), imagCpow = C_power.imag();
    C_power.real(C.real() * realCpow - C.imag() * imagCpow);
    C_power.imag(C.imag() * realCpow + C.real() * imagCpow);
    
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
  results[results < 5.55e-17] = 0;
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  //results = results / sum(results);
  norm_dpb(results);
  
  // return final results
  return results[obs];
}

// [[Rcpp::export]]
NumericVector ppb_dftcf(IntegerVector obs, NumericVector probs){
  // probability masses
  NumericVector pmf = dpb_dftcf(IntegerVector(Range(0, max(obs))), probs);
  
  // return final results
  return ppb_generic(obs, probs, pmf);
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
NumericVector dpb_gmba(IntegerVector obs, NumericVector probs, bool anti){
  // cumulative probabilities
  NumericVector cdf = ppb_gmba(IntegerVector(Range(0, max(obs))), probs, anti);
  
  // compute and return results
  return dpb_generic(obs, probs, cdf);
}

// [[Rcpp::export]]
NumericVector ppb_gmba(IntegerVector obs, NumericVector probs, bool anti){
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
  
  // compute and return results
  return results;
}

// Normal Approximations (NA, RNA)
// [[Rcpp::export]]
NumericVector dpb_na(IntegerVector obs, NumericVector probs, bool refined){
  // cumulative probabilities
  NumericVector cdf = ppb_na(IntegerVector(Range(0, max(obs))), probs, refined);
  
  // compute and return results
  return dpb_generic(obs, probs, cdf);
}

// [[Rcpp::export]]
NumericVector ppb_na(IntegerVector obs, NumericVector probs, bool refined){
  // maximum observed value
  int max_q = max(obs);
  // mu
  double mu = sum(probs);
  // p * q
  NumericVector pq = probs * (1 - probs);
  // sigma
  double sigma = std::sqrt(sum(pq));
  // standardised observations with continuity correction
  NumericVector obs_std = (NumericVector(IntegerVector(Range(0, max_q))) + 0.5 - mu)/sigma;
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
  
  // return final results
  return results[obs];
}
