// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>


// [[Rcpp::export]]
SEXP CorScoreCalcArmaAltLoop_t(arma::mat gem, Rcpp::NumericVector next_seed, Rcpp::List next_seed_sq, Rcpp::List next_seed_pair, 
                               int samp_num, arma::rowvec unary_sum, arma::rowvec sq_sum, arma::mat pair_sum){
  
  Rcpp::NumericVector cor_scores(next_seed.size());
  
  // Increase number of samples by 1
  samp_num += 1;
  
  for (int seed = 0; seed < next_seed.size(); seed++) {
    arma::rowvec curr_seed_sq = Rcpp::as<arma::rowvec>(next_seed_sq[next_seed[seed]-1]);
    arma::mat curr_seed_pair = Rcpp::as<arma::mat>(next_seed_pair[next_seed[seed]-1]);
    
    // Sum of samples update
    arma::rowvec new_unary_sum = unary_sum + gem.row(next_seed[seed]-1);
    // Sum of samples*2 update
    arma::rowvec new_sq_sum = sq_sum + curr_seed_sq;
    // Matrix of sample - sample paris
    arma::mat new_pair_sum = pair_sum + curr_seed_pair;
    
    arma::mat cov_ij = samp_num * new_pair_sum - new_unary_sum.t() * new_unary_sum;
    arma::mat cov_ii = samp_num * new_sq_sum - new_unary_sum % new_unary_sum;
    
    arma::mat a = abs(cov_ij/sqrt(cov_ii.t() * cov_ii));
    
    cor_scores[seed] = accu(a)/a.n_elem;
  }
  
  return Rcpp::wrap(cor_scores);
}

// [[Rcpp::export]]
SEXP CorScoreCalcEigenAltLoop_t(Eigen::MatrixXd gem, Rcpp::NumericVector next_seed, Rcpp::List next_seed_sq, Rcpp::List next_seed_pair,
                                int samp_num, Eigen::RowVectorXd unary_sum, Eigen::RowVectorXd sq_sum, Eigen::MatrixXd pair_sum){
  
  Rcpp::NumericVector cor_scores(next_seed.size());
  
  samp_num += 1;
  
  for (int seed = 0; seed < next_seed.size(); seed++) {
    Eigen::RowVectorXd curr_seed_sq = Rcpp::as<Eigen::RowVectorXd>(next_seed_sq[next_seed[seed]-1]);
    Eigen::MatrixXd curr_seed_pair = Rcpp::as<Eigen::MatrixXd>(next_seed_pair[next_seed[seed]-1]);
    
    // Sum of samples update
    Eigen::MatrixXd new_unary_sum = unary_sum.matrix() + gem.row(next_seed[seed]-1);
    // Sum of samples*2 update
    Eigen::MatrixXd new_sq_sum = sq_sum.matrix() + curr_seed_sq.matrix();
    // Matrix of sample - sample paris
    Eigen::MatrixXd new_pair_sum = pair_sum + curr_seed_pair;
    
    Eigen::MatrixXd cov_ij = samp_num * new_pair_sum - new_unary_sum.transpose() * new_unary_sum;
    Eigen::MatrixXd cov_ii = samp_num * new_sq_sum - new_unary_sum.cwiseProduct(new_unary_sum);
    
    Eigen::MatrixXd a = (cov_ij.cwiseQuotient((cov_ii.transpose() * cov_ii).cwiseSqrt())).cwiseAbs();
    
    cor_scores[seed] = a.sum()/a.size();
  }
  
  return Rcpp::wrap(cor_scores);
}