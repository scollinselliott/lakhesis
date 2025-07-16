#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @useDynLib lakhesis
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::NumericVector rss_rotation(arma::mat& x_r, arma::mat& ref_r, double ref_mid_rad) {
    int nr = x_r.n_rows;
    Rcpp::NumericVector rss (nr);
    arma::mat rotate (2,2);
    double theta = 0;
    double rss_ = 0;
    Rcpp::NumericVector k0 (nr);
    Rcpp::NumericVector k1 (nr);
    Rcpp::NumericVector r0 (nr);
    Rcpp::NumericVector r1 (nr);
    Rcpp::NumericVector v_ (nr);
    Rcpp::NumericVector min_ (nr);

    r0 = ref_r.col(0);
    r1 = ref_r.col(1);

    for (int i = 0; i < nr; i++) {
        rss_ = 0;
        theta = atan2(x_r(i,1), x_r(i,0)) - ref_mid_rad;
        rotate(0,0) = cos(theta);
        rotate(1,0) = sin(theta);
        rotate(0,1) = -sin(theta);
        rotate(1,1) = cos(theta);
        arma::mat x_r_rot = x_r * rotate;

        for (int k = 0; k < nr; k++) {
            k0.fill(x_r_rot(k , 0));
            k1.fill(x_r_rot(k , 1));
            v_ = ((k0 - r0) * (k0 - r0)) + ((k1 - r1) * (k1 - r1));
            Rcpp::NumericVector::iterator it = std::min_element(v_.begin(), v_.end());
            rss_ += *it;
        }

    rss(i) = rss_;

    }
    
    return rss;
}

    // rss <- c()
    // for (j in 1:nrow(x.r)) {

    //     theta <- atan2(x.r[j,2], x.r[j,1]) - ref.mid.rad
    //     rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
    //     x.r.rot <- x.r %*% rotate
    //     #x.c.rot <- x.c %*% rotate

    //     rss1 <- 0
    //     for (k in 1:nrow(x.r)) {
    //     rss1 <- rss1 + min( rowSums(( t(matrix(x.r.rot[k,],  nrow = 2, ncol = nrow(x.r.rot)))  - ref.r)^2) )
    //     }

    //     rss <- c(rss, rss1)
    // }


// NumericMatrix mmult(NumericMatrix m, NumericMatrix v)
// {
//   Environment base("package:base");
//   Function mat_Mult = base["%*%"];
//   return(mat_Mult(m, v));
// }

