// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// Squared Euclidean distance: Taking the square root is not necessary
// because EEMS uses the distances to find closest points. For example,
//   pairwise_distance(X,Y).col(0).minCoeff( &closest )
// finds the row/point in X that is closest to the first row/point in Y
Eigen::MatrixXd euclidean_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y) {
    return(  X.rowwise().squaredNorm().eval().replicate(1,Y.rows())
             + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(),1)
             - 2.0 * X * Y.transpose() );
}

//' @title compute contour vals for each pixel
//'
//' @description compute one contour, by filling in each of the pixels/marks
//'
//' @param marks
//' @param now_rates
//' @param now_seeds
//'
//' @return values vectoring storing the value for each pixel
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd compute_contour_vals(const Eigen::MatrixXd &marks,
                                     const Eigen::VectorXd &now_rates, const Eigen::MatrixXd &now_seeds) {
    int nmrks = marks.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Zero(nmrks);
    Eigen::MatrixXd distances = euclidean_dist(marks, now_seeds);
    for ( int row = 0, closest = 0 ; row < nmrks ; row++ ) {
        distances.row(row).minCoeff( &closest );
        values(row) = now_rates(closest);
        
    }
    
    return values;
}

