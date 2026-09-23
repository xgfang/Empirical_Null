// file: src/utils.h
#ifndef UTILS_H
#define UTILS_H

#include <cmath>

inline double pnorm_cpp(double x) {
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}


inline arma::vec pnorm_vec(const arma::vec& x) {
  arma::vec output = x;
  output.transform([](double val) { return pnorm_cpp(val); });
  return output;
}


#endif