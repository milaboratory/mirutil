#include <Rcpp.h>
using namespace Rcpp;

void safeNormalize(NumericVector p, bool apc) {
  double Sp = 0;

  for (int i = 0; i < p.size(); i++) {
    if (R_IsNA(p[i])) {
      p[i] = 0;
    }
    if (apc) {
      p[i] += 1;
    }
    Sp += p[i];
  }

  for (int i = 0; i < p.size(); i++) {
    p[i] /= Sp;
  }
}

// Jensen-Shannon divergence
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double jsdDistSafe(NumericVector p,
                   NumericVector q,
                   LogicalVector addPseudocount) {
  bool apc = addPseudocount[0] == TRUE;
  safeNormalize(p, apc);
  safeNormalize(q, apc);

  double res = 0;

  for (int i = 0; i < p.size(); i++) {
    auto Pi = p[i], Qi = q[i], m = (Pi + Qi) / 2.0;
    if(Pi != 0) {
      res += Pi * log(Pi / m);
    }
    if(Qi != 0) {
      res += Qi * log(Qi / m);
    }
  }

  return sqrt(res / 2.0 / log(2));
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
jsdDistSafe(c(0.1, 0.5, 0.4), c(0.1, 0.5, 0.4), F)
jsdDistSafe(c(1.0, 0.0, 0.0), c(0.0, 0.5, 0.5), F)
jsdDistSafe(c(1.0, 0.0, NA), c(0.0, 0.5, 0.5), F)
jsdDistSafe(c(0.1 * 10, 0.5 * 10, 0.4 * 10), c(0.1, 0.5, 0.4), F)
jsdDistSafe(c(0.1 * 1000, 0.5 * 1000, 0.4 * 1000),
            c(0.1 * 100, 0.5 * 100, 0.4 * 100), T)
*/
