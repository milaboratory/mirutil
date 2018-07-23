#include <Rcpp.h>
using namespace Rcpp;

int getId(char c) {
  switch(c) {
  case 'A':
    return 0;
  case 'T':
    return 1;
  case 'G':
    return 2;
  case 'C':
    return 3;
  }
  return -1;
}

int getId(char c1, char c2) {
  return 4 * getId(c1) + getId(c2);
}

CharacterVector nt1 = CharacterVector::create('A', 'A', 'A', 'A',
                                              'T', 'T', 'T', 'T',
                                              'G', 'G', 'G', 'G',
                                              'C', 'C', 'C', 'C');
CharacterVector nt2 = CharacterVector::create('A', 'T', 'G', 'C',
                                              'A', 'T', 'G', 'C',
                                              'A', 'T', 'G', 'C',
                                              'A', 'T', 'G', 'C');

// Get unweighted and weighted dinucleotide frequencies.
// Both forward and reverse sequences
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
DataFrame getDiNtFreq(StringVector seqs,
                      NumericVector wts,
                      LogicalVector reverse) {
  NumericVector f(16), f_wt(16);

  if (reverse[0] == TRUE) { // 3'->5'
    for (int i = 0; i < seqs.size(); ++i) {
      auto chars = seqs[i];
      auto wt = wts[i];

      int nchar = seqs[i].size();

      // last base with unknown previous
      int idLast = getId(chars[nchar-1]);
      for (int k = 0; k < 4; ++k) {
        int id = 4 * k + idLast;
        f[id] += 0.25;
        f_wt[id] += 0.25 * wt;
      }

      // remaining dinucleotides (reversed)
      for (int j = 0; j < nchar - 1; ++j) {
        int id = getId(chars[j + 1], chars[j]);
        f[id]++;
        f_wt[id] += wt;
      }
    }
  } else { // 5'->3'
    for (int i = 0; i < seqs.size(); ++i) {
      auto chars = seqs[i];
      auto wt = wts[i];

      // first base with unknown previous
      int id0 = getId(chars[0]);
      for (int k = 0; k < 4; ++k) {
        int id = 4 * k + id0;
        f[id] += 0.25;
        f_wt[id] += 0.25 * wt;
      }

      // remaining dinucleotides
      for (int j = 1; j < seqs[i].size(); ++j) {
        int id = getId(chars[j - 1], chars[j]);
        f[id]++;
        f_wt[id] += wt;
      }
    }
  }

  return DataFrame::create(Named("nt.1") = nt1,
                           Named("nt.2") = nt2,
                           Named("count.clonotypes") = f,
                           Named("count.reads") = f_wt);
}

/*** R
# Dont forget devtools::document()
getDiNtFreq("ATGC", 2, F)
getDiNtFreq(c("ATGC", "CGTA"), c(2, 5), F)
getDiNtFreq("CGTA", 2, T)
*/
