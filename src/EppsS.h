#ifndef EPPSS_H
#define EPPSS_H

#include <vector>

/* Epps-Singleton statistic for sorted values v with counts cx and cy of the
   two samples. If rank is not NULL it receives the rank of the estimated
   covariance matrix (the degrees of freedom of the asymptotic null distribution). */
double eppss_core(const std::vector<double>& v,
                  const std::vector<double>& cx,
                  const std::vector<double>& cy,
                  int* rank);

#endif
