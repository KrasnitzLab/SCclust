#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::export]]
int setBits(int n)
{
	return  __builtin_popcount (n);
}
