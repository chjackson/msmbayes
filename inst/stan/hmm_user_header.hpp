// Header passed to cmdstanr as `user_header` when hmm.stan is compiled
// with cmdstanr (see cmdstanr_fit()), supplying the C++ definition of
// matrix_exp_ad() declared in hmm.stan. The rstan models compiled into
// the package get the same definition via inst/include/stan_meta_header.hpp.
#include <stan/math.hpp>
namespace hmm_model_namespace {
#include "../include/matrix_exp_ad.hpp"
}
