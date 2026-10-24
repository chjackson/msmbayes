// Forward declaration of matrix_exp_ad() (defined in matrix_exp_ad.hpp),
// force-included (-include, see configure / configure.win) at the top of
// every C++ file compiled in src/.
//
// This is needed because the stanc bundled with rstan emits no declaration
// for functions declared without a body in the Stan program, and rstantools
// includes stan_meta_header.hpp (and hence the definition) only after the
// Stan-generated functions that call matrix_exp_ad(). Those calls are in
// templates, but the arguments' types are Eigen/Stan types, so
// argument-dependent lookup at instantiation would not find a definition in
// the model's namespace; the name must already be declared.
#ifndef MSMBAYES_MATRIX_EXP_AD_FWD_HPP
#define MSMBAYES_MATRIX_EXP_AD_FWD_HPP
#include <iosfwd>
namespace model_hmm_namespace {
template <typename T_A>
auto matrix_exp_ad(const T_A& A, std::ostream* pstream__);
}
namespace model_msm_namespace {
template <typename T_A>
auto matrix_exp_ad(const T_A& A, std::ostream* pstream__);
}
#endif
