// Matrix exponential with an analytic reverse-mode derivative, for use in
// the Stan models as
//   matrix matrix_exp_ad(matrix A);
// declared (without a body) in the Stan functions block.
//
// Stan's own matrix_exp() has no reverse-mode specialisation, so for a
// matrix of parameters it records every scalar operation of the Pade
// approximation and of each scaling-and-squaring step on the autodiff
// tape. That costs hundreds of KB per call, and grows with the number of
// squarings, i.e. with log2 of the norm of A -- which can be huge at the
// extreme parameter values visited by an optimizer's line search.
//
// Here exp(A) is computed in doubles, and only its K*K entries are placed
// on the tape. In the reverse pass, the adjoint of A is obtained from the
// adjoint G of exp(A) using the Frechet derivative L(X, E) of the matrix
// exponential, whose adjoint satisfies <G, L(A, E)> = <L(A^T, G), E>, so
//   adj(A) += L(A^T, G) = top-right K x K block of exp([[A^T, G], [0, A^T]])
// (Najfeld and Havel 1995; Al-Mohy and Higham 2009). This is exact, and
// costs one 2K x 2K matrix exponential in doubles.
//
// This file is included inside the model's namespace: by
// stan_meta_header.hpp for the rstan models compiled into the package, and
// by inst/stan/hmm_user_header.hpp for cmdstanr. It therefore contains no
// #include directives or namespace of its own. See matrix_exp_ad_fwd.hpp
// for why the rstan build also needs a forward declaration.

// The template signature must match the forward declaration in
// matrix_exp_ad_fwd.hpp.
template <typename T_A>
auto matrix_exp_ad(const T_A& A, std::ostream* pstream__) {
  using stan::math::var;
  if constexpr (!stan::is_var<stan::value_type_t<T_A>>::value) {
    return Eigen::MatrixXd(stan::math::matrix_exp(A));
  } else {
    stan::math::check_square("matrix_exp_ad", "A", A);
    const int K = A.rows();
    stan::arena_t<Eigen::Matrix<var, -1, -1>> A_arena = A;
    stan::arena_t<Eigen::MatrixXd> A_val = A_arena.val();
    stan::arena_t<Eigen::Matrix<var, -1, -1>> res
        = stan::math::matrix_exp(Eigen::MatrixXd(A_val));
    stan::math::reverse_pass_callback([A_arena, A_val, res, K]() mutable {
      Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2 * K, 2 * K);
      B.topLeftCorner(K, K) = A_val.transpose();
      B.bottomRightCorner(K, K) = A_val.transpose();
      B.topRightCorner(K, K) = res.adj();
      A_arena.adj() += stan::math::matrix_exp(B).topRightCorner(K, K);
    });
    return Eigen::Matrix<var, -1, -1>(res);
  }
}
