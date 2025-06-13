#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;

  // Data
  DATA_STRING(dist);
  DATA_VECTOR(c_i);  // counts for observation i
  DATA_VECTOR(depth_i_obs);  // depth at observations

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  // DATA_STRUCT(spde, spde_t);

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  // DATA_SPARSE_MATRIX(A_gs);

  // Parameters
  PARAMETER(beta0);
  PARAMETER_VECTOR(beta_j);
  PARAMETER(ln_tau);
  PARAMETER(ln_kappa);

  // Random effects
  PARAMETER_VECTOR(omega_s);

  // Objective funcction
  Type jnll = 0.0;

  // Derived quantities
  Eigen::SparseMatrix<Type> Q = (exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)* M1 + M2) * exp(2*ln_tau);
  jnll += GMRF(Q)(omega_s);

  // Eigen::SparseMatrix<Type> Q = R_inla::Q_spde(spde, exp(ln_kappa));
  // jnll += SCALE(GMRF(Q), 1. / exp(ln_tau))(omega_s);

  // Probability of random effects

  vector<Type> depth_i = depth_i_obs;

  // vector<Type> depth_g = A_gs * depth_s;
  // vector<Type> omega_g = A_gs * omega_s;
  // REPORT(depth_s);
  // REPORT(depth_i);
  // REPORT(depth_g);
  vector<Type> omega_i = A_is * omega_s;
  vector<Type> pdepth_i = depth_i * beta_j(0);
  vector<Type> mu_i = exp(beta0 + omega_i + pdepth_i);

  // Probability of data conditional on random effects
  if (dist=="Poisson"){
    for(int i=0; i<c_i.size(); i++){
      jnll -= dpois(c_i(i), mu_i(i), true);
      // SIMULATE{c_i(i) = rpois(mu_i(i));}
    }
  }
  if (dist=="Tweedie"){
    PARAMETER(ln_phi);
    PARAMETER(finv_power);
    for(int i=0; i<c_i.size(); i++){
      jnll -= dtweedie(c_i(i), mu_i(i), exp(ln_phi), 1.0 + invlogit(finv_power), true);
      // SIMULATE{c_i(i) = rtweedie(mu_i(i), exp(ln_phi), 1.0 + invlogit(finv_power));}
    }
  }
  if (dist=="LNP"){
    PARAMETER(ln_sigma_eta);
    PARAMETER_VECTOR(eta_i);
    for(int i=0; i<c_i.size(); i++){
      jnll -= dnorm(eta_i(i), Type(0.0), exp(ln_sigma_eta), true);
      jnll -= dpois(c_i(i), mu_i(i) * exp(eta_i(i)), true);
      // SIMULATE{
      //   eta_i(i) = rnorm(Type(0.0), exp(ln_sigma_eta));
      //   c_i(i) = rpois(mu_i(i) * exp(eta_i(i)));
      // }
    }
  }
  // vector<Type> mu_g = exp(beta0 + omega_g + pdepth_g);

  // Reporting
  REPORT(Q);
  //REPORT(omega_s);
  // REPORT(omega_g);
  REPORT(mu_i);
  // REPORT(mu_g);
  // REPORT(pdepth_i);
  // REPORT(pdepth_g);

  // SIMULATE {
  //   REPORT(c_i);
  // }

  return jnll;
}
