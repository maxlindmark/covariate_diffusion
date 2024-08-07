
#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_STRING(method);
  DATA_STRING(dist);
  DATA_VECTOR(c_i);  // counts for observation i
  DATA_VECTOR(depth_s);  // counts for observation i
  DATA_INTEGER(sim_gmrf) // simulate GMRFs?

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_SPARSE_MATRIX(invsqrtM0);
  DATA_SPARSE_MATRIX(invM0);

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);

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
  Eigen::SparseMatrix<Type> Q = (exp(4*ln_kappa)*M0 +
    Type(2.0)*exp(2*ln_kappa)* M1 + M2) * exp(2*ln_tau);
  jnll += GMRF(Q)(omega_s);

  if (sim_gmrf) SIMULATE {GMRF(Q).simulate(omega_s);}

  // Probability of random effects
  if(method == "diffusion"){
    PARAMETER(ln_kappa2);
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu;
    Eigen::SparseMatrix<Type> invD = invM0 * (M0 + exp(2*ln_kappa2)*M0 + M1) /
      (1.0 + exp(2*ln_kappa2)); // Denominator corrects for proportionality constant
    lu.compute(invD);
    REPORT(invD);

    // Solve
    matrix<Type> depth_sz = lu.solve(depth_s.matrix());;

    depth_s = depth_sz.col(0); // turn one column matrix to vector

    //Type SigmaD = 1 / sqrt(4 * M_PI * exp(2*ln_tau2) * exp(2*ln_kappa2));
    Type Range2 = sqrt(8) / exp(ln_kappa2);
    REPORT(Range2);
    //REPORT(SigmaD);
    REPORT(ln_kappa2);
  }
  vector<Type> depth_i = A_is * depth_s;
  vector<Type> depth_g = A_gs * depth_s;
  REPORT(depth_s);
  REPORT(depth_i);
  REPORT(depth_g);
  vector<Type> omega_i = A_is * omega_s;
  vector<Type> omega_g = A_gs * omega_s;
  vector<Type> pdepth_i = depth_i*beta_j(0) + pow(depth_i,2)*beta_j(1);
  vector<Type> pdepth_g = depth_g*beta_j(0) + pow(depth_g,2)*beta_j(1);

  // Probability of data conditional on random effects
  vector<Type> mu_i = exp(beta0 + omega_i + pdepth_i);
  if(dist=="Poisson"){
    for(int i=0; i<c_i.size(); i++){
      jnll -= dpois(c_i(i), mu_i(i), true);
      SIMULATE{c_i(i) = rpois(mu_i(i));}
    }
  }
  if(dist=="Tweedie"){
    PARAMETER(ln_phi);
    PARAMETER(finv_power);
    for(int i=0; i<c_i.size(); i++){
      jnll -= dtweedie(c_i(i), mu_i(i), exp(ln_phi), 1.0 + invlogit(finv_power), true);
      SIMULATE{c_i(i) = rtweedie(mu_i(i), exp(ln_phi), 1.0 + invlogit(finv_power));}
    }
  }
  if(dist=="LNP"){
    PARAMETER(ln_sigma_eta);
    PARAMETER_VECTOR(eta_i);
    for(int i=0; i<c_i.size(); i++){
      jnll -= dnorm(eta_i(i), Type(0.0), exp(ln_sigma_eta), true);
      jnll -= dpois(c_i(i), mu_i(i) * exp(eta_i(i)), true);
      SIMULATE{
        eta_i(i) = rnorm(Type(0.0), exp(ln_sigma_eta));
        c_i(i) = rpois(mu_i(i) * exp(eta_i(i)));
      }
    }
  }
  vector<Type> mu_g = exp(beta0 + omega_g + pdepth_g);

  // Reporting
  REPORT(Q);
  REPORT(omega_g); //omega_g ??
  REPORT(mu_i);
  REPORT(mu_g);
  REPORT(pdepth_i);
  REPORT(pdepth_g);

  SIMULATE {
    REPORT(c_i);
  }

  return jnll;
}
