
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
  DATA_VECTOR(weights_i); // new: weights for calculating cAIC
  DATA_VECTOR(pop_dens_s);  // counts for observation i
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

  // Objective function
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
    //matrix<Type> depth_sz = lu.solve(depth_s.matrix());;
    matrix<Type> pop_dens_sz = lu.solve(pop_dens_s.matrix());;

    //depth_s = depth_sz.col(0); // turn one column matrix to vector
    pop_dens_s = pop_dens_sz.col(0); // turn one column matrix to vector

    //Type SigmaD = 1 / sqrt(4 * M_PI * exp(2*ln_tau2) * exp(2*ln_kappa2));
    Type Range2 = sqrt(8) / exp(ln_kappa2);
    REPORT(Range2);
    //REPORT(SigmaD);
    REPORT(ln_kappa2);
  }
  vector<Type> pop_dens_i = A_is * pop_dens_s;
  vector<Type> pop_dens_g = A_gs * pop_dens_s;
  REPORT(pop_dens_s);
  REPORT(pop_dens_i);
  REPORT(pop_dens_g);
  vector<Type> omega_i = A_is * omega_s;
  vector<Type> omega_g = A_gs * omega_s;
  vector<Type> ppop_dens_i = pop_dens_i*beta_j(0);
  vector<Type> ppop_dens_g = pop_dens_g*beta_j(0);

  // Probability of data conditional on random effects
  // vector<Type> mu_i = exp(beta0 + omega_i + ppop_dens_i);
  // if(dist=="Poisson"){
  //   for(int i=0; i<c_i.size(); i++){
  //     jnll -= dpois(c_i(i), mu_i(i), true);
  //     SIMULATE{c_i(i) = rpois(mu_i(i));}
  //   }
  // }
  // I need to add weights here so that I can set them to 0 when calculating conditional AIC
  vector<Type> mu_i = exp(beta0 + omega_i + ppop_dens_i);
  if(dist=="Poisson"){
  for(int i=0; i<c_i.size(); i++){
    if (weights_i(i) > Type(0.0)) {
      jnll -= dpois(c_i(i), mu_i(i), true);
      }
    SIMULATE {
      if (weights_i(i) > Type(0.0)) {
        c_i(i) = rpois(mu_i(i));
        }
      }
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
  // if(dist=="LNP"){
  //   PARAMETER(ln_sigma_eta);
  //   PARAMETER_VECTOR(eta_i);
  //   for(int i=0; i<c_i.size(); i++){
  //     jnll -= dnorm(eta_i(i), Type(0.0), exp(ln_sigma_eta), true);
  //     jnll -= dpois(c_i(i), mu_i(i) * exp(eta_i(i)), true);
  //     SIMULATE{
  //       eta_i(i) = rnorm(Type(0.0), exp(ln_sigma_eta));
  //       c_i(i) = rpois(mu_i(i) * exp(eta_i(i)));
  //     }
  //   }
  // }
  if(dist=="LNP"){
    PARAMETER(ln_sigma_eta);
    PARAMETER_VECTOR(eta_i);
    for(int i=0; i<c_i.size(); i++){
      if (weights_i(i) > Type(0.0)) {
        jnll -= dnorm(eta_i(i), Type(0.0), exp(ln_sigma_eta), true);
        jnll -= dpois(c_i(i), mu_i(i) * exp(eta_i(i)), true);
      }
      SIMULATE {
        if (weights_i(i) > Type(0.0)) {
          eta_i(i) = rnorm(Type(0.0), exp(ln_sigma_eta));
          c_i(i) = rpois(mu_i(i) * exp(eta_i(i)));
        }
      }
    }
  }
  vector<Type> mu_g = exp(beta0 + omega_g + ppop_dens_g);

  // Reporting
  REPORT(Q);
  REPORT(omega_g); //omega_s
  REPORT(mu_i);
  REPORT(mu_g);
  REPORT(ppop_dens_i);
  REPORT(ppop_dens_g);

  SIMULATE {
    REPORT(c_i);
  }

  return jnll;
}
