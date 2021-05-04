
functions {
  real[] system(real time,
                real[] y,
                real[] theta,
                real[] x_r, int[] x_i) {
    real ka = theta[1];
    real V = theta[2];
    real Vm = theta[3];
    real Km = theta[4];
    real C = y[2] / V;
    real dydt[2];

    dydt[1] = - ka * y[1];
    dydt[2] = ka * y[1] - Vm * C / (Km + C);
    return dydt;              
  }
}

data {
  int N;
  vector[N] y;
  real t[N];
  int stiff_solver;
  real y0[2];  // initial drug mass in the gut.

  int n_patients;
  int<lower = 1> start[n_patients];
  int<lower = 1> end[n_patients];
}

transformed data {
  int n_obs = N / n_patients;
  int n_parm = 4;
  int n_cmt = 2;
  real t0 = 0;
  real x_r[0];
  int x_i[0];
  // vector<lower = 0>[n_patients] omega = rep_vector(0.25, 5);
  real rel_tol = 1e-6;
  real abs_tol = 1e-6;
  int max_num_steps = 10000;  // 1e4
}

parameters {
  real<lower = 0> ka_pop;
  real<lower = 0> V_pop;
  real<lower = 0> Vm_pop;
  real<lower = 0> Km_pop;
  vector<lower = 0>[n_parm] omega;
  real<lower = 0> sigma;
  real<lower = 0> theta[n_parm, n_patients];
  // real eta[n_parm, n_patients];
}

transformed parameters {
  vector[N] concentration;

  for (i in 1:n_patients) {
    real x[n_obs, n_cmt];

    if (!stiff_solver) {
      x = integrate_ode_rk45(system, y0, t0, t[start[i]:end[i]], 
                             theta[, i], x_r, x_i,
                             rel_tol, abs_tol, max_num_steps);
    } else {
      x = integrate_ode_bdf(system, y0, t0, t[start[i]:end[i]], 
                            theta[, i], x_r, x_i,
                            rel_tol, abs_tol, max_num_steps);
    }

    concentration[start[i]:end[i]] = to_vector(x[, 2]) / theta[2, i];
  }
}

model {
  // priors -- CHECK if these are reasonable values
  ka_pop ~ lognormal(log(2.5), 3);  // old model uses 3 for sd.
  V_pop ~ lognormal(log(35), 0.5);
  Vm_pop ~ lognormal(log(10), 0.5);
  Km_pop ~ lognormal(log(2.5), 3);
  sigma ~ normal(0, 1);

  omega ~ lognormal(0.25, 0.1);  // CHECK

  theta[1, ] ~ lognormal(log(ka_pop), omega[1]);
  theta[2, ] ~ lognormal(log(V_pop), omega[2]);
  theta[3, ] ~ lognormal(log(Vm_pop), omega[3]);
  theta[4, ] ~ lognormal(log(Km_pop), omega[4]);

  // likelihood
  y ~ lognormal(log(concentration), sigma);
}

generated quantities {
  real y_pred[N] = lognormal_rng(log(concentration), sigma);
}
