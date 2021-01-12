
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
}

transformed data {
  int n_parm = 4;
  int n_cmt = 2;
  real t0 = 0;
  real x_r[0];
  int x_i[0];
  // real rel_tol = 1e-6;
  // real abs_tol = 1e-6;
  // int max_num_steps = 10000;  // 1e4
}

parameters {
  real<lower = 0> ka;
  real<lower = 0> V;
  real<lower = 0> Vm;
  real<lower = 0> Km;
  real<lower = 0> sigma;
}

transformed parameters {
  real theta[4] = {ka, V, Vm, Km};
  vector[N] concentration;

  {
    real x[N, n_cmt];
    if (!stiff_solver) {
      x = integrate_ode_rk45(system, y0, t0, t, theta, x_r, x_i);
                             // rel_tol, abs_tol, max_num_steps);
    } else {
      x = integrate_ode_bdf(system, y0, t0, t, theta, x_r, x_i);
                            // rel_tol, abs_tol, max_num_steps);
    }

    concentration = to_vector(x[, 2]) / V;
  }
}

model {
  // priors -- CHECK if these are reasonable values
  ka ~ lognormal(log(2.5), 3);  // 1.85
  V ~ lognormal(log(35), 0.5);
  Vm ~ lognormal(log(10), 0.5);
  Km ~ lognormal(log(2.5), 3);  // 2
  sigma ~ normal(0, 1);

  // likelihood
  y ~ normal(concentration, sigma);
}

generated quantities {
  real y_pred[N] = normal_rng(concentration, sigma);
}
