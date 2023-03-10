data {
    int<lower=0> n_genes;
    int<lower=0> n_samples;
    int <lower=0>counts[n_genes, n_samples];
    int <lower=0>psi[n_samples];
    real <lower=0>mean_counts[n_genes];
}

parameters {
    real <lower=0, upper=1> pi;
    real <lower=0> sigma;
    real <lower=0> phi;
    real <lower=0> tau;

    vector[n_samples] outlier;
    vector[n_samples] outlier_effect;
    vector<lower=0>[n_genes] phi_genes;
}

transformed parameters {
    vector[n_samples] log_outlier_effect;

    for (n in 1:n_samples) {
        log_outlier_effect[n] = 2^(psi[n] * outlier_effect[n] + (1- psi[n]) * 0);
    }
}

model {
    pi ~ beta(0.3, 1);
    sigma ~ gamma(1, 0.1);
    phi ~ gamma (1, 1);
    tau ~ gamma(1, 0.1);

    for (n in 1:n_samples) {
        psi[n] ~ binomial(1, pi);
        outlier[n] ~ normal(0, sigma);
        outlier_effect[n] ~ normal(outlier[n], tau);
    }
    for (g in 1:n_genes) {
        phi_genes[g] ~ gamma(phi, 1);
        for (n in 1:n_samples) {
            counts[g,n] ~ neg_binomial_2(mean_counts[g] * log_outlier_effect[n], phi_genes[g]);
        }
    }
}