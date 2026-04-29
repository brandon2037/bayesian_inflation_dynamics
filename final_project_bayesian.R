# Data Cleaning 
library(rstan) 
library(dplyr) 
library(ggplot2) 
library(bayesplot) 


data <- read.csv( 
  "pce_data.csv", 
  skip = 5,        # skips title + revised date + year row 
  header = TRUE, 
  fill = TRUE, 
  stringsAsFactors = FALSE 
) 

raw <- read.csv("pce data - Sheet1.csv",     # raw data 
                header = FALSE, 
                stringsAsFactors = FALSE, 
                fill = TRUE) 
years  <- as.numeric(as.character(unlist(raw[1, -c(1,2)]))) # all year values 

year_table <-table(years) #freq table for years 

months <- as.character(unlist(raw[2, -c(1,2)])) # all month values 
month_table <- table(months) #freq table for months 

month_num <- match(months, 
                   c("JAN","FEB","MAR","APR","MAY","JUN",     # match values with label 
                     "JUL","AUG","SEP","OCT","NOV","DEC")) 

dates <- as.Date(paste(years, month_num, "01", sep = "-")) # change to data Format 

pce  <- as.numeric(unlist(raw[3, -c(1,2)])) #extract pce values 
core <- as.numeric(unlist(raw[4, -c(1,2)])) #extract core values 

inflation <- data.frame( #final complete model 
  date = dates, 
  year = as.numeric(format(dates, "%Y")), 
  month = as.numeric(format(dates, "%m")), 
  PCE = pce, 
  CorePCE = core 
) 

if ("pi" %in% names(inflation)) { 
  inflation <- inflation %>% mutate(pi = as.numeric(pi)) 
} else if ("PCE" %in% names(inflation)) { 
  inflation <- inflation %>% mutate(pi = as.numeric(PCE)) 
} 

inflation <- inflation %>% filter(!is.na(pi), is.finite(pi)) 

rstan_options(auto_write = TRUE) 
options(mc.cores = 2)  # 2 for now. More if runtime isn't high 

set.seed(77) 


covid_date <- as.Date("2020-03-01") 

# dates for PPC 
pre_start  <- as.Date("2010-01-01") 
pre_end    <- as.Date("2019-12-01") 
post_start <- as.Date("2021-01-01") 
post_end   <- as.Date("2024-12-01") 

delta_phi  <- 0.05     # 5 
ratio_sig  <- 1.20     #20% 


inflation <- inflation %>% 
  mutate(date = as.Date(date)) %>% 
  arrange(date) 

# parameter pi 
pi <- as.numeric(inflation$pi) 
dates_pi <- inflation$date 
T <- length(pi) 


# STAN CODE 
stan_code <- " 
data { 
  int<lower=2> T; 
  vector[T] pi; 
} 
parameters { 
  real alpha1; 
  real zphi1; 
  real h1; 
 
  vector[T-1] alpha_innov_raw; 
  vector[T-1] zphi_innov_raw; 
  vector[T-1] h_innov_raw; 
 
  real<lower=0> sigma_alpha; 
  real<lower=0> sigma_zphi; 
  real<lower=0> sigma_h; 
} 
transformed parameters { 
  vector[T] alpha; 
  vector[T] zphi; 
  vector[T] phi; 
  vector[T] h; 
 
  alpha[1] = alpha1; 
  zphi[1]  = zphi1; 
  h[1]     = h1; 
 
  for (t in 2:T) { 
    alpha[t] = alpha[t-1] + sigma_alpha * alpha_innov_raw[t-1]; 
    zphi[t]  = zphi[t-1]  + sigma_zphi  * zphi_innov_raw[t-1]; 
    h[t]     = h[t-1]     + sigma_h     * h_innov_raw[t-1]; 
  } 
 
  phi = tanh(zphi); // ensures |phi_t|<1 
} 
model { 
  // Weakly-informative priors (tuned for monthly percent-change scale) 
  alpha1 ~ normal(0, 2); 
  zphi1  ~ normal(0, 1); 
  h1     ~ normal(0, 2); 
 
  sigma_alpha ~ normal(0, 0.2); 
  sigma_zphi  ~ normal(0, 0.2); 
  sigma_h     ~ normal(0, 0.2); 
 
  alpha_innov_raw ~ std_normal(); 
  zphi_innov_raw  ~ std_normal(); 
  h_innov_raw     ~ std_normal(); 
 
  // Likelihood (include t=1 weakly) 
  pi[1] ~ normal(alpha[1], exp(0.5 * h[1])); 
  for (t in 2:T) { 
    pi[t] ~ normal(alpha[t] + phi[t] * pi[t-1], exp(0.5 * h[t])); 
  } 
} 
generated quantities { 
  vector[T] pi_rep; 
 
  // Posterior predictive replicate 
  pi_rep[1] = normal_rng(alpha[1], exp(0.5 * h[1])); 
  for (t in 2:T) { 
    pi_rep[t] = normal_rng(alpha[t] + phi[t] * pi[t-1], exp(0.5 * h[t])); 
  } 
} 
" 

sm <- stan_model(model_code = stan_code) 

# OLS done to get around error in fit 
ols_df <- data.frame(y = pi[-1], ylag = pi[-T]) 
ols_fit <- lm(y ~ ylag, data = ols_df) 

alpha_init <- as.numeric(coef(ols_fit)[1]) 
phi_init   <- as.numeric(coef(ols_fit)[2]) 
phi_init   <- max(min(phi_init, 0.95), -0.95)   # clamp for stability 
zphi_init  <- atanh(phi_init) 

h_init <- log(var(pi, na.rm = TRUE) + 1e-6)     # log variance start 

init_fun <- function() list( 
  alpha1 = alpha_init, 
  zphi1  = zphi_init, 
  h1     = h_init, 
  sigma_alpha = 0.05, 
  sigma_zphi  = 0.05, 
  sigma_h     = 0.10, 
  alpha_innov_raw = rep(0, T-1), 
  zphi_innov_raw  = rep(0, T-1), 
  h_innov_raw     = rep(0, T-1) 
) 

# 1200 iter, 600 warmup, chains 2 were the most I could set for my computer. 
fit <- sampling( 
  sm, 
  data = list(T = T, pi = pi), 
  seed = 77, 
  chains = 2, 
  iter = 1200, 
  warmup = 600, 
  init = init_fun, 
  control = list(adapt_delta = 0.995, max_treedepth = 14), 
  refresh = 200 
) 

# diagnostic check 
print(fit, pars = c("sigma_alpha","sigma_zphi","sigma_h","alpha1","zphi1","h1")) 
# r hats are all generall ~1. n_eff are good. 

summ <- summary(fit)$summary 
diag_tbl <- data.frame( 
  max_Rhat = max(summ[, "Rhat"], na.rm = TRUE), 
  min_n_eff = min(summ[, "n_eff"], na.rm = TRUE) 
) 


ex <- rstan::extract(fit) 
phi_mat <- ex$phi            
h_mat   <- ex$h               
pi_rep  <- ex$pi_rep          

sigma_mat <- exp(0.5 * h_mat) 

# Summaries 
phi_mean <- colMeans(phi_mat) 
phi_lo   <- apply(phi_mat, 2, quantile, 0.025) 
phi_hi   <- apply(phi_mat, 2, quantile, 0.975) 

sig_mean <- colMeans(sigma_mat) 
sig_lo   <- apply(sigma_mat, 2, quantile, 0.025) 
sig_hi   <- apply(sigma_mat, 2, quantile, 0.975) 

results_m1 <- data.frame( 
  date = dates_pi, 
  phi = phi_mean, phi_lo = phi_lo, phi_hi = phi_hi, 
  sigma = sig_mean, sigma_lo = sig_lo, sigma_hi = sig_hi 
) 

# Figure 1 
p_phi <- ggplot(results_m1, aes(x = date)) + 
  geom_ribbon(aes(ymin = phi_lo, ymax = phi_hi), alpha = 0.2) + 
  geom_line(aes(y = phi)) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = covid_date, linetype = "dotted") + 
  labs(title = "Model 1 (TVP-AR(1)-SV): Time-varying persistence", 
       subtitle = "Posterior mean and 95% credible interval", 
       x = "Date", y = expression(phi[t])) + 
  theme_bw() 

print(p_phi) 

pre_band <- results_m1 %>% filter(date >= pre_start, date <= pre_end) 
band_lo <- quantile(pre_band$sigma, 0.05, na.rm = TRUE) 
band_hi <- quantile(pre_band$sigma, 0.95, na.rm = TRUE) 

# Figure 2 
p_sig <- ggplot(results_m1, aes(x = date)) + 
  annotate("rect", 
           xmin = min(results_m1$date), xmax = max(results_m1$date), 
           ymin = band_lo, ymax = band_hi, alpha = 0.12) + 
  geom_ribbon(aes(ymin = sigma_lo, ymax = sigma_hi), alpha = 0.2) + 
  geom_line(aes(y = sigma), linewidth = 0.8) + 
  geom_vline(xintercept = covid_date, linetype = "dotted") + 
  labs(title = "Model 1 (TVP-AR(1)-SV): Time-varying volatility", 
       subtitle = "pre-COVID baseline (2010–2019) 5–95% range of posterior mean", 
       x = "Date", y = expression(sigma[t])) + 
  theme_bw() 

print(p_sig) 


pre_idx  <- which(dates_pi >= pre_start  & dates_pi <= pre_end) 
post_idx <- which(dates_pi >= post_start & dates_pi <= post_end) 


phi_pre_draw  <- rowMeans(phi_mat[, pre_idx, drop = FALSE]) 
phi_post_draw <- rowMeans(phi_mat[, post_idx, drop = FALSE]) 

sig_pre_draw  <- rowMeans(sigma_mat[, pre_idx, drop = FALSE]) 
sig_post_draw <- rowMeans(sigma_mat[, post_idx, drop = FALSE]) 

evidence_tbl <- data.frame( 
  quantity = c( 
    "Pr(phi_post > phi_pre)", 
    paste0("Pr(phi_post - phi_pre > ", delta_phi, ")"), 
    "Pr(sigma_post > sigma_pre)", 
    paste0("Pr(sigma_post / sigma_pre > ", ratio_sig, ")") 
  ), 
  value = c( 
    mean(phi_post_draw > phi_pre_draw), 
    mean(phi_post_draw - phi_pre_draw > delta_phi), 
    mean(sig_post_draw > sig_pre_draw), 
    mean(sig_post_draw / sig_pre_draw > ratio_sig) 
  ) 
) 

# Figure 1 
print(transform(evidence_tbl, value = round(value, 3))) 


# Figure 3 
ppc_dens_overlay(y = pi, yrep = pi_rep[sample(1:nrow(pi_rep), 200), ]) 

draw_ids <- sample(1:nrow(pi_rep), 40) 
ppc_ts_df <- data.frame(date = dates_pi, y = pi) 

rep_long <- data.frame( 
  date = rep(dates_pi, times = length(draw_ids)), 
  draw = rep(seq_along(draw_ids), each = T), 
  yrep = as.numeric(t(pi_rep[draw_ids, ])) 
) 

p_ts <- ggplot() + 
  geom_line(data = rep_long, aes(x = date, y = yrep, group = draw), alpha = 0.15) + 
  geom_line(data = ppc_ts_df, aes(x = date, y = y), linewidth = 0.7) + 
  geom_vline(xintercept = covid_date, linetype = "dotted") + 
  labs(title = "Model 1 PPC: Posterior predictive TS Graph", 
       x = "Date", y = "Inflation (percent change)") + 
  theme_bw() 

# Figure 4 
print(p_ts) 



################################ MODEL 2 
set.seed(77) 

# standardized pi 
y_raw  <- pi 
y_mean <- mean(y_raw) 
y_sd   <- sd(y_raw) 
y      <- (y_raw - y_mean) / y_sd 


K <- 2 
diag_strength <- 12.0   

stan_code_msar <- sprintf(" 
functions { 
  real obs_logpdf(real yt, real ylag, real alpha, real phi, real sigma) { 
    return normal_lpdf(yt | alpha + phi * ylag, sigma); 
  } 
} 
data { 
  int<lower=2> T; 
  vector[T] y; 
  int<lower=2> K; 
} 
parameters { 
  ordered[K] alpha;             // identify regimes by mean level 
  vector[K] zphi;               // unconstrained 
  positive_ordered[K] sigma;    // identify by volatility: sigma[1] < sigma[2] 
  simplex[K] P[K]; 
  simplex[K] pi0; 
} 
transformed parameters { 
  vector[K] phi = tanh(zphi);   // constrain |phi|<1 smoothly 
} 
model { 
  // Priors (weakly informative) 
  alpha ~ normal(0, 1.5); 
  zphi  ~ normal(0, 1.0); 
  sigma ~ lognormal(-0.3, 0.5); 
 
  // Transition priors: encourage persistence but allow switching 
  for (j in 1:K) { 
    vector[K] a = rep_vector(1.0, K); 
    a[j] = %f; 
    P[j] ~ dirichlet(a); 
  } 
  pi0 ~ dirichlet(rep_vector(2.0, K)); 
 
  // Forward algorithm (marginalize discrete states) 
  { 
    vector[K] log_fwd; 
    vector[K] log_fwd_next; 
 
    // t=1 (no lag) 
    for (k in 1:K) 
      log_fwd[k] = log(pi0[k]) + normal_lpdf(y[1] | alpha[k], sigma[k]); 
 
    for (t in 2:T) { 
      for (k in 1:K) { 
        vector[K] tmp; 
        for (j in 1:K) 
          tmp[j] = log_fwd[j] + log(P[j, k]); 
 
        log_fwd_next[k] = log_sum_exp(tmp) + 
          obs_logpdf(y[t], y[t-1], alpha[k], phi[k], sigma[k]); 
      } 
      log_fwd = log_fwd_next; 
    } 
    target += log_sum_exp(log_fwd); 
  } 
} 
generated quantities { 
  matrix[T, K] filt_prob; 
  vector[T] y_rep; 
 
  // --- Filtered probabilities --- 
  { 
    vector[K] log_fwd; 
    vector[K] log_fwd_next; 
 
    for (k in 1:K) 
      log_fwd[k] = log(pi0[k]) + normal_lpdf(y[1] | alpha[k], sigma[k]); 
 
    { 
      real c = log_sum_exp(log_fwd); 
      for (k in 1:K) filt_prob[1, k] = exp(log_fwd[k] - c); 
    } 
 
    for (t in 2:T) { 
      for (k in 1:K) { 
        vector[K] tmp; 
        for (j in 1:K) 
          tmp[j] = log_fwd[j] + log(P[j, k]); 
 
        log_fwd_next[k] = log_sum_exp(tmp) + 
          normal_lpdf(y[t] | alpha[k] + phi[k] * y[t-1], sigma[k]); 
      } 
      log_fwd = log_fwd_next; 
 
      { 
        real c = log_sum_exp(log_fwd); 
        for (k in 1:K) filt_prob[t, k] = exp(log_fwd[k] - c); 
      } 
    } 
  } 
 
  // --- Lightweight PPC: simulate one replicated series per draw --- 
  { 
    int s; 
    s = categorical_rng(pi0); 
    y_rep[1] = normal_rng(alpha[s], sigma[s]); 
 
    for (t in 2:T) { 
      s = categorical_rng(P[s]); 
      y_rep[t] = normal_rng(alpha[s] + phi[s] * y_rep[t-1], sigma[s]); 
    } 
  } 
} 
", diag_strength) 

sm_msar2 <- stan_model(model_code = stan_code_msar) 

fit2 <- sampling( 
  sm_msar2, 
  data = list(T = T, y = y, K = K), 
  seed = 77, 
  chains = 2,             
  iter = 1200, 
  warmup = 600, 
  init = 0, 
  control = list(adapt_delta = 0.995, max_treedepth = 15) 
) 
# Diagnostics are fine. Rhat near 1 and no extreme n_eff values  
print(fit2, pars = c("alpha","phi","sigma","P[1,1]","P[1,2]","P[2,1]","P[2,2]")) 
summ2 <- summary(fit2)$summary 
diag_tbl2 <- data.frame( 
  max_Rhat  = max(summ2[, "Rhat"], na.rm = TRUE), 
  min_n_eff = min(summ2[, "n_eff"], na.rm = TRUE) 
) 
ex2 <- rstan::extract(fit2) 
alpha_draw <- ex2$alpha     
phi_draw   <- ex2$phi 
sigma_draw <- ex2$sigma 
P_draw     <- ex2$P        
pi0_draw   <- ex2$pi0 
filt_prob  <- ex2$filt_prob 
y_rep      
<- ex2$y_rep 
\post_mean <- function(x) apply(x, 2, mean) 
post_q <- function(x, p) apply(x, 2, quantile, p) 
alpha_m <- post_mean(alpha_draw); alpha_lo <- post_q(alpha_draw, 0.025); alpha_hi <- 
  post_q(alpha_draw, 0.975) 
phi_m   <- post_mean(phi_draw);   phi_lo   <- post_q(phi_draw, 0.025);   phi_hi   <- post_q(phi_draw, 
                                                                                            0.975) 
sig_m   <- post_mean(sigma_draw); sig_lo   <- post_q(sigma_draw, 0.025); sig_hi   <- 
  post_q(sigma_draw, 0.975) 
# Long-run mean (scaled) mu_k = alpha_k / (1 - phik) 
mu_draw <- alpha_draw / (1 - phi_draw) 
mu_m  <- post_mean(mu_draw); mu_lo <- post_q(mu_draw, 0.025); mu_hi <- post_q(mu_draw, 0.975) 
param_tbl_scaled <- data.frame( 
  regime = 1:K, 
  alpha_mean = alpha_m, alpha_lo = alpha_lo, alpha_hi = alpha_hi, 
  phi_mean   = phi_m,   phi_lo   = phi_lo,   phi_hi   = phi_hi, 
  sigma_mean = sig_m,   sigma_lo = sig_lo,   sigma_hi = sig_hi, 
  lrmean_mean = mu_m,   lrmean_lo = mu_lo,   lrmean_hi = mu_hi 
) 
# Back-transform key items to raw inflation units (approx): 
# y_raw = y_mean + y_sd * y  => alpha_raw = y_mean + y_sd*alpha ; sigma_raw = y_sd*sigma 
param_tbl_raw <- data.frame( 
  regime = 1:K, 
  alpha_raw_mean   = y_mean + y_sd * alpha_m, 
  sigma_raw_mean   = y_sd * sig_m, 
  lrmean_raw_mean  = y_mean + y_sd * mu_m 
) 
# Figure 5 
print(round(param_tbl_scaled, 3)) 


P_hat <- apply(P_draw, c(2,3), mean) 
dur_hat <- 1 / (1 - diag(P_hat))  # expected duration in months 

# Figure 6 
print(round(P_hat, 3)) # P 
print(dur_hat) #Expected duration (months) 

# high-vol regime is regime 2 because sigma[2] > sigma[1] 
p2_mean <- apply(filt_prob[, , 2], 2, mean) 
p2_lo   <- apply(filt_prob[, , 2], 2, quantile, 0.05) 
p2_hi   <- apply(filt_prob[, , 2], 2, quantile, 0.95) 

reg_df <- data.frame(date = dates_pi, p_highvol = p2_mean, lo = p2_lo, hi = p2_hi) 

p_filt2 <- ggplot(reg_df, aes(x = date)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) + 
  geom_line(aes(y = p_highvol), linewidth = 0.8) + 
  geom_vline(xintercept = covid_date, linetype = "dotted") + 
  labs(title = "Model 2 (MS-AR(1)): Filtered probability of high-volatility regime (regime 2)", 
       x = "Date", y = "Pr(high-vol regime | y_{1:t})") + 
  theme_minimal() 
print(p_filt2) # Figure 7 

# Figure 8 
idx <- sample(1:nrow(y_rep), min(50, nrow(y_rep))) 
matplot(dates_pi, t(y_rep[idx, ]), type = "l", lty = 1, col = "gray", 
        xlab = "Date", ylab = "Scaled inflation", main = "Model 2 PPC: TS graph") 
lines(dates_pi, y, lwd = 2) 
abline(v = covid_date, lty = 3) 
legend("topright", legend = c("Replicated draws","Observed"), 
       lwd = c(1,2), bty = "n")