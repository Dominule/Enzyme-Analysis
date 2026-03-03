# Enzyme Kinetics non-linear regression
library(ggplot2)

# Loading data
data <- read.table("enzyme_kinetics.tsv", header = FALSE, col.names = c("S", "I", "v"))
head(data)

data_control <- subset(data, I == 0)
data_inhibitor <- subset(data, I == 1)

# Fitting non-linear regression
fit_control <- nls(v ~ (Vmax * S) / (Km + S), data = data_control, start = list(Vmax = 2, Km = 1))
fit_inhibitor <- nls(v ~ (Vmax * S) / (Km + S), data = data_inhibitor, start = list(Vmax = 0.5, Km = 0.5))

cf_control <- coef(fit_control)
cf_inhibitor <- coef(fit_inhibitor)

# Plotting
ggplot(data, aes(x = S, y = v, color = factor(I))) +
  geom_point(size = 3) +
  geom_smooth(method = "nls", 
              formula = y ~ (Vmax * x) / (Km + x), 
              method.args = list(start = list(Vmax = 2, Km = 1)), 
              se = FALSE) +
  expand_limits(x = 0, y = 0) +
  labs(title = "Michaelis-Menten", x = "S [mmol/l]", y = "v [µmol/l/min]", color = "Inhibitor") +
  theme_minimal()

ggplot(data, aes(x = 1/S, y = 1/v, color = factor(I))) +
  geom_point(size = 3) +
  geom_abline(intercept = 1/cf_control["Vmax"], slope = cf_control["Km"]/cf_control["Vmax"], color = "#F8766D") +
  geom_abline(intercept = 1/cf_inhibitor["Vmax"], slope = cf_inhibitor["Km"]/cf_inhibitor["Vmax"], color = "#00BFC4") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(limits = c(-2.5, 1.1)) + 
  scale_y_continuous(limits = c(0, 4)) +
  labs(title = "Lineweaver-Burk Plot (extrapolovaný)",
       x = "1/S [l/mmol]", 
       y = "1/v [l.min/µmol]",
       color = "Inhibitor") +
  theme_minimal()

# Showing type of inhibition (uncompetitive)
v_ctrl <- coef(fit_control)["Vmax"]
k_ctrl <- coef(fit_control)["Km"]
v_inh  <- coef(fit_inhibitor)["Vmax"]
k_inh  <- coef(fit_inhibitor)["Km"]

alpha_v_prime <- v_ctrl / v_inh
alpha_k_prime <- k_ctrl / k_inh

cat("Alpha' calculated from Vmax:", round(alpha_v_prime, 3), "\n")
cat("Alpha' calculated from Km:  ", round(alpha_k_prime, 3), "\n")
cat("Difference:                ", round(abs(alpha_v_prime - alpha_k_prime), 3), "\n")

avg_alpha <- (alpha_v_prime + alpha_k_prime) / 2

# Inhibition constant Ki
Ki_final <- 1.0 / (avg_alpha - 1) # [I] = 1.0

cat("\nInhibition Constant (Ki):", round(Ki_final, 3), "µmol/l\n")