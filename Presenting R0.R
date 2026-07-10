
dist.mat <- as.matrix(dist(cds))

K <- (b / (2 * pi * theta^2 * gamma(2 / b))) *
  exp(-(dist.mat / theta)^b)

diag(K) <- 0

row.sum <- rowSums(K)

Kbar <- mean(row.sum)

Kbar

failed.epidemics <- history.all %>%
  group_by(sim, omega, Nvectors) %>%
  summarise(
    max.I = max(I),
    .groups = "drop"
  ) %>%
  group_by(omega, Nvectors) %>%
  summarise(
    failed = sum(max.I <= 5),
    total = n(),
    proportion.failed = failed / total,
    .groups = "drop"
  )

failed.epidemics$R0 <- (failed.epidemics$omega*failed.epidemics$Nvectors*Kbar)/(nrow(cds)*removal.rate)

failed.epidemics$proportion.success <-
  1 - failed.epidemics$proportion.failed

library(ggplot2)

ggplot(failed.epidemics,
       aes(R0, proportion.success)) +
  
  geom_point(size = 4) +
  
  geom_smooth(method = "glm",
              method.args = list(family = binomial),
              se = FALSE) +
  
  labs(
    x = expression(R[0]),
    y = "Probability of epidemic establishment"
  ) +
  
  theme_bw()