library(tidyverse)

set.seed(1)

sensitivity_pcr <- 1
specificity_pcr <- 1

# https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.16.2000508
sensitivity_symptoms <- 0.91
specificity_symptoms <- 0.55

prob_pcr <- 0.2

prob_symptomatic_covid <- 0.8
prob_symptomatic_noncovid <- 0.05

prob_infection_base <- 0.15
prob_infection_treatment <- 0.1125

compute_test <- function (x) {
  if (x == 1) {
    if (rbernoulli(1, p = prob_symptomatic_covid)) {
      if (rbernoulli(1, p = prob_pcr)) {
        return (as.integer(rbernoulli(1, p = sensitivity_pcr)))
      } else {
        return (as.integer(rbernoulli(1, p = sensitivity_symptoms)))
      }
    } else {
      return (0L)
    }
  } else {
    if (rbernoulli(1, p = prob_symptomatic_noncovid)) {
      if (rbernoulli(1, p = prob_pcr)) {
        return (as.integer(rbernoulli(1, p = 1 - specificity_pcr)))
      } else {
        return (as.integer(rbernoulli(1, p = 1 - specificity_symptoms)))
      }
    } else {
      return (0L)
    }
  }
}

n_simulations <- 10000

results <- tibble(
  significant = rep(0, n_simulations),
  effect = rep(0, n_simulations)
)

for (i in 1:n_simulations) {
  true_positive_treatment <- rbinom(n = 414, size = 1, prob = prob_infection_treatment)
  true_positive_control <- rbinom(n = 407, size = 1, prob = prob_infection_base)
  
  test_positive_treatment <- map_int(true_positive_treatment, compute_test)
  test_positive_control<- map_int(true_positive_control, compute_test)
  
  contingency_table <- rbind(table(test_positive_treatment), table(test_positive_control))
  
  results$significant[i] <- fisher.test(contingency_table)$p.value < 0.05
  results$effect[i] <- sum(test_positive_control) / length(test_positive_control) -
    sum(test_positive_treatment) / length(test_positive_treatment)
}

observed_effect <- 2.5

mean_effect <- mean(results$effect)

min_significant <- min(filter(results, significant == TRUE & effect > 0)$effect)

ggplot(results, aes(x = effect * 100)) + 
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  geom_vline(aes(xintercept = observed_effect, linetype = "oe"), color = "blue", size = 1) +
  geom_vline(aes(xintercept = mean_effect * 100, linetype = "moe"), color = "red", size = 1) +
  geom_vline(aes(xintercept = (prob_infection_base - prob_infection_treatment) * 100, linetype = "re"), color = "purple", size = 1) +
  theme_minimal() +
  ylab("Density") +
  xlab("Observed effect size") +
  scale_linetype_manual(
    name = "lines",
    values = c(
      "oe" = 1,
      "moe" = 1,
      "re" = 1
    ),
    labels = c(
      "Observed effect",
      "Mean observed effect",
      "Real effect"
    ),
    guide = guide_legend(
      title = "",
      override.aes = list(color = c("blue", "red", "purple"))
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  ggsave(
    "Distribution of observed effect sizes in simulation of Boulware et al.'s RCT of hydroxychloroquine as post-exposure prophylactic treatment.png",
    width = 12,
    height = 6
  )