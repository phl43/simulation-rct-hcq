library(tidyverse)

set.seed(1)

sensitivity_pcr <- 1
specificity_pcr <- 1

# https://pubmed.ncbi.nlm.nih.gov/32347200/
sensitivity_symptoms <- 0.91
specificity_symptoms <- 0.55

prob_pcr <- 0.25

prob_symptomatic_covid <- 0.7
prob_symptomatic_noncovid <- 0.1

prob_infection_base <- 0.15
prob_infection_treatment <- 0.1

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

significant <- rep(0, 10000)

for (i in 1:10000) {
  true_positive_treatment <- rbinom(n = 414, size = 1, prob = prob_infection_treatment)
  true_positive_control <- rbinom(n = 407, size = 1, prob = prob_infection_base)
  
  test_positive_treatment <- map_int(true_positive_treatment, compute_test)
  test_positive_control<- map_int(true_positive_control, compute_test)
  
  contingency_table <- rbind(table(test_positive_treatment), table(test_positive_control))
  
  significant[i] <- fisher.test(contingency_table)$p.value < 0.05
}

sum(significant) / length(significant) * 100