library(tidyverse)
source("summarize_beta_logistic.R")
# Need: betaXY

# Coverage
df %>% filter(Type == "Coverage") %>% ggplot(aes(x = p00, y = Values, col = Parameters)) +
  geom_line() + geom_point(aes(shape = Parameters)) + geom_hline(yintercept = 0.95, linetype = "dashed") +
  ylab("Coverage Probability") + xlab(expression(p[11])) + scale_color_discrete(labels = c(expression(beta[0]), expression(beta[1]))) +
  scale_shape_discrete(labels = c(expression(beta[0]), expression(beta[1])))+ylim(c(0.9, 1))

# Estimate
df %>% filter(Type == "Estimate") %>% ggplot(aes(x = p00, y = Values, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) + geom_hline(aes(yintercept = ifelse(Parameters == "beta[0]", betaXY[1], betaXY[2])), linetype = "dashed") +
  facet_wrap( ~ Parameters, scales = "free", labeller = label_parsed) + ylab("Estimated Values") +
  xlab(expression(p[11]))

# SE
df %>% filter(Type %in% c("SE", "SD")) %>% ggplot(aes(x = p00, y = Values, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) +
  facet_wrap( ~ Parameters, scales = "free", labeller = label_parsed) + ylab("Standard Deviation") +
  xlab(expression(p[11]))

# Efficiency
df %>% filter(Type == "Efficiency") %>% ggplot(aes(x = p00, y = Values)) +
  geom_line() + geom_point() +
  facet_wrap( ~ Parameters, labeller = label_parsed) + ylab("Effciency") +
  xlab(expression(p[11]))+geom_hline(yintercept = 1, linetype = "dashed")
