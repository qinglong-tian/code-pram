library(tidyverse)
dir_2_dat <- "dat1/"
file_list <- list.files(dir_2_dat)

trueVal <- c(-1, 1.5)
df <- NULL

for (file in file_list)
{
  tmp <- readRDS(paste(dir_2_dat, file, sep = ""))
  df <- rbind(df, tmp)
}

nVec <- seq(from = 1000, to = 2000, by = 200)
pVec <- seq(from = 0.75, to = 0.95, by = 0.1)


npList <- expand.grid(pVec, nVec)
df %>% mutate(p00 = npList[p00, 1], p11 = npList[p11, 1]) -> df
df %>% mutate(
  p00Label = ifelse(
    p00 == 0.75,
    'p[1][1]*"=0.75"',
    ifelse(p00 == 0.85, 'p[1][1]*"=0.85"', 'p[1][1]*"=0.95"')
  ),
  Method = factor(
    Method,
    levels = c("Efficient", "Naive", "Oracle", "U1", "U2", "Est. Efficient"),
    labels = c("Efficient", "Naive", "Oracle", "C1", "C2", "Est. Eff.")
  )
) -> df

# Coverage

df %>% filter(Type == "Coverage") %>% ggplot(aes(x = n, y = Values)) +
  geom_line() + geom_point() + geom_hline(yintercept = 0.95, linetype = "dashed") +
  ylab("Coverage Probability") + xlab("n") +
  ylim(c(0.9, 1)) + facet_grid(rows = vars(Parameters),
                               cols = vars(p00Label),
                               labeller = label_parsed) + theme(axis.text.x = element_text(angle =
                                                                                             90, hjust = 1))

# Estimate
df %>% filter(Type == "Estimate") %>% mutate(
  Bias = ifelse(Parameters == "beta[0]",
                Values + 1,
                Values - 1.5)
  ) %>% ggplot(aes(x = n, y = Bias, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) + geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Parameters ~ p00Label, scales = "free", labeller = label_parsed) + ylab("Bias") +
  xlab("n") + theme(axis.text.x = element_text(angle =
                                                 90, hjust = 1))

# SE
df %>% filter(Type %in% c("SE", "SD"))%>% ggplot(
  aes(x = n, y = Values, col = Method)
) + geom_line(
  aes(linetype = Method)
) + geom_point(
  aes(shape = Method)
) + facet_grid(Parameters ~ p00Label, scales = "free", labeller = label_parsed) + ylab("Standard Deviation") +
  xlab("n") + theme(axis.text.x = element_text(angle =
                                                 90, hjust = 1))

# Efficiency
df %>% filter(Type == "Efficiency") %>% ggplot(aes(x = n, y = Values)) +
  geom_line() + geom_point() +
  facet_grid(Parameters ~ p00Label, labeller = label_parsed) + ylab("Effciency") +
  xlab(expression(p[11]))+geom_hline(yintercept = 1, linetype = "dashed")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# For plotting the heatmap
dir2heat <- "dat3/"
files <- list.files(dir2heat)
df_all <- NULL
for (file in files)
{
  tmp <- readRDS(paste(dir2heat, file, sep = ""))
  df_all <- rbind(df_all, tmp)
}

df_all %>% filter(Type == "Efficiency") %>% mutate(Efficiency = Values) %>%  ggplot(aes(x = p00, y = p11, fill = Efficiency)) +
  geom_tile() + xlab(expression(p[11])) + ylab(expression(p[22])) + facet_grid(cols = vars(Parameters), labeller = label_parsed) +
  scale_fill_gradient(low = "blue", high = "white")

                                                                         