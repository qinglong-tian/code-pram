######################
# Logistic
######################

library(tidyverse)
dir_2_dat <- "dat10/"
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
    labels = c("Proposed", "Naive", "Oracle", "Model1", "Model2", "Est. Eff.")
  ),
  Parameters = factor(Parameters, levels = c("beta[0]", "beta[1]"), labels = c("beta[1]", "beta[2]"))
) -> df

# Bias
df %>% filter(Type == "Estimate") %>% mutate(
  Bias = ifelse(Parameters == "beta[1]",
                Values + 1,
                Values - 1.5)
) %>% ggplot(aes(x = n, y = Bias, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) + geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Parameters ~ p00Label, scales = "free", labeller = label_parsed) + ylab("Empirical Bias") +
  xlab("n")

# MSE
dir_2_dat <- "dat11/"
file_list <- list.files(dir_2_dat)
trueVal <- c(-1, 1.5)
df <- NULL

for (file in file_list)
{
  tmp <- readRDS(paste(dir_2_dat, file, sep = ""))
  df <- rbind(df, Read_in_Data_MSE(dir_2_dat, file, trueVal))
}

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
    labels = c("Proposed", "Naive", "Oracle", "Model1", "Model2", "Est. Eff.")
  ),
  Parameters = factor(Parameters, levels = c("beta[0]", "beta[1]"), labels = c("beta[1]", "beta[2]"))
) -> df

df %>% filter(Type == "MSE") %>% ggplot(aes(x = n, y = Values, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) +
  facet_grid(Parameters ~ p00Label, scales = "free", labeller = label_parsed) +
  ylab("MSE")

# Efficiency
df %>% filter(Type == "Efficiency") %>% ggplot(aes(x = n, y = Values)) +
  geom_line() + geom_point() +
  facet_grid(Parameters ~ p00Label, labeller = label_parsed) + ylab("Relative Efficiency") +
  xlab("n")+geom_hline(yintercept = 1, linetype = "dashed")

# For plotting the heatmap
dir2heat <- "dat3/"
files <- list.files(dir2heat)
df_all <- NULL
for (file in files)
{
  tmp <- readRDS(paste(dir2heat, file, sep = ""))
  df_all <- rbind(df_all, tmp)
}

df_all %>% filter(Type == "Efficiency") %>% mutate(RE = Values) %>%
  mutate(Parameters = factor(
    Parameters,
    levels = c("beta[0]", "beta[1]"),
    labels = c("beta[1]", "beta[2]")
  )) %>%  ggplot(aes(x = p00, y = p11, fill = RE)) +
  geom_tile() + xlab(expression(p[0][0])) + ylab(expression(p[1][1])) + facet_grid(cols = vars(Parameters), labeller = label_parsed) +
  scale_fill_gradient(low = "blue", high = "white")

# Table

df %>% filter(Type == "Estimate", Method == "Proposed") %>% mutate(
  Bias = ifelse(Parameters == "beta[1]",
                Values + 1,
                Values - 1.5)) -> df_table
bias1 <- df_table %>% filter(Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Bias)
bias2 <- df_table %>% filter(Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Bias)

df %>% filter(Type == "SE", Method == "Proposed") -> df_table
empsd1 <- df_table %>% filter(Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Values)
empsd2 <- df_table %>% filter(Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Values)

se1 <- df %>% filter(Type == "SD", Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Values)
se2 <- df %>% filter(Type == "SD", Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Values)

cp1 <- df %>% filter(Type == "Coverage", Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Values)
cp2 <- df %>% filter(Type == "Coverage", Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Values)

nn <- rep(seq(1000, 2000, by = 200), 3)

printDF <- data.frame(nn, bias1, bias2, empsd1, empsd2, se1, se2, cp1, cp2)
xtable(printDF[1:18,], digits = 3)

######################
# Linear
######################

dir2dat <- "dat20/"
files <- list.files(dir2dat)
df_all <- NULL
for(file in files)
{
  tmp <- readRDS(paste(dir2dat, file, sep = ""))
  df_all <- rbind(df_all, tmp)
}

df_all %>% mutate(
  p00Label = ifelse(
    p00 == 0.75,
    'p[1][1]*"=0.75"',
    ifelse(p00 == 0.85, 'p[1][1]*"=0.85"', 'p[1][1]*"=0.95"')
  ),
  Method = factor(
    Method,
    levels = c("Efficient", "Naive", "Oracle", "U1", "U2", "Est. Efficient"),
    labels = c("Proposed", "Naive", "Oracle", "Model1", "Model2", "Est. Eff.")
  ),
  Parameters = factor(Parameters, levels = c("beta[0]", "beta[1]"), labels = c("beta[1]", "beta[2]"))
) -> df_all

# Empirical Bias
df_all %>% filter(Type == "Estimate") %>% mutate(
  Bias = ifelse(Parameters == "beta[1]",
                Values + 1,
                Values - 1)
) %>% ggplot(aes(x = n, y = Bias, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) + geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Parameters ~ p00Label, scales = "free", labeller = label_parsed) + ylab("Empirical Bias") +
  xlab("n")

aesMatch1 <- c("Proposed" = "solid", "Oracle" = "dashed", "Model1" = "dashed")
aesMatch2 <- c("Proposed" = "circle", "Oracle" = "square", "Model1" = "plus")

df_all %>% filter(Type == "Estimate", Method %in% c("Proposed", "Oracle", "Model1")) %>% mutate(
  Bias = ifelse(Parameters == "beta[1]",
                Values + 1,
                Values - 1)
) %>% ggplot(aes(x = n, y = Bias, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) + geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Parameters ~ p00Label, scales = "free", labeller = label_parsed) + ylab("Empirical Bias") +
  xlab("n")+scale_linetype_manual(name = "Method", values = aesMatch1)+scale_shape_manual(name = "Method", values = aesMatch2)

# MSE

dir_2_dat <- "dat21/"
file_list <- list.files(dir_2_dat)
trueVal <- c(-1, 1)
df <- NULL

for (file in file_list)
{
  tmp <- readRDS(paste(dir_2_dat, file, sep = ""))
  df <- rbind(df, Read_in_Data_MSE(dir_2_dat, file, trueVal))
}

df %>% mutate(
  p00Label = ifelse(
    p00 == 0.75,
    'p[1][1]*"=0.75"',
    ifelse(p00 == 0.85, 'p[1][1]*"=0.85"', 'p[1][1]*"=0.95"')
  ),
  Method = factor(
    Method,
    levels = c("Efficient", "Naive", "Oracle", "U1", "U2", "Est. Efficient"),
    labels = c("Proposed", "Naive", "Oracle", "Model1", "Model2", "Est. Eff.")
  ),
  Parameters = factor(Parameters, levels = c("beta[0]", "beta[1]"), labels = c("beta[1]", "beta[2]"))
) -> df

df %>% filter(Type == "MSE", Method %in% c("Proposed", "Oracle", "Model1")) %>% ggplot(aes(x = n, y = Values, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) +
  facet_grid(Parameters ~ p00Label, scales = "free", labeller = label_parsed) +
  ylab("MSE")+scale_linetype_manual(name = "Method", values = aesMatch1)+scale_shape_manual(name = "Method", values = aesMatch2)

# Effciency
df %>% filter(Type == "Efficiency") %>% ggplot(aes(x = n, y = Values)) +
  geom_line() + geom_point() +
  facet_grid(Parameters ~ p00Label, labeller = label_parsed) + ylab("Relative Efficiency") +
  xlab("n")+geom_hline(yintercept = 1, linetype = "dashed")

dir2heat <- "dat4/"
files <- list.files(dir2heat)
df_all <- NULL
for (file in files)
{
  tmp <- readRDS(paste(dir2heat, file, sep = ""))
  df_all <- rbind(df_all, tmp)
}

df_all %>% filter(Type == "Efficiency") %>% mutate(RE = Values) %>%
  mutate(Parameters = factor(
    Parameters,
    levels = c("beta[0]", "beta[1]"),
    labels = c("beta[1]", "beta[2]")
  )) %>%  ggplot(aes(x = p00, y = p11, fill = RE)) +
  geom_tile() + xlab(expression(p[0][0])) + ylab(expression(p[1][1])) + facet_grid(cols = vars(Parameters), labeller = label_parsed) +
  scale_fill_gradient(low = "blue", high = "white")

# Summary

df %>% filter(Type == "Estimate", Method == "Proposed") %>% mutate(
  Bias = ifelse(Parameters == "beta[1]",
                Values + 1,
                Values - 1)) -> df_table
bias1 <- df_table %>% filter(Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Bias)
bias2 <- df_table %>% filter(Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Bias)

df %>% filter(Type == "SE", Method == "Proposed") -> df_table
empsd1 <- df_table %>% filter(Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Values)
empsd2 <- df_table %>% filter(Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Values)

se1 <- df %>% filter(Type == "SD", Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Values)
se2 <- df %>% filter(Type == "SD", Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Values)

cp1 <- df %>% filter(Type == "Coverage", Parameters == "beta[1]") %>% arrange(p00, n) %>% select(Values)
cp2 <- df %>% filter(Type == "Coverage", Parameters == "beta[2]") %>% arrange(p00, n) %>% select(Values)

nn <- rep(seq(1000, 2000, by = 200), 3)

printDF <- data.frame(nn, bias1, bias2, empsd1, empsd2, se1, se2, cp1, cp2)
xtable(printDF[1:18,], digits = 4)

