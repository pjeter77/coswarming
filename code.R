library("brms")
library("dplyr")
library("bayesplot")
library("igraph")
library("reshape2")
library("hrbrthemes")
library("clipr")
library("ggplot2")
library("rstan")
library("cmdstanr")
library("broom")
library("loo")
library("scales")
library("FSA")
library("rcompanion")

################################
# Negative binomial brms model #
################################

data <- read.csv("brms_data.csv", header = T, sep = ",")
head(data, n = 25)

fit_nb <- brm(
  coswarming ~ sexage * relatedness + offset(log(max_co_occurrence * days + 1)) +
    (1 | mm(bat1, bat2)),
  data = data,
  family = negbinomial(),
  save_pars = save_pars(latent = TRUE, all = FALSE),
  chains = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

saveRDS(fit_nb, "fit_nb.rds", compress = "xz")
fit_nb <- readRDS("fit_nb.rds")

summary(fit_nb)
p <- plot(fit_nb)

pp_check(fit_nb, type)
pp_check(fit_nb, type = "stat", bins = 50, stat = "max")

ce <- conditional_effects(fit_nb, effects = "sexage:relatedness")
plot(ce, points = FALSE)

conds <- data.frame(relatedness = c(0.00, 0.25, 0.50))
ce <- conditional_effects(fit_nb, effects = "sexage", conditions = conds)
class(p)
p <- plot(ce, points = FALSE)
p[[1]] + 
  scale_y_log10() +
  hrbrthemes::theme_ipsum()

######################################################
# Zero-inflated negative binomial (ZINB) brms model  #
######################################################

data <- read.csv("brms_data.csv", header = T, sep = ",")

head(data, n = 25)

formula_zinb <- bf(
  coswarming ~ sexage * relatedness + offset(log(max_co_occurrence * days + 1)) +
    (1 | mm(bat1, bat2)),
  zi ~ 1
)

priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(-8, 3), class = "Intercept"),
  prior(student_t(3, 0, 1), class = "sd"),
  prior(gamma(2, 0.2), class = "shape"),
  prior(normal(0, 1), class = "Intercept", dpar = "zi")
)

fit_zinb <- brm(
  formula = formula_zinb,
  data = data,
  family = zero_inflated_negbinomial(),
  save_pars = save_pars(latent = TRUE, all = FALSE),
  prior = priors,
  chains = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  seed = 12345,
  cores = parallel::detectCores()
)

saveRDS(fit_zinb, "fit_zinb.rds", compress = "xz")
fit_zinb <- readRDS("fit_zinb.rds")

summary(fit_zinb)
plot(fit_zinb)
pp_check(fit_zinb)

ce <- conditional_effects(fit_zinb, effects = "sexage:relatedness")
plot(ce, points = FALSE)

conds <- data.frame(relatedness = c(0.00, 0.25, 0.50))
ce <- conditional_effects(fit_zinb, effects = "sexage", conditions = conds)
class(p)
p <- plot(ce, points = FALSE)
p[[1]] + 
  scale_y_log10() +
  hrbrthemes::theme_ipsum()

####################
# Comparing models #
####################

fit_nb <- readRDS("fit_nb.rds")
fit_zinb <- readRDS("fit_zinb.rds")

loo_nb  <- loo(fit_nb, reloo = FALSE)
loo_zinb <- loo(fit_zinb, reloo = FALSE)
loo_compare(loo_nb, loo_zinb)

##################
# Swarming index #
##################

individuals <- read.csv("index_data.csv")

individuals$index <- ((((individuals$events)/(individuals$cooccurrences5s+1))-
                         (sum(individuals$events)/sum(individuals$cooccurrences5s))))

groupwiseMedian(index ~ sexage, data = individuals, conf = 0.95, R = 1000, 
                percentile = TRUE, bca = FALSE, digits = 3)

kruskal.test(index ~ sexage, data = individuals)

DT = dunnTest(index ~ sexage, method = "bonferroni", data = individuals)
DT
PT = DT$res
cldList(P.adj ~ Comparison, data = PT, threshold = 0.05)

ggplot(individuals, aes(x=sexage, y=index)) +
  geom_jitter(aes(color = days), width = 0.2, size=1, alpha=0.5) +
  scale_colour_gradient2(low="#ffbfcaff", midpoint = 60, mid = "#ff0000ff", high="#8a0000ff") +
  theme_ipsum() + scale_y_continuous(breaks=c(-5,0,5,10,15,20)) +
  xlab("category") + ylab("swarming index") +
  scale_x_discrete() + stat_summary(fun=median, geom="point", size=2, color="black") +
  geom_segment(x = 1, xend = 1, y = 0.61, yend = 1.54) +
  geom_segment(x = 2, xend = 2, y = -0.55, yend = 0.96) +
  geom_segment(x = 3, xend = 3, y = -1.60, yend = -0.10) +
  coord_cartesian(ylim = c(-4, 11))


