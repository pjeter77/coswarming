library("readxl")
library("ggplot2")
library("patchwork")
library("hrbrthemes")
library("scales")
library("FSA")
library("rcompanion")
library("rstatix")
library("dplyr") 
library("STRAND")
library("network")
library("sna")

##################################################################
# Number of interactions of individuals                          #
##################################################################

individuals <- read.csv("indexdata.csv")

individuals$index <- ((((individuals$events)/(individuals$interactions5s+1))-
    (sum(individuals$events)/sum(individuals$interactions5s))))*2

identify_outliers(as.data.frame(individuals$index))

individuals <- subset(individuals, index < 22.8)

nrow(individuals)

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
  geom_segment(x = 1, xend = 1, y = 0.58, yend = 2.52) +
  geom_segment(x = 2, xend = 2, y = -1.20, yend = 1.66) +
  geom_segment(x = 3, xend = 3, y = -3.86, yend = -0.50) +
  coord_cartesian(ylim = c(-7, 22))

##################################################################
# Generative model of social network data                        #
##################################################################

swarming <- read.csv("swarming5.csv",header=T,row.names=1,stringsAsFactors=FALSE)
swarming <- as.matrix(swarming)
isSymmetric(swarming)

days <- read.csv("days.csv",header=T,row.names=1,stringsAsFactors=FALSE)
days <- as.matrix(days)
isSymmetric(days)

swarming = round(swarming/days*100,0)
swarming[is.nan(swarming)] = 0
isSymmetric(swarming)

noopportunity <- read.csv("noopportunity.csv",header=T,row.names=1,stringsAsFactors=FALSE)
noopportunity <- as.matrix(noopportunity)
isSymmetric(noopportunity)

relatedness <- read.csv("relatedness.csv",header=T,row.names=1,stringsAsFactors=FALSE)
relatedness <- as.matrix(relatedness)
isSymmetric(relatedness)

sexage <- read.csv("sexage.csv",header=F)
sexage <- sexage$V1

data <- list(swarming, noopportunity, relatedness, sexage)
names(data) <- c("swarming", "noopportunity", "relatedness", "sexage")

data

nets = list(swarm=data$swarming)
dyad = list(relatedness=data$relatedness, noopportunity=data$noopportunity)
groups = data.frame(sexage=as.factor(data$sexage))

dat = make_strand_data(
  outcome = nets,
  block_covariates = groups,
  individual_covariates = NULL,
  dyadic_covariates = dyad,
  outcome_mode = "poisson"
)

dat

fit =
  fit_block_plus_social_relations_model(
    data = dat,
    block_regression = ~ sexage,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ relatedness*noopportunity,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500,
      iter_sampling = 1500)
  )

res = summarize_strand_results(fit)

vis1 = strand_caterpillar_plot(res, normalized=T,  only_slopes=T) + 
  theme_ipsum()


vis2 = strand_caterpillar_plot(res, normalized=T, only_technicals = T, only_slopes=F) +
  theme_ipsum()

vis1 / vis2 + 
  plot_layout(heights = c(2, 0.75))
