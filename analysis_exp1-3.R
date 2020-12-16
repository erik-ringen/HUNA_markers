library(rethinking)
library(tidyverse)
library(ggridges)
library(patchwork)
library(wesanderson)

d <- read_csv("data_exp.csv")
d <- d[!is.na(d$y),] # drop rows where the individual did not get asked

###### Data Dictionary ######
### site = village 
## 0 = Bevondrorano
## 2 = Antaimbalabo
## 4 = Beangolo
## 6 = Ankililoake
## 7 = Andravitsazo
## 8 = Faramasay

### expcond = experiment num
## 1-4 indicates which experiment

### id = participant id

### sex = female (0), male (1)

### agecat1 = researcher-defined age category
## 0 = young adult
## 1 = adult
## 2 = old adult

### ethn1 = ethnicity (based on residence)
## 1 = Masikoro
## 2 = Mikea
## 3 = Vezo

### magic = person is an ambiasa, tromba, mpanompotromba, or mpitankazomanga
## 0 = none
## 1 = one or more

### christian
## 0 = no
## 1 = yes

### edu = education (years)
## 0 - 12 indicates years of formal education

### alter = id of individual in photo

### y = response to photo on given condition
## 1-3 responses indicate ethnicity guesses (exp1-3), 4 indicates trust (exp4), 5 indicates mistrust (exp4), 10 indicates unsure (exp1-3)

### ethn_alter = ethnicity of alter in photo
## 1 = Masikoro, 2 = Mikea, 3 = Vezo

#############################
##### Experiment One #######
d <- d %>% filter(expcond == 1) %>% 
  mutate(correct = ifelse(y == ethn_alter, 1, 0))

table(d$correct)
#############################
# NOTE: only 6 individuals who appear in both experiments 1 and 4
#############################

##### Prepare data ##########
# We need to make integer indices for grouping variables
village <- match(d$site, unique(d$site))
id <- match(d$id, unique(d$id))
alter_id <- match(d$alter, unique(d$alter))

# Not many older adults (19 individuals), so lets collapse into binary
d$age <- ifelse(d$agecat1 > 0, 1, 0)

#############################
N_obs <- nrow(d)
N_id <- max(id) 
N_eth <- max(d$ethn1)
N_village <- max(village)
N_alter <- max(alter_id)

#############################
data_list <- list(
  N_obs = N_obs,
  N_id = N_id,
  N_eth = N_eth,
  N_village = N_village,
  N_alter = N_alter,
  correct = d$correct,
  id = id,
  alter_id = alter_id,
  eth = d$ethn1,
  eth_alter = d$ethn_alter,
  village = village,
  edu = d$edu,
  sex = d$sex,
  sex_alter = d$sex_alter,
  age = d$age,
  magic = d$magic,
  christian = d$christian
)  

### Fit model in Stan ##################################
fit_1 <- stan( file="exp_model.stan", data=data_list, iter=5000, chains=4, cores=4, init="0", control = list(adapt_delta=0.96) )
saveRDS(fit_1, "fit_1.rds")

# After you've fit for the first time, un-comment the next line so that you can read it in without re-fitting
#fit_1 <- readRDS("fit_1.rds")

post <- extract.samples(fit_1) # Extract posterior samples
n_samps <- length(post$lp__)

##### ggridges plot for fixed effects ##################
b_long <- as.data.frame(post$b)
names(b_long) <- c("Education", "Male (judge)", "Male (alter)", "Male (judge) * Male (alter)", "Age", "Magic", "Christian")
b_long$samp <- 1:n_samps

b_long <- b_long %>% pivot_longer(-samp)

exp1_fixed <- ggplot(b_long, aes(x=value, y=fct_reorder(name,abs(value)))) +
  geom_hline(yintercept = c(1:7), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_vline(aes(xintercept=0), lty="dashed") + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab("Log Odds Ratio") + 
  annotate("blank", x = 0, y=8)

##### R^2 for random effects and fixed effects #########
mo <- function(scale, edu) {
  if (edu > 0) {
    b_edu <- c()
    for (i in 1:n_samps) b_edu[i] <- sum(scale[i,1:edu])
    return(b_edu)
  }
  else return(0);
}

id_mu <- matrix( NA, nrow=n_samps, ncol=N_obs )
eth_mu <- id_mu; village_mu <- id_mu; alter_mu <- id_mu; f_mu <- id_mu

for (n in 1:N_obs) {
  id_mu[,n] = post$id_v[,data_list$id[n],1] + post$id_v[,data_list$id[n],data_list$eth_alter[n] + 1] 
  eth_mu[,n] = post$eth_v[,data_list$eth[n],1] + post$eth_v[,data_list$eth_alter[n],2] + post$eth_v[,data_list$eth[n],data_list$eth_alter[n] + 2]
    village_mu[,n] = post$village_v[,data_list$village[n],1] + post$village_v[,data_list$village[n],data_list$eth_alter[n] + 1] 
  alter_mu[,n] = post$alter_v[,data_list$alter_id[n]]
  f_mu[,n] = post$b[,1]*mo(post$scale_edu, data_list$edu[n]) + post$b[,2]*data_list$sex[n] + post$b[,3]*data_list$sex_alter[n] + post$b[,4]*data_list$sex_alter[n] + post$b[,5]*data_list$age[n] + post$b[,6]*data_list$magic[n] + post$b[,7]*data_list$christian[n]
}

id_var <- apply(id_mu, 1, var)
eth_var <- apply(eth_mu, 1, var)
village_var <- apply(village_mu, 1, var)
alter_var <- apply(alter_mu, 1, var)
fixed_var <- apply(f_mu, 1, var)

# pi^2/3 is the constant variance induced by the logit link function
tot_var <- id_var + eth_var + village_var + alter_var + fixed_var + pi^2/3
# Organize variances
var_df <- data.frame(
  est = c(id_var, eth_var, village_var, alter_var, fixed_var)/tot_var,
  var = rep(c("ID (judge)", "Ethnicity", "Village", "ID (alter)", "Fixed Effects"), each=n_samps)
)

exp1_r2 <- ggplot(var_df, aes(x=est, y=fct_reorder(var,est))) +
  geom_hline(yintercept = c(1:5), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05), limits=c(0,1)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab(expression(R^2)) + 
  ggtitle("Experiment 1:\nEverday Photos") +
  annotate("blank", x = 0, y=6)

##### Raw data for comparison ##########################
d$eth <- factor(d$ethn1, labels=c("Masikoro", "Mikea", "Vezo"))

# Calculate how often each individual correctly identified
exp1_summary <- d %>% group_by(id) %>%
  summarise(
    prob_correct = mean(correct),
    eth = unique(eth)
  )

exp1_ethn <- d %>% group_by(eth) %>%
  summarise(
    eth_prob = mean(correct)
  )

exp1_ethn <- bind_rows(exp1_ethn, data.frame(eth="Random Guessing", eth_prob=1/3))
exp1_ethn$eth <- fct_relevel(exp1_ethn$eth, "Random Guessing", after=Inf) # make sure guessing is the last factor level

# Unify the factor levels in both dataframes
exp1_ethn$eth <-  fct_unify( list(exp1_ethn$eth, exp1_summary$eth) )[[1]]
exp1_summary$eth <-  fct_unify( list(exp1_ethn$eth, exp1_summary$eth) )[[2]]

# Plot distributions for each ethnicity
exp1_overall <- ggplot(exp1_summary, aes(x=prob_correct, fill=eth)) + geom_vline(data=exp1_ethn, aes(xintercept=eth_prob, color=eth, linetype=eth),lwd=1) + geom_density(data=exp1_summary, alpha=0.7) +
  scale_y_continuous(expand=c(0,0),limits=c(0,5.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1), breaks=c(0, 1/3, 2/3, 1), labels=c("0", "1/3", "2/3", "1")) +
  scale_color_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_fill_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_linetype_manual(values=c( rep("solid",3), "dashed"),drop=F) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank()) +
  xlab("Pr(Correct Classification)") +
  ylab("Density") + 
  ggtitle("Experiment 1: Everday Photos")


# Now, break up into ethnicity * ethnicity
exp1_summary2 <- d %>% group_by(id, ethn_alter) %>%
  summarise(
    prob_correct = sum(correct)/n(),
    eth = unique(eth)
  )

exp1_ethn2 <- d %>% group_by(eth, ethn_alter) %>%
  summarise(
    eth_prob = mean(correct)
  )

exp1_ethn2 <- bind_rows(exp1_ethn2, data.frame(eth=rep("Random Guessing",3), eth_prob=rep(1/3,3), ethn_alter=1:3))
exp1_ethn2$eth <- fct_relevel(exp1_ethn2$eth, "Random Guessing", after=Inf) # make sure guessing is the last factor level

# Unify the factor levels in both dataframes
exp1_ethn2$eth <-  fct_unify( list(exp1_ethn2$eth, exp1_summary2$eth) )[[1]]
exp1_summary2$eth <-  fct_unify( list(exp1_ethn2$eth, exp1_summary2$eth) )[[2]]

# Label the alter ethnicity facets
exp1_summary2$ethn_alter <- factor(exp1_summary2$ethn_alter, labels=c("Alter = Masikoro", "Alter = Mikea", "Alter = Vezo"))
exp1_ethn2$ethn_alter <- factor(exp1_ethn2$ethn_alter, labels=c("Alter = Masikoro", "Alter = Mikea", "Alter = Vezo"))


exp1_alter <- ggplot(exp1_summary2, aes(x=prob_correct, fill=eth)) + facet_wrap(~ethn_alter) + geom_vline(data=exp1_ethn2, aes(xintercept=eth_prob, color=eth, linetype=eth),lwd=1) + geom_density(data=exp1_summary2, alpha=0.7) +
  scale_y_continuous(expand=c(0,0), limits=c(0,4.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1),breaks=c(0, 1/3, 2/3, 1), labels=c("0", "1/3", "2/3", "1")) +
  scale_color_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_fill_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_linetype_manual(values=c( rep("solid",3), "dashed")) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(),panel.spacing = unit(1.5, "lines")) +
  xlab("Pr(Correct Classification)") +
  ylab("Density") +
  ggtitle("Experiment 1: Everday Photos")

svg( file="exp1_descriptive.svg", width=7, height=5, pointsize=12 )
exp1_overall / exp1_alter + plot_layout(guides="collect")
dev.off()

########################################################
########################################################
##### Experiment Two ###################################
d <- read_csv("data_exp.csv")
d <- d[!is.na(d$y),] # drop rows where the individual did not get asked

d <- d %>% filter(expcond == 2) %>% 
  mutate(correct = ifelse(y == ethn_alter, 1, 0))

table(d$correct)

##### Prepare data ##########
# We need to make integer indices for grouping variables
village <- match(d$site, unique(d$site))
id <- match(d$id, unique(d$id))
alter_id <- match(d$alter, unique(d$alter))

# Not many older adults (19 individuals), so lets collapse into binary
age <- ifelse(d$agecat1 > 0, 1, 0)

#############################
N_obs <- nrow(d)
N_id <- max(id) 
N_eth <- max(d$ethn1)
N_village <- max(village)
N_alter <- max(alter_id)

#############################
data_list <- list(
  N_obs = N_obs,
  N_id = N_id,
  N_eth = N_eth,
  N_village = N_village,
  N_alter = N_alter,
  correct = d$correct,
  id = id,
  alter_id = alter_id,
  eth = d$ethn1,
  eth_alter = d$ethn_alter,
  village = village,
  edu = d$edu,
  sex = d$sex,
  sex_alter = d$sex_alter,
  age = age,
  magic = d$magic,
  christian = d$christian
)  

### Fit model in Stan ##################################
fit_2 <- stan( file="exp_model.stan", data=data_list, iter=5000, chains=4, cores=4, init="0", control = list(adapt_delta=0.96) )
saveRDS(fit_2, "fit_2.rds")

# After you've fit for the first time, un-comment the next line so that you can read it in without re-fitting
#fit_2 <- readRDS("fit_2.rds")

post <- extract.samples(fit_2) # Extract posterior samples
n_samps <- length(post$lp__)

##### ggridges plot for fixed effects ##################
b_long <- as.data.frame(post$b)
names(b_long) <- c("Education", "Male (judge)", "Male (alter)", "Male (judge) * Male (alter)", "Age", "Magic", "Christian")
b_long$samp <- 1:n_samps

b_long <- b_long %>% pivot_longer(-samp)

exp2_fixed <- ggplot(b_long, aes(x=value, y=fct_reorder(name,abs(value)))) +
  geom_hline(yintercept = c(1:7), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_vline(aes(xintercept=0), lty="dashed") + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab("Log Odds Ratio") + 
  annotate("blank", x = 0, y=8)

##### R^2 for random effects and fixed effects #########
mo <- function(scale, edu) {
  if (edu > 0) {
    b_edu <- c()
    for (i in 1:n_samps) b_edu[i] <- sum(scale[i,1:edu])
    return(b_edu)
  }
  else return(0);
}

id_mu <- matrix( NA, nrow=n_samps, ncol=N_obs )
eth_mu <- id_mu; village_mu <- id_mu; alter_mu <- id_mu; f_mu <- id_mu

for (n in 1:N_obs) {
  id_mu[,n] = post$id_v[,data_list$id[n],1] + post$id_v[,data_list$id[n],data_list$eth_alter[n] + 1] 
  eth_mu[,n] = post$eth_v[,data_list$eth[n],1] + post$eth_v[,data_list$eth_alter[n],2] + post$eth_v[,data_list$eth[n],data_list$eth_alter[n] + 2]
    village_mu[,n] = post$village_v[,data_list$village[n],1] + post$village_v[,data_list$village[n],data_list$eth_alter[n] + 1] 
  alter_mu[,n] = post$alter_v[,data_list$alter_id[n]]
  f_mu[,n] = post$b[,1]*mo(post$scale_edu, data_list$edu[n]) + post$b[,2]*data_list$sex[n] + post$b[,3]*data_list$sex_alter[n] + post$b[,4]*data_list$sex_alter[n] + post$b[,5]*data_list$age[n] + post$b[,6]*data_list$magic[n] + post$b[,7]*data_list$christian[n]
}

id_var <- apply(id_mu, 1, var)
eth_var <- apply(eth_mu, 1, var)
village_var <- apply(village_mu, 1, var)
alter_var <- apply(alter_mu, 1, var)
fixed_var <- apply(f_mu, 1, var)

tot_var <- id_var + eth_var + village_var + alter_var + fixed_var + pi^2/3
# Organize variances
var_df <- data.frame(
  est = c(id_var, eth_var, village_var, alter_var, fixed_var)/tot_var,
  var = rep(c("ID (judge)", "Ethnicity", "Village", "ID (alter)", "Fixed Effects"), each=n_samps)
)

exp2_r2 <- ggplot(var_df, aes(x=est, y=fct_reorder(var,est))) +
  geom_hline(yintercept = c(1:5), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05), limits=c(0,1)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  ggtitle("Experiment 2:\nPurposeful Marking") +
  xlab(expression(R^2)) + 
  annotate("blank", x = 0, y=6)

svg(file="exp2_effects.svg", width=7, height=5, pointsize=12)

exp2_fixed + exp2_r2

dev.off()

##### Raw data for comparison ##########################
d$eth <- factor(d$ethn1, labels=c("Masikoro", "Mikea", "Vezo"))

# Calculate how often each individual correctly identified
exp2_summary <- d %>% group_by(id) %>%
  summarise(
    prob_correct = mean(correct),
    eth = unique(eth)
  )

exp2_ethn <- d %>% group_by(eth) %>%
  summarise(
    eth_prob = mean(correct)
  )

exp2_ethn <- bind_rows(exp2_ethn, data.frame(eth="Random Guessing", eth_prob=1/3))
exp2_ethn$eth <- fct_relevel(exp2_ethn$eth, "Random Guessing", after=Inf) # make sure guessing is the last factor level

# Unify the factor levels in both dataframes
exp2_ethn$eth <-  fct_unify( list(exp2_ethn$eth, exp2_summary$eth) )[[1]]
exp2_summary$eth <-  fct_unify( list(exp2_ethn$eth, exp2_summary$eth) )[[2]]

# Plot distributions for each ethnicity
exp2_overall <- ggplot(exp2_summary, aes(x=prob_correct, fill=eth)) + geom_vline(data=exp2_ethn, aes(xintercept=eth_prob, color=eth, linetype=eth),lwd=1) + geom_density(data=exp2_summary, alpha=0.7) +
  scale_y_continuous(expand=c(0,0),limits=c(0,5.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1), breaks=c(0, 1/3, 2/3, 1), labels=c("0", "1/3", "2/3", "1")) +
  scale_color_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_fill_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_linetype_manual(values=c( rep("solid",3), "dashed"),drop=F) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank()) +
  xlab("Pr(Correct Classification)") +
  ylab("Density") + 
  ggtitle("Experiment 2: Purposefully Mark")


# Now, break up into ethnicity * ethnicity
exp2_summary2 <- d %>% group_by(id, ethn_alter) %>%
  summarise(
    prob_correct = sum(correct)/n(),
    eth = unique(eth)
  )

exp2_ethn2 <- d %>% group_by(eth, ethn_alter) %>%
  summarise(
    eth_prob = mean(correct)
  )

exp2_ethn2 <- bind_rows(exp2_ethn2, data.frame(eth=rep("Random Guessing",3), eth_prob=rep(1/3,3), ethn_alter=1:3))
exp2_ethn2$eth <- fct_relevel(exp2_ethn2$eth, "Random Guessing", after=Inf) # make sure guessing is the last factor level

# Unify the factor levels in both dataframes
exp2_ethn2$eth <-  fct_unify( list(exp2_ethn2$eth, exp2_summary2$eth) )[[1]]
exp2_summary2$eth <-  fct_unify( list(exp2_ethn2$eth, exp2_summary2$eth) )[[2]]

# Label the alter ethnicity facets
exp2_summary2$ethn_alter <- factor(exp2_summary2$ethn_alter, labels=c("Alter = Masikoro", "Alter = Mikea", "Alter = Vezo"))
exp2_ethn2$ethn_alter <- factor(exp2_ethn2$ethn_alter, labels=c("Alter = Masikoro", "Alter = Mikea", "Alter = Vezo"))


exp2_alter <- ggplot(exp2_summary2, aes(x=prob_correct, fill=eth)) + facet_wrap(~ethn_alter) + geom_vline(data=exp2_ethn2, aes(xintercept=eth_prob, color=eth, linetype=eth),lwd=1) + geom_density(data=exp2_summary2, alpha=0.7) +
  scale_y_continuous(expand=c(0,0), limits=c(0,4.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1),breaks=c(0, 1/3, 2/3, 1), labels=c("0", "1/3", "2/3", "1")) +
  scale_color_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_fill_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_linetype_manual(values=c( rep("solid",3), "dashed")) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(),panel.spacing = unit(1.5, "lines")) +
  xlab("Pr(Correct Classification)") +
  ylab("Density") +
  ggtitle("Experiment 2: Purposeful Marking")

svg( file="exp2_descriptive.svg", width=7, height=5, pointsize=12 )
exp2_overall / exp2_alter + plot_layout(guides="collect")
dev.off()


########################################################
########################################################
##### Experiment Three #################################
d <- read_csv("data_exp.csv")
d <- d[!is.na(d$y),] # drop rows where the individual did not get asked

d <- d %>% filter(expcond == 3) %>% 
  mutate(correct = ifelse(y == ethn_alter, 1, 0))

table(d$correct)

##### Prepare data ##########
# We need to make integer indices for grouping variables
village <- match(d$site, unique(d$site))
id <- match(d$id, unique(d$id))
alter_id <- match(d$alter, unique(d$alter))

# Not many older adults (19 individuals), so lets collapse into binary
age <- ifelse(d$agecat1 > 0, 1, 0)

#############################
N_obs <- nrow(d)
N_id <- max(id) 
N_eth <- max(d$ethn1)
N_village <- max(village)
N_alter <- max(alter_id)

#############################
data_list <- list(
  N_obs = N_obs,
  N_id = N_id,
  N_eth = N_eth,
  N_village = N_village,
  N_alter = N_alter,
  correct = d$correct,
  id = id,
  alter_id = alter_id,
  eth = d$ethn1,
  eth_alter = d$ethn_alter,
  village = village,
  edu = d$edu,
  sex = d$sex,
  sex_alter = d$sex_alter,
  age = age,
  magic = d$magic,
  christian = d$christian
)  

### Fit model in Stan ##################################
fit_3 <- stan( file="exp_model.stan", data=data_list, iter=5000, chains=4, cores=4, init="0", control = list(adapt_delta=0.96) )
saveRDS(fit_3, "fit_3.rds")

# After you've fit for the first time, un-comment the next line so that you can read it in without re-fitting
#fit_3 <- readRDS("fit_3.rds")

post <- extract.samples(fit_3) # Extract posterior samples
n_samps <- length(post$lp__)

##### ggridges plot for fixed effects ##################
b_long <- as.data.frame(post$b)
names(b_long) <- c("Education", "Male (judge)", "Male (alter)", "Male (judge) * Male (alter)", "Age", "Magic", "Christian")
b_long$samp <- 1:n_samps

b_long <- b_long %>% pivot_longer(-samp)

exp3_fixed <- ggplot(b_long, aes(x=value, y=fct_reorder(name,abs(value)))) +
  geom_hline(yintercept = c(1:7), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_vline(aes(xintercept=0), lty="dashed") + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab("Log Odds Ratio") + 
  annotate("blank", x = 0, y=8)

##### R^2 for random effects and fixed effects #########
mo <- function(scale, edu) {
  if (edu > 0) {
    b_edu <- c()
    for (i in 1:n_samps) b_edu[i] <- sum(scale[i,1:edu])
    return(b_edu)
  }
  else return(0);
}

id_mu <- matrix( NA, nrow=n_samps, ncol=N_obs )
eth_mu <- id_mu; village_mu <- id_mu; alter_mu <- id_mu; f_mu <- id_mu

for (n in 1:N_obs) {
  id_mu[,n] = post$id_v[,data_list$id[n],1] + post$id_v[,data_list$id[n],data_list$eth_alter[n] + 1] 
  eth_mu[,n] = post$eth_v[,data_list$eth[n],1] + post$eth_v[,data_list$eth_alter[n],2] + post$eth_v[,data_list$eth[n],data_list$eth_alter[n] + 2]
  village_mu[,n] = post$village_v[,data_list$village[n],1] + post$village_v[,data_list$village[n],data_list$eth_alter[n] + 1] 
  alter_mu[,n] = post$alter_v[,data_list$alter_id[n]]
  f_mu[,n] = post$b[,1]*mo(post$scale_edu, data_list$edu[n]) + post$b[,2]*data_list$sex[n] + post$b[,3]*data_list$sex_alter[n] + post$b[,4]*data_list$sex_alter[n] + post$b[,5]*data_list$age[n] + post$b[,6]*data_list$magic[n] + post$b[,7]*data_list$christian[n]
}

id_var <- apply(id_mu, 1, var)
eth_var <- apply(eth_mu, 1, var)
village_var <- apply(village_mu, 1, var)
alter_var <- apply(alter_mu, 1, var)
fixed_var <- apply(f_mu, 1, var)

tot_var <- id_var + eth_var + village_var + alter_var + fixed_var + pi^2/3
# Organize variances
var_df <- data.frame(
  est = c(id_var, eth_var, village_var, alter_var, fixed_var)/tot_var,
  var = rep(c("ID (judge)", "Ethnicity", "Village", "ID (alter)", "Fixed Effects"), each=n_samps)
)

exp3_r2 <- ggplot(var_df, aes(x=est, y=fct_reorder(var,est))) +
  geom_hline(yintercept = c(1:5), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05), limits=c(0,1)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab(expression(R^2)) + 
  ggtitle("Experiment 3:\nMarket") +
  annotate("blank", x = 0, y=6)

svg(file="exp3_effects.svg", width=7, height=5, pointsize=12)

exp3_fixed + exp3_r2

dev.off()

#### Now, raw data ###########################################
d$eth <- factor(d$ethn1, labels=c("Masikoro", "Mikea", "Vezo"))

# Calculate how often each individual correctly identified
exp3_summary <- d %>% group_by(id) %>%
  summarise(
    prob_correct = mean(correct),
    eth = unique(eth)
  )

exp3_ethn <- d %>% group_by(eth) %>%
  summarise(
    eth_prob = mean(correct)
  )

exp3_ethn <- bind_rows(exp3_ethn, data.frame(eth="Random Guessing", eth_prob=1/2))
exp3_ethn$eth <- fct_relevel(exp3_ethn$eth, "Random Guessing", after=Inf) # make sure guessing is the last factor level

# Unify the factor levels in both dataframes
exp3_ethn$eth <-  fct_unify( list(exp3_ethn$eth, exp3_summary$eth) )[[1]]
exp3_summary$eth <-  fct_unify( list(exp3_ethn$eth, exp3_summary$eth) )[[2]]

# Plot distributions for each ethnicity
exp3_overall <- ggplot(exp3_summary, aes(x=prob_correct, fill=eth)) + geom_vline(data=exp3_ethn, aes(xintercept=eth_prob, color=eth, linetype=eth),lwd=1) + geom_density(data=exp3_summary, alpha=0.7) +
  scale_y_continuous(expand=c(0,0),limits=c(0,3.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1), breaks=c(0, 1/2, 1), labels=c("0", "1/2", "1")) +
  scale_color_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_fill_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_linetype_manual(values=c( rep("solid",3), "dashed"),drop=F) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank()) +
  xlab("Pr(Correct Classification)") +
  ylab("Density") + 
  ggtitle("Experiment 3: Market")


# Now, break up into ethnicity * ethnicity
exp3_summary2 <- d %>% group_by(id, ethn_alter) %>%
  summarise(
    prob_correct = sum(correct)/n(),
    eth = unique(eth)
  )

exp3_ethn2 <- d %>% group_by(eth, ethn_alter) %>%
  summarise(
    eth_prob = mean(correct)
  )

exp3_ethn2 <- bind_rows(exp3_ethn2, data.frame(eth=rep("Random Guessing",3), eth_prob=rep(1/2,3), ethn_alter=1:3))
exp3_ethn2$eth <- fct_relevel(exp3_ethn2$eth, "Random Guessing", after=Inf) # make sure guessing is the last factor level

# Unify the factor levels in both dataframes
exp3_ethn2$eth <-  fct_unify( list(exp3_ethn2$eth, exp3_summary2$eth) )[[1]]
exp3_summary2$eth <-  fct_unify( list(exp3_ethn2$eth, exp3_summary2$eth) )[[2]]

# Label the alter ethnicity facets
exp3_summary2$ethn_alter <- factor(exp3_summary2$ethn_alter, labels=c("Alter = Masikoro", "Alter = Mikea", "Alter = Vezo"))
exp3_ethn2$ethn_alter <- factor(exp3_ethn2$ethn_alter, labels=c("Alter = Masikoro", "Alter = Mikea", "Alter = Vezo"))


exp3_alter <- ggplot(exp3_summary2, aes(x=prob_correct, fill=eth)) + facet_wrap(~ethn_alter) + geom_vline(data=exp3_ethn2, aes(xintercept=eth_prob, color=eth, linetype=eth),lwd=1) + geom_density(data=exp3_summary2, alpha=0.7) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1),breaks=c(0, 1/2, 1), labels=c("0", "1/2", "1")) +
  scale_color_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_fill_manual(values=c("goldenrod", "seagreen","#00AFBB", "darkgrey"), drop=F) +
  scale_linetype_manual(values=c( rep("solid",3), "dashed")) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(),panel.spacing = unit(1.5, "lines")) +
  xlab("Pr(Correct Classification)") +
  ylab("Density") + 
  ggtitle("Experiment 3: Market")

svg( file="exp3_descriptive.svg", width=7, height=5, pointsize=12 )
exp3_overall / exp3_alter + plot_layout(guides="collect")
dev.off()

##################################################################
###### Hypothesis testing ########################################

#### (1) In experiment 1, do participants classify subjects by ethnicity at better than chance rate?
fit_1 <- readRDS("fit_1.rds")
post_1 <- extract.samples(fit_1)
n_samps <- length(post_1$lp__)

# Averaging over id, ethnicity, alter, assuming no education (the modal case), averaging over sex differences in ego and alter, averaging over age, assuming not magic and not christian (the modal case).
exp1_odds <- logistic( post_1$b0 + (post_1$b[,2] + post_1$b[,3] + post_1$b[,4] + post_1$b[,5])/2 ) / (1/3)
h1 <- length(exp1_odds[exp1_odds > 1])/length(exp1_odds) # posterior probability that performance is better than chance


#### (2) In experiment 2, do participants classify subjects by ethnicity at better than chance rate?
fit_2 <- readRDS("fit_2.rds")
post_2 <- extract.samples(fit_2)

# Averaging over id, ethnicity, alter, assuming no education (the modal case), averaging over sex differences in ego and alter, averaging over age, assuming not magic and not christian (the modal case).
exp2_odds <- logistic( post_2$b0 + (post_2$b[,2] + post_2$b[,3] + post_2$b[,4] + post_2$b[,5])/2 ) / (1/3)
h2 <- length(exp2_odds[exp2_odds > 1])/length(exp2_odds) # posterior probability that performance is better than chance


#### (3) Do the judges in experiment 2 have better success rates than experiment 1?
odds2_odds1 <- exp2_odds / exp1_odds
h3 <- length(odds2_odds1[odds2_odds1 > 1])/length(odds2_odds1)


#### (4) In experiment 3, do participants classify subjects by ethnicity at better than chance rate?
fit_3 <- readRDS("fit_3.rds")
post_3 <- extract.samples(fit_3)

# Averaging over id, ethnicity, alter, assuming no education (the modal case), averaging over sex differences in ego and alter, averaging over age, assuming not magic and not christian (the modal case).
exp3_odds <- logistic( post_3$b0 + (post_3$b[,2] + post_3$b[,3] + post_3$b[,4] + post_3$b[,5])/2 ) / (1/2)
h4 <- length(exp3_odds[exp3_odds > 1])/length(exp3_odds) # posterior probability that performance is better than chance


#### (5) Do the judges in experiment 3 have better success rates than experiment 2?
odds3_odds2 <- exp3_odds / exp2_odds
h5 <- length(odds3_odds2[odds3_odds2 > 1])/length(odds3_odds2)


##### Bring it all together to plot ####
odds_df <- data.frame(
  est = c(exp1_odds, exp2_odds, odds2_odds1, exp3_odds, odds3_odds2),
  hypothesis = rep(c("Experiment 1", "Experiment 2", "2/1", "Experiment 3", "3/2"), each=n_samps)
)

odds_df_long <- odds_df %>% pivot_longer(-hypothesis)
odds_df_long$type <- ifelse(odds_df_long$hypothesis %in% c("2/1", "3/2"), "CR Ratio", "Classification Ratio (Pr Correct:Random Guessing)")


cols <- wes_palette("BottleRocket2")
cols[1] <- "tomato4"
cols[2] <- "deeppink4"

svg(file="exp1-3_hypotheses.svg", width=8, height=5, pointsize=12)

exp123_hyp <- ggplot(filter(odds_df_long), aes(x=value, group=hypothesis, fill=hypothesis)) +
  facet_wrap(~type, scales="free_y") + 
  geom_density(alpha=0.8) +
  geom_vline(aes(xintercept=1), linetype="dashed") +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.75)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,4)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(),panel.spacing = unit(1.5, "lines"), legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) + xlab("") + ylab("")

dev.off()

####################################################################
##### Model estimates for each experiment ##########################
svg(file="exp1-3_effects.svg", width=12, height=7, pointsize=12)

(exp1_r2 + exp2_r2 + exp3_r2) / (exp1_fixed + exp2_fixed + exp3_fixed)

dev.off()
####################################################################
##### Save all three experiments raw data ##########################
svg(file="exp1-3_raw.svg", width=6, height=7, pointsize=12)

exp1_alter / exp2_alter / exp3_alter + plot_layout(guides="collect")

dev.off()

