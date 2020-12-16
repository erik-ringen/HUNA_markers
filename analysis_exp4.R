library(rethinking)
library(tidyverse)
library(igraph)
library(ggridges)
library(patchwork)

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
## 1-3 responses indicate ethnicity guesses (exp1-3), 4 indicates trust (exp4), 5 indicates mistrust (exp4), 10 indicates unsure

### ethn_alter = ethnicity of alter in photo
## 1 = Masikoro, 2 = Mikea, 3 = Vezo

#############################
##### Experiment Four #######
d <- d %>% filter(d$expcond == 4)

##### Now, subset to just the co-ethnic data ###########
d <- d %>% filter(ethn1 == ethn_alter)

##### Prepare data ##########
# We need to make integer indices for grouping variables
village <- match(d$site, unique(d$site))
id <- match(d$id, unique(d$id))
alter_id <- match(d$alter, unique(d$alter))

# Not many older adults (19 individuals), so lets collapse into binary
d$age <- ifelse(d$agecat1 > 0, 1, 0)

# Set NA to -99 for Stan, we'll impute within the model. We're missing some education responses and some age responses (and sometimes both)
d[is.na(d)] <- -99

# Recode trust response
d$trust <- ifelse( d$y == 4, 1, NA ) # 1 indicates trust
d$trust <- ifelse( d$y == 5, 2, d$trust ) # 2 indicates mistrust
d$trust <- ifelse( d$y == 10, 3, d$trust ) # 3 indicates neither trust nor mistrust

#############################
N_obs <- nrow(d)
N_id <- max(id) 
N_eth <- max(d$ethn1)
N_village <- max(village)
N_alter <- max(alter_id)

d_id <- d %>% group_by(id) %>% summarise(age_id = median(age), edu_id = median(edu))

#############################
data_list <- list(
  N_obs = N_obs,
  N_id = N_id,
  N_eth = N_eth,
  N_village = N_village,
  N_alter = N_alter,
  trust = d$trust,
  id = id,
  alter_id = alter_id,
  eth = d$ethn1,
  village = village,
  edu = d$edu,
  sex = d$sex,
  sex_alter = d$sex_alter,
  age = d$age,
  magic = d$magic,
  christian = d$christian,
  age_id = d_id$age_id,
  edu_id = d_id$edu_id
)  

### Fit model in Stan ##################################
fit <- stan( file="trust_model.stan", data=data_list, iter=4000, chains=4, cores=4, init="0", control = list(adapt_delta=0.96) )
saveRDS(fit, "fit_trust.rds")

# After you've fit for the first time, un-comment the next line so that you can read it in without re-fitting
#fit <- readRDS("fit_trust.rds")

post <- extract.samples(fit) # Extract posterior samples
n_samps <- length(post$lp__)

########################################################
##### Plot fixed and random effects ####################
## Trust effects
b_trust <- as.data.frame(cbind(post$b0[,1], post$b[,,1]))
names(b_trust) <- c("Intercept", "Education", "Male (judge)", "Male (alter)", "Male (judge) * Male (alter)", "Age", "Magic", "Christian")
b_trust$samp <- 1:n_samps

b_trust <- b_trust %>% pivot_longer(-samp)

trust_fixed <- ggplot(b_trust, aes(x=value, y=fct_reorder(name,abs(value)))) +
  geom_hline(yintercept = c(1:8), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_vline(aes(xintercept=0), lty="dashed") + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab(expression(beta)) + 
  ggtitle("Trust") + 
  annotate("blank", x = 0, y=9)

# random effects
sd_trust <- as.data.frame(cbind( post$sigma_id[,1], post$sigma_alter[,1], post$sigma_eth[,1], post$sigma_village[,1] ))
names(sd_trust) <- c("ID (judge)", "ID (alter)", "Ethnicity", "Village")
sd_trust$samp <- 1:n_samps

sd_trust <- sd_trust %>% pivot_longer(-samp)

trust_re <- ggplot(sd_trust, aes(x=value, y=fct_reorder(name,abs(value)))) +
  geom_hline(yintercept = c(1:4), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab(expression(sigma)) + 
  ggtitle("Trust") + 
  annotate("blank", x = 0, y=5)


## Now, mistrust
b_mistrust <- as.data.frame(cbind(post$b0[,2], post$b[,,2]))
names(b_mistrust) <- c("Intercept", "Education", "Male (judge)", "Male (alter)", "Male (judge) * Male (alter)", "Age", "Magic", "Christian")
b_mistrust$samp <- 1:n_samps

b_mistrust <- b_mistrust %>% pivot_longer(-samp)

mistrust_fixed <- ggplot(b_mistrust, aes(x=value, y=fct_reorder(name,abs(value)))) +
  geom_hline(yintercept = c(1:8), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_vline(aes(xintercept=0), lty="dashed") + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab(expression(beta)) + 
  ggtitle("Mistrust") + 
  annotate("blank", x = 0, y=9)

# random effects
sd_mistrust <- as.data.frame(cbind( post$sigma_id[,2], post$sigma_alter[,2], post$sigma_eth[,2], post$sigma_village[,2] ))
names(sd_mistrust) <- c("ID (judge)", "ID (alter)", "Ethnicity", "Village")
sd_mistrust$samp <- 1:n_samps

sd_mistrust <- sd_mistrust %>% pivot_longer(-samp)

mistrust_re <- ggplot(sd_mistrust, aes(x=value, y=fct_reorder(name,abs(value)))) +
  geom_hline(yintercept = c(1:4), alpha=0.6) + 
  geom_density_ridges(scale=0.9) + 
  scale_x_continuous(expand = c(0.00, 0.05)) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank()) + ylab("") + 
  xlab(expression(sigma)) + 
  ggtitle("Mistrust") + 
  annotate("blank", x = 0, y=5)

## Put em together
#svg( "trust_effects.svg", height=8, width=8, pointsize=12 )
#(trust_fixed + trust_re) / (mistrust_fixed + mistrust_re)
#dev.off()

##################################################################
##### Test hypothesis: preferential coethnic trust ###############

##### Softmax functions #####
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax2 <- function (x) {
  exp(x - logsumexp(x))
}

# Prob trust depending on sex interaction. Averaging over age (0.5)
f_f <- cbind( post$b0[,1] + post$b[,5,1]*0.5, post$b0[,2] + post$b[,5,2]*0.5, 0 )
f_m <- cbind( post$b0[,1] + post$b[,3,1] + post$b[,5,1]*0.5, post$b0[,2] + post$b[,3,2] + post$b[,5,2]*0.5, 0 )
m_f <- cbind( post$b0[,1] + post$b[,2,1] + post$b[,5,1]*0.5, post$b0[,2] + post$b[,2,2] + post$b[,5,2]*0.5, 0 )
m_m <- cbind( post$b0[,1] + post$b[,2,1] + post$b[,3,1] + post$b[,4,1] + post$b[,5,1]*0.5, post$b0[,2] + post$b[,2,2] + post$b[,3,2] + post$b[,4,2] + post$b[,5,2]*0.5, 0 )

# Converting to prob scale
for (i in 1:nrow(f_f)) {
  f_f[i,] <- softmax2(f_f[i,])
  f_m[i,] <- softmax2(f_m[i,])
  m_f[i,] <- softmax2(m_f[i,])
  m_m[i,] <- softmax2(m_m[i,])
}

round( median(f_f[,1]), 3 )
round( HPDI(f_f[,1], prob=0.9), 3 )

round( median(f_m[,1]), 3 )
round( HPDI(f_m[,1], prob=0.9), 3 )


### Organize predictions for plotting
pred_df <- as.data.frame(rbind(f_f, f_m, m_f, m_m))
names(pred_df) <- c("Trust", "Mistrust", "Neutral")
pred_df$diff <- pred_df$Trust - pred_df$Mistrust

pred_df$condition <- rep( c("Female sorting Female", "Female sorting Male", "Male sorting Female", "Male sorting Male"), each=n_samps )

pred_long <- pred_df %>% pivot_longer(-condition)

pred_long$name <- factor(pred_long$name, labels=c("Delta(Trust - Mistrust)", "Pr(Mistrust Coethnic)", "Neutral", "Pr(Trust Coethnic)"))
pred_long <- filter(pred_long, name != "Neutral")

cols <- c("darkgrey", "darkred", "cornflowerblue")

#### Evaluate overall difference
round(median(pred_long$value[pred_long$name == "Delta(Trust - Mistrust)"]),2)
# PP
round( length( pred_long$value[pred_long$name == "Delta(Trust - Mistrust)"][ pred_long$value[pred_long$name == "Delta(Trust - Mistrust)"] > 0])/length(pred_long$value[pred_long$name == "Delta(Trust - Mistrust)"]), 2 )

## Condition specific PP
# FF
round( length( pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Female sorting Female"][ pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Female sorting Female"] > 0])/length(pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Female sorting Female"]), 2 )

# FM
round( length( pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Female sorting Male"][ pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Female sorting Male"] > 0])/length(pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Female sorting Male"]), 2 )

# MF
round( length( pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Male sorting Female"][ pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Male sorting Female"] > 0])/length(pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Male sorting Female"]), 2 )

# MM
round( length( pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Male sorting Male"][ pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Male sorting Male"] > 0])/length(pred_long$value[pred_long$name == "Delta(Trust - Mistrust)" & pred_long$condition == "Male sorting Male"]), 2 )

## Save Plot
exp4_hyp <- ggplot(pred_long, aes(x=value)) + 
  facet_wrap(~condition) + 
  geom_density(aes(fill=name), alpha=0.7) + 
  geom_vline(xintercept = 0, linetype="dashed") + 
  scale_fill_manual(values=cols, label=c(expression(Delta(Trust - Mistrust)), "Pr(Mistrust Coethnie)", "Pr(Trust Coethnie)")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,7.5)) +
  xlab("Probability") + 
  ylab("") +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(),panel.spacing = unit(1.5, "lines"),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(
    label.position = "left",
    label.hjust = 1)
  )

### Plotting everything together
svg("experiment4_hyp.svg", width=12, height=6, pointsize=12)

exp4_hyp + ((trust_re + mistrust_re) / (trust_fixed + mistrust_fixed)) 

dev.off()


