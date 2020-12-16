library(tidyverse)
library(fastDummies)
library(ggridges)
library(patchwork)
library(rethinking)

d <- read_csv("data_norms.csv")

###### Data Dictionary ######
### site = village 
## 0 = Bevondrorano
## 1 = Namonte
## 2 = Antaimbalabo
## 3 = Tsiloakarivo
## 4 = Beangolo
## 5 = Ampasilava

### id = participant id

### sex = female (0), male (1)

### agecat1 = researcher-defined age category
## 0 = young adult
## 1 = adult
## 2 = old adult

### birthplace = type of environment individual was born in


### ethn1 = ethnicity (based on residence)
## 1 = Masikoro
## 2 = Mikea
## 3 = Vezo

### christian
## 0 = no
## 1 = yes

### edu = education (years)
## 0 - 12 indicates years of formal education

### magic = person is an ambiasa, tromba, mpanompotromba, or mpitankazomanga
## 0 = none
## 1 = one or more

### y = response to norm/knowledge questions

### category = domain of norm/knowledge
#############################

d <- d[!is.na(d$y),] # removing cases where individual was not asked the question
d$y <- ifelse(d$y == 9, NA, d$y) # changing "don't know" responses to NA 

d <- filter(d, name != "b2_A") # filtering out one question with effectively 0 variance

### Individual ethnic and village membership
d_ID <- d %>% group_by(id) %>%
  summarise(ethn1 = unique(ethn1), site=unique(site), birthplace=unique(birthplace)[1])

### Separate binary responses
d_bin <- d[!(d$category %in% c("Foraging Knowledge", "Maritime Knowledge","Agricultural Knowledge")),]

### And categorical
d_cat <- d[(d$category %in% c("Foraging Knowledge", "Maritime Knowledge","Agricultural Knowledge")),]

## Categorical vers A
d_catA <- filter(d_cat, ver == 1) %>% select(name, y, id) %>% pivot_wider(names_from=name, values_from=y)

# Collapsing response categories with 4 or fewer answers into '99'
d_catA$c6_A <- ifelse(d_catA$c6_A %in% c(2,3.5,6,10), 99, d_catA$c6_A)
d_catA$c7_A <- ifelse(d_catA$c7_A %in% c(3,0), 99, d_catA$c7_A)
d_catA$c8_A <- ifelse(d_catA$c8_A %in% c(4,5,6,12,13,14,15), 99, d_catA$c8_A)
d_catA$c5_A <- ifelse(d_catA$c5_A %in% c(0,4,5), 99, d_catA$c5_A)

# Making the categories "one-hot" encoded
onehot_A <- dummy_cols(d_catA, select_columns = c("c6_A", "c7_A", "c8_A", "c5_A"), ignore_na = T, remove_selected_columns = T)

## Categorical vers B
d_catB <- filter(d_cat, ver == 2) %>% select(name, y, id) %>% pivot_wider(names_from=name, values_from=y)

# Collapsing response categories with 2 or fewer answers into '99'
d_catB$c2_B <- ifelse(d_catB$c2_B %in% c(3,4), 99, d_catB$c2_B)
d_catB$c3_B <- ifelse(d_catB$c3_B %in% c(1,3.5,4.5,6,8,10,15,21,90), 99, d_catB$c3_B)
d_catB$c4_B <- ifelse(d_catB$c4_B %in% c(2,3,5,5.5,6.5,7.5,8.5,10), 99, d_catB$c4_B)
d_catB$c6_B <- ifelse(d_catB$c6_B %in% c(4,7,8,10,11,12,15,16,17,18), 99, d_catB$c6_B)
d_catB$c7_B <- ifelse(d_catB$c7_B %in% c(6,10,11,12,14,16,17,18,19), 99, d_catB$c7_B)
d_catB$c10_B <- ifelse(d_catB$c10_B %in% c(3,5), 99, d_catB$c10_B)

# Making the categories "one-hot" encoded
onehot_B <- dummy_cols(d_catB, select_columns = c("c2_B", "c3_B", "c4_B", "c5_B", "c6_B", "c7_B", "c10_B"), ignore_na = T, remove_selected_columns = T)

#### Bring back in individual covariates, and integrate the categorical variables
d_catA2 <- left_join( onehot_A, d_ID ) %>%
  pivot_longer(-c(id,ethn1,site,birthplace), values_to="y")

d_catA2$category <- "Subsistence Knowledge"

d_catB2 <- left_join( onehot_B, d_ID ) %>%
  pivot_longer(-c(id,ethn1,site,birthplace), values_to="y")

d_catB2$category <- "Subsistence Knowledge"

d_A_comb <- bind_rows(filter(d_bin, ver == 1), d_catA2)
d_B_comb <- bind_rows(filter(d_bin, ver == 2), d_catB2)

#########################################################
#### Calculating Cultural FST, Between Ethnic Groups ####
# Average frequency among all ethnicities
pki_A_eth <- d_A_comb %>% select(name, y, id, ethn1) %>%
  group_by(name, ethn1) %>%
  summarise(pk = mean(y,na.rm=T), n_eth=sum(!is.na(y)))

# Overall frequency of norm
pk_A_eth <- d_A_comb %>% select(name, y, id) %>% group_by(name) %>%
  summarise( pki = mean(y,na.rm=T), n_resp=sum(!is.na(y)) )

# Between-village CFST
CFST_A_eth <- left_join(pki_A_eth, pk_A_eth) %>%
  mutate(sq_diff = (n_eth/n_resp)*(pki - pk)^2) %>%
  ungroup() %>%
  group_by(name) %>%
  summarise( CFST_eth = sum(sq_diff) / (mean(pk)*(1-mean(pk))) )

# Bring back categories
CFST_A_eth <- left_join( CFST_A_eth, d_A_comb %>% group_by(name) %>% summarise(category = unique(category)) )

##### Version B #########
# Group frequency of norm
pki_B_eth <- d_B_comb %>% select(name, y, id, ethn1) %>%
  group_by(name, ethn1) %>%
  summarise(pk = mean(y,na.rm=T), n_eth=sum(!is.na(y)))

# Overall frequency of norm
pk_B_eth <- d_B_comb %>% select(name, y, id) %>% group_by(name) %>%
  summarise( pki = mean(y,na.rm=T), n_resp=sum(!is.na(y)) )

# Between-village CFST
CFST_B_eth <- left_join(pki_B_eth, pk_B_eth) %>%
  mutate(sq_diff = (n_eth/n_resp)*(pki - pk)^2) %>%
  ungroup() %>%
  group_by(name) %>%
  summarise( CFST_eth = sum(sq_diff) / (mean(pk)*(1-mean(pk))) )

# Bring back categories
CFST_B_eth <- left_join( CFST_B_eth, d_B_comb %>% group_by(name) %>% summarise(category = unique(category)) )
#########
# Combine
CFST_eth_both <- bind_rows(CFST_A_eth, CFST_B_eth)
#########################################################

#### Calculating Cultural FST, Between birth zones ####
# Average frequency among all zones
pki_A_zone <- d_A_comb %>%  filter(birthplace != "far") %>%
  select(name, y, id, birthplace) %>%
  group_by(name, birthplace) %>%
  summarise(pk = mean(y,na.rm=T), n_zone=sum(!is.na(y)))

# Overall frequency of norm
pk_A_zone <- d_A_comb %>% filter(birthplace != "far") %>% 
  select(name, y, id) %>% group_by(name) %>%
  summarise( pki = mean(y,na.rm=T), n_resp=sum(!is.na(y)) )

# Between-village CFST
CFST_A_zone <- left_join(pki_A_zone, pk_A_zone) %>%
  mutate(sq_diff = (n_zone/n_resp)*(pki - pk)^2) %>%
  ungroup() %>%
  group_by(name) %>%
  summarise( CFST_zone = sum(sq_diff) / (mean(pk)*(1-mean(pk))) )

# Bring back categories
CFST_A_zone <- left_join( CFST_A_zone, d_A_comb %>% group_by(name) %>% summarise(category = unique(category)) )

##### Version B #########
# Group frequency of norm
pki_B_zone <- d_B_comb %>% filter(birthplace != "far") %>%
  select(name, y, id, birthplace) %>%
  group_by(name, birthplace) %>%
  summarise(pk = mean(y,na.rm=T), n_zone=sum(!is.na(y)))

# Overall frequency of norm
pk_B_zone <- d_B_comb %>% filter(birthplace != "far") %>%
  select(name, y, id) %>% group_by(name) %>%
  summarise( pki = mean(y,na.rm=T), n_resp=sum(!is.na(y)) )

# Between-village CFST
CFST_B_zone <- left_join(pki_B_zone, pk_B_zone) %>%
  mutate(sq_diff = (n_zone/n_resp)*(pki - pk)^2) %>%
  ungroup() %>%
  group_by(name) %>%
  summarise( CFST_zone = sum(sq_diff) / (mean(pk)*(1-mean(pk))) )

# Bring back categories
CFST_B_zone <- left_join( CFST_B_zone, d_B_comb %>% group_by(name) %>% summarise(category = unique(category)) )
#########
# Combine
CFST_zone_both <- bind_rows(CFST_A_zone, CFST_B_zone)
#########################################################

#########################################################
#### Calculating Cultural FST, Between Villages #########
### Version A ##
# Group frequency of norm at village level
pki_A_v <- d_A_comb %>% select(name, y, id, ethn1, site) %>% group_by(name, site,ethn1) %>%
  summarise( pki = mean(y,na.rm=T), n_village=sum(!is.na(y)) )

# Average frequency among all villages within eth
pk_A_v <- d_A_comb %>% select(name, y, id, ethn1, site) %>%
  group_by(name, ethn1) %>%
  summarise(pk = mean(y,na.rm=T), n_eth=sum(!is.na(y)))

# Between-village CFST
CFST_A_village <- left_join(pki_A_v, pk_A_v) %>%
  mutate(sq_diff = (n_village/n_eth)*(pki - pk)^2) %>%
  ungroup() %>%
  group_by(name, ethn1) %>%
  summarise( CFST = sum(sq_diff) / (mean(pk)*(1-mean(pk))) ) %>%
  mutate( CFST = ifelse(is.nan(CFST), NA, CFST) ) %>%
  summarise( CFST_village = mean(CFST, na.rm = T) )

# Bring back categories
CFST_A_village <- left_join( CFST_A_village, d_A_comb %>% group_by(name) %>% summarise(category = unique(category)) )

##### Version B #########
# Group frequency of norm
pki_B_v <- d_B_comb %>% select(name, y, id, ethn1, site) %>% group_by(name, site,ethn1) %>%
  summarise( pki = mean(y,na.rm=T), n_village=sum(!is.na(y)) )

# Average frequency among all villages within eth
pk_B_v <- d_B_comb %>% select(name, y, id, ethn1, site) %>%
  group_by(name, ethn1) %>%
  summarise(pk = mean(y,na.rm=T), n_eth=sum(!is.na(y)))

# Between-village CFST
CFST_B_village <- left_join(pki_B_v, pk_B_v) %>%
  mutate(sq_diff = (n_village/n_eth)*(pki - pk)^2) %>%
  ungroup() %>%
  group_by(name, ethn1) %>%
  summarise( CFST = sum(sq_diff) / (mean(pk)*(1-mean(pk))) ) %>%
  mutate( CFST = ifelse(is.nan(CFST), NA, CFST) ) %>%
  summarise( CFST_village = mean(CFST, na.rm = T) )

# Bring back categories
CFST_B_village <- left_join( CFST_B_village, d_B_comb %>% group_by(name) %>% summarise(category = unique(category)) )
#########
# Combine
CFST_v_both <- bind_rows(CFST_A_village, CFST_B_village)

###############################################
### Now, integrate the ethnic, zone, and village CFST
CFST_all <- left_join(CFST_eth_both, CFST_v_both)
CFST_all <- left_join(CFST_all, CFST_zone_both)

CFST_all$diff_EV <- CFST_all$CFST_eth - CFST_all$CFST_village
CFST_all$diff_EZ <- CFST_all$CFST_eth - CFST_all$CFST_zone

# We only want to compare CFST where we have data at all level
CFST_long <- CFST_all[complete.cases(CFST_all),] %>% 
  group_by(name) %>%
  pivot_longer(-c(category,name), names_to="level", values_to="CFST")

CFST_long$level <- factor(CFST_long$level, labels=c("Ethnicity", "Village", "Natal\nEnvironment", "Difference1", "Difference2"))

CFST_long$grouping <- ifelse(CFST_long$level %in% c("Ethnicity", "Village", "Natal\nEnvironment"), "Cultural FST", NA)
CFST_long$grouping <- ifelse(CFST_long$level == "Difference1", "Diff(Ethnicity FST - Village FST)", CFST_long$grouping)
CFST_long$grouping <- ifelse(CFST_long$level == "Difference2", "Diff(Ethnicity FST - Natal Environment FST)", CFST_long$grouping)

CFST_long$category <- factor(CFST_long$category)
CFST_long$category <- fct_relevel(CFST_long$category, "Subsistence Knowledge", "Food Taboos", "Social Organization & Gender", "Supernatural")

#####################################################
##### Organize data for Stan ########################
CFST <- filter(CFST_long, grouping == "Cultural FST")

norm <- match(CFST$name, unique(CFST$name))
cat <- match(CFST$category, unique(CFST$category))
level <- match(CFST$level, unique(CFST$level))

data_list <- list(
  N = nrow(CFST),
  CFST = CFST$CFST,
  norm = norm,
  cat = cat,
  level = level
)

fit <- stan( file="cfst_model.stan", data=data_list, iter=5000, chains=4, cores=4, init="0", control=list(adapt_delta=0.99) )
#fit <- readRDS("fit_cfst.rds")

post <- extract.samples(fit)
n_samps <- length(post$lp__)


##### CFST for subsistence knowledge ##############
## Level 1 = eth, 2 = village, 3= natal environment
## Cat 1 = Social Org & Gender, 2 = Subsistence Knowledge, 3 = Supernatural, 4 = Food Taboos

# Ethnicity
eth_sub <- logistic(post$a + post$level_v[,1] + post$level_cat_v[,2,1] + post$level_cat_v[,2,(1+1)])
round( median(eth_sub), 3 )
round( HPDI(eth_sub, prob=0.9), 3 )

# Village
village_sub <- logistic(post$a + post$level_v[,2] + post$level_cat_v[,2,1] + post$level_cat_v[,2,(2+1)])
round( median(village_sub), 3 )
round( HPDI(village_sub, prob=0.9), 3 )

# Natal environment
natal_sub <- logistic(post$a + post$level_v[,3] + post$level_cat_v[,2,1] + post$level_cat_v[,2,(3+1)])
round( median(natal_sub), 3 )
round( HPDI(natal_sub, prob=0.9), 3 )

# Diff env - ethnicity
round( median(natal_sub - eth_sub), 3 )
round( HPDI(natal_sub - eth_sub, prob=0.9), 3 )
length((natal_sub - eth_sub)[(natal_sub - eth_sub) > 0]) / length((natal_sub - eth_sub)) # PP

# Diff env - village
round( median(natal_sub - village_sub), 3 )
round( HPDI(natal_sub - village_sub, prob=0.9), 3 )
length((natal_sub - village_sub)[(natal_sub - village_sub) > 0]) / length((natal_sub - village_sub)) # PP

##################################################

##### Make posterior predictions for each category*level


# Expected values
cfst_mu <- array( NA, dim = c(n_samps, 4, 3), dimnames=list( samp=1:n_samps, category = c("Social Organization & Gender","Subsistence Knowledge", "Supernatural", "Food Taboos"), level = c("Ethnicity", "Village", "Natal\nEnvironment")) )
# Predictions with dispersion
cfst_pred <- cfst_mu

for ( cat in 1:4 )
  for ( l in 1:3 ) {
    cfst_mu[,cat,l] <- logistic( post$a + post$level_v[,l] + post$level_cat_v[,cat,1] + post$level_cat_v[,cat,(l+1)] ) 
    
    for (i in 1:n_samps)  cfst_pred[i,cat,l] <- rbeta( 1, logistic( post$a[i] + post$level_v[i,l] + post$level_cat_v[i,cat,1] + post$level_cat_v[i,cat,(l+1)] )*(1/post$phi[i]),  (1-logistic( post$a[i] + post$level_v[i,l] + post$level_cat_v[i,cat,1] + post$level_cat_v[i,cat,(l+1)] ))*(1/post$phi[i]) )
  }

## Convert to long form and put together
cfst_pred_long <- cfst_pred %>% 
  as.tbl_cube(met_name = "CFST") %>% 
  as_tibble

cfst_mu_long <- cfst_mu %>% 
  as.tbl_cube(met_name = "CFST") %>% 
  as_tibble

## Get contrasts
cfst_contrasts <- cfst_mu_long %>% pivot_wider(names_from=level, values_from=CFST) %>%
  mutate(  eth_vill = Ethnicity - Village, env_eth = `Natal\nEnvironment` - Ethnicity) %>%
  select(c(eth_vill, env_eth, samp, category))

## Plot
cols <- c("darkorchid4", "forestgreen", "indianred")

CFST_main <- ggplot(cfst_mu_long, aes(x=CFST, y=fct_rev(category), color=level, fill=level)) +
  geom_hline(yintercept = c(1:4)) +
  geom_density_ridges(scale=0.9, alpha=0.6, color="black") + 
  geom_jitter(data=CFST, aes(x=CFST, y=fct_rev(category), color=level, fill=level), alpha=0.8, height=0.025, width=0.01) +
  scale_x_continuous(expand = c(0.00, 0.01), breaks=c(0,0.25,0.5,0.75,1), labels=c("0",".25",".5",".75","1")) + 
  scale_y_discrete(expand = c(0.00, 0)) + 
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(), legend.position = "top") +
  guides(color=F) +
  xlab("Cultural FST") + 
  ylab("")  + 
  annotate("blank", x = 0, y=5)

## PP annotations
contrast1_ann <- data.frame(
  category = c("Supernatural","Subsistence Knowledge", "Social Organization & Gender", "Food Taboos"),
  y = 1:4  +  0.7,
  PP = c(
    sum(cfst_contrasts$eth_vill[cfst_contrasts$category == "Supernatural"] >= 0)/length(cfst_contrasts$eth_vill[cfst_contrasts$category == "Supernatural"]),
    sum(cfst_contrasts$eth_vill[cfst_contrasts$category == "Subsistence Knowledge"] >= 0)/length(cfst_contrasts$eth_vill[cfst_contrasts$category == "Subsistence Knowledge"]),
    sum(cfst_contrasts$eth_vill[cfst_contrasts$category == "Social Organization & Gender"] >= 0)/length(cfst_contrasts$eth_vill[cfst_contrasts$category == "Social Organization & Gender"]),
    sum(cfst_contrasts$eth_vill[cfst_contrasts$category == "Food Taboos"] >= 0)/length(cfst_contrasts$eth_vill[cfst_contrasts$category == "Food Taboos"])
  ),
  x = -0.25
)

contrast1_ann$PP <- paste("PP =", round(contrast1_ann$PP,2))


contrast1 <- ggplot(cfst_contrasts, aes(x=eth_vill, y=fct_rev(category))) +
  geom_hline(yintercept = c(1:4)) +
  geom_density_ridges(scale=0.9, alpha=0.6, color="black") + 
  scale_x_continuous(expand = c(0.00, 0.01), limits=c(-0.4, 0.25)) +
  scale_y_discrete(expand = c(0.00, 0)) + 
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(), legend.position = "top") +
  xlab(expression(paste(Delta, "CFST(Ethnicity - Village)"))) + 
  ylab("")  + 
  annotate("text", x=contrast1_ann$x, y=contrast1_ann$y, label=contrast1_ann$PP, size=4) +
  annotate("blank", x = 0, y=5)


## PP annotations
contrast2_ann <- data.frame(
  category = c("Supernatural","Subsistence Knowledge", "Social Organization & Gender", "Food Taboos"),
  y = 1:4  +  0.7,
  PP = c(
    sum(cfst_contrasts$env_eth[cfst_contrasts$category == "Supernatural"] >= 0)/length(cfst_contrasts$env_eth[cfst_contrasts$category == "Supernatural"]),
    sum(cfst_contrasts$env_eth[cfst_contrasts$category == "Subsistence Knowledge"] >= 0)/length(cfst_contrasts$env_eth[cfst_contrasts$category == "Subsistence Knowledge"]),
    sum(cfst_contrasts$env_eth[cfst_contrasts$category == "Social Organization & Gender"] >= 0)/length(cfst_contrasts$env_eth[cfst_contrasts$category == "Social Organization & Gender"]),
    sum(cfst_contrasts$env_eth[cfst_contrasts$category == "Food Taboos"] >= 0)/length(cfst_contrasts$env_eth[cfst_contrasts$category == "Food Taboos"])
  ),
  x = -0.25
)

contrast2_ann$PP <- paste("PP =", round(contrast2_ann$PP,2))
contrast2_ann$PP[2] <- "PP = ~1" # asymptotically close to 1, but technically the posterior prob here couldn't be exactly 1

contrast2 <- ggplot(cfst_contrasts, aes(x=env_eth, y=fct_rev(category))) +
  geom_hline(yintercept = c(1:4)) +
  geom_density_ridges(scale=0.9, alpha=0.6, color="black") + 
  scale_x_continuous(expand = c(0.00, 0.01), limits=c(-0.4, 0.25)) +
  scale_y_discrete(expand = c(0.00, 0)) + 
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  theme_bw(base_size=12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(), legend.position = "top") +
  xlab(expression(paste(Delta, "CFST(Natal Env. - Ethnicity)"))) + 
  ylab("")  + 
  annotate("text", x=contrast2_ann$x, y=contrast2_ann$y, label=contrast2_ann$PP, size=4) +
  annotate("blank", x = 0, y=5)


#### Plot and save the distribution of CFST #################
svg(filename="CFST_data.svg", width=8.5, height=7, pointsize=12)

CFST_main / (contrast1 + contrast2)

dev.off()


#### Export summaries as tables ####
CFST_med <- cfst_mu_long %>% group_by(level, category) %>% summarise(median_CFST = median(CFST), lower_90 = HPDI(CFST, prob=0.9)[1], upper_90 = HPDI(CFST, prob=0.9)[2])

eth_village <- cfst_contrasts %>% group_by(category) %>% summarise(diff = median(eth_vill), lower_90 = HPDI(eth_vill, prob=0.9)[1], upper_90 = HPDI(eth_vill, prob=0.9)[2])

env_eth <- cfst_contrasts %>% group_by(category) %>% summarise(diff = median(env_eth), lower_90 = HPDI(env_eth, prob=0.9)[1], upper_90 = HPDI(env_eth, prob=0.9)[2])

write_csv(CFST_med, "CFST_estimates.csv")
write_csv(eth_village, "eth_village_contrasts.csv")
write_csv(env_eth, "env_eth_contrasts.csv")



