# replace "~/Desktop/Shared Code for Soon" by "folder-location"

# load libraries 
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
n_chains <- 4


dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix <<- prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {


loo_m05 <- readRDS(paste(pfix, parameters["loo_m05", 2], sep="/"))
loo_pop_2_entceae <- readRDS(paste(pfix, parameters["loo_pop_2_entceae", 2], sep="/"))
loo_pop_2_entles <- readRDS(paste(pfix, parameters["loo_pop_2_entles", 2], sep="/"))
loo_pop_2_e4 <- readRDS(paste(pfix, parameters["loo_pop_2_e4", 2], sep="/"))
loo_pop_1a2_entceae <- readRDS(paste(pfix, parameters["loo_pop_1a2_entceae", 2], sep="/"))
loo_pop_1a2_entles <- readRDS(paste(pfix, parameters["loo_pop_1a2_entles", 2], sep="/"))
loo_pop_1a2_e4 <- readRDS(paste(pfix, parameters["loo_pop_1a2_e4", 2], sep="/"))
loo_pop_2_notax <- readRDS(paste(pfix, parameters["loo_pop_2_notax", 2], sep="/"))
loo_pop_1a1_notax <- readRDS(paste(pfix, parameters["loo_pop_1a1_notax", 2], sep="/"))
pop_1a2_e4 <- readRDS(paste(pfix, parameters["pop_1a2_e4", 2], sep="/"))
n_atbs <- readRDS(paste(pfix, parameters["n_atbs", 2], sep="/"))

# Model comparison --------------------------------------------------------


wts2 <- loo_model_weights(
                list(loo_m05,
                     loo_pop_2_entceae,
                     loo_pop_2_entles,
                     loo_pop_2_e4,
                     loo_pop_1a2_entceae,
                     loo_pop_1a2_entles,
                     loo_pop_1a2_e4,
                     loo_pop_2_notax,
                     loo_pop_1a1_notax),
                method = "pseudobma",
                optim_control = list(reltol=1e-10)
)

df_comp <- compare( loo_m05,
                    loo_pop_2_entceae,
                    loo_pop_2_entles,
                    loo_pop_2_e4,
                    loo_pop_1a2_entceae,
                    loo_pop_1a2_entles,
                    loo_pop_1a2_e4,
                    loo_pop_2_notax,
                    loo_pop_1a1_notax )
df_comp <- df_comp %>% as_tibble(rownames = "model") %>% 
                select(model,elpd_loo,elpd_diff) %>% 
                left_join( tibble(
                                model=c("loo_pop_1a2_e4",
                                        "loo_pop_1a2_entceae",
                                        "loo_pop_2_e4",
                                        "loo_pop_2_entceae",
                                        "loo_pop_1a2_entles",
                                        "loo_pop_2_entles",
                                        "loo_pop_1a1_notax",
                                        "loo_pop_2_notax",
                                        "loo_m05"), 
                                model_nname=c("ALL_e4",
                                              "ALL_entceae",
                                              "DEF_e4",
                                              "DEF_entceae",
                                              "ALL_entles",
                                              "DEF_entles",
                                              "ALL_notax",
                                              "DEF_notax",
                                              "baseline")
                ), by="model") %>% 
                mutate( model=model_nname ) %>% select(-model_nname) %>% 
                rename( Model=model,
                        loo_prediction=elpd_loo,
                        loo_diff_to_bestmodel=elpd_diff)

wts2 %>% enframe() %>% 
                mutate(name=c("baseline",
                              "DEF_entceae",
                              "DEF_entles",
                              "DEF_e4",
                              "ALL_entceae",
                              "ALL_entles",
                              "ALL_e4",
                              "DEF_notax",
                              "ALL_notax")) %>% 
                mutate(value=round(value,3)) %>% 
                rename(BMA_weight=value,
                       Model=name) %>% 
                left_join( df_comp, by="Model" ) %>% arrange(desc(BMA_weight)) %>% 
                mutate(loo_prediction=round(loo_prediction,2),
                       loo_diff_to_bestmodel=round(loo_diff_to_bestmodel,2)) 



# Plot model predictions with data ----------------------------------------

mrel <- 1
n_sett <- 3
n_genes <- n_atbs
# model with metagenomic information
fit_post <- rstan::extract(pop_1a2_e4)
post_pred <- fit_post$pred_count_sr
pred_lower <- matrix(NA,nrow=n_sett,ncol=n_genes)
pred_upper <- matrix(NA,nrow=n_sett,ncol=n_genes)
pred_mean <- matrix(NA,nrow=n_sett,ncol=n_genes)
for (i in 1:n_sett) {
                for (g in 1:n_genes){
                                post_PI <- bayestestR::hdi(post_pred[,i,g],ci=0.95 ) %>% unlist() 
                                pred_lower[i,g] <- post_PI[2]
                                pred_upper[i,g] <- post_PI[3]
                                pred_mean[i,g] <- mean(post_pred[,i,g])
                }
}

# wrrangle real data
count_s <- readRDS(paste(pfix, parameters["count_s", 2], sep="/"))
count_sr <- readRDS(paste(pfix, parameters["count_sr", 2], sep="/"))
count_s_for_pred <- readRDS(paste(pfix, parameters["count_s_for_pred", 2], sep="/"))
setting_v <- readRDS(paste(pfix, parameters["setting_v", 2], sep="/"))
df1 <- count_sr %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=real,-Setting)
tot <- count_s_for_pred %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=tot,-Setting) %>% .$tot
tot_withNA <- count_s %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=tot,-Setting) %>% .$tot

# wrangle predictions
pred_mean_v <- pred_mean %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=pred_mean,-Setting) %>% .$pred_mean
pred_lower_v <- pred_lower %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=pred_lower,-Setting) %>% .$pred_lower
pred_upper_v <- pred_upper %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=pred_upper,-Setting) %>% .$pred_upper

df2 <- df1 %>% mutate( pred_lower=pred_lower_v,
                       pred_upper=pred_upper_v,
                       pred_mean=pred_mean_v,
                       tot=tot,
                       tot_withNA=tot_withNA,
                       close_brack=")",
                       open_brack=" (") %>% 
                mutate(Setting=replace(Setting,Setting=="CAPOP","Cambodia"),
                       Setting=replace(Setting,Setting=="KEPOP","Kenya"),
                       Setting=replace(Setting,Setting=="UKPOP","UK"))

sett_df <- tibble( old_sett=c("CAMBODIA_POP_POOL",
                              "KENYA_POP_POOL",
                              "UK_POP_POOL"),
                   new_sett=c("Cambodia",
                              "Kenya",
                              "UK"))

(p1 <- df2 %>% arrange(Setting) %>% unite(col = Sett_Anti,Setting,Antibiotic,sep="-",remove=F) %>% 
                                unite(col = Antibiotic_count,Antibiotic,open_brack,tot_withNA,close_brack,sep="",remove=F) %>% 
                                # for color of bars
                                mutate( col_bars=ifelse(is.na(tot_withNA), "no_clin_cases" , Setting) ) %>% 
                                ggplot( aes(y=Antibiotic_count) ) + 
                                geom_point( aes(x=pred_mean, col=col_bars ),alpha=0.7,size=4,pch="|" ,show.legend = F) +
                                geom_segment( aes(yend=Antibiotic_count,x=pred_lower,xend=pred_upper,
                                                  col=col_bars),alpha=0.50,size=2.3,show.legend = F ) + # 1model
                                scale_color_manual(values=c('CAMBODIA_POP_POOL'='#efc750',
                                                            'KENYA_POP_POOL'='#a6cdd9',
                                                            'UK_POP_POOL'='#b7b079',
                                                            'no_clin_cases'='darkgrey') ) +
                                geom_point( aes(x=real), size=2) + 
                                facet_wrap(~Setting,ncol=1,scales="free") +
                                labs(x="Resistance Proportion",y="",title="")) 

}
