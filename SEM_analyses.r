# Works in C:\Program Files\R\R-4.0.2\bin\x64

## Install piecewiseSEM from github
#library(devtools)
#install_github("jslefche/piecewiseSEM@devel", build_vignette = TRUE) # requires R version > 3.5

# Load library
library(piecewiseSEM)

##### load data
setwd('~/StabilityScaleUp')
varpart_All_prod <- read.csv('varpart_All_prod.csv')
head(varpart_All_prod)

#######################
#### Piecewise SEM ####
#######################

coefs_df = data.frame() 
for(i in unique(varpart_All_prod$year_trt)) {
  for(j in unique(varpart_All_prod$trt)) {
     df <- varpart_All_prod[varpart_All_prod$year_trt==i&varpart_All_prod$trt==j,]
        model = psem(
        lm(log_gamma_stab_prod ~ log_spatial_asynch_prod + log_alpha_stab_prod, data=df),
        lm(log_alpha_stab_prod ~ log_spp_stab_cover + log_spp_asynch_cover, data=df),
        lm(log_spp_stab_cover ~ avg_richness2, data=df),
        lm(log_spp_asynch_cover ~ avg_richness2, data=df),
        lm(log_spatial_asynch_prod ~ beta_invsimpsons, data=df),
        avg_richness2 %~~% beta_invsimpsons,
        log_spp_stab_cover %~~% log_spp_asynch_cover,
        log_spp_stab_cover %~~% log_spatial_asynch_prod,
        log_spp_asynch_cover %~~% log_spatial_asynch_prod
        )
        coefs <- summary(model, .progressBar = F)$coefficients
        coefs$year <- i
        coefs$trt <- j
        names(coefs)[names(coefs) == ''] <- 'Sig'
        coefs_df = plyr::rbind.fill(coefs_df, coefs)
  }
}
# rename
names(coefs_df)[names(coefs_df) == 'Std.Estimate'] <- 'TE'; names(coefs_df)[names(coefs_df) == 'Std.Error'] <- 'seTE'
head(coefs_df) 

# replace SE of 0 with 0.01
coefs_df$seTE[coefs_df$seTE==0] <- 0.01
# replace SE of NA with 0.01
coefs_df$seTE[is.na(coefs_df$seTE)] <- 0.01

# replace by numerical values
coefs_df$seTE <- as.numeric(as.character(coefs_df$seTE))

# replace NA with 0.01
coefs_df[is.na(coefs_df)] <- 0.01

#######################
#### Meta analysis ####
#######################

library(meta)

p1 <- coefs_df[coefs_df$Response=='log_gamma_stab_prod'&coefs_df$Predictor=='log_spatial_asynch_prod',]
p2 <- coefs_df[coefs_df$Response=='log_gamma_stab_prod'&coefs_df$Predictor=='log_alpha_stab_prod',]
p3 <- coefs_df[coefs_df$Response=='log_alpha_stab_prod'&coefs_df$Predictor=='log_spp_stab_cover',]
p4 <- coefs_df[coefs_df$Response=='log_alpha_stab_prod'&coefs_df$Predictor=='log_spp_asynch_cover',]
p5 <- coefs_df[coefs_df$Response=='log_spp_stab_cover'&coefs_df$Predictor=='avg_richness2',]
p6 <- coefs_df[coefs_df$Response=='log_spp_asynch_cover'&coefs_df$Predictor=='avg_richness2',]
p8 <- coefs_df[coefs_df$Response=='log_spatial_asynch_prod'&coefs_df$Predictor=='beta_invsimpsons',]
p9 <- coefs_df[coefs_df$Response=='~~avg_richness2'&coefs_df$Predictor=='~~beta_invsimpsons',]
p10 <- coefs_df[coefs_df$Response=='~~log_spp_stab_cover'&coefs_df$Predictor=='~~log_spp_asynch_cover',]
p11 <- coefs_df[coefs_df$Response=='~~log_spp_stab_cover'&coefs_df$Predictor=='~~log_spatial_asynch_prod',]
p12 <- coefs_df[coefs_df$Response=='~~log_spp_asynch_cover'&coefs_df$Predictor=='~~log_spatial_asynch_prod',]

# check whether estimates vary by treatment
out <- do.call(rbind, lapply(list(p1, p2, p3, p4, p5, p6, p8, p9, p10, p11, p12), function(i) {
  r <- i$Response[1]
  p <- i$Predictor[1]
  
  # Fit model with interaction
  model <- metagen(TE, seTE, studlab = paste(Response), comb.fixed = TRUE, comb.random = TRUE, prediction = TRUE, sm = "SMD", byvar = trt,
    data = subset(coefs_df, Response == r & Predictor == p))

  out <- data.frame(
    Response = r,
    Predictor = p,
    estimate_C = model$TE.fixed.w[1],
    ci.lb_C = model$lower.fixed.w[1], 
    ci.ub_C = model$upper.fixed.w[1],
    estimate_F = model$TE.fixed.w[2],
    ci.lb_F = model$lower.fixed.w[2], 
    ci.ub_F = model$upper.fixed.w[2]
  )

  return(out)
  
} ) )
out

