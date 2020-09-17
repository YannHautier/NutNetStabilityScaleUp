# Works in C:\Program Files\R\R-4.0.2\bin\x64

##### Load packages                     
library(tidyverse) 
library(dplyr)     
library(ggplot2)     
library(nlme)
library(gridExtra)
library(viridis)

##### load data
setwd('~/StabilityScaleUp')
varpart_All_prod <- read.csv('varpart_All_prod.csv')
varpart_All_prod$year_trt <- factor(varpart_All_prod$year_trt)
varpart_All_prod$trt <- factor(varpart_All_prod$trt)
head(varpart_All_prod)

varpart_All_LRR_prod <- read.csv('varpart_All_LRR_prod.csv')
varpart_All_LRR_prod$year_trt <- factor(varpart_All_LRR_prod$year_trt)
head(varpart_All_LRR_prod)

###########################################
# correlation between richness and other local diversity indices #
###########################################

cor.test(varpart_All_prod$avg_richness2, varpart_All_prod$avg_evenness_Sh)
cor.test(varpart_All_prod$avg_richness2, varpart_All_prod$avg_invsimpsons2)
cor.test(varpart_All_prod$avg_richness2, varpart_All_prod$avg_shannons)

###########################################
###### figure MAT MAP              ########
###########################################

jpeg('FigureS1.jpg', width = 4, height = 4, units='in', res=300)
qplot(TEMP_VAR_v2/100, MAP_VAR_v2, data = varpart_All_prod, xlab='Temperature seasonality (°C)', ylab='Precipitation seasonality (mm)', alpha = I(0)) + theme_bw() + geom_point(size = 2)
dev.off()
system('open FigureS1.jpg')

###########################################
###### fertilization impact on diversity  ########
###########################################

# richness
gls1 <- gls(avg_richness2 ~ trt * year_trt, method='ML', data = varpart_All_prod)
gls2 <- gls(avg_richness2 ~ trt * year_trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

#beta
gls1 <- gls(beta_invsimpsons ~ trt * year_trt, method='ML', data = varpart_All_prod)
gls2 <- gls(beta_invsimpsons ~ trt * year_trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

###########################################
###### fertilization impact on stability  ########
###########################################

# alpha stability
gls1 <- gls(log_alpha_stab_prod ~ trt * year_trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_alpha_stab_prod ~ trt * year_trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

# gamma stability
gls1 <- gls(log_gamma_stab_prod ~ trt * year_trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_gamma_stab_prod ~ trt * year_trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

###########################################
###### figure treatment effect on diversity and stability           ########
###########################################

pdf('FigureS2.pdf', width = 8, height = 8, useDingbats=FALSE)
grid.arrange(
qplot(year_trt, avg_richness2, data = varpart_All_prod, xlab='Post treatment year', ylab='Species richness (m-2)', alpha = I(0)) + theme_bw() + geom_boxplot(aes(fill=trt)) + scale_fill_manual(values=c('#ffffff', '#999999')) + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')),
qplot(year_trt, beta_invsimpsons, data = varpart_All_prod, xlab='Post treatment year', ylab='Beta diversity (abundance-based)', alpha = I(0)) + theme_bw() + geom_boxplot(aes(fill=trt)) + scale_fill_manual(values=c('#ffffff', '#999999')),
qplot(year_trt, log_alpha_stab_prod, data = varpart_All_prod, xlab='Post treatment year', ylab='Log (alpha stability)', alpha = I(0)) + theme_bw() + geom_boxplot(aes(fill=trt)) + scale_fill_manual(values=c('#ffffff', '#999999')) + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')),
qplot(year_trt, log_gamma_stab_prod, data = varpart_All_prod, xlab='Post treatment year', ylab='Log (gamma stability)', alpha = I(0)) + theme_bw() + geom_boxplot(aes(fill=trt)) + scale_fill_manual(values=c('#ffffff', '#999999')),
nrow=2)
dev.off()
system('open FigureS2.pdf')

###########################################
###### figure mean/sd           ########
###########################################

pdf('FigureS3.pdf', width = 8, height = 3, useDingbats=FALSE)
grid.arrange(
qplot(avg_richness2, log(temp_mean), col=year_trt, data = varpart_All_prod, xlab=expression('Species richness' ~ ('m'^-2)), ylab='Temporal mean of productivity', alpha = I(0)) + theme_bw() + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_richness2, log(temp_sd), col=year_trt, data = varpart_All_prod, xlab=expression('Species richness' ~ ('m'^-2)), ylab='Temporal sd of productivity', alpha = I(0)) + theme_bw() + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
ncol=2)
dev.off()
system('open FigureS3.pdf')

gls1 <- gls(log(temp_mean) ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log(temp_mean) ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(log(temp_sd) ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log(temp_sd) ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

###########################################
###### analyses and figures LRR    ########
###########################################

########################################################################################################################################################################
# prod
########################################################################################################################################################################

pdf('FigureS9.pdf', width = 4, height = 6, useDingbats=FALSE)
grid.arrange(
qplot(LRR_avg_richness, LRR_alpha_stab_prod, col=year_trt, data = varpart_All_LRR_prod, xlab=expression('Changes in species richness' ~ ('m'^-2)), ylab='Changes in alpha stability', alpha = I(0)) + theme_bw() + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + scale_color_viridis(discrete=TRUE),
qplot(LRR_beta_invsimpsons, LRR_gamma_stab_prod, col=year_trt, data = varpart_All_LRR_prod, xlab='Changes in beta diversity (abundance-based)', ylab='Changes in gamma stability', alpha = I(0)) + theme_bw() + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + scale_color_viridis(discrete=TRUE),
ncol=1)
dev.off()
system('open FigureS9.pdf')

# richness
gls1 <- gls(LRR_alpha_stab_prod ~ LRR_avg_richness * year_trt, method='ML', data = varpart_All_LRR_prod)
gls2 <- gls(LRR_alpha_stab_prod ~ LRR_avg_richness * year_trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_LRR_prod)
anova(gls1, gls2)
anova(gls2)

# beta
gls1 <- gls(LRR_gamma_stab_prod ~ LRR_beta_invsimpsons * year_trt, method='ML', data = varpart_All_LRR_prod)
gls2 <- gls(LRR_gamma_stab_prod ~ LRR_beta_invsimpsons * year_trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_LRR_prod)
anova(gls1, gls2)
anova(gls2)

###########################################
###### analyses and figures control vs NPK      ########
###########################################

########################################################################################################################################################################
# prod
########################################################################################################################################################################

pdf('Figure2.pdf', width = 8, height = 6, useDingbats=FALSE)
grid.arrange(
qplot(avg_richness2, log_alpha_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Species richness (m-2)', ylab='Log (alpha stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(beta_invsimpsons, log_spatial_asynch_prod, col=year_trt, data = varpart_All_prod, xlab='Beta diversity (abundance-based)', ylab='Log (spatial asynchrony)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_richness2, log_gamma_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Species richness (m-2)', ylab='Log (gamma stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(beta_invsimpsons, log_gamma_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Beta diversity (abundance-based)', ylab='Log (gamma stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + facet_wrap('trt') + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + scale_color_viridis(discrete=TRUE),
ncol=2)
dev.off()
system('open Figure2.pdf')

pdf('FigureS4.pdf', width = 8, height = 6, useDingbats=FALSE)
grid.arrange(
qplot(avg_richness2, res_log_alpha_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Species richness (m-2)', ylab='Log (alpha stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(beta_invsimpsons, res_log_spatial_asynch_prod, col=year_trt, data = varpart_All_prod, xlab='Beta diversity (abundance-based)', ylab='Log (spatial asynchrony)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_richness2, res_log_gamma_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Species richness (m-2)', ylab='Log (gamma stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(beta_invsimpsons, res_log_gamma_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Beta diversity (abundance-based)', ylab='Log (gamma stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + facet_wrap('trt') + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + scale_color_viridis(discrete=TRUE),
ncol=2)
dev.off()
system('open FigureS4.pdf')

pdf('FigureS5.pdf', width = 12, height = 6, useDingbats=FALSE)
grid.arrange(
qplot(avg_invsimpsons2, log_alpha_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Inverse Simpson', ylab='Log (alpha stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_shannons, log_alpha_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Shannon', ylab='Log (alpha stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_evenness_Sh, log_alpha_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Species evenness', ylab='Log (alpha stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_invsimpsons2, res_log_gamma_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Inverse Simpson', ylab='Log (gamma stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_shannons, log_gamma_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Shannon', ylab='Log (gamma stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_evenness_Sh, log_gamma_stab_prod, col=year_trt, data = varpart_All_prod, xlab='Species evenness', ylab='Log (alpha stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
ncol=3)
dev.off()
system('open FigureS5.pdf')

pdf('FigureS8.pdf', width = 5, height = 6, useDingbats=FALSE)
grid.arrange(
qplot(avg_richness2, log_spp_stab_cover, col=year_trt, data = varpart_All_prod, xlab=expression('Species richness' ~ ('m'^-2)), ylab='Log (species stability)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + scale_x_continuous(breaks=c(0,1,2,3,4,5), labels=c('1','2','4','8','16','32')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
qplot(avg_richness2, log_spp_asynch_cover, col=year_trt, data = varpart_All_prod, xlab=expression('Species richness' ~ ('m'^-2)), ylab='Log (species asynchrony)', alpha = I(0)) + theme_bw() + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c('0','1','2','3','4','5')) + geom_point(size = 1, alpha = I(0.3), shape=16) + geom_smooth(method=lm, se=F) + facet_wrap('trt') + scale_color_viridis(discrete=TRUE),
ncol=1)
dev.off()
system('open FigureS8.pdf')

########################################
# stats
########################################

# richness
gls1 <- gls(log_alpha_stab_prod ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_alpha_stab_prod ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(log_spp_stab_cover ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_spp_stab_cover ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(log_spp_asynch_cover ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_spp_asynch_cover ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(log_gamma_stab_prod ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_gamma_stab_prod ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

# abundance-based multiplicative beta
gls1 <- gls(log_spatial_asynch_prod ~ beta_invsimpsons * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_spatial_asynch_prod ~ beta_invsimpsons * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(log_gamma_stab_prod ~ beta_invsimpsons * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(log_gamma_stab_prod ~ beta_invsimpsons * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

########################################
# stats with env
########################################

# richness
gls1 <- gls(res_log_alpha_stab_prod ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_alpha_stab_prod ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(res_log_gamma_stab_prod ~ avg_richness2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_gamma_stab_prod ~ avg_richness2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

# abundance-based multiplicative beta
gls1 <- gls(res_log_spatial_asynch_prod ~ beta_invsimpsons * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_spatial_asynch_prod ~ beta_invsimpsons * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(res_log_gamma_stab_prod ~ beta_invsimpsons * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_gamma_stab_prod ~ beta_invsimpsons * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

########################################
# stats for other local diversity indices
########################################

# Simpson
gls1 <- gls(res_log_alpha_stab_prod ~ avg_invsimpsons2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_alpha_stab_prod ~ avg_invsimpsons2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(res_log_gamma_stab_prod ~ avg_invsimpsons2 * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_gamma_stab_prod ~ avg_invsimpsons2 * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

# Shannon
gls1 <- gls(res_log_alpha_stab_prod ~ avg_shannons * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_alpha_stab_prod ~ avg_shannons * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(res_log_gamma_stab_prod ~ avg_shannons * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_gamma_stab_prod ~ avg_shannons * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

# evenness
gls1 <- gls(res_log_alpha_stab_prod ~ avg_evenness_Sh * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_alpha_stab_prod ~ avg_evenness_Sh * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

gls1 <- gls(res_log_gamma_stab_prod ~ avg_evenness_Sh * year_trt * trt, method='ML', data = varpart_All_prod)
gls2 <- gls(res_log_gamma_stab_prod ~ avg_evenness_Sh * year_trt * trt, correlation = corAR1(form=~1), method='ML', data = varpart_All_prod)
anova(gls1, gls2)
anova(gls2)

