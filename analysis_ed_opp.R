
# analysis

library(openxlsx)
library(brms)
library(ggplot2)
library(lme4)
library(lmerTest)
library(sjPlot)

# ==== read in data, make sure everything matches ====

setwd('/Users/adamkaminski/Desktop/DRCMR/stimproj_followup/')
adhd_dat <- read.csv('stimulant_data_wide.csv')
edu_dat <- read.xlsx('ABCD_education.xlsx')

# merge dfs based on sub id
merge_dat <- merge(adhd_dat, edu_dat, by = "src_subject_id", all.x = TRUE)
# clean anything you must
merge_dat$sex <- as.factor(merge_dat$sex)
merge_dat$nonstim_med_bi <- as.factor(merge_dat$nonstim_med_bi)
merge_dat$stim_group_bi <- as.factor(merge_dat$stim_group_bi)

# matches data Hannah sent?
adhd_hannah_dat <- read.csv('merged_adhd_IDs_ABCD_slope_intercept.csv')
all(adhd_hannah_dat$src_subject_id==merge_dat$src_subject_id)
all(adhd_hannah_dat$eventname==merge_dat$eventname,na.rm = TRUE)
all(round(adhd_hannah_dat$led_sch_seda_s_mn_avg_eb,4)==round(merge_dat$led_sch_seda_s_mn_avg_eb,4),na.rm = TRUE)
all(round(adhd_hannah_dat$led_sch_seda_s_mn_coh_eb,4)==round(merge_dat$led_sch_seda_s_mn_coh_eb,4),na.rm = TRUE)
# yes it does!

# ==== exploring variables ====

ggplot(merge_dat, aes(x = led_sch_seda_s_mn_avg_eb, fill = stim_group_bi)) +
  scale_fill_manual(name = "Stimulant exposed", 
                    values = c("red", "blue"), 
                    labels = c("no", "yes")) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  xlab('Educational opportunity (intercept)')

t.test(led_sch_seda_s_mn_avg_eb ~ stim_group_bi, data = merge_dat)

cor.test(merge_dat$led_sch_seda_s_mn_avg_eb,merge_dat$ADHD.Problems.2)

merge_dat$ADHD.Problems.Change <- merge_dat$ADHD.Problems.2 - merge_dat$ADHD.Problems.0

tab_model(lmer(ADHD.Problems.2 ~ stim_group_bi*led_sch_seda_s_mn_avg_eb + ADHD.Problems.0 +
             sex.0 + interview_age.0 + ses + nonstim_med_bi + comorbs + 
               nihtbx_picvocab_agecorrected.0 + (1|site.0), data=merge_dat))

# ==== ANALYTIC STEP 1 ====
# ==== L caudate ====
slope.prior.all <- c(
  prior(normal(0,10), class=b, coef=led_sch_seda_s_mn_avg_eb),
  prior(normal(0,10), class=b, coef=meanFD.0),
  prior(normal(0,10), class=b, coef=meanFD.2),
  prior(normal(0,10), class=b, coef=FD.under.0.2.0_sqrt),
  prior(normal(0,10), class=b, coef=FD.under.0.2.2_sqrt),
  prior(normal(0,10), class=b, coef=ses),
  prior(normal(0,10), class=b, coef=sex1),
  #prior(normal(0,10), class=b, coef=stim_group_bi1),
  prior(normal(0,10), class=b, coef=nonstim_med_bi1),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.0),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.2),
  prior(normal(0,10), class=b, coef=comorbs),
  prior(normal(0,10), class=b, coef=ADHD.Problems.0),
  prior(normal(0,10), class=b, coef=ADHD.Problems.2),
  prior(normal(0,10), class=b, coef=au_cdelh),
  prior(normal(0,10), class=b, coef=cercsa_cdelh),
  prior(normal(0,10), class=b, coef=copa_cdelh),
  prior(normal(0,10), class=b, coef=df_cdelh),
  prior(normal(0,10), class=b, coef=dsa_cdelh),
  prior(normal(0,10), class=b, coef=fopa_cdelh),
  prior(normal(0,10), class=b, coef=rst_cdelh),
  prior(normal(0,10), class=b, coef=smhsmm_cdelh),
  prior(normal(0,10), class=b, coef=vta_cdelh),
  prior(normal(0,10), class=b, coef=vs_cdelh)
)

cdelh_step_1<-brm(as.numeric(stim_group_bi) ~ led_sch_seda_s_mn_avg_eb + meanFD.0 + meanFD.2 +
              FD.under.0.2.0_sqrt + FD.under.0.2.2_sqrt + ses + sex + nonstim_med_bi + nihtbx_picvocab_agecorrected.0 +
              nihtbx_picvocab_agecorrected.2 + ADHD.Problems.0 + ADHD.Problems.2 + comorbs + au_cdelh + cercsa_cdelh + 
              copa_cdelh + df_cdelh + dsa_cdelh + fopa_cdelh + rst_cdelh + smhsmm_cdelh + vta_cdelh + vs_cdelh +
              (1|site.0), data = merge_dat, family="bernoulli", prior = slope.prior.all, silent = FALSE,
              warmup = 2000, iter = 10000)

# ==== R caudate ====
slope.prior.all <- c(
  prior(normal(0,10), class=b, coef=led_sch_seda_s_mn_avg_eb),
  prior(normal(0,10), class=b, coef=meanFD.0),
  prior(normal(0,10), class=b, coef=meanFD.2),
  prior(normal(0,10), class=b, coef=FD.under.0.2.0_sqrt),
  prior(normal(0,10), class=b, coef=FD.under.0.2.2_sqrt),
  prior(normal(0,10), class=b, coef=ses),
  prior(normal(0,10), class=b, coef=sex1),
  #prior(normal(0,10), class=b, coef=stim_group_bi1),
  prior(normal(0,10), class=b, coef=nonstim_med_bi1),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.0),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.2),
  prior(normal(0,10), class=b, coef=comorbs),
  prior(normal(0,10), class=b, coef=ADHD.Problems.0),
  prior(normal(0,10), class=b, coef=ADHD.Problems.2),
  prior(normal(0,10), class=b, coef=au_cderh),
  prior(normal(0,10), class=b, coef=cercsa_cderh),
  prior(normal(0,10), class=b, coef=copa_cderh),
  prior(normal(0,10), class=b, coef=df_cderh),
  prior(normal(0,10), class=b, coef=dsa_cderh),
  prior(normal(0,10), class=b, coef=fopa_cderh),
  prior(normal(0,10), class=b, coef=rst_cderh),
  prior(normal(0,10), class=b, coef=smhsmm_cderh),
  prior(normal(0,10), class=b, coef=vta_cderh),
  prior(normal(0,10), class=b, coef=vs_cderh)
)

cderh_step_1<-brm(as.numeric(stim_group_bi) ~ led_sch_seda_s_mn_avg_eb + meanFD.0 + meanFD.2 +
                    FD.under.0.2.0_sqrt + FD.under.0.2.2_sqrt + ses + sex + nonstim_med_bi + nihtbx_picvocab_agecorrected.0 +
                    nihtbx_picvocab_agecorrected.2 + ADHD.Problems.0 + ADHD.Problems.2 + comorbs + au_cderh + cercsa_cderh + 
                    copa_cderh + df_cderh + dsa_cderh + fopa_cderh + rst_cderh + smhsmm_cderh + vta_cderh + vs_cderh +
                    (1|site.0), data = merge_dat, family="bernoulli", prior = slope.prior.all, silent = FALSE,
                  warmup = 2000, iter = 10000)
# ==== L putamen ====
slope.prior.all <- c(
  prior(normal(0,10), class=b, coef=led_sch_seda_s_mn_avg_eb),
  prior(normal(0,10), class=b, coef=meanFD.0),
  prior(normal(0,10), class=b, coef=meanFD.2),
  prior(normal(0,10), class=b, coef=FD.under.0.2.0_sqrt),
  prior(normal(0,10), class=b, coef=FD.under.0.2.2_sqrt),
  prior(normal(0,10), class=b, coef=ses),
  prior(normal(0,10), class=b, coef=sex1),
  #prior(normal(0,10), class=b, coef=stim_group_bi1),
  prior(normal(0,10), class=b, coef=nonstim_med_bi1),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.0),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.2),
  prior(normal(0,10), class=b, coef=comorbs),
  prior(normal(0,10), class=b, coef=ADHD.Problems.0),
  prior(normal(0,10), class=b, coef=ADHD.Problems.2),
  prior(normal(0,10), class=b, coef=au_ptlh),
  prior(normal(0,10), class=b, coef=cercsa_ptlh),
  prior(normal(0,10), class=b, coef=copa_ptlh),
  prior(normal(0,10), class=b, coef=df_ptlh),
  prior(normal(0,10), class=b, coef=dsa_ptlh),
  prior(normal(0,10), class=b, coef=fopa_ptlh),
  prior(normal(0,10), class=b, coef=rst_ptlh),
  prior(normal(0,10), class=b, coef=smhsmm_ptlh),
  prior(normal(0,10), class=b, coef=vta_ptlh),
  prior(normal(0,10), class=b, coef=vs_ptlh)
)

ptlh_step_1<-brm(as.numeric(stim_group_bi) ~ led_sch_seda_s_mn_avg_eb + meanFD.0 + meanFD.2 +
                    FD.under.0.2.0_sqrt + FD.under.0.2.2_sqrt + ses + sex + nonstim_med_bi + nihtbx_picvocab_agecorrected.0 +
                    nihtbx_picvocab_agecorrected.2 + ADHD.Problems.0 + ADHD.Problems.2 + comorbs + au_ptlh + cercsa_ptlh + 
                    copa_ptlh + df_ptlh + dsa_ptlh + fopa_ptlh + rst_ptlh + smhsmm_ptlh + vta_ptlh + vs_ptlh +
                    (1|site.0), data = merge_dat, family="bernoulli", prior = slope.prior.all, silent = FALSE,
                  warmup = 2000, iter = 10000)
# ==== R putamen ====
slope.prior.all <- c(
  prior(normal(0,10), class=b, coef=led_sch_seda_s_mn_avg_eb),
  prior(normal(0,10), class=b, coef=meanFD.0),
  prior(normal(0,10), class=b, coef=meanFD.2),
  prior(normal(0,10), class=b, coef=FD.under.0.2.0_sqrt),
  prior(normal(0,10), class=b, coef=FD.under.0.2.2_sqrt),
  prior(normal(0,10), class=b, coef=ses),
  prior(normal(0,10), class=b, coef=sex1),
  #prior(normal(0,10), class=b, coef=stim_group_bi1),
  prior(normal(0,10), class=b, coef=nonstim_med_bi1),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.0),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.2),
  prior(normal(0,10), class=b, coef=comorbs),
  prior(normal(0,10), class=b, coef=ADHD.Problems.0),
  prior(normal(0,10), class=b, coef=ADHD.Problems.2),
  prior(normal(0,10), class=b, coef=au_ptrh),
  prior(normal(0,10), class=b, coef=cercsa_ptrh),
  prior(normal(0,10), class=b, coef=copa_ptrh),
  prior(normal(0,10), class=b, coef=df_ptrh),
  prior(normal(0,10), class=b, coef=dsa_ptrh),
  prior(normal(0,10), class=b, coef=fopa_ptrh),
  prior(normal(0,10), class=b, coef=rst_ptrh),
  prior(normal(0,10), class=b, coef=smhsmm_ptrh),
  prior(normal(0,10), class=b, coef=vta_ptrh),
  prior(normal(0,10), class=b, coef=vs_ptrh)
)

ptrh_step_1<-brm(as.numeric(stim_group_bi) ~ led_sch_seda_s_mn_avg_eb + meanFD.0 + meanFD.2 +
                    FD.under.0.2.0_sqrt + FD.under.0.2.2_sqrt + ses + sex + nonstim_med_bi + nihtbx_picvocab_agecorrected.0 +
                    nihtbx_picvocab_agecorrected.2 + ADHD.Problems.0 + ADHD.Problems.2 + comorbs + au_ptrh + cercsa_ptrh + 
                    copa_ptrh + df_ptrh + dsa_ptrh + fopa_ptrh + rst_ptrh + smhsmm_ptrh + vta_ptrh + vs_ptrh +
                    (1|site.0), data = merge_dat, family="bernoulli", prior = slope.prior.all, silent = FALSE,
                  warmup = 2000, iter = 10000)
# ==== L NAc ====
slope.prior.all <- c(
  prior(normal(0,10), class=b, coef=led_sch_seda_s_mn_avg_eb),
  prior(normal(0,10), class=b, coef=meanFD.0),
  prior(normal(0,10), class=b, coef=meanFD.2),
  prior(normal(0,10), class=b, coef=FD.under.0.2.0_sqrt),
  prior(normal(0,10), class=b, coef=FD.under.0.2.2_sqrt),
  prior(normal(0,10), class=b, coef=ses),
  prior(normal(0,10), class=b, coef=sex1),
  #prior(normal(0,10), class=b, coef=stim_group_bi1),
  prior(normal(0,10), class=b, coef=nonstim_med_bi1),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.0),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.2),
  prior(normal(0,10), class=b, coef=comorbs),
  prior(normal(0,10), class=b, coef=ADHD.Problems.0),
  prior(normal(0,10), class=b, coef=ADHD.Problems.2),
  prior(normal(0,10), class=b, coef=au_aalh),
  prior(normal(0,10), class=b, coef=cercsa_aalh),
  prior(normal(0,10), class=b, coef=copa_aalh),
  prior(normal(0,10), class=b, coef=df_aalh),
  prior(normal(0,10), class=b, coef=dsa_aalh),
  prior(normal(0,10), class=b, coef=fopa_aalh),
  prior(normal(0,10), class=b, coef=rst_aalh),
  prior(normal(0,10), class=b, coef=smhsmm_aalh),
  prior(normal(0,10), class=b, coef=vta_aalh),
  prior(normal(0,10), class=b, coef=vs_aalh)
)

aalh_step_1<-brm(as.numeric(stim_group_bi) ~ led_sch_seda_s_mn_avg_eb + meanFD.0 + meanFD.2 +
                   FD.under.0.2.0_sqrt + FD.under.0.2.2_sqrt + ses + sex + nonstim_med_bi + nihtbx_picvocab_agecorrected.0 +
                   nihtbx_picvocab_agecorrected.2 + ADHD.Problems.0 + ADHD.Problems.2 + comorbs + au_aalh + cercsa_aalh + 
                   copa_aalh + df_aalh + dsa_aalh + fopa_aalh + rst_aalh + smhsmm_aalh + vta_aalh + vs_aalh +
                   (1|site.0), data = merge_dat, family="bernoulli", prior = slope.prior.all, silent = FALSE,
                 warmup = 2000, iter = 10000)

# ==== R NAc ====
slope.prior.all <- c(
  prior(normal(0,10), class=b, coef=led_sch_seda_s_mn_avg_eb),
  prior(normal(0,10), class=b, coef=meanFD.0),
  prior(normal(0,10), class=b, coef=meanFD.2),
  prior(normal(0,10), class=b, coef=FD.under.0.2.0_sqrt),
  prior(normal(0,10), class=b, coef=FD.under.0.2.2_sqrt),
  prior(normal(0,10), class=b, coef=ses),
  prior(normal(0,10), class=b, coef=sex1),
  #prior(normal(0,10), class=b, coef=stim_group_bi1),
  prior(normal(0,10), class=b, coef=nonstim_med_bi1),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.0),
  prior(normal(0,10), class=b, coef=nihtbx_picvocab_agecorrected.2),
  prior(normal(0,10), class=b, coef=comorbs),
  prior(normal(0,10), class=b, coef=ADHD.Problems.0),
  prior(normal(0,10), class=b, coef=ADHD.Problems.2),
  prior(normal(0,10), class=b, coef=au_aarh),
  prior(normal(0,10), class=b, coef=cercsa_aarh),
  prior(normal(0,10), class=b, coef=copa_aarh),
  prior(normal(0,10), class=b, coef=df_aarh),
  prior(normal(0,10), class=b, coef=dsa_aarh),
  prior(normal(0,10), class=b, coef=fopa_aarh),
  prior(normal(0,10), class=b, coef=rst_aarh),
  prior(normal(0,10), class=b, coef=smhsmm_aarh),
  prior(normal(0,10), class=b, coef=vta_aarh),
  prior(normal(0,10), class=b, coef=vs_aarh)
)

aarh_step_1<-brm(as.numeric(stim_group_bi) ~ led_sch_seda_s_mn_avg_eb + meanFD.0 + meanFD.2 +
                   FD.under.0.2.0_sqrt + FD.under.0.2.2_sqrt + ses + sex + nonstim_med_bi + nihtbx_picvocab_agecorrected.0 +
                   nihtbx_picvocab_agecorrected.2 + ADHD.Problems.0 + ADHD.Problems.2 + comorbs + au_aarh + cercsa_aarh + 
                   copa_aarh + df_aarh + dsa_aarh + fopa_aarh + rst_aarh + smhsmm_aarh + vta_aarh + vs_aarh +
                   (1|site.0), data = merge_dat, family="bernoulli", prior = slope.prior.all, silent = FALSE,
                 warmup = 2000, iter = 10000)

# ==== analytic step 2 ====
# ==== analytic step 3 ====

