#Supplementary Code for
#Fear of the Human “Super Predator” Restricts Elephant Ecosystem Engineering and Crop Damage      
#Authors: Michael B. Kowalski, Anderson Larpei, Michael Clinchy, Justin Suraci, Liana Y. Zanette, Christopher C. Wilmers
#Corresponding author: mbkowals@ucsc.edu

##########
#Libraries
##########
library(tidyverse); library(activity); library(R2jags); library(abind);
library(glmmTMB); library(emmeans); library(patchwork); library(DHARMa); 
library(performance); library(see)

##########
#Data
##########
detections <- read.csv("DataS2.csv")
trees <- read.csv("DataS1.csv")
conflict <- read.csv("DataS3.csv", colClasses = c(Date = "character"))
farm.effort <- read.csv("DataS4.csv")
anthro <- read.csv("DataS5.csv")
load(file = "DataS7.rda", verbose = T) 
load("DataS8.rda")
set.seed(1)

###########################  MRC Data Analyses  ###############################
##########
#Human Presence in Experimental Grids
##########
#Mann-Whitney U-Test
anthro_camera <- anthro %>%
  group_by(grid_id, camera_id) %>%
  summarise(total_occurrences = n(), .groups = "drop") %>%
  mutate( camera_nights = 10 * 7,
          rate_per_night = total_occurrences / camera_nights)

wilcox.test(
  rate_per_night ~ grid_id,
  data  = anthro_camera,
  exact = FALSE)

##########
#Bayesian Hierarchical Model
##########
#Elephant presence with spatial Gaussian process
sink("MFP_Elephant_GausProc.jags")
cat("model{

    # Prior for gamma distribution (rho)
    theta <- exp(logtheta)
    logtheta ~ dunif(-5,5)
    
    ############
    # Spatial Parameters for Gaussian Process
    # Range parameters
    rho.psi ~ dunif(0.1, 3)     
    rho.lambda ~ dunif(0.1, 3)  
            
    # Marginal variance parameters
    sigma2.psi ~ dunif(0.01, 3)     
    sigma2.lambda ~ dunif(0.01, 3)  
    
    ############
    # Priors for probability of occurrence parameters (psi)
    a0 ~ dnorm(0, 0.001)      #intercept
    a1 ~ dnorm(0, 0.001)      #site effect
    a2 ~ dnorm(0, 0.001)      #treatment effect
    a3 ~ dnorm(0, 0.001)      #site x treatment interaction
    
    ############
    # Priors for intensity of site use parameters (lambda)
    b0 ~ dnorm(0, 0.001)      #intercept
    b1 ~ dnorm(0, 0.001)      #site effect
    b2 ~ dnorm(0, 0.001)      #treatment effect
    b3 ~ dnorm(0, 0.001)      #site x treatment interaction
    
    ############
    # Covariance matrices for spatial effects
    for(j in 1:M){
        for(k in 1:M){
            Sigma.psi[j,k] <- sigma2.psi *
                              exp(-pow(D[j,k],2)/(2*pow(rho.psi,2))) +
                              1.0E-6 * equals(j,k)
            Sigma.lambda[j,k] <- sigma2.lambda *
                                exp(-pow(D[j,k],2)/(2*pow(rho.lambda,2)))+
                                1.0E-6 * equals(j,k)
        }
    }
    Omega.psi[1:M,1:M] <- inverse(Sigma.psi[,])
    Omega.lambda[1:M,1:M] <- inverse(Sigma.lambda[,])
    
    ############
    # Spatial random effects
    phi.psi[1:M] ~ dmnorm(zeros[], Omega.psi[,])
    phi.lambda[1:M] ~ dmnorm(zeros[], Omega.lambda[,])
    
    # Vector of zeros for MVN mean
    for(i in 1:M){
        zeros[i] <- 0
    }
    
    ############
    # Couples each k survey period to the appropriate treatment period
    t.track<-c(1,1,1,1,1,2,2,2,2,2)  
    
    ############
    # Logistic regression submodel for probability of occurence (psi)
    for(j in 1:M){                #loop through camera traps (1 to 24)         
        for(t in 1:2){            #loop through treatment periods (1 to 2)
            z[j,t] ~ dbern(psi[j,t])
            
            logit(psi[j,t]) <- b0 + b1*Trt[j,t] + b2*S.psi[j,t] + 
            b3*Trt[j,t]*S.psi[j,t] + phi.psi[j]
                
            loglik.psi[j,t] <- z[j,t] * log(psi[j,t] + 1.0E-10) + 
                               (1 - z[j,t]) * log(1 - psi[j,t] + 1.0E-10)
        }
    }
    
    ############
    # Negative binomial submodel for intensity of site use (lambda)
    for(j in 1:M){                  #loop through camera traps (1 to 24)
        for(k in 1:K){              #loop through survey periods (1 to 10)
            log(lambda[j,k]) <- a0 + a1*Trt[j,k] + a2*S[j,k] + 
            a3*Trt[j,k]*S[j,k] + phi.lambda[j]
            
            rho[j,k] ~ dgamma(theta, theta)
            
            mu[j,k] <- rho[j,k]*lambda[j,k]
            
            y[j,k] ~ dpois(z[j,t.track[k]] * mu[j,k])
            
            loglik.lambda[j,k] <- logdensity.pois(y[j,k], 
            z[j,t.track[k]] * mu[j,k])
    
            # Create new data for calculating Bayesian p-values
            y_new[j,k] ~ dpois(z[j,t.track[k]] * mu[j,k])
            
            # Calculate statistics for Bayesian p-values
            eval[j,k] <- z[j,t.track[k]] * mu[j,k]
            
            # Freeman-Tukey Residual
            Terr[j,k] <- pow(pow(y[j,k],.5) - pow(eval[j,k],.5),2)
            Terrnew[j,k] <- pow(pow(y_new[j,k],.5) - pow(eval[j,k],.5),2)
            
            # Chi-squared stat
            ch.err[j,k] <- pow((y[j,k] - eval[j,k]),2)/ (eval[j,k] + 0.5)
            ch.errnew[j,k] <- pow((y_new[j,k] - eval[j,k]),2)/ (eval[j,k]         
            + 0.5)
        }
    }
    
    # Posterior predictive checks
    Tobs <- sum(Terr[,])
    Tnew <- sum(Terrnew[,])
    Chisq.obs <- sum(ch.err[,])
    Chisq.new <- sum(ch.errnew[,])
           
    ############
    # Derived parameters
    
    # Sum within experimental sites for both treatments
    for(t in 1:2){
        z.sum.S2[t] <- sum(z[1:12,t])
        z.sum.N2[t] <- sum(z[13:24,t])
    }
    
    # Calculate average across sites for each treatment
    z.sum.C <- mean(c(z.sum.S2[2], z.sum.N2[1]))
    z.sum.H <- mean(c(z.sum.S2[1], z.sum.N2[2]))
    
    # Get average value of rho
    rho.ave <- mean(rho[,])
    
    # Estimate mean detection rate (lambda)
    det.N2.c <- mean(lambda[13:24,1:5])
    det.N2.h <- mean(lambda[13:24,6:10])
    det.S2.h <- mean(lambda[1:12,1:5])
    det.S2.c <- mean(lambda[1:12,6:10])

    # Calculate average across sites for each treatment
    det.C <- mean(c(det.N2.c, det.S2.c))
    det.H <- mean(c(det.N2.h, det.S2.h))
    
    # Derived parameters for spatial correlation assessment
    eff.range.psi <- rho.psi * sqrt(-2 * log(0.05))
    eff.range.lambda <- rho.lambda * sqrt(-2 * log(0.05))
}",fill = TRUE)
sink()

# Array dimensions
M <- nrow(ele.det)  # Number of camera sites (24)
K <- ncol(ele.det)  # Number of survey periods (10)

# Parameters to monitor
parameters <- c("a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3",
                "Tobs", "Tnew", "Chisq.obs", "Chisq.new",
                "det.N2.c", "det.N2.h", "det.S2.h", "det.S2.c",
                "z.sum.C", "z.sum.H", "det.C", "det.H",
                "rho.psi", "rho.lambda", "sigma2.psi", "sigma2.lambda",
                "eff.range.psi", "eff.range.lambda")

# MCMC settings
ni <- 300000 # iterations
nb <- 250000 # burn-in
nthin <- 50 # thinning
nc <- 3 # number of chains

# Data list
data <- list(
  y    = ele.det,
  K    = K,
  M    = M,
  Trt  = T.hc,
  S    = S.hc,
  S.psi = S.psi,
  D    = D
)

# Initial values for z (occupancy states)
z1  <- apply(ele.det[, 1:5],  1, sum); z1 <- ifelse(z1 > 0, 1, 0)
z2  <- apply(ele.det[, 6:10], 1, sum); z2 <- ifelse(z2 > 0, 1, 0)
zst <- cbind(z1, z2)

inits <- list(
  list(a0 = runif(1,-1,1), a1 = runif(1,-1,1), a2 = runif(1,-1,1), a3 = runif(1,-1,1),
       b0 = runif(1,-0.5,0.5), b1 = runif(1,-0.5,0.5), b2 = runif(1,-0.5,0.5), b3 = runif(1,-0.5,0.5),
       rho.psi = runif(1,0.5,1.2), rho.lambda = runif(1,0.5,1.2),
       sigma2.psi = runif(1,0.4,1.5), sigma2.lambda = runif(1,0.4,1.5),
       logtheta = runif(1,-2,2), z = zst),
  list(a0 = runif(1,-1,1), a1 = runif(1,-1,1), a2 = runif(1,-1,1), a3 = runif(1,-1,1),
       b0 = runif(1,-0.5,0.5), b1 = runif(1,-0.5,0.5), b2 = runif(1,-0.5,0.5), b3 = runif(1,-0.5,0.5),
       rho.psi = runif(1,0.5,1.2), rho.lambda = runif(1,0.5,1.2),
       sigma2.psi = runif(1,0.4,1.5), sigma2.lambda = runif(1,0.4,1.5),
       logtheta = runif(1,-2,2), z = zst),
  list(a0 = runif(1,-1,1), a1 = runif(1,-1,1), a2 = runif(1,-1,1), a3 = runif(1,-1,1),
       b0 = runif(1,-0.5,0.5), b1 = runif(1,-0.5,0.5), b2 = runif(1,-0.5,0.5), b3 = runif(1,-0.5,0.5),
       rho.psi = runif(1,0.5,1.2), rho.lambda = runif(1,0.5,1.2),
       sigma2.psi = runif(1,0.4,1.5), sigma2.lambda = runif(1,0.4,1.5),
       logtheta = runif(1,-2,2), z = zst)
)

# Run model
det.out.wZT6.ele <- jags(data = data, inits = inits, parameters.to.save = parameters,
                         model.file = "MFP_Elephant_GausProc.jags",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nthin,
                         progress.bar = "text")

#To extract posterior summaries for derived parameters
det.out.wZT6.ele$BUGSoutput$summary[
  c("det.C", "det.H", "z.sum.C", "z.sum.H",
    "eff.range.psi", "eff.range.lambda",
    "Tobs", "Tnew", "Chisq.obs", "Chisq.new"),]

##########
#Diel Shift Analysis
##########
#Probability of nocturnal detections on experimental grids
elediurnal <- detections %>%
  mutate(treatment  = factor(treatment, levels = c("control","human")),
         daylight   = as.integer(daylight),
         camera_num = factor(camera_num))

m_elenocturnal <- glmmTMB(I(1 - daylight) ~ treatment + grid_id + (1 | camera_num),
                          family = binomial(link = "logit"),
                          data   = elediurnal)

summary(m_elenocturnal)
emmeans(m_elenocturnal, ~ treatment, type = "response")

##########
#Elephant Tree Damage Analyses
##########
#tree damage frequency in experimental grids
trees_f <- trees %>% mutate(
  treatment = factor(treatment),
  grid      = factor(grid),
  transect  = factor(transect)) %>%
  count(treatment, grid, transect, name = "count") %>%
  complete(treatment, grid, transect, fill = list(count = 0)) %>%
  filter(substr(transect, 1, 2) == grid)

m_forage_f <- glmmTMB(
  count ~ treatment + (1 | transect),
  family = nbinom2(),
  data = trees_f)

summary(m_forage_f)
emmeans(m_forage_f, ~ treatment, type = "response")

#Tree damage percent in experimental grids
trees_d <- trees %>% filter(percent>=5)
trees_d <- trees_d %>%
  group_by(grid, transect, treatment) %>%
  summarise(
    mean_percent = mean(percent),
    n_trees = n(),
    .groups = "drop")

m_forage_d <- glmmTMB(
  log(mean_percent) ~ treatment * grid + (1 | transect),
  family = gaussian(),
  data = trees_d)

summary(m_forage_d)
emmeans(m_forage_d, ~ treatment, type = "response")

# Hierarchical bootstrap: tree damage proportionality
# Tests whether elephant tree damage declined more than expected from the 53%
# drop in elephant detections under the human treatment in experimental grids.
# D̂ < 0 → tree damage suppressed beyond proportional expectation.
#
# D̂ = (1/G) Σ_g [ log(D_H,g / D_C,g) - log(N_H,g / N_C,g) ]
#   where D = tree counts/transect, N = detections/camera-week
#
# Weekly elephant detections per camera
ele_dw <- detections %>% mutate(
    grid      = factor(grid_id),
    treatment = factor(treatment, levels = c("control", "human")),
    week      = factor(week_id),
    camera    = factor(camera_id)) %>%
  count(grid, treatment, camera, week, name = "det_week") %>%
  complete(grid, treatment, camera, week, fill = list(det_week = 0))

# Elephant-damaged tree counts per transect 
trees_tr <- trees %>% mutate(
    grid      = factor(grid),
    treatment = factor(treatment, levels = c("control", "human")),
    transect  = factor(transect)) %>%
  count(grid, treatment, transect, name = "n_tree") %>%
  complete(grid, treatment, transect, fill = list(n_tree = 0)) %>%
  filter(substr(as.character(transect), 1, 2) == as.character(grid))

#Point estimate: D̂ on observed data
compute_D <- function(eps = 0.5) {
  det_sum <- ele_dw %>%
    group_by(grid, treatment) %>%
    summarise(N = sum(det_week), .groups = "drop")
  tree_sum <- trees_tr %>%
    group_by(grid, treatment) %>%
    summarise(D = sum(n_tree), .groups = "drop")
  dat <- full_join(det_sum, tree_sum, by = c("grid", "treatment")) %>%
    mutate(N = N + eps, D = D + eps) %>%
    pivot_wider(names_from = treatment, values_from = c(N, D))
  D_grid <- with(dat, log(D_human / D_control) - log(N_human / N_control))
  mean(D_grid, na.rm = TRUE)}

#Hierarchical bootstrap
boot_D <- function(eps = 0.5) {
  # Resample cameras within each grid
  cams_b <- ele_dw %>%
    distinct(grid, camera) %>%
    group_by(grid) %>%
    reframe(camera = sample(camera, size = n(), replace = TRUE))
  
  det_sum <- ele_dw %>%
    inner_join(
      cams_b %>% count(grid, camera, name = "w"),
      by = c("grid", "camera")) %>%
    group_by(grid, treatment) %>%
    summarise(N = sum(det_week * w), .groups = "drop")
  
  # Resample transects within each grid
  trs_b <- trees_tr %>%
    distinct(grid, transect) %>%
    group_by(grid) %>%
    reframe(transect = sample(transect, size = n(), replace = TRUE))
  
  tree_sum <- trees_tr %>%
    inner_join(
      trs_b %>% count(grid, transect, name = "w"),
      by = c("grid", "transect")) %>%
    group_by(grid, treatment) %>%
    summarise(D = sum(n_tree * w), .groups = "drop")
  
  dat <- full_join(det_sum, tree_sum, by = c("grid", "treatment")) %>%
    mutate(N = N + eps, D = D + eps) %>%
    pivot_wider(names_from = treatment, values_from = c(N, D))
  
  D_grid <- with(dat, log(D_human / D_control) - log(N_human / N_control))
  mean(D_grid, na.rm = TRUE)
}

D_boot <- replicate(10000, boot_D())

#Bootstrap Outputs
D_hat <- compute_D(); D_hat
quantile(D_boot, c(0.025, 0.975)) 
mean(D_boot >= 0)                  

#############
#Abatement
#############
#Week-over-week change in elephant detections in experimental grids
mrc_ele_det <- detections %>%
  group_by(treatment, grid_id,week_id) %>% summarise(detections = n(), .groups = "drop")
full_weeks <- expand.grid(treatment = c("control", "human"),grid_id   = c("N2", "S2"), week_id   = 1:10)
mrc_ele_det <- full_weeks %>%
  left_join(mrc_ele_det, by = c("treatment", "grid_id", "week_id")) %>%
  mutate( detections = tidyr::replace_na(detections, 0))

mrc_ele_abt <- glmmTMB(
  detections ~ treatment + week_id + grid_id, family = nbinom2, data = mrc_ele_det)

summary(mrc_ele_abt)

###########################  Farm Data Analyses  ###############################
##########
#Generalized Linear Mixed-Effects Models
##########
#Elephant presence on farms
eledetect <- conflict %>%
  filter(Ele.present == "Yes") %>%
  count(Farm, Treatment, Cycle, name = "n_conflicts") %>%
  left_join(farm.effort, by = c("Farm" = "FarmID")) %>%
  mutate(
    nights = case_when(
      Treatment == "Humans"   & Cycle == 1 ~ night_human_c1,
      Treatment == "Humans"   & Cycle == 2 ~ night_human_c2,
      Treatment == "Crickets" & Cycle == 1 ~ night_control_c1,
      Treatment == "Crickets" & Cycle == 2 ~ night_control_c2
    ),
    Cycle = factor(Cycle, levels = c(1, 2))
  )

ele_rate <- glmmTMB(
  n_conflicts ~ Treatment + Cycle + offset(log(nights)) + (1 | Farm),
  family = poisson,
  data = eledetect
)

summary(ele_rate)
emmeans(ele_rate, ~ Treatment, type = "response", offset = log(30))

#Elephant proximity to speakers
eleprox <- glmmTMB(
  log(Proximity.m.) ~ Treatment + Cycle + (1 | Farm),
  family = gaussian(),
  data   = conflict %>% filter(!is.na(Proximity.m.)))

summary(eleprox)
emmeans(eleprox, ~ Treatment, type = "response")

#Elephant crop raiding probability
eleraid <- conflict %>%
  filter(Ele.present == "Yes") %>%
  mutate(conflict = as.integer(Cropdamaged == "Yes"),
         Cycle   = factor(Cycle, levels = c(1, 2)))

eleraidprob <- glmmTMB(
  conflict ~ Treatment + Cycle + (1 | Farm),
  family = binomial(link = "logit"),
  data   = eleraid)

summary(eleraidprob)
emmeans(eleraidprob, ~ Treatment, type = "response")

#Total area of damaged crops per farm
eledamage <- conflict %>%
  group_by(Farm, Treatment, Cycle) %>%
  summarise(sum_damage = sum(Area.m2., na.rm = TRUE), .groups = "drop") %>%
  left_join(farm.effort, by = c("Farm" = "FarmID")) %>%
  mutate(nights = case_when(
    Treatment == "Humans"   & Cycle == 1 ~ night_human_c1,
    Treatment == "Humans"   & Cycle == 2 ~ night_human_c2,
    Treatment == "Crickets" & Cycle == 1 ~ night_control_c1,
    Treatment == "Crickets" & Cycle == 2 ~ night_control_c2),
    Cycle = factor(Cycle, levels = c(1, 2)))

eledamrate <- glmmTMB(
  sum_damage ~ Treatment + Cycle + offset(log(nights)) + (1 | Farm),
  family = nbinom2(link = "log"),
  data   = eledamage)

summary(eledamrate)
emm_eledam <- emmeans(eledamrate, ~ Treatment, type = "response", offset = log(30))
emm_eledam 

#############
#Abatement
#############
#Night-over-night changes in elephant detection rate
agc_abtdetect <- conflict %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%y"), Farm = factor(Farm), Treatment = factor(Treatment), Cycle = factor(Cycle)) %>%
  arrange(Farm, Date)
agc_abtdetect <- agc_abtdetect%>% group_by(Farm) %>%
  mutate(exposure_night = row_number()) %>% ungroup()
agc_abtdetect <- agc_abtdetect %>%
  mutate(ele_present_bin = ifelse(Ele.present == "Yes", 1, 0))

m_abtdetect <- glmmTMB(ele_present_bin ~ exposure_night + Treatment + Cycle + (1 | Farm),
                       family = binomial(link = "logit"), data = agc_abtdetect)

summary(m_abtdetect)

#Night-over-night changes in elephant damage rate
agc_abtdamage <- agc_abtdetect %>%
  mutate(Area_m2 = as.numeric(Area.m2.)) %>% filter(Area_m2 > 0)

m_abtdamage <- glmmTMB(
  Area_m2 ~ exposure_night + Treatment + Cycle + (1 | Farm),
  family = Gamma(link = "log"),  data = agc_abtdamage)

summary(m_abtdamage)

################################## Plots ######################################
make_plot <- function(df, ylab, show_x = FALSE) {
  ggplot(df, aes(x = xlab, y = mean, fill = xlab)) +
    geom_col(width = 0.7, color = "black", linewidth = 1.5) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 2.5) +
    scale_fill_manual(values = fill_vals) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 52, base_family = "Palatino") +
    theme(
      legend.position   = "none",
      axis.text.x       = if (show_x) element_text(size = 46, lineheight = 0.9) else element_blank(),
      axis.ticks.x      = if (show_x) element_line(linewidth = 1.5) else element_blank(),
      axis.text.y       = element_text(size = 46),
      axis.title.y      = element_text(size = 52, margin = margin(r = 20)),
      axis.line         = element_line(linewidth = 1.5),
      axis.ticks        = element_line(linewidth = 1.5),
      axis.ticks.length = unit(0.3, "cm"))
}

text_block <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             family = "Palatino", size = 18, fontface = "bold") +
    theme_void()
}

damage_f      <- tibble(xlab = factor(lvls, lvls),
                        mean = c(26.80,  5.41),        lo = c(16.24,  3.09),        hi = c(44.24,   9.46))
ele_cropdam   <- tibble(xlab = factor(lvls, lvls),
                        mean = c(268.6,  49.5),         lo = c(136.1,  25.5),         hi = c(530.3,  96.1))
det_grid      <- tibble(xlab = factor(lvls, lvls),
                        mean = c(1.431, 0.664) * 12,   lo = c(1.110, 0.486) * 12,   hi = c(1.808, 0.898) * 12)
eledetect_emm <- tibble(xlab = factor(lvls, lvls),
                        mean = c(7.29,   1.79),         lo = c(5.19,   1.12),         hi = c(10.24,   2.86))
damage_d      <- tibble(xlab = factor(lvls, lvls),
                        mean = c(19.97,  6.96),         lo = c(15.54,  5.02),         hi = c(25.65,   9.64))
ele_conv_emm  <- tibble(xlab = factor(lvls, lvls),
                        mean = c(65.0,  21.8),          lo = c(55.1,  13.0),          hi = c(73.8,   34.3))
elenocturnal  <- tibble(xlab = factor(lvls, lvls),
                        mean = c(27.9,  72.0),          lo = c(20.0,  59.3),          hi = c(37.3,   82.0))
eleprox_emm   <- tibble(xlab = factor(lvls, lvls),
                        mean = c(8.2,   78.2),          lo = c(5.73,  49.81),         hi = c(11.7,  122.6))

############################ Manuscript Figures ###############################
fill_vals <- c("Control" = "#2C7BB6", "Human" = "#D7191C")
lvls      <- c("Control", "Human")

#Figure 1
F1A <- make_plot(damage_f,      "Trees damaged\nper transect on grids")
F1B <- make_plot(ele_cropdam,   "Monthly crop damage\nper farm (m²)")
F1C <- make_plot(det_grid,      "Elephant detections\nper week on grids",  show_x = TRUE)
F1D <- make_plot(eledetect_emm, "Elephant detections\nper month on farms", show_x = TRUE)

F1 <- (text_block("Reserve") | text_block("Farms")) /
      (F1A | F1B) /
      (F1C | F1D) /
      text_block("Playback Treatment") +
      plot_layout(heights = c(0.12, 1, 1, 0.12))
F1

#Figure 2
F2A <- make_plot(damage_d,     "Damage per tree\non grids (percent)",              show_x = TRUE)
F2B <- make_plot(ele_conv_emm, "Crop raided if elephant\non farm (% probability)", show_x = TRUE)

F2 <- (text_block("Reserve") | text_block("Farms")) /
  (F2A | F2B) /
  text_block("Playback Treatment") +
  plot_layout(heights = c(0.12, 1, 0.12))
F2

#Figure 3
F3A <- make_plot(elenocturnal, "Nighttime detections\non grids (percent)", show_x = TRUE)
F3B <- make_plot(eleprox_emm,  "Distance to speaker\non farms (m)",        show_x = TRUE)

F3 <- (text_block("Reserve") | text_block("Farms")) /
  (F3A | F3B) /
  text_block("Playback Treatment") +
  plot_layout(heights = c(0.12, 1, 0.12))
F3

###################### Supplemental Material Plots #############################
#Figure S4 (Mpala Experimental Grids)
fill_vals <- c(
  "Wk 1–5\nNorth Grid\nControl"  = "#2C7BB6",
  "Wk 6–10\nNorth Grid\nHuman"   = "#D7191C",
  "Wk 1–5\nSouth Grid\nHuman"    = "#D7191C",
  "Wk 6–10\nSouth Grid\nControl" = "#2C7BB6")

inner_levels <- names(fill_vals)

tree_freq    <- tibble(xlab = factor(inner_levels, inner_levels),
                       mean = c(23.54,  4.75,  6.37, 31.58),
                       lo   = c(13.31,  2.53,  3.38, 18.50),
                       hi   = c(41.64,  8.91, 11.99, 53.88))
weekly_det   <- tibble(xlab = factor(inner_levels, inner_levels),
                       mean = c(19.1,  10.1,  5.85, 15.3),
                       lo   = c(13.4,   6.46,  3.56, 10.6),
                       hi   = c(26.5,  15.0,   8.87, 21.1))
tree_dam     <- tibble(xlab = factor(inner_levels, inner_levels),
                       mean = c(13.59,  4.74, 10.23, 29.33),
                       lo   = c( 9.59,  3.08,  6.98, 20.70),
                       hi   = c(19.26,  7.28, 14.98, 41.56))
nocturnality <- tibble(xlab = factor(inner_levels, inner_levels),
                       mean = 100 - c(72.8, 28.6, 27.3, 71.4),
                       lo   = 100 - c(82.4, 43.9, 42.8, 81.7),
                       hi   = 100 - c(60.5, 17.0, 15.8, 58.3))

S4A <- make_plot(tree_freq,    "Trees damaged\nper transect on grids")
S4B <- make_plot(weekly_det,   "Elephant detections\nper week on grids")
S4C <- make_plot(tree_dam,     "Damage per tree\non grids (percent)")
S4D <- make_plot(nocturnality, "Nighttime detections\non grids (percent)", show_x = TRUE)

S4 <- S4A / S4B / S4C / S4D / text_block("Playback Treatment") +
  plot_layout(heights = c(1, 1, 1, 1, 0.12))
S4

#Figure S5 (Farms)
fill_vals <- c(
  "Wk 1–4\nHuman"     = "#D7191C",
  "Wk 5–8\nControl"   = "#2C7BB6",
  "Wk 9–12\nHuman"    = "#D7191C",
  "Wk 13–16\nControl" = "#2C7BB6")

week_levels <- names(fill_vals)

ele_dam  <- tibble(xlab = factor(week_levels, week_levels),
                   mean = c( 55.5, 405.6,  44.0, 167.9),
                   lo   = c( 22.4, 164.4,  16.8,  62.4),
                   hi   = c(137.0, 1001,  115.0, 452.0))
ele_rate <- tibble(xlab = factor(week_levels, week_levels),
                   mean = c(1.88,  7.72,  4.09, 10.87),
                   lo   = c(1.25,  5.89,  2.94,  8.40),
                   hi   = c(2.82, 10.12,  5.70, 14.07))
ele_conv <- tibble(xlab = factor(week_levels, week_levels),
                   mean = c(29.2, 80.4, 15.2, 47.0),
                   lo   = c(14.8, 70.3,  7.7, 36.2),
                   hi   = c(49.6, 87.7, 27.7, 58.1))
ele_prox <- tibble(xlab = factor(week_levels, week_levels),
                   mean = c( 60.8,   6.38, 100.48,  10.53),
                   lo   = c( 37.09,  4.28,  63.5,    7.19),
                   hi   = c( 99.67,  9.51, 159.0,   15.44))

S5A <- make_plot(ele_dam,  "Crop damage\nper farm (m²)")
S5B <- make_plot(ele_rate, "Elephant detections\nper month on farms")
S5C <- make_plot(ele_conv, "Crop raided if elephant\non farm (% probability)")
S5D <- make_plot(ele_prox, "Distance to speaker\non farm (m)", show_x = TRUE)

S5 <- S5A / S5B / S5C / S5D / text_block("Playback Treatment") +
  plot_layout(heights = c(1, 1, 1, 1, 0.12))
S5
