#' Load packages
library(glmmTMB)
library(parallel)
library(ggplot2)

#' Clear objects from environment
rm(list = ls())

#' Functions to report messages (e.g. progress) from inside mclapply
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}


#' Define the design  

#' n.family breeding pairs of fish produce n.family families of eggs  
n.family <- 10

#' each family is split evenly among n.tank tanks
n.tank <- 4

#' n.fish.per.family fish per family are sampled
n.fish.per.family <- 2

#' each fish is measured n.measure times
n.measure <- 3

#' Put all of these design choices into a template data set, ready for recording the response
dat <- expand.grid(measure = 1:n.measure, fish = 1:n.fish.per.family, family = 1:n.family, tank = 1:n.tank)
n <- nrow(dat)
n

#' What about treatment? There are two temperature treatments, let's assume the odd-numbered tanks
#' get the low temperature and the even numbered tanks get the high temperature
dat$hi.temp <- as.integer(dat$tank > (n.tank + 1)/2)
table(dat$tank, dat$hi.temp)

#' Turn each ID variable into a factor
dat$tank <- factor(paste0(c("LOW-", "HIGH-")[dat$hi.temp + 1], "tank", dat$tank))
dat$family <- factor(paste0("fam", dat$family))
dat$fish <- factor(paste0(dat$tank, "-", dat$family, "-fish", dat$fish))
dat$measure <- factor(paste0(dat$fish, "-measure", dat$measure))
table(dat$fish, dat$tank)
table(dat$fish)
table(dat$measure)

#' Check that there are the right number of unique levels for each design variable
cbind(sapply(dat, function(x) length(unique(x))))

#' Make assumptions about means, effects, SDs
#' in units of mg of oxygen per hour
lo.temp.mean <- 0.8 
hi.temp.effect <- 0.2
family.sd <- 0.05
tank.sd <- 0.2
fish.sd <- 0.05
resid.sd <- 0.15

head(dat)

#' I could simulate this from a function (e.g. sim.glmm from my GLMMmisc package on github)
#' but to illustrate I'll build it up from scratch.
#' Create all the terms of the LMM as columns in the data frame:
dat$fe.intercept <- lo.temp.mean
dat$fe.effect <- hi.temp.effect * dat$hi.temp
dat$re.family <- rnorm(nlevels(dat$family), sd = family.sd)[dat$family]
dat$re.tank <- rnorm(nlevels(dat$tank), sd = tank.sd)[dat$tank]
dat$re.fish <- rnorm(nlevels(dat$fish), sd = fish.sd)[dat$fish]
dat$re.resid <- rnorm(n, sd = resid.sd)

#' Add them together to give the simulated response
dat$response <- 
  dat$fe.intercept + dat$fe.effect + dat$re.family + dat$re.tank + dat$re.fish + dat$re.resid 

hist(dat$response)


ggplot(data = dat, mapping = aes(x = tank, y = response, colour = family, shape = factor(hi.temp))) +
  geom_jitter(width = 0.2) 
diff(tapply(dat$response, dat$hi.temp, mean))

#' Fit a LMM and test for the hi temp effect, deliberately excluding important random terms.  
#' 
#' assuming 
#' (1) replication at the level of the individual alone, 
#' (2) then with +family, then 
#' (3) with +tank effects

control <- lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))

fit.fishonly <- lmer(response ~ hi.temp + (1 | fish), data = dat, control = control)
coef(summary(fit.fishonly))

fit.fishfamily <- lmer(response ~ hi.temp + (1 | fish) + (1 | family), data = dat, control = control)
coef(summary(fit.fishfamily))

fit.right <- lmer(response ~ hi.temp + (1 | family) + (1 | tank) + (1 | fish), data = dat, control = control)
coef(summary(fit.right))


n.sim <- 100
sim.results <-
  replicate(n.sim, {
    
    #' Create all the terms of the LMM as columns in the data frame:
    dat$fe.intercept <- lo.temp.mean
    dat$fe.effect <- hi.temp.effect * dat$hi.temp
    dat$re.family <- rnorm(nlevels(dat$family), sd = family.sd)[dat$family]
    dat$re.tank <- rnorm(nlevels(dat$tank), sd = tank.sd)[dat$tank]
    dat$re.fish <- rnorm(nlevels(dat$fish), sd = fish.sd)[dat$fish]
    dat$re.resid <- rnorm(n, sd = resid.sd)
    
    #' Add them together to give the simulated response
    dat$response <- 
      dat$fe.intercept + dat$fe.effect + dat$re.family + dat$re.tank + dat$re.fish + dat$re.resid 
    
    #' Fit LMMs
    #' 
    #' assuming 
    #' (1) replication at the level of the individual alone, 
    #' (2) then with +family, then 
    #' (3) with +tank effects
    
    
    fit.fishonly <- lmer(response ~ hi.temp + (1 | fish), data = dat, control = control)
    coef.tab.fishonly <- coef(summary(fit.fishonly))
    p.value.fishonly <- coef.tab.fishonly["hi.temp", "Pr(>|t|)"]
    
    fit.fishfamily <- lmer(response ~ hi.temp + (1 | fish) + (1 | family), data = dat, control = control)
    coef.tab.fishfamily <- coef(summary(fit.fishfamily))
    p.value.fishfamily <- coef.tab.fishfamily["hi.temp", "Pr(>|t|)"]
    
    fit.right <- lmer(response ~ hi.temp + (1 | family) + (1 | tank) + (1 | fish), data = dat, control = control)
    coef.tab.right <- coef(summary(fit.right))
    p.value.right <- coef.tab.right["hi.temp", "Pr(>|t|)"]
    
    c(p.value.fishonly = p.value.fishonly, 
      p.value.fishfamily = p.value.fishfamily, 
      p.value.right = p.value.right)
    
  })


apply(sim.results < 0.05, 1, mean)



rm(dat)
#' Investigate the effect of differing balance in within vs between variation
#' on the numbers of fish required, 
par.tab <-
  expand.grid(
    n.sim = 500,
    n.family = 10,
    n.tank = c(4, 10, 20),
    n.fish.per.family = c(2, 5, 10),
    n.measure = 3,
    lo.temp.mean = 0.8,
    hi.temp.effect = c(0.15, 0.5),
    family.sd = 0.05,
    tank.sd = c(0.05, 0.2),
    resid.sd = 0.1)

par.tab$total.n.fish <- par.tab$n.fish.per.family * par.tab$n.family *  par.tab$n.tank

par.tab$fish.sd <- round(sqrt(sum(unique(par.tab$tank.sd)^2) - (par.tab$tank.sd)^2), 5)

table(fish.sd = par.tab$fish.sd, tank.sd = par.tab$tank.sd)

start.time <- Sys.time()
sim.results.list <-
  mclapply(1:nrow(par.tab), function(j) {
    
    # j=36
    
    if(j/nrow(par.tab) == round(j/nrow(par.tab), 2)) 
      message_parallel(paste0(100 * (j/nrow(par.tab)), "%"))

    dat <- 
      expand.grid(measure = 1:par.tab$n.measure[j], 
                  fish = 1:par.tab$n.fish.per.family[j], 
                  family = 1:par.tab$n.family[j], 
                  tank = 1:par.tab$n.tank[j])
    n <- nrow(dat)
    dat$hi.temp <- as.integer(dat$tank > (n.tank + 1)/2)
    dat$tank <- factor(paste0(c("LOW-", "HIGH-")[dat$hi.temp + 1], "tank", dat$tank))
    dat$family <- factor(paste0("fam", dat$family))
    dat$fish <- factor(paste0(dat$tank, "-", dat$family, "-fish", dat$fish))
    dat$measure <- factor(paste0(dat$fish, "-measure", dat$measure))
    
    
    sim.p.values <-
      replicate(par.tab$n.sim[j], {
        dat$fe.intercept <- par.tab$lo.temp.mean[j]
        dat$fe.effect <- par.tab$hi.temp.effect[j] * dat$hi.temp
        dat$re.family <- rnorm(nlevels(dat$family), sd = par.tab$family.sd[j])[dat$family]
        dat$re.tank <- rnorm(nlevels(dat$tank), sd = par.tab$tank.sd[j])[dat$tank]
        dat$re.fish <- rnorm(nlevels(dat$fish), sd = par.tab$fish.sd[j])[dat$fish]
        dat$re.resid <- rnorm(n, sd = par.tab$resid.sd[j])
        
        #' Add them together to give the simulated response
        dat$response <- 
          dat$fe.intercept + dat$fe.effect + dat$re.family + dat$re.tank + dat$re.fish + dat$re.resid 
        fit <- lmer(response ~ hi.temp + (1 | family) + (1 | tank) + (1 | fish), data = dat)
        coef.tab <- coef(summary(fit))
        p.value <- coef.tab["hi.temp", "Pr(>|t|)"]
        
        c(p.value = p.value)
      })
    sim.p.values
  }, mc.cores = detectCores())
finish.time <- Sys.time()
finish.time - start.time

par.tab$power <- sapply(sim.results.list, function(x) mean(x < 0.05))

hist(par.tab$power)

apply(par.tab, 2, unique)

par.tab$n.fish.per.family <- factor(par.tab$n.fish.per.family)
par.tab$tank.sd <- factor(par.tab$tank.sd)
par.tab$fish.sd <- factor(par.tab$fish.sd)
ggplot(data = droplevels(par.tab[par.tab$hi.temp.effect != 0, ]), 
       mapping = aes(x = total.n.fish, y = power, colour = n.fish.per.family, shape = n.fish.per.family)) +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Paired") +
  facet_wrap(~ tank.sd + fish.sd + hi.temp.effect, labeller = "label_both")
  
