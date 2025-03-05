
#### Load the appropriate packages

library(brms)
library(rstudioapi)


###Open the data
dw = read.csv(file ="./For github/Data/Data_tracto.csv", header= T)



### Insure that stg/mtg, hemipsheres, sex and ID are factors
dw$stg_mtg = as.factor(as.character(dw$stg_mtg))
dw$hemisphere = as.factor(as.character(dw$hemisphere))
dw$ID = as.factor(as.character(dw$ID))
dw$Sex = as.factor(as.character(dw$Sex))

### create centered mtg_stg and left_right to put in random slopes

dw$hemisphere.c = as.numeric(as.factor(dw$hemisphere)) - mean(as.numeric(as.factor(dw$hemisphere)))

dw$stg_mtg.c = as.numeric(as.factor(dw$stg_mtg)) - mean(as.numeric(as.factor(dw$stg_mtg)))


## scale age

dw$age.z = as.vector(as.numeric(scale(dw$Age)))

#### center sex, more convenient for the plot (Does not impact the results)

dw$sex.centered = as.numeric(as.factor(dw$Sex)) - mean(as.numeric(as.factor(dw$Sex)))

#### Transforming the response

dw$log.w = log(dw$weight)


#### Model the data (log.w = log transformed tract weight) ####
### We use default prior from brms


m.w = brm(log.w ~  hemisphere*stg_mtg*species + sex.centered*species + age.z*species +
            (1+hemisphere.c+stg_mtg.c|ID), data = dw, chains =4, 
          cores = 4, warmup = 2000, 
          iter = 4000)


#Posterior predictive check

pp_check(m.w, ndraws = 100) 

#Check model results
summary(m.w)


#### Prepare conditional effect posterior data for plotting

ce = conditional_effects(m.w, effects = "hemisphere:stg_mtg", conditions = make_conditions(m.w, vars = c("species")))

data.frame(ce["hemisphere:stg_mtg"])


### compile the posterior distribution

post.mw = posterior_samples(m.w)[, c( 'b_Intercept', 'b_hemisphereR',
                                      'b_stg_mtgstg','b_speciesxchimp', 'b_hemisphereR:stg_mtgstg',
                                      'b_hemisphereR:speciesxchimp', 'b_stg_mtgstg:speciesxchimp', 
                                      'b_hemisphereR:stg_mtgstg:speciesxchimp')]

# Calculate the posterior support for each effect

round(sum(post.mw[,"b_stg_mtgstg"] < 0.00)/
        length(post.mw[,"b_stg_mtgstg"]),3)  #100 % of the posterior supports 
# more weight in MTG than STG in HUMANS

round(sum(post.mw[,"b_hemisphereR"] < 0.00)/
        length(post.mw[,"b_hemisphereR"]),3) 

## 99.3 % support for left lateralisation in humans 

round(sum(post.mw[,"b_hemisphereR:stg_mtgstg:speciesxchimp"] > 0.00)/
        length(post.mw[,"b_hemisphereR:stg_mtgstg:speciesxchimp"]),3)

# 56.6 % so no support for the three way interaction (no species difference 
  # in difference of lateralisation between STG and MTG)


round(sum(post.mw[,"b_hemisphereR:stg_mtgstg"] > 0.00)/
        length(post.mw[,"b_hemisphereR:stg_mtgstg"]),3) 

# 93.6% support: Some support for STG to be less lateralised than MTG in humans (and same in chimps).

round(sum(post.mw[,"b_hemisphereR:speciesxchimp"] > 0.00)/
        length(post.mw[,"b_hemisphereR:speciesxchimp"]),3) 

### only 51.2% support for difference in lateralisation between chimps and humans

round(sum(post.mw[,"b_stg_mtgstg:speciesxchimp"] > 0.00)/
        length(post.mw[,"b_stg_mtgstg:speciesxchimp"]),3) 

# 100% posterior support for STG being stronger than MTG in chimps compared to the opposite in humans. 



### effect size of the model

### Effect size

#cond.r2 ### proportion of variance explained by both fixed and random effects
bayes_R2(m.w, re_formula = NULL, summary = T)

# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.8695214 0.02180214 0.8247918 0.9086023


#marginal.r2 ... proportion of variance explained only by the fixed effect

bayes_R2(m.w, re_formula = NA, summary = T)
# 
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.6020323 0.04996803 0.4846003 0.6783597




#### Extract predictive dataframe

pred.data = data.frame(expand.grid(
  stg_mtg = c('mtg','stg'),
  hemisphere = c('L','R'), species = c('human','xchimp'), sex.centered = 0, age.z =0))

pred = predict(m.w, newdata = pred.data, re_formula = NA, summary = F)

 ### fill in pred data in log space with mean and median and 95% CI
# 
for(i in 1:8) (pred.data$mean[i] = mean(pred[,i]))
for(i in 1:8) (pred.data$median[i] = median(pred[,i]))
for(i in 1:8) (pred.data$l95ci[i] = quantile(pred[,i], prob =0.025))
for(i in 1:8) (pred.data$h95ci[i] = quantile(pred[,i], prob =0.975))
#

## plot the data and the model log scale just for Chimpanzees (Figure XX) ####


### define position for STG and MTG in the plot


xx = c(1,2,3,4,5)
jitter.fac = 6
alpha = 0.4
cex = 1.5
ylim = c(-3, 15)
xlim = c(0.5,5.5)


windows(10,7)
### STG left


par (oma = c(2,2,0.5,0.5))

boxplot(dw$log.w[dw$stg_mtg =='stg'&dw$hemisphere =='L'&dw$species =='xchimp'], xlab="", ylab= "",
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[1], width = 0.8, lwd = 2)

par(new = T)

plot( x = jitter(rep(xx[1],length(dw$log.w[dw$stg_mtg =='stg'&dw$hemisphere =='L'&dw$species =='xchimp'])),factor = jitter.fac), 
      y= dw$log.w[dw$stg_mtg =='stg'&dw$hemisphere =='L'&dw$species =='xchimp'], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('blueviolet', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


## STG Right

par (new =T)
boxplot(dw$log.w[dw$stg_mtg =='stg'&dw$hemisphere =='R'&dw$species =='xchimp'], xlab="", 
        ylab= "", col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[2], width = 0.8, lwd = 2)

par(new = T)
plot(x = jitter(rep(xx[2],length(dw$log.w[dw$stg_mtg =='stg'&dw$hemisphere =='R'&dw$species =='xchimp'])),factor = jitter.fac/2), 
     y= dw$log.w[dw$stg_mtg =='stg'&dw$hemisphere =='R'&dw$species =='xchimp'], 
     xlim = xlim, axes = F, pch = 16, col = adjustcolor('blue', alpha = alpha),
     ylab ='', xlab ='', cex = cex, ylim = ylim)


#MTG left

par (new =T)
boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$hemisphere =='L'&dw$species =='xchimp'], xlab="", ylab= "", col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[4], width = 0.8, lwd = 2)

par(new = T)
plot( x = jitter(rep(xx[4],length(dw$log.w[dw$stg_mtg =='mtg'&dw$hemisphere =='L'&dw$species =='xchimp'])),factor = jitter.fac/4), 
      y= dw$log.w[dw$stg_mtg =='mtg'&dw$hemisphere =='L'&dw$species =='xchimp'], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('darkred', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


#MTG right

par (new =T)
boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$hemisphere =='R'&dw$species =='xchimp'], xlab="", ylab= "", col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[5], width = 0.8, lwd = 2)

par(new = T)
plot(x = jitter(rep(xx[5],length(dw$log.w[dw$stg_mtg =='mtg'&dw$hemisphere =='R'&dw$species =='xchimp'])),factor = jitter.fac/5), 
     y= dw$log.w[dw$stg_mtg =='mtg'&dw$hemisphere =='R'&dw$species =='xchimp'], 
     xlim = xlim, axes = F, pch = 16, col = adjustcolor('orange', alpha = alpha),
     ylab ='', xlab ='', cex = cex, ylim = ylim)


### add the model lines and 95% CI

#model line
sqt=0.25


lines(x=c(1-sqt,1+sqt), y= c(pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'L'&pred.data$species =='xchimp'], 
                             pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'L'&pred.data$species =='xchimp']), 
      col="red", lwd=5)

lines(x=c(2-sqt,2+sqt), y= c(pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'R'&pred.data$species =='xchimp'], 
                             pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'R'&pred.data$species =='xchimp']), 
      col="red", lwd=5)


lines(x=c(4-sqt,4+sqt), y= c(pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&pred.data$species =='xchimp'], 
                             pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&pred.data$species =='xchimp']), 
      col="red", lwd=5)

lines(x=c(5-sqt,5+sqt), y= c(pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&pred.data$species =='xchimp'], 
                             pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&pred.data$species =='xchimp']), 
      col="red", lwd=5)


### 95% CI

arrows(x0 = 1, x1 = 1, y0 = pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'L'&
                                               pred.data$species =='xchimp'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'L'&
                               pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 1, x1 = 1, y0 = pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'L'&
                                             pred.data$species =='xchimp'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'L'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')


arrows(x0 = 2, x1 = 2, y0 = pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'R'&
                                               pred.data$species =='xchimp'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'R'&
                               pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 2, x1 = 2, y0 = pred.data$mean[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'R'&
                                               pred.data$species =='xchimp'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='stg'& pred.data$hemisphere == 'R'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')



arrows(x0 = 4, x1 = 4, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                                               pred.data$species =='xchimp'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                               pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 4, x1 = 4, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                                               pred.data$species =='xchimp'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')



arrows(x0 = 5, x1 = 5, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                                               pred.data$species =='xchimp'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                               pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 5, x1 = 5, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                                               pred.data$species =='xchimp'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')


axis(2, labels = c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
     at = log(c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)), las =1)
axis(2, labels = c('',''), at = c(-100, 100))


mtext(c('left', 'right', 'left', 'right'), at =c(1,2,4,5), side =1, cex = 1.2)
mtext(c('STG', 'MTG'), at =c(1.5,4.5), line = 1.5, cex =1.4, side =1)
mtext('Waytotal (log scale)', side =2, line =4, cex =1.3)



####### Plot STG MTG difference in both species (Figure XX)####


#extract datasets from the conditinal effect fonction

ce = conditional_effects(m.w)


## Get the estimates and CI for stg:mtg:species

names(ce)

zz = data.frame(ce["stg_mtg:species"])

names(zz)


### define position for chimps and humans in the plot


xx = c(1,2,3,4,5)
jitter.fac = 6
alpha = 0.4
cex = 1.5
ylim = c(-3, 16)
xlim = c(0.5,5.5)


windows(10,7)

par (oma = c(2,2,0.5,0.5))

### STG chimps

boxplot(dw$log.w[dw$stg_mtg =='stg'&dw$species =='xchimp'], xlab="", ylab= "",
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[1], width = 0.8, lwd = 2)

par(new = T)

plot( x = jitter(rep(xx[1],length(dw$log.w[dw$stg_mtg =='stg'&dw$species =='xchimp'])),factor = jitter.fac), 
      y= dw$log.w[dw$stg_mtg =='stg'&dw$species =='xchimp'], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('blueviolet', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


#MTG chimps

par (new =T)
boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'], xlab="", ylab= "", 
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[2], width = 0.8, lwd = 2)

par(new = T)
plot( x = jitter(rep(xx[2],length(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'])),
                 factor = jitter.fac/4), 
      y= dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('orange', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


### STG humans

par(new = T)
boxplot(dw$log.w[dw$stg_mtg =='stg'&dw$species =='human'], xlab="", ylab= "",
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[4], width = 0.8, lwd = 2)

par(new = T)

plot( x = jitter(rep(xx[4],length(dw$log.w[dw$stg_mtg =='stg'&dw$species =='human'])),
                 factor = jitter.fac/4), 
      y= dw$log.w[dw$stg_mtg =='stg'&dw$species =='human'], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('blueviolet', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


#MTG humans

par (new =T)
boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'], xlab="", ylab= "", 
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[5], width = 0.8, lwd = 2)

par(new = T)
plot( x = jitter(rep(xx[5],length(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'])),
                 factor = jitter.fac/4), 
      y= dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('orange', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


### add the model lines and 95% CI
sqt=0.25



lines(x=c(1-sqt,1+sqt), y= c(zz$stg_mtg.species.estimate__[4], 
                             zz$stg_mtg.species.estimate__[4]), 
      col="red", lwd=5)

lines(x=c(2-sqt,2+sqt), y= c(zz$stg_mtg.species.estimate__[2], 
                             zz$stg_mtg.species.estimate__[2]), 
      col="red", lwd=5)


lines(x=c(4-sqt,4+sqt), y= c(zz$stg_mtg.species.estimate__[3], 
                             zz$stg_mtg.species.estimate__[3]), 
      col="red", lwd=5)

lines(x=c(5-sqt,5+sqt), y= c(zz$stg_mtg.species.estimate__[1], 
                             zz$stg_mtg.species.estimate__[1]), 
      col="red", lwd=5)


### add the CI

arrows(x0 = 1, x1 = 1, y0 = zz$stg_mtg.species.lower__[4], y1 = zz$stg_mtg.species.upper__[4],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 2, x1 = 2, y0 = zz$stg_mtg.species.lower__[2], y1 = zz$stg_mtg.species.upper__[2],
       angle=90, length=0.1, lwd=2, col ='red')


arrows(x0 = 4, x1 = 4, y0 = zz$stg_mtg.species.lower__[3], y1 = zz$stg_mtg.species.upper__[3],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 5, x1 = 5, y0 = zz$stg_mtg.species.lower__[1], y1 = zz$stg_mtg.species.upper__[1],
       angle=90, length=0.1, lwd=2, col ='red')


arrows(x0 = 1, x1 = 1, y1 = zz$stg_mtg.species.lower__[4], y0 = zz$stg_mtg.species.upper__[4],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 2, x1 = 2, y1 = zz$stg_mtg.species.lower__[2], y0 = zz$stg_mtg.species.upper__[2],
       angle=90, length=0.1, lwd=2, col ='red')


arrows(x0 = 4, x1 = 4, y1 = zz$stg_mtg.species.lower__[3], y0 = zz$stg_mtg.species.upper__[3],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 5, x1 = 5, y1 = zz$stg_mtg.species.lower__[1], y0 = zz$stg_mtg.species.upper__[1],
       angle=90, length=0.1, lwd=2, col ='red')



axis(2, labels = c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), 
     at = log(c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)), las =1)
axis(2, labels = c('',''), at = c(-100, 100))


mtext(c('STG', 'MTG', 'STG', 'MTG'), at =c(1,2,4,5), side =1, cex = 1.2)
mtext(c('Chimpanzee', 'Human'), at =c(1.5,4.5), line = 1.5, cex =1.4, side =1)
mtext('Waytotal (log scale)', side =2, line =4, cex =1.3)



##### Figure S8 comparison lateralisation MTG human vs. chimps ####


xx = c(1,2,3,4,5)
jitter.fac = 6
alpha = 0.4
cex = 1.5
ylim = c(-3, 17)
xlim = c(0.5,5.5)


windows(10,7)
par (oma = c(2,2,0.5,0.5))

### MTG chimps left

boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'&dw$hemisphere=="L"], 
        xlab="", ylab= "",
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[1], width = 0.8, lwd = 2)

par(new = T)

plot( x = jitter(rep(xx[1],
                     length(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'&dw$hemisphere=="L"])),factor = jitter.fac), 
      y= dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'&dw$hemisphere=="L"], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('darkred', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


#MTG chimps right

par (new =T)

boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'&dw$hemisphere=="R"], 
        xlab="", ylab= "",
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[2], width = 0.8, lwd = 2)

par(new = T)

plot( x = jitter(rep(xx[2],
                     length(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'&dw$hemisphere=="R"])),
                 factor = jitter.fac/xx[2]), 
      y= dw$log.w[dw$stg_mtg =='mtg'&dw$species =='xchimp'&dw$hemisphere=="R"], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('orange', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


### MTG humans left

par (new =T)

boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'&dw$hemisphere=="L"], 
        xlab="", ylab= "",
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[4], width = 0.8, lwd = 2)

par(new = T)

plot( x = jitter(rep(xx[4],
                     length(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'&dw$hemisphere=="L"])),
                 factor = jitter.fac/xx[4]), 
      y= dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'&dw$hemisphere=="L"], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('darkred', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


#MTG humans

par (new =T)

boxplot(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'&dw$hemisphere=="R"], 
        xlab="", ylab= "",
        col="white", outline=F, whisklty = 0, staplelty = 0, 
        ylim = ylim, xlim = xlim, axes=F, at =xx[5], width = 0.8, lwd = 2)

par(new = T)

plot( x = jitter(rep(xx[5],
                     length(dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'&dw$hemisphere=="R"])),
                 factor = jitter.fac/xx[5]), 
      y= dw$log.w[dw$stg_mtg =='mtg'&dw$species =='human'&dw$hemisphere=="R"], 
      xlim = xlim, axes = F, pch = 16, col = adjustcolor('orange', alpha = alpha),
      ylab ='', xlab ='', cex = cex, ylim = ylim)


### add the model lines and 95% CI
sqt=0.25

lines(x=c(1-sqt,1+sqt), y= rep(pred.data$mean[pred.data$stg_mtg=="mtg"& 
                                                  pred.data$hemisphere=="L"&
                                                  pred.data$species=="xchimp"],2), 
      col="red", lwd=5)

lines(x=c(2-sqt,2+sqt), y= rep(pred.data$mean[pred.data$stg_mtg=="mtg"& 
                                                  pred.data$hemisphere=="R"&
                                                  pred.data$species=="xchimp"],2), 
      col="red", lwd=5)


lines(x=c(4-sqt,4+sqt), y= rep(pred.data$mean[pred.data$stg_mtg=="mtg"& 
                                                  pred.data$hemisphere=="L"&
                                                  pred.data$species=="human"],2), 
      col="red", lwd=5)

lines(x=c(5-sqt,5+sqt), y= rep(pred.data$mean[pred.data$stg_mtg=="mtg"& 
                                                  pred.data$hemisphere=="R"&
                                                  pred.data$species=="human"],2), 
      col="red", lwd=5)


### add the CI

arrows(x0 = 1, x1 = 1, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                                             pred.data$species =='xchimp'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 1, x1 = 1, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                                             pred.data$species =='xchimp'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')



arrows(x0 = 2, x1 = 2, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                                             pred.data$species =='xchimp'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 2, x1 = 2, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                                             pred.data$species =='xchimp'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                              pred.data$species =='xchimp'],
       angle=90, length=0.1, lwd=2, col ='red')


arrows(x0 = 4, x1 = 4, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                                             pred.data$species =='human'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                              pred.data$species =='human'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 4, x1 = 4, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                                             pred.data$species =='human'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'L'&
                              pred.data$species =='human'],
       angle=90, length=0.1, lwd=2, col ='red')



arrows(x0 = 5, x1 = 5, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                                             pred.data$species =='human'],
       y1 = pred.data$l95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                              pred.data$species =='human'],
       angle=90, length=0.1, lwd=2, col ='red')

arrows(x0 = 5, x1 = 5, y0 = pred.data$mean[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                                             pred.data$species =='human'],
       y1 = pred.data$h95ci[pred.data$stg_mtg=='mtg'& pred.data$hemisphere == 'R'&
                              pred.data$species =='human'],
       angle=90, length=0.1, lwd=2, col ='red')


axis(2, labels = c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), 
     at = log(c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)), las =1)
axis(2, labels = c('',''), at = c(-100, 100))

mtext(c('Left', 'Right', 'Left', 'Right'), at =c(1,2,4,5), side =1, cex = 1.2)
mtext(c('Chimpanzee', 'Human'), at =c(1.5,4.5), line = 1.5, cex =1.4, side =1)
mtext('Waytotal (log scale)', side =2, line =4, cex =1.3)









