library(lme4)
library(ggplot2)
library(Rmisc)
library(lmerTest)
library(leaps)
library(ggthemes)
library(car)

##FUNCTIONS + other info ####

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

emplog = function(p, totalN) {
  y = p*totalN
  return(log( (y+.5)/ (totalN-y + .5) ))
}
emplogweight = function(p, totalN) {
  y = p*totalN
  return(1 / ( 1 / (y+.5) + 1 / (totalN-y + .5) ))
}

meanVerbOn = 12
meanNounOn = 30.10815



### ANALYSIS ####

## Import the data [change this to where you saved the GAMMData.csv files]

setwd("~/Downloads/Raw Data")

d<- read.csv('GAMMData.csv') #or however you prefer to import your data

View(d) #let's walk through the headers together. What do they mean, what has already 
#been done to the data?





#Let's look at the data together:

Overall <- summarySE(d, measurevar = c("targetHit"), groupvars = c("bin", "Condition"))
Overall <- Overall[which(Overall$Condition != "-1"),]
Overall$time <- Overall$bin*50

ggplot(Overall, aes(time, targetHit, color = Condition)) + geom_point(size = 2) + 
  geom_ribbon(alpha = 0.5, aes(fill = Condition, ymin = targetHit -se, ymax = targetHit + se), colour = NA) + 
  geom_vline(xintercept = meanVerbOn*50) + geom_vline(xintercept = meanNounOn*50) + 
  theme_classic(base_family = "serif", base_size = 14) +
  labs(x = "Time (ms)", y = "Proportion of Target Fixations")


#Let's zoom in on the critical window:

CriticalWindow <- summarySE(d, measurevar = c("targetHit"), groupvars = c("bin2","Condition"))
CriticalWindow <- CriticalWindow[which(CriticalWindow$Condition != "-1"),]

ggplot(CriticalWindow, aes(bin2, targetHit, color = Condition)) + geom_point(size = 2) + 
  geom_ribbon(alpha = 0.5, aes(fill = Condition, ymin = targetHit -se, ymax = targetHit + se), colour = NA) + 
  theme_classic(base_family = "serif", base_size = 14) +
  labs(x = "# bins", y = "Proportion of Target Fixations")




#Using lme4

Model <- summarySE(d, measurevar = c("targetHit"), groupvars = 
                     c("bin2","Condition","SubjectNumber","CurrentSentence"))

Model <- Model[which(Model$Condition != "-1"),]

Model$Condition <- as.factor(Model$Condition)
Model$SubjectNumber <- as.factor(Model$SubjectNumber)
Model$CurrentSentence <- as.factor(Model$CurrentSentence)

contrasts(Model$Condition) <- contr.sum(2)

model1 <- lmer(emplog(targetHit, 3) ~ bin2*Condition + (1|SubjectNumber) + 
                 (1|CurrentSentence), data = Model, weight = emplogweight(targetHit, 3))

summary(model1)

ggplot(CriticalWindow,aes(bin2, targetHit, color = Condition)) + geom_point(size = 2) + 
  geom_smooth(method= 'lm') + theme_classic(base_family = "serif", base_size = 14) +
  labs(x = "# bins", y = "Proportion of Target Fixations")






# Let's GAMM-ify!
library(mgcv)
library(itsadug)




model1 <- lmer(emplog(targetHit, 3) ~ bin2*Condition + (1|SubjectNumber) + 
                 (1|CurrentSentence), data = Model, weight = emplogweight(targetHit, 3))






model2 <- bam(emplog(targetHit,3)~ Condition + s(bin2, by = Condition) + 
                s(bin2, SubjectNumber, bs = "fs", m = 1) + 
                s(bin2, CurrentSentence, bs = "fs", m = 1), 
              data = Model, weights = emplogweight(targetHit, 3), rho = -0.03781602)







#Let's run some models

model3 <- bam(emplog(targetHit,3)~ Condition + s(bin2, by = Condition) + 
                s(bin2, SubjectNumber, bs = "fs", m = 1) + 
                s(bin2, CurrentSentence, bs = "fs", m = 1), 
              data = Model, weights = emplogweight(targetHit, 3))

#For the rho value, we run the model first, then apply the function, add the value, and run it again
start_value_rho(model3)

summary(model3) #let's look at the output







#But wait, how do I know if they are different from each other?

#Method #1:

#Make the factor ordered:

Model$OFCondition <-as.ordered(Model$Condition)

contrasts(Model$OFCondition) <- contr.sum(2)

model4 <- bam(emplog(targetHit,3)~ OFCondition + s(bin2) + s(bin2, by = OFCondition) + 
                s(bin2, SubjectNumber, bs = "fs", m = 1) + 
                s(bin2, CurrentSentence, bs = "fs", m = 1), 
              data = Model, weights = emplogweight(targetHit, 3), rho = -0.08901097)


#Method #2:

plot_diff(model3,view = 'bin2', comp=list(Condition = c('res','nonres')))



#Model Diagnostics

#Whatever you use for glm you can use here too.

gam.check(model4)
?choose.k

#Model comparisons

#Whatever you use for glm you can use here too.

compareML
anova(model3,model4)

















