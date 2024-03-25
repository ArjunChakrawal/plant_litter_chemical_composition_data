

library("readxl")
library(modelsummary)
library(lmerTest)
library(lme4)
library(MuMIn)
library(dplyr)
library(lmtest)
library(nlme)
library(effects)

corrdata = read_excel("../collated data/corrdata.xlsx")
corrdata$Csource= as.factor(corrdata$Csource)
names(corrdata)[names(corrdata) == "Csource"] <- "Litter_category"

corrdata$Decomp = ""
corrdata$Decomp[corrdata$timeDay == 0] ="Fresh"
# corrdata$Decomp[corrdata$timeDay>0 & corrdata$timeDay < 9999] ="D1"
corrdata$Decomp[corrdata$timeDay > 0] ="Decomposed"
corrdata$Decomp=as.factor(corrdata$Decomp)
# corrdata$Decomp <- relevel(corrdata$Decomp, ref = 'Low')

Litter_category  = unique(corrdata$Litter_category)

unique(corrdata$Decomp)

#Model selection for Carbohydrates as a function of AS, AIS and Decomp---- 
Carbohydrates.1 = lm(Carbohydrates ~ (AS+ AIS + Decomp)^2, data = corrdata)
Carbohydrates.2 = lme(Carbohydrates ~ (AS+ AIS + Decomp)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Carbohydrates.3 = lme(Carbohydrates ~ (AS+ AIS + Decomp)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Carbohydrates.4 = lme(Carbohydrates ~ (AS+ AIS + Decomp)^2 ,random =~ AIS|Litter_category,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

# Carbohydrates.5 = lme(Carbohydrates ~ (AS+ AIS + Litter_category)^2 ,random =~ 1|Decomp,data = corrdata, method = "ML",
#                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbohydrates.6 = lme(Carbohydrates ~ (AS+ AIS + Litter_category)^2 ,random =~ AS|Decomp,data = corrdata, method = "ML",
#                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbohydrates.7 = lme(Carbohydrates ~ (AS+ AIS + Litter_category)^2 ,random =~ AIS|Decomp,data = corrdata, method = "ML",
#                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

Carbohydrates.8 = lme(Carbohydrates ~ (AS+ AIS)^2 ,random =~ 1|Litter_category/Decomp,data = corrdata, method = "ML")
Carbohydrates.9 = lme(Carbohydrates ~ (AS+ AIS)^2 ,random =~ AS|Litter_category/Decomp,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbohydrates.10 = lme(Carbohydrates ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category/Decomp,data = corrdata, method = "ML",
#                  control = lmeControl(maxIter = 15000, msMaxIter = 15000,msMaxEval = 15000))

Carbohydrates.11 = lme(Carbohydrates ~ (AS+ AIS)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Carbohydrates.12 = lme(Carbohydrates ~ (AS+ AIS)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Carbohydrates.13 = lme(Carbohydrates ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category,data = corrdata, method = "ML",
                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

AIC(Carbohydrates.1, Carbohydrates.2, Carbohydrates.3, Carbohydrates.4,Carbohydrates.8,Carbohydrates.9,
    Carbohydrates.11,Carbohydrates.12,Carbohydrates.13)
summary(Carbohydrates.2)



Carbohydrates.2A <- update(Carbohydrates.2, .~. -AIS:Decomp)
anova(Carbohydrates.2, Carbohydrates.2A)
summary(Carbohydrates.2A)

Carbohydrates.2B <- update(Carbohydrates.2A, .~. -Decomp)
anova(Carbohydrates.2A, Carbohydrates.2B)
summary(Carbohydrates.2B)

Carbohydrates.2B <- update(Carbohydrates.2A, .~. -AIS)
anova(Carbohydrates.2A, Carbohydrates.2B)
summary(Carbohydrates.2B)

Carbohydrates = lme(Carbohydrates ~ AS + Decomp + AS:AIS + AS:Decomp ,random =~ 1|Litter_category,data = corrdata, method = "REML")
summary(Carbohydrates)
plot(allEffects(Carbohydrates, partial.residuals =F))
r.squaredGLMM(Carbohydrates)

E2 <- resid(Carbohydrates, type = "normalized")
F2 <- fitted(Carbohydrates) 
op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 2)) 
MyYlab <- "Residuals" 
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
plot(E2 ~ AS, data = corrdata, main = "AS", ylab = MyYlab)
plot(E2 ~ AIS, data = corrdata, main = "AS", ylab = MyYlab)
boxplot(E2 ~ Litter_category, data = corrdata, main = "Litter_category", ylab = MyYlab)
boxplot(E2 ~ Decomp, data = corrdata, main = "Decomp", ylab = MyYlab)
par(op)

#Model selection for Proteins as a function of AS, AIS and Decomp---------
Proteins.1 = lm(Proteins ~ (AS+ AIS + Decomp)^2, data = corrdata)
Proteins.2 = lme(Proteins ~ (AS+ AIS + Decomp)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Proteins.3 = lme(Proteins ~ (AS+ AIS + Decomp)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Proteins.4 = lme(Proteins ~ (AS+ AIS + Decomp)^2 ,random =~ AIS|Litter_category,data = corrdata, method = "ML",
#                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

# Proteins.5 = lme(Proteins ~ (AS+ AIS + Litter_category)^2 ,random =~ 1|Decomp,data = corrdata, method = "ML",
#                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Proteins.6 = lme(Proteins ~ (AS+ AIS + Litter_category)^2 ,random =~ AS|Decomp,data = corrdata, method = "ML",
#                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Proteins.7 = lme(Proteins ~ (AS+ AIS + Litter_category)^2 ,random =~ AIS|Decomp,data = corrdata, method = "ML",
#                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

Proteins.8 = lme(Proteins ~ (AS+ AIS)^2 ,random =~ 1|Litter_category/Decomp,data = corrdata, method = "ML")
# Proteins.9 = lme(Proteins ~ (AS+ AIS)^2 ,random =~ AS|Litter_category/Decomp,data = corrdata, method = "ML",
#                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Proteins.10 = lme(Proteins ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category/Decomp,data = corrdata, method = "ML",
#                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

Proteins.11 = lme(Proteins ~ (AS+ AIS)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Proteins.12 = lme(Proteins ~ (AS+ AIS)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Proteins.13 = lme(Proteins ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category,data = corrdata, method = "ML",
#                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

AIC(Proteins.1, Proteins.2, Proteins.3, Proteins.8,Proteins.11,Proteins.12)
summary(Proteins.2)



Proteins.2A <- update(Proteins.2, .~. -Decomp)
anova(Proteins.2, Proteins.2A)
summary(Proteins.2A)

# Proteins.2B <- update(Proteins.2A, .~. -AS:Decomp)
# anova(Proteins.2A, Proteins.2B)
# summary(Proteins.2B)

Proteins = lme(Proteins ~ AS + AIS + AS:AIS + AS:Decomp + Decomp:AIS  ,random =~ 1|Litter_category,data = corrdata, method = "REML")
summary(Proteins)

plot(allEffects(Proteins, partial.residuals =F))


r.squaredGLMM(Proteins)
E2 <- resid(Proteins, type = "normalized")
F2 <- fitted(Proteins) 
op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 2)) 
MyYlab <- "Residuals" 
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
plot(E2 ~ AS, data = corrdata, main = "AS", ylab = MyYlab)
plot(E2 ~ AIS, data = corrdata, main = "AS", ylab = MyYlab)
boxplot(E2 ~ Decomp, data = corrdata, main = "Decomp", ylab = MyYlab)
boxplot(E2 ~ Litter_category, data = corrdata, main = "Litter_category", ylab = MyYlab)
par(op)

#Model selection for Lignins as a function of AS, AIS and Decomp---------

Lignins.1 = lm(Lignins ~ (AS+ AIS + Decomp)^2, data = corrdata)
Lignins.2 = lme(Lignins ~ (AS+ AIS + Decomp)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Lignins.3 = lme(Lignins ~ (AS+ AIS + Decomp)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lignins.4 = lme(Lignins ~ (AS+ AIS + Decomp)^2 ,random =~AIS|Litter_category,data = corrdata, method = "ML",
#                  control = lmeControl(maxIter = 10000, msMaxIter = 10000,msMaxEval = 10000))

# Lignins.5 = lme(Lignins ~ (AS+ AIS + Litter_category)^2 ,random =~ 1|Decomp,data = corrdata, method = "ML",
#                  control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lignins.6 = lme(Lignins ~ (AS+ AIS + Litter_category)^2 ,random =~ AS|Decomp,data = corrdata, method = "ML",
#                  control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lignins.7 = lme(Lignins ~ (AS+ AIS + Litter_category)^2 ,random =~ AIS|Decomp,data = corrdata, method = "ML",
#                  control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

Lignins.8 = lme(Lignins ~ (AS+ AIS)^2 ,random =~ 1|Litter_category/Decomp,data = corrdata, method = "ML")
Lignins.9 = lme(Lignins ~ (AS+ AIS)^2 ,random =~ AS|Litter_category/Decomp,data = corrdata, method = "ML",
                 control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lignins.10 = lme(Lignins ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category/Decomp,data = corrdata, method = "ML",
#                   control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

Lignins.11 = lme(Lignins ~ (AS+ AIS)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Lignins.12 = lme(Lignins ~ (AS+ AIS)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Lignins.13 = lme(Lignins ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category,data = corrdata, method = "ML",
                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

AIC(Lignins.1, Lignins.2,Lignins.3,Lignins.4, Lignins.8,Lignins.9,Lignins.11,Lignins.12,Lignins.13)
summary(Lignins.8)


Lignins.8A <- update(Lignins.8, .~. -AS:AIS)
anova(Lignins.8, Lignins.8A)
summary(Lignins.8A)

Lignins.8B <- update(Lignins.8A, .~. -AS)
anova(Lignins.8A, Lignins.8B)
summary(Lignins.8B)


Lignins = lme(Lignins ~ AIS, random =~ 1|Litter_category/Decomp,data = corrdata, method = "REML",
             control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

summary(Lignins)
plot(allEffects(Lignins, partial.residuals =T))
plot(Lignins)
qqnorm(residuals(Lignins))
qqline(residuals(Lignins))
r.squaredGLMM(Lignins)

E2 <- resid(Lignins, type = "normalized")
F2 <- fitted(Lignins) 
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals" 
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
plot(E2 ~ AS, data = corrdata, main = "AS", ylab = MyYlab)
boxplot(E2 ~ Litter_category, data = corrdata, main = "Litter_category", ylab = MyYlab)
boxplot(E2 ~ Decomp, data = corrdata, main = "Decomp", ylab = MyYlab)
par(op)


# fixed_coef <- fixef(Lignins)
# rand_coef <- ranef(Lignins)$Litter_category[, 1]
# 
# op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
# # Plot 1: AIS vs. Carbohydrates
# plot(corrdata$AIS, corrdata$Lignins, main = "AIS vs. Carbohydrates", xlab = "AIS", ylab = "Carbohydrates")
# for (i in 1:length(Litter_category)) {
#   abline(a = fixed_coef[1] + ranef(Lignins)$Litter_category[i, 1], b = fixed_coef[2]+ ranef(Lignins)$Litter_category[i, 1], col = i)
# }
# # Plot 2: AS vs. Carbohydrates
# plot(corrdata$AS, corrdata$Lignins, main = "AS vs. Carbohydrates", xlab = "AS", ylab = "Carbohydrates")
# for (i in 1:length(Litter_category)) {
#   abline(a = fixed_coef[1] + ranef(Lignins)$Litter_category[i, 1], b = fixed_coef[3]+ ranef(Lignins)$Litter_category[i, 1], col = i)
# }
# par(op)

#Model selection for Lipids as a function of AS, AIS and Decomp---------

Lipids.1 = lm(Lipids ~ (AS+ AIS + Decomp)^2, data = corrdata)
Lipids.2 = lme(Lipids ~ (AS+ AIS + Decomp)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lipids.3 = lme(Lipids ~ (AS+ AIS + Decomp)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lipids.4 = lme(Lipids ~ (AS+ AIS + Decomp)^2 ,random =~ AIS|Litter_category,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lipids.5 = lme(Lipids ~ (AS+ AIS + Litter_category)^2 ,random =~ 1|Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lipids.6 = lme(Lipids ~ (AS+ AIS + Litter_category)^2 ,random =~ AS|Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lipids.7 = lme(Lipids ~ (AS+ AIS + Litter_category)^2 ,random =~ AIS|Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Lipids.8 = lme(Lipids ~ (AS+ AIS)^2 ,random =~ 1|Litter_category/Decomp,data = corrdata, method = "ML")
# Lipids.9 = lme(Lipids ~ (AS+ AIS)^2 ,random =~ AS|Litter_category/Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Lipids.10 = lme(Lipids ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category/Decomp,data = corrdata, method = "ML",
#                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

AIC(Lipids.1, Lipids.2,Lipids.8)
summary(Lipids.2)


Lipids.2A <- update(Lipids.2, .~. -AS:Decomp)
anova(Lipids.2, Lipids.2A)
summary(Lipids.2A)

Lipids.2B <- update(Lipids.2A, .~. -AIS)
anova(Lipids.2A, Lipids.2B)
summary(Lipids.2B)

# Lipids.2C <- update(Lipids.2B, .~. -AS:AIS)
# anova(Lipids.2B, Lipids.2C)
# summary(Lipids.2C)

Lipids = lme(Lipids ~ AS + Decomp + AS:AIS + AIS:Decomp , random =~ 1|Litter_category,data = corrdata, method = "REML",
             control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
summary(Lipids)
r.squaredGLMM(Lipids)

E2 <- resid(Lipids, type = "normalized")
F2 <- fitted(Lipids) 
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals" 
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
plot(E2 ~ AS, data = corrdata, main = "AS", ylab = MyYlab)
plot(E2 ~ AIS, data = corrdata, main = "AS", ylab = MyYlab)
boxplot(E2 ~ Decomp, data = corrdata, main = "Decomp", ylab = MyYlab)
par(op)


#Model selection for Carbonyls as a function of AS, AIS and Decomp---------

Carbonyls.1 = lm(Carbonyls ~ (AS+ AIS + Decomp)^2, data = corrdata)
Carbonyls.2 = lme(Carbonyls ~ (AS+ AIS + Decomp)^2 ,random =~ 1|Litter_category,data = corrdata, method = "ML",
              control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbonyls.3 = lme(Carbonyls ~ (AS+ AIS + Decomp)^2 ,random =~ AS|Litter_category,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbonyls.4 = lme(Carbonyls ~ (AS+ AIS + Decomp)^2 ,random =~ AIS|Litter_category,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbonyls.5 = lme(Carbonyls ~ (AS+ AIS + Litter_category)^2 ,random =~ 1|Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbonyls.6 = lme(Carbonyls ~ (AS+ AIS + Litter_category)^2 ,random =~ AS|Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbonyls.7 = lme(Carbonyls ~ (AS+ AIS + Litter_category)^2 ,random =~ AIS|Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
Carbonyls.8 = lme(Carbonyls ~ (AS+ AIS)^2 ,random =~ 1|Litter_category/Decomp,data = corrdata, method = "ML")
# Carbonyls.9 = lme(Carbonyls ~ (AS+ AIS)^2 ,random =~ AS|Litter_category/Decomp,data = corrdata, method = "ML",
#               control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
# Carbonyls.10 = lme(Carbonyls ~ (AS+ AIS)^2 ,random =~ AIS|Litter_category/Decomp,data = corrdata, method = "ML",
#                control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))

AIC(Carbonyls.1, Carbonyls.2,Carbonyls.8)
summary(Carbonyls.2)


Carbonyls.2A <- update(Carbonyls.2, .~. -AS)
anova(Carbonyls.2, Carbonyls.2A)
summary(Carbonyls.2A)

Carbonyls.2B <- update(Carbonyls.2A, .~. -AS:Decomp)
anova(Carbonyls.2A, Carbonyls.2B)
summary(Carbonyls.2B)

Carbonyls = lme(Carbonyls ~ AS + Decomp + AS:AIS + AIS:Decomp , random =~ 1|Litter_category,data = corrdata, method = "REML",
            control = lmeControl(maxIter = 5000, msMaxIter = 5000,msMaxEval = 5000))
summary(Carbonyls)
r.squaredGLMM(Carbonyls)

E2 <- resid(Carbonyls, type = "normalized")
F2 <- fitted(Carbonyls) 
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals" 
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
plot(E2 ~ AS, data = corrdata, main = "AS", ylab = MyYlab)
plot(E2 ~ AIS, data = corrdata, main = "AS", ylab = MyYlab)
boxplot(E2 ~ Litter_category, data = corrdata, main = "Decomp", ylab = MyYlab)
par(op)



# -----
models <- list("Carbohydrates" = Carbohydrates,"Proteins"=Proteins,"Lignins"=Lignins, "Lipids"=Lipids, "Carbonyls"=Carbonyls)
bbmle::AICtab(models,base=TRUE,logLik=TRUE)
custom_notes <- list('Significance levels: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1, + p > 0.1')
model_formulas <- lapply(models, function(model) {
  formula(model)
})
model_formulas_char <- sapply(model_formulas, as.character)

modelsummary(models,fmt = 2,
             estimate  = "{estimate} ({std.error}){stars}",
             statistic = NULL,
             notes = custom_notes,
             gof_omit = 'ICC|RMSE|BIC',
             addrows = list("Model Formula" = model_formulas_char),
             output = "tableNMR.docx"
)

# Model selection for AS as a function of Carbohydrates + Proteins + Lignins + Lipids + Decomp---------


AS.1 = lme(AS ~ (Carbohydrates + Proteins + Lignins + Lipids + Decomp)^2,
           random =~ 1|Litter_category, data = corrdata,method = "ML")
summary(AS.1)

AS.2 = lme(AS ~ (Carbohydrates + Proteins + Lignins + Lipids )^2,
         random =~ 1|Litter_category, data = corrdata,method = "ML")
anova(AS.1, AS.2)
summary(AS.2)

AS.2A <- update(AS.2, .~. -Carbohydrates:Lignins - Carbohydrates:Lipids -
                  Proteins:Lignins -Proteins:Lipids)

anova(AS.2, AS.2A)
summary(AS.2A)

AS.2B <- update(AS.2A, .~. -Lignins-Lipids)
anova(AS.2A, AS.2B)
summary(AS.2B)

AS.2C <- update(AS.2B, .~. -Proteins)
anova(AS.2B, AS.2C)
summary(AS.2C)

AS.2D <- update(AS.2C, .~. -Lignins:Lipids)
anova(AS.2C, AS.2D)
summary(AS.2D)

r.squaredGLMM(AS.2D)
randcoefI = ranef(AS.2D)[,1]
fixedcoef = fixef(AS.2D)

plot(corrdata$Carbohydrates, corrdata$AS, main = "Carbohydrates vs. AS", 
     ylab = "AS", xlab = "Carbohydrates")
for(i in 1:length(randcoefI)){
  abline(a = fixedcoef[1] + randcoefI[i] , b = fixedcoef[2], col = i)
}

plot(AS.2D)
qqnorm(residuals(AS.2D))
qqline(residuals(AS.2D))
x = ranef(AS.2D)
qqnorm(x$`(Intercept)`)
qqline(x$`(Intercept)`)

E2 <- resid(AS.2D, type = "normalized")
F2 <- fitted(AS.2D) 
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals" 
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
plot(E2 ~ Carbohydrates, data = corrdata, main = "AS", ylab = MyYlab)
plot(E2 ~ Proteins, data = corrdata, main = "AS", ylab = MyYlab)
boxplot(E2 ~ Litter_category, data = corrdata, main = "Decomp", ylab = MyYlab)
par(op)


# Model selection for AIS as a function of Carbohydrates + Proteins + Lignins + Lipids + Decomp---------
AIS.1 = lme(AIS ~ (Carbohydrates + Proteins + Lignins + Lipids + Decomp)^2,
           random =~ 1|Litter_category, data = corrdata,method = "ML")
summary(AIS.1)

AIS.2 = lme(AIS ~ (Carbohydrates + Proteins + Lignins + Lipids )^2,
           random =~ 1|Litter_category, data = corrdata,method = "ML")
anova(AIS.1, AIS.2)
summary(AIS.2)

AIS.3 = lme(AIS ~ (Carbohydrates +Proteins+ Lignins  +Decomp)^2,
            random =~ 1|Litter_category, data = corrdata,method = "ML")
anova(AIS.1, AIS.3)
summary(AIS.3)

AIS.1A <- update(AIS.1, .~. -Lipids:Decomp-Lignins:Decomp -Carbohydrates:Decomp- 
                   Lignins:Lipids -Proteins:Lignins -Carbohydrates:Lipids -Carbohydrates:Lignins)
anova(AIS.1, AIS.1A)
summary(AIS.1A)

AIS.1B <- update(AIS.1A, .~. -Proteins)
anova(AIS.1A, AIS.1B)
summary(AIS.1B)

r.squaredGLMM(AIS.1B)
randcoefI = ranef(AIS.1B)[,1]
fixedcoef = fixef(AIS.1B)

plot(corrdata$Carbohydrates, corrdata$AIS, main = "Carbohydrates vs. AIS", 
     ylab = "AIS", xlab = "Carbohydrates")
for(i in 1:length(randcoefI)){
  abline(a = fixedcoef[1] + randcoefI[i] , b = fixedcoef[2], col = i)
}

plot(corrdata$Lignins, corrdata$AIS, main = "Lignins vs. AIS", 
     ylab = "AIS", xlab = "Lignins")
for(i in 1:length(randcoefI)){
  abline(a = fixedcoef[1] + randcoefI[i] , b = fixedcoef[3], col = i)
}

plot(corrdata$Lipids, corrdata$AIS, main = "Lipids vs. AIS", 
     ylab = "AIS", xlab = "Lipids")
for(i in 1:length(randcoefI)){
  abline(a = fixedcoef[1] + randcoefI[i] , b = fixedcoef[4], col = i)
}

plot(AIS.1B)
qqnorm(residuals(AIS.1B))
qqline(residuals(AIS.1B))
x = ranef(AIS.1B)
qqnorm(x$`(Intercept)`)
qqline(x$`(Intercept)`)

E2 <- resid(AIS.1B, type = "normalized")
F2 <- fitted(AIS.1B) 
op <- par(mfrow = c(4, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals" 
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab,panel.first = grid())
plot(E2 ~ Carbohydrates, data = corrdata, main = "AIS", ylab = MyYlab,panel.first = grid())
plot(E2 ~ Proteins, data = corrdata, main = "AIS", ylab = MyYlab,panel.first = grid())
plot(E2 ~ Lignins, data = corrdata, main = "AIS", ylab = MyYlab,panel.first = grid())
plot(E2 ~ Lipids, data = corrdata, main = "AIS", ylab = MyYlab,panel.first = grid())
boxplot(E2 ~ Decomp, data = corrdata, main = "Decomp", ylab = MyYlab, panel.first = grid())
boxplot(E2 ~ Litter_category, data = corrdata, main = "Litter_category", ylab = MyYlab,panel.first = grid())
par(op)



models <- list("AS" = AS.2D,"AIS"=AIS.1B)
for (model in models) {
  print(r.squaredGLMM(model))
  print(AIC(model))
}
r.squaredGLMM(AS.2D)
r.squaredGLMM(AIS.1B)

custom_notes <- list('Significance levels: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1, + p > 0.1')

modelsummary(models,fmt = 2,
             estimate  = "{estimate} ({std.error}){stars}",
             statistic = NULL,
             notes = custom_notes,
             gof_omit = 'ICC|RMSE|BIC',
             output = "tablePA.docx"
)


