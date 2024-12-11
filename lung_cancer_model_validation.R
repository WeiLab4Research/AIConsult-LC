library(survival)
load('lung_cancer_model_validation.RData')

# The algorithms of the ReducedHUNT model, the JPHC model, and the Shanghai-LCM
# are not included because the authors did not give permission.

# risk calculating fuction for each model------------

#' Model: Bach
#' Reference：Bach PB, Kattan MW, Thornquist MD, et al. Variations in lung cancer risk among smokers. J Natl Cancer Inst. 2003;95(6):470-478. doi:10.1093/jnci/95.6.470
#' Published years: 2003
#' Country: US
#' Targeted population: Ever smokers
#' Prediction time range: 1-10 years
risk.bach <- function(data, nyears) {
  bach.formula <- function(outcome) {
    # age = continuous
    # female = binary indicator
    # cpd = cigarettes per day
    # smkyears = years smoked
    # qtyears = years quit
    # asb = asbestos exposure binary indicator
    
    # Formula from appendix of Bach's paper minus coefficients
    formula.string <- "female+
      asb+
      age+
      I(I((age-53.459001)^3)*I(age>53))+
      I(I((age-61.954825)^3)*I(age>61))+
      I(I((age-70.910335)^3)*I(age>70))+
      qtyears+
      I((qtyears)^3)+
      I(I((qtyears-0.50513347)^3)*I(qtyears>0))+
      I(I((qtyears-12.295688)^3)*I(qtyears>12))+
      smkyears+
      I(I((smkyears-27.6577)^3)*I(smkyears>27))+
      I(I((smkyears-40)^3)*I(smkyears>40))+
      I(I((smkyears-50.910335)^3)*I(smkyears>50))+
      cpd+
      I(I((cpd-15)^3)*I(cpd>15))+
      I(I((cpd-20.185718)^3)*I(cpd>20))+
      I(I((cpd-40)^3)*I(cpd>40))"
    
    formula(paste(outcome, "~", formula.string, collapse = ""))
  }
  # Coefficients from diagnosis of lung cancer appendix
  bach.coef <- c(
    intercept = -9.7960571,
    female = -0.05827261,
    asb = 0.2153936,
    age = 0.070322812,
    age2 = -0.00009382122,
    age3 = 0.00018282661,
    age4 = -0.000089005389,
    qtyears = -0.085684793,
    qtyears2 = 0.0065499693,
    qtyears3 = -0.0068305845,
    qtyears4 = 0.00028061519,
    smkyears = 0.11425297,
    smkyears1 = -0.000080091477,
    smkyears2 = 0.00017069483,
    smkyears3 = -0.000090603358,
    cpd = 0.060818386,
    cpd1 = -0.00014652216,
    cpd2 = 0.00018486938,
    cpd3 = -0.000038347226
  )
  
  # Coefficients from probability of death from something other than lung cancer
  bach.coef.mort <- c(
    Intercept = -7.2036219, female = -0.49042298, asb = 0.06084611,
    age = 0.099168033, age1 = 6.2174577e-06, age2 = -1.2115774e-05, age3 = 5.8983164e-06,
    qtyears = -0.023358962, qtyears2 = 0.0019208669, qtyears3 = -0.0020031611, qtyears4 = 8.2294194e-05,
    smkyears = 0.020041889, smkyears1 = 6.5443781e-06, smkyears2 = -1.3947696e-05, smkyears3 = 7.4033175e-06,
    cpd = 0.015490665, cpd1 = -1.737645e-05, cpd2 = 2.1924149e-05, cpd3 = -4.5476985e-06
  )
  
  vars <- c(
    "female",
    "age",
    "cpd",
    "smkyears",
    "qtyears",
    "asb"
  )
  
  var.not.in <- names(data)[!(names(data) %in% vars) & sapply(data, function(x) sum(is.na(x))) == 0][1] # SELECT ARBITRARY VAR FOR OUTCOME
  
  P.lung <- matrix(0, nrow(data), ncol = 10)
  P.death <- matrix(0, nrow(data), ncol = 10)
  
  for (i in 1:max(nyears)) {
    # Calculate RRs and probabilities
    rr.lung <- model.matrix(bach.formula(var.not.in), data) %*% bach.coef
    rr.death <- model.matrix(bach.formula(var.not.in), data) %*% bach.coef.mort
    P.lung[complete.cases(subset(data, select = vars)), i] <- .99629^exp(rr.lung)
    P.death[complete.cases(subset(data, select = vars)), i] <- 0.9917663^exp(rr.death)
    # Update age, smoking and quit years for next cycle
    data$age <- ifelse(is.na(data$age), NA, data$age + 1)
    data$smkyears <- ifelse(is.na(data$smkyears), NA, data$smkyears + ifelse(data$qtyears == 0, 1, 0)) # Add smoking year for current smokers
    data$qtyears <- ifelse(is.na(data$qtyears), NA, data$qtyears + ifelse(data$qtyears == 0, 0, 1)) # Add quit year for former smokers
  }
  
  # CUMULATIVE PROBABILITIES OF NO CANCER AND NO DEATH
  
  C.lung <- P.lung
  # Delete last column, put 1s in first column and shift other columns over by one
  C.lung[, 2:10] <- C.lung[, 1:9]
  C.lung[, 1] <- 1
  C.lung <- t(apply(C.lung, 1, cumprod))
  
  C.mort <- P.death
  C.mort[, 2:10] <- C.mort[, 1:9]
  C.mort[, 1] <- 1
  C.mort <- t(apply(C.mort, 1, cumprod))
  
  # The final risk calculation
  # ifelse(rowSums((1-P.lung)*P.death*C.lung*C.mort)==0,NA,rowSums((1-P.lung)*P.death*C.lung*C.mort))
  tmp1 <- (1 - P.lung) * P.death * C.lung * C.mort
  ret <- data.frame(number = 1:nrow(data))
  for (i in nyears) {
    ret.ny <- rowSums(tmp1[, 1:i, drop = FALSE])
    ret.ny <- ifelse(data$smkstatus == 1 | data$smkstatus == 2, ret.ny, NA)
    tmp2 <- ret.ny %in% 0
    if (any(tmp2)) ret.ny[tmp2] <- NA
    ret <- cbind(ret, ret.ny)
  }
  ret <- data.frame(ret[, -1])
  colnames(ret) <- paste0("bach.", nyears, "y")
  ret
}

#' Model: PLCOm2012
#' Reference：Tammemägi MC, Katki HA, Hocking WG, et al. Selection criteria for lung-cancer screening [published correction appears in N Engl J Med. 2013 Jul 25;369(4):394]. N Engl J Med. 2013;368(8):728-736. doi:10.1056/NEJMoa1211776
#' Published years: 2013
#' Country: US
#' Targeted population: Ever smokers
#' Prediction time range: 6 years
risk.plcom2012 <- function(data) {
  plcom2012.formula <- function(outcome) {
    # age = continuous
    # race = categorical
    # edu6 = ordinal
    # bmi = continuous
    # fam.lung.cancer = binary
    # prior.cancer = binary
    # copd = binary
    # cpd = cigarettes per day
    # current = binary
    # qtyears = years quit
    # smkyears = years smoked
    
    # Formula for PLCOm2012 lung cancer risk prediction model
    formula.string <- "I(age-62)+
                    black+hispanic+asian+indian+islander+
                    I(edu6-4)+
                    I(bmi-27)+
                    fam.lung.cancer+
                    prior.cancer+
                    copd+
                    current+
                    I(10/cpd-.4021541613)+
                    I(qtyears-10)+
                    I(smkyears-27)"
    
    formula(paste(outcome, "~", formula.string, collapse = ""))
  }
  
  # Coefficients for the logistic regression model based on the PLCOm2012 risk model
  plcom2012.coef <- c(
    intercept = -4.532506,
    age = 0.0778868,
    black = 0.3944778,
    hispanic = -0.7434744,
    asian = -0.466585,
    indian = 1.027152,
    islander = 0,
    edu6 = -0.0812744,
    bmi = -0.0274194,
    fam.lung.cancer = 0.587185,
    prior.cancer = 0.4589971,
    copd = 0.3553063,
    current = 0.02597431,
    cpd = -1.822606,
    qtyears = -0.0308572,
    smkyears = 0.0317321
  )
  
  # Create binary variables for race and smoking history based on input data
  data$black <- 1 * (data$race == 1)
  data$hispanic <- 1 * (data$race == 2)
  data$fam.lung.cancer <- 1 * (data$fam.lung.trend > 0)
  data$current <- 1 * (data$smkstatus == 2)
  data$edu6 <- ifelse(data$edu == 0, 1, data$edu)
  vars <- c(
    "age",
    "black", "hispanic", "asian", "islander", "indian",
    "edu6",
    "bmi",
    "fam.lung.cancer",
    "prior.cancer",
    "copd",
    "cpd",
    "current",
    "smkyears",
    "qtyears"
  )
  
  var.not.in <- names(data)[!(names(data) %in% vars) & sapply(data, function(x) sum(is.na(x))) == 0][1] # SELECT ARBITRARY VAR FOR OUTCOME
  
  # Generate the model matrix for logistic regression
  X <- model.matrix(plcom2012.formula(var.not.in), data)
  
  # Linear predictor (lp): log-odds of lung cancer, calculated for smokers
  lp <- ifelse(data$smkstatus == 1 | data$smkstatus == 2, ifelse(complete.cases(subset(data, select = vars)), X %*% plcom2012.coef, NA), NA)
  
  # Calculate the probability of lung cancer risk using logistic function
  r.plcom2012 <- data.frame(exp(lp) / (1 + exp(lp)))
  colnames(r.plcom2012) <- "plcom2012.6y"
  r.plcom2012
}

#' Model: Korean Men
#' Reference：Park S, Nam BH, Yang HR, et al. Individualized risk prediction model for lung cancer in Korean men. PLoS One. 2013;8(2):e54823. doi:10.1371/journal.pone.0054823
#' Published years: 2013
#' Country: Korea
#' Targeted population: General population
#' Prediction time range: 8 years
risk.koreanmen <- function(data) {
  # Calculate the linear predictor (S)
  S <- 0.1668 * (data$age - 44.66) -
    0.0020 * ((data$age - 44.66)^2 - 107.3007) +
    0.4180 * (I(data$smkstatus == 1) - 0.1511) +
    0.4444 * (I(data$smkstatus == 2 & data$cpd < 10) - 0.0905) +
    0.9414 * (I(data$smkstatus == 2 & data$cpd >= 10 & data$cpd < 20) - 0.3330) +
    1.3889 * (I(data$smkstatus == 2 & data$cpd >= 20) - 0.1390) +
    0.3306 * (I(data$bmi < 18.5) - 0.0239) -
    0.2468 * (I(data$bmi >= 23 & data$bmi < 25) - 0.2837) -
    0.3386 * (I(data$bmi >= 25) - 0.2851) -
    0.0909 * (I(data$activity == 1) - 0.1590) -
    0.1412 * (I(data$activity == 2) - 0.2952) -
    0.0521 * (I(data$activity == 3) - 0.0676) +
    0.0792 * (I(data$glucose >= 126) - 0.0605)
  
  S <- ifelse(data$smkstatus==0,S,S+0.2194 * (I(data$smoke.age.start >= 30 & data$smoke.age.start < 40) - 0.1441) +
                0.2809 * (I(data$smoke.age.start >= 19 & data$smoke.age.start < 30) - 0.3570) +
                0.5249 * (I(data$smoke.age.start >= 16 & data$smoke.age.start < 19) - 0.0164) +
                0.7120 * (I(data$smoke.age.start < 16) - 0.0064) )
  # Exclude females from the calculation as this model is for men only
  S <- ifelse(data$female == 1, NA, S)
  # Calculate the 8-year lung cancer risk
  r.koreanmen <- data.frame(1 - 0.9983894078**exp(S))
  colnames(r.koreanmen) <- "koreanmen.8y"
  r.koreanmen
}

#' Model: PLCOall2014
#' Reference：Tammemägi MC, Church TR, Hocking WG, et al. Evaluation of the lung cancer risks at which to screen ever- and never-smokers: screening rules applied to the PLCO and NLST cohorts [published correction appears in PLoS Med. 2015 Jan 28;12(1):e1001787. doi: 10.1371/journal.pmed.1001787]. PLoS Med. 2014;11(12):e1001764. Published 2014 Dec 2. doi:10.1371/journal.pmed.1001764
#' Published years: 2014
#' Country: US
#' Targeted population: General population
#' Prediction time range: 6 years
risk.plcoall2014 <- function(data) {
  plcoall2014.never.formula <- function(outcome) {
    # age = continuous
    # race = categorical
    # edu6 = ordinal
    # bmi = continuous
    # fam.lung.cancer = binary
    # prior.cancer = binary
    # copd = binary
    
    # Formula for never smokers
    formula.string <- "I(age-62)+
                    black+hispanic+asian+indian+islander+
                    I(edu6-4)+
                    I(bmi-27)+
                    fam.lung.cancer+
                    prior.cancer+
                    copd"
    
    formula(paste(outcome, "~", formula.string, collapse = ""))
  }
  
  # Formula for former smokers
  plcoall2014.former.formula <- function(outcome) {
    # age = continuous
    # race = categorical
    # edu6 = ordinal
    # bmi = continuous
    # fam.lung.cancer = binary
    # prior.cancer = binary
    # copd = binary
    # former = binary
    # cpd = cigarettes per day
    # qtyears = years quit
    # smkyears = years smoked
    formula.string <- "I(age-62)+
                    black+hispanic+asian+indian+islander+
                    I(edu6-4)+
                    I(bmi-27)+
                    fam.lung.cancer+
                    prior.cancer+
                    copd+
                    former+
                    I(100/cpd-4.021541613)+
                    I(qtyears-8.593417626)+
                    I(smkyears-27)"
    
    formula(paste(outcome, "~", formula.string, collapse = ""))
  }
  
  # Formula for current smokers
  plcoall2014.current.formula <- function(outcome) {
    # age = continuous
    # race = categorical
    # edu6 = ordinal
    # bmi = continuous
    # fam.lung.cancer = binary
    # prior.cancer = binary
    # copd = binary
    # current = binary
    # cpd = cigarettes per day
    # smkyears = years smoked
    formula.string <- "I(age-62)+
                    black+hispanic+asian+indian+islander+
                    I(edu6-4)+
                    I(bmi-27)+
                    fam.lung.cancer+
                    prior.cancer+
                    copd+
                    current+
                    I(100/cpd-4.021541613)+
                    I(smkyears-27)"
    
    formula(paste(outcome, "~", formula.string, collapse = ""))
  }
  
  # Coefficients for never smokers
  plcoall2014.never.coef <- c(
    intercept = -7.02198,
    age = 0.079597,
    black = 0.3211605,
    hispanic = -0.8203332,
    asian = -0.5241286,
    indian = 0.952699,
    islander = -1.364379,
    edu6 = -0.0879289,
    bmi = -0.028948,
    fam.lung.cancer = 0.5856777,
    prior.cancer = 0.4845208,
    copd = 0.3457265
  )
  
  # Coefficients for former smokers
  plcoall2014.former.coef <- c(
    intercept = -7.02198,
    age = 0.079597,
    black = 0.3211605,
    hispanic = -0.8203332,
    asian = -0.5241286,
    indian = 0.952699,
    islander = -1.364379,
    edu6 = -0.0879289,
    bmi = -0.028948,
    fam.lung.cancer = 0.5856777,
    prior.cancer = 0.4845208,
    copd = 0.3457265,
    former = 2.542472,
    cpd = -0.1815486,
    qtyears = -0.0321362,
    smkyears = 0.0305566
  )
  
  # Coefficients for current smokers
  plcoall2014.current.coef <- c(
    intercept = -7.02198,
    age = 0.079597,
    black = 0.3211605,
    hispanic = -0.8203332,
    asian = -0.5241286,
    indian = 0.952699,
    islander = -1.364379,
    edu6 = -0.0879289,
    bmi = -0.028948,
    fam.lung.cancer = 0.5856777,
    prior.cancer = 0.4845208,
    copd = -0.1815486,
    current = 2.799727,
    cpd = -0.1815486,
    smkyears = 0.0305566
  )
  
  # Create binary variables for race and smoking history based on input data
  data$black <- 1 * (data$race == 1)
  data$hispanic <- 1 * (data$race == 2)
  data$fam.lung.cancer <- 1 * (data$fam.lung.trend > 0)
  data$former <- 1 * (data$smkstatus == 1)
  data$current <- 1 * (data$smkstatus == 2)
  data$edu6 <- ifelse(data$edu == 0, 1, data$edu)
  vars <- c(
    "age",
    "black", "hispanic", "asian", "islander", "indian",
    "edu6",
    "bmi",
    "fam.lung.cancer",
    "prior.cancer",
    "copd",
    "cpd",
    "former",
    "current",
    "smkyears",
    "qtyears"
  )
  
  var.not.in <- names(data)[!(names(data) %in% vars) & sapply(data, function(x) sum(is.na(x))) == 0][1] # SELECT ARBITRARY VAR FOR OUTCOME
  
  # Linear predictor for never smokers, calculated only for complete cases
  lp.never <- ifelse(complete.cases(subset(data, select = c(
    "age", "black", "hispanic", "asian", "islander", "indian", "edu6", "bmi", "fam.lung.cancer",
    "prior.cancer", "copd"
  ))), model.matrix(plcoall2014.never.formula(var.not.in), data) %*% plcoall2014.never.coef, NA)
  # Linear predictor for former smokers, calculated only for complete cases
  lp.former <- ifelse(complete.cases(subset(data, select = c(
    "age", "black", "hispanic", "asian", "islander", "indian", "edu6", "bmi", "fam.lung.cancer",
    "prior.cancer", "copd", "cpd", "former", "smkyears", "qtyears"
  ))), model.matrix(plcoall2014.former.formula(var.not.in), data) %*% plcoall2014.former.coef, NA)
  # Linear predictor for current smokers, calculated only for complete cases
  lp.current <- ifelse(complete.cases(subset(data, select = c(
    "age", "black", "hispanic", "asian", "islander", "indian", "edu6", "bmi", "fam.lung.cancer",
    "prior.cancer", "copd", "cpd", "current", "smkyears"
  ))), model.matrix(plcoall2014.current.formula(var.not.in), data) %*% plcoall2014.current.coef, NA)
  
  # Combine the linear predictors based on smoking status: never (0), former (1), or current (2)
  lp <- ifelse(data$smkstatus == 0, lp.never, ifelse(data$smkstatus == 1, lp.former, ifelse(data$smkstatus == 2, lp.current, NA)))
  
  # Calculate the probability of lung cancer risk based on the linear predictor
  r.plcoall2014 <- data.frame(exp(lp) / (1 + exp(lp)))
  colnames(r.plcoall2014) <- "plcoall2014.6y"
  r.plcoall2014
}

#' Model: Pittsburgh Predictor
#' Reference：Wilson DO, Weissfeld J. A simple model for predicting lung cancer occurrence in a lung cancer screening program: The Pittsburgh Predictor. Lung Cancer. 2015;89(1):31-37. doi:10.1016/j.lungcan.2015.03.021
#' Published years: 2015
#' Country: US
#' Targeted population: Ever smokers
#' Prediction time range: 6 years
risk.pittsburgh <- function(data){
  S <- -10*I(data$smkyears<30)+8*I(data$smkyears>=40&data$smkyears<50)+14*I(data$smkyears>=50)+
    4*I(data$smkyears<30&data$age>=57)+4*I(data$smkyears>=30&data$smkyears<40&data$age>=59)+
    4*I(data$smkyears>=40&data$smkyears<50&data$age>=61)+4*I(data$smkyears>=50&data$age>=68)-
    3*I(data$qtyears>=1)-4*I(data$cpd<20)+2*I(data$cpd>=30&data$cpd<40)+5*I(data$cpd>=40)
  
  r.pittsburgh <- ifelse(data$smkstatus==1|data$smkstatus==2,as.numeric(1/(1+exp(-1*(-4.2195+0.10*S)))),NA)
  r.pittsburgh <- data.frame(r.pittsburgh)
  colnames(r.pittsburgh) <- "pittsburgh.6y"
  r.pittsburgh
}

#' Model: LLPi
#' Reference：Marcus MW, Chen Y, Raji OY, Duffy SW, Field JK. LLPi: Liverpool Lung Project Risk Prediction Model for Lung Cancer Incidence. Cancer Prev Res (Phila). 2015;8(6):570-575. doi:10.1158/1940-6207.CAPR-14-0438
#' Published years: 2015
#' Country: UK
#' Targeted population: General population
#' Prediction time range: 8.7 years
risk.llpi <- function(data) {
  data$male <- as.numeric(data$female == 0)
  # family history of lung cancer treated as late onset
  lp <- 0.036 * data$age + 0.391 * data$male + 0.043 * data$smkyears + 0.890 * data$copd + 1.044 * data$prior.cancer + 0.521 * I(data$fam.cancer.onset == 1) + 0.071 * I(data$fam.cancer.onset == 2)
  
  # 0.9728386 is the baseline survival at 8.7 years in the liverpool data
  # 3.556 is the mean linear predictor in the liverpool data
  r.llpi <- data.frame(1 - 0.9728386**(exp(lp - 3.556)))
  colnames(r.llpi) <- "llpi.8.7y"
  r.llpi
}

#' Model: LCRAT
#' Reference：Katki HA, Kovalchik SA, Berg CD, Cheung LC, Chaturvedi AK. Development and Validation of Risk Models to Select Ever-Smokers for CT Lung Cancer Screening. JAMA. 2016;315(21):2300-2311. doi:10.1001/jama.2016.6255
#' Published years: 2016
#' Country: US
#' Targeted population: Ever smokers
#' Prediction time range: 1-5 years
risk.lcrat <- function(data, nyears) {
  risk.kovalchik <- function(begin, end, newdata, coxph1, ...) {
    # Default value for end is 5
    if (!length(end)) end <- 5
    
    c.coxph.risk <- function(coxph, ...) {
      models <- list(coxph)
      if (!missing(..1)) {
        models <- c(models, list(...))
      }
      models
    }
    
    cuts <- function(data, npieces, simple.labels = TRUE, ...) {
      if (length(npieces) == 0 | any(is.na(npieces)) | any(!is.numeric(npieces))) {
        stop("npieces must be a numeric scalar or a vector of breakpoints")
      }
      
      if (length(npieces) == 1) {
        # Just break into quantiles.  Use quantile labelling: Q1,Q2,Q3,etc. if you want
        if (simple.labels) {
          cutdata <- cut(data,
                         breaks = quantile(data, seq(0, 1, 1 / npieces), na.rm = T, ...),
                         include.lowest = T, labels = paste("Q", 1:npieces, sep = ""), ...
          )
        } else {
          cutdata <- cut(data,
                         breaks = quantile(data, seq(0, 1, 1 / npieces), na.rm = T, ...),
                         include.lowest = T, ...
          )
        }
      } else
        # Break into pieces as specified by the breakpoints
        if (simple.labels) {
          cutdata <- cut(data,
                         breaks = npieces, include.lowest = T,
                         labels = paste("Q", 1:(length(npieces) - 1), sep = ""), ...
          )
        } else {
          cutdata <- cut(data, breaks = npieces, include.lowest = T, ...)
        }
      
      return(cutdata)
    }
    
    
    coxph.relrisk.uncentered <- function(coxph.object, newdata) {
      center <- coxph.object$means %*% coef(coxph.object)
      if (!missing(newdata)) {
        lp <- predict(coxph.object, newdata, type = "lp")
      } else {
        lp <- coxph.object$linear.predictor
      }
      exp(lp + rep(center, length(lp)))
    }
    
    projection.relrisk <- function(object, data) {
      if (is.numeric(object)) {
        return(object)
      } else if (class(object) == "coxph") {
        if (missing(data)) {
          return(coxph.relrisk.uncentered(object))
        } else {
          return(coxph.relrisk.uncentered(object, data))
        }
      } else {
        stop(cat("No method for class", class(object)))
      }
    }
    
    # calculate absolute risk given hazards/survival and relative risks
    risk.fixed.interval <- function(H, RR) {
      if (!is.list(H)) {
        0
      } else {
        absrisk <- H[[1]]$surv^RR[1] * H[[1]]$haz * RR[1]
        for (i in 2:length(H)) {
          absrisk <- absrisk * (H[[i]]$surv^RR[i])
        }
        sum(absrisk)
      }
    }
    
    models <- c.coxph.risk(coxph1, ...)
    
    # calculate relative risk for each subject and store as list
    rr <- sapply(models, projection.relrisk, data = newdata)
    
    if (is.matrix(rr)) {
      rr.list <- lapply(1:nrow(rr), function(x) rr[x, ])
    } else {
      rr.list <- list(rr)
    }
    
    # estimate risk
    if (length(begin) == 1) {
      AllVars <- unique(unlist(sapply(models, function(x) all.vars(x$formula))))
      in.interval <- function(x, begin, end) x >= begin & x <= end
      if ("lung.cancer.death" %in% AllVars) {
        H <- list(models[[1]]$basehaz[in.interval(models[[1]]$basehaz$time, begin, end), ], models[[2]]$basehaz_LCDRAT[in.interval(models[[1]]$basehaz$time, begin, end), ])
      } else if ("case" %in% AllVars) {
        H <- list(models[[1]]$basehaz[in.interval(models[[1]]$basehaz$time, begin, end), ], models[[2]]$basehaz_LCRAT[in.interval(models[[1]]$basehaz$time, begin, end), ])
      }
      # calculate absolute risk given hazards/survival for average covariate values and relative risks
      risks <- mapply(risk.fixed.interval, RR = rr.list, MoreArgs = list(H = H))
    } else {
      stop("No code for length(begin) != 1")
      # risks <- mapply(risk, begin = begin, end = end, RR = rr.list, MoreArgs = list(models = models))
    }
    
    ifelse(complete.cases(newdata), risks, NA)
  }
  data.lcrat <- data[, 1:11]
  
  colnames(data.lcrat)[10] <- c("edu6")
  
  data.lcrat$edu6 <- ifelse(data.lcrat$edu6 == 0, 1, data.lcrat$edu6)
  
  r.lcrat <- data.frame(number = 1:nrow(data))
  for (i in nyears) {
    r.LCDRAT <- risk.kovalchik(0, i, data.lcrat, LCDRAT, cox.death)
    r.LCRAT <- risk.kovalchik(0, i, data.lcrat, LCRAT, cox.death)
    r.lcrat.ny <- pmax(r.LCRAT, r.LCDRAT)
    r.lcrat <- cbind(r.lcrat, r.lcrat.ny)
  }
  r.lcrat <- data.frame(r.lcrat[, -1])
  colnames(r.lcrat) <- paste0("lcrat.", nyears, "y")
  r.lcrat
}
#' function "cuts" is needed for "risk.lcrat"
cuts <- function(data, npieces, simple.labels = TRUE, ...) {
  if (length(npieces) == 0 | any(is.na(npieces)) | any(!is.numeric(npieces))) {
    stop("npieces must be a numeric scalar or a vector of breakpoints")
  }
  
  if (length(npieces) == 1) {
    # Just break into quantiles.  Use quantile labelling: Q1,Q2,Q3,etc. if you want
    if (simple.labels) {
      cutdata <- cut(data,
                     breaks = quantile(data, seq(0, 1, 1 / npieces), na.rm = T, ...),
                     include.lowest = T, labels = paste("Q", 1:npieces, sep = ""), ...
      )
    } else {
      cutdata <- cut(data,
                     breaks = quantile(data, seq(0, 1, 1 / npieces), na.rm = T, ...),
                     include.lowest = T, ...
      )
    }
  } else
    # Break into pieces as specified by the breakpoints
    if (simple.labels) {
      cutdata <- cut(data,
                     breaks = npieces, include.lowest = T,
                     labels = paste("Q", 1:(length(npieces) - 1), sep = ""), ...
      )
    } else {
      cutdata <- cut(data, breaks = npieces, include.lowest = T, ...)
    }
  
  return(cutdata)
}

#' Model: HUNT
#' Reference：Markaki M, Tsamardinos I, Langhammer A, Lagani V, Hveem K, Røe OD. A Validated Clinical Risk Prediction Model for Lung Cancer in Smokers of All Ages and Exposure Types: A HUNT Study [published correction appears in EBioMedicine. 2022 Aug;82:104187. doi: 10.1016/j.ebiom.2022.104187]. EBioMedicine. 2018;31:36-46. doi:10.1016/j.ebiom.2018.03.027
#' Published years: 2018
#' Country: Norway
#' Targeted population: Ever smokers
#' Prediction time range: 6 years
risk.hunt <- function(data) {
  # Create a binary variable 'male', where 1 indicates male and 0 indicates female
  data$male <- 1 * (data$female == 0)
  # Calculate the linear predictor (S)
  S <- 1.18203062 + 0.31573217 * data$male - 1.98496138 * 100 / data$age + 1.11994217 * log(data$pkyr.cat + 1) - 0.04002877 * data$cpd -
    0.24019955 * log(data$qtyears + 1) - 1.70238304 * log(data$bmi + 1) + 0.0807242 * log(data$smkexp + 1) + 0.49212668 * (data$cough.daily)
  # Calculate the 6-year lung cancer risk using the logistic function
  # If the person is a smoker (current or former), calculate the risk; otherwise, set it to NA
  r.hunt <- data.frame(ifelse(data$smkstatus == 1 | data$smkstatus == 2, 1 / (1 + exp(-S)), NA))
  colnames(r.hunt) <- "hunt.6y"
  r.hunt
}

#' Model: LLPv3
#' Reference：Field JK, Vulkan D, Davies MPA, Duffy SW, Gabe R. Liverpool Lung Project lung cancer risk stratification model: calibration and prospective validation. Thorax. 2021;76(2):161-168. doi:10.1136/thoraxjnl-2020-215158
#' Published years: 2021
#' Country: UK
#' Targeted population: General population
#' Prediction time range: 5 years
risk.llpv3 <- function(data){
  llpv3.formula <- function(outcome){
    # pneu = any diagnosis of pneumonia
    # asb = asbestos exposure binary indicator
    # prior.cancer = any prior cancer
    # fam.cancer.onset = 1 or more first degree ; early onset (3 categories)
    # smkyears.cat = 5 categories
    formula.string <- "pneu+asb+prior.cancer+fam.cancer.onset+smkyears.cat"
    formula(paste(outcome,"~",formula.string,collapse=""))
  }
  
  llpv3.coef <- c(intercept =1,
                  pneu=0.6025,
                  asb=0.6343,
                  prior.cancer=0.6754,
                  fam.cancer.onset1=0.7034,
                  fam.cancer.onset2=0.1677,
                  smkyears.cat1=0.7692,
                  smkyears.cat2=1.4516,
                  smkyears.cat3=2.5072,
                  smkyears.cat4=2.7243)
  
  # To make the probabilities interpretable, they used Liverpool data to estimate the beta0, intercept here is from model
  llpv3.intercepts <- data.frame(
    intercept = -c(9.84, 8.94, 8.09, 7.41, 6.75, 6.34, 6.09, 5.61, 5.46,
                   10.37, 8.53, 7.93, 6.97, 6.69, 6.46, 5.96, 5.70, 5.89),
    age = rep(1:9,2),
    female = rep(c(0,1),each=9))
  
  llpv3.intercept <- function(age, female){
    #Assign age groups
    age.group <- function(age){
      if(age<45) 1
      else if(age<50) 2
      else if(age<55) 3
      else if(age<60) 4
      else if(age<65) 5
      else if(age<70) 6
      else if(age<75) 7
      else if(age<80) 8
      else 9
    }
    
    value <- function(age, age.cat, female){
      if(age.cat==1)
        llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                     llpv3.intercepts$age==age.cat]
      else if(age.cat==2)
        sum(c(50-age-0.5, 0.5+age-45)*llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                                                   llpv3.intercepts$age>=2&llpv3.intercepts$age<=3])/5
      else if(age.cat==3)
        sum(c(55-age-0.5, 0.5+age-50)*llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                                                   llpv3.intercepts$age>=3&llpv3.intercepts$age<=4])/5
      else if(age.cat==4)
        sum(c(60-age-0.5, 0.5+age-55)*llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                                                   llpv3.intercepts$age>=4&llpv3.intercepts$age<=5])/5
      else if(age.cat==5)
        sum(c(65-age-0.5, 0.5+age-60)*llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                                                   llpv3.intercepts$age>=5&llpv3.intercepts$age<=6])/5
      else if(age.cat==6)
        sum(c(70-age-0.5, 0.5+age-65)*llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                                                   llpv3.intercepts$age>=6&llpv3.intercepts$age<=7])/5
      else if(age.cat==7)
        sum(c(75-age-0.5, 0.5+age-70)*llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                                                   llpv3.intercepts$age>=7&llpv3.intercepts$age<=8])/5
      else if(age.cat==8)
        sum(c(80-age-0.5, 0.5+age-75)*llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                                                   llpv3.intercepts$age>=8&llpv3.intercepts$age<=9])/5
      else
        llpv3.intercepts$intercept[llpv3.intercepts$female==female&
                                     llpv3.intercepts$age==age.cat]
    }
    
    value(age, age.group(age), female)
  }
  
  vars <- c("female",
            "age",
            "asb",
            "smkyears",
            "pneu",
            "prior.cancer",
            "fam.cancer.onset")
  
  var.not.in <- names(data)[!(names(data)%in%vars) & sapply(data, function(x) sum(is.na(x)))==0][1] # SELECT ARBITRARY VAR FOR OUTCOME
  expit <- function(x) exp(x)/(1+exp(x))
  data$smkyears.cat <- cut(data$smkyears, c(-1,1,20,40,60,max(data$smkyears,na.rm=TRUE)+1), right=FALSE, lab=1:5)
  #Manually construct the intercept
  intercept <- mapply(llpv3.intercept, age=data$age, female=data$female)
  X <- model.matrix(llpv3.formula(var.not.in), data)
  X[,1] <- intercept[complete.cases(subset(data,select= c("asb",
                                                          "smkyears",
                                                          "pneu",
                                                          "prior.cancer",
                                                          "fam.cancer.onset")))]
  
  r.llpv3 <- ifelse(complete.cases(subset(data,select=c("asb",
                                                        "smkyears",
                                                        "pneu",
                                                        "prior.cancer",
                                                        "fam.cancer.onset"))),expit(X%*%llpv3.coef),NA)
  r.llpv3 <- data.frame(r.llpv3)
  colnames(r.llpv3) <- "llpv3.5y"
  r.llpv3
}

#' Model: LCRS
#' Reference：Ma Z, Lv J, Zhu M, et al. Lung cancer risk score for ever and never smokers in China. Cancer Commun (Lond). 2023;43(8):877-895. doi:10.1002/cac2.12463
#' Published years: 2023
#' Country: China
#' Targeted population: General population
#' Prediction time range: 3/5/6/10 years
risk.lcrs <- function(data, nyears) {
  # Calculate LCRS for ever smokers (current or former smokers)
  lcrs.ever <- ifelse(data$smkstatus == 1 | data$smkstatus == 2,
                      0.559 * I(data$age >= 40 & data$age < 45) +
                        1.023 * I(data$age >= 45 & data$age < 50) +
                        1.565 * I(data$age >= 50 & data$age < 55) +
                        1.906 * I(data$age >= 55 & data$age < 60) +
                        2.182 * I(data$age >= 60 & data$age < 65) +
                        2.525 * I(data$age >= 65 & data$age < 70) +
                        2.659 * I(data$age >= 70) +
                        0.417 * data$urban +
                        0.393 * I(data$edu >= 2 & data$edu < 5) +
                        0.455 * I(data$edu <= 1) +
                        0.022 * I(data$height >= 160 & data$height < 165) +
                        0.128 * I(data$height >= 165 & data$height < 170) +
                        0.135 * I(data$height >= 170) +
                        0.201 * I(data$bmi >= 18.5 & data$bmi < 24) +
                        0.552 * I(data$bmi < 18.5) +
                        0.282 * data$frequent.cough +
                        0.332 * data$copd +
                        0.589 * data$prior.cancer +
                        0.300 * I(data$fam.lung.trend >= 2) +
                        0.122 * I(data$cpd >= 10 & data$cpd < 15) +
                        0.232 * I(data$cpd >= 15 & data$cpd < 20) +
                        0.402 * I(data$cpd >= 20 & data$cpd < 25) +
                        0.546 * I(data$cpd >= 25 & data$cpd < 30) +
                        0.582 * I(data$cpd >= 30 & data$cpd < 35) +
                        0.636 * I(data$cpd >= 35) +
                        0.345 * I(data$smkyears >= 10 & data$smkyears < 20) +
                        0.495 * I(data$smkyears >= 20 & data$smkyears < 30) +
                        0.773 * I(data$smkyears >= 30 & data$smkyears < 40) +
                        1.046 * I(data$smkyears >= 40 & data$smkyears < 50) +
                        1.199 * I(data$smkyears >= 50) +
                        0.220 * data$smk.into.lung +
                        0.188 * I(data$qtyears <= 5), NA
  )
  
  # Calculate LCRS for never smokers
  lcrs.never <- ifelse(data$smkstatus == 0,
                       0.802 * I(data$age >= 40 & data$age < 45) +
                         1.184 * I(data$age >= 45 & data$age < 50) +
                         1.596 * I(data$age >= 50 & data$age < 55) +
                         1.890 * I(data$age >= 55 & data$age < 60) +
                         2.364 * I(data$age >= 60 & data$age < 65) +
                         2.570 * I(data$age >= 65 & data$age < 70) +
                         2.615 * I(data$age >= 70) +
                         0.174 * data$urban +
                         0.128 * I(data$height >= 150 & data$height < 155) +
                         0.142 * I(data$height >= 155 & data$height < 160) +
                         0.229 * I(data$height >= 160) +
                         0.159 * I(data$bmi >= 18.5 & data$bmi < 24) +
                         0.371 * I(data$bmi < 18.5) +
                         0.122 * I(data$low.activity) +
                         0.186 * data$frequent.cough +
                         0.395 * data$copd +
                         0.792 * data$prior.cancer +
                         0.277 * I(data$fam.lung.trend >= 2), NA
  )
  
  # Combine lcrs for ever smokers and never smokers based on smoking status
  lcrs <- ifelse(data$smkstatus == 0, lcrs.never,
                 ifelse(data$smkstatus == 1 | data$smkstatus == 2, lcrs.ever, NA)
  )
  
  # Calculate the hazard ratio for lung cancer risk based on lcrs
  hr.lcrs <- ifelse(data$smkstatus == 0, exp(lcrs - 1.795),
                    ifelse(data$smkstatus == 1 | data$smkstatus == 2, exp(lcrs - 3.460), NA)
  )
  
  # Calculate n-year risk
  r.lcrs <- data.frame(number = 1:nrow(data))
  if (3 %in% nyears) {
    r.lcrs.3y <- ifelse(data$smkstatus == 0, 1 - 0.998941**hr.lcrs,
                        ifelse(data$smkstatus == 1 | data$smkstatus == 2, 1 - 0.997501**hr.lcrs, NA)
    )
    r.lcrs <- cbind(r.lcrs, lcrs.3y = r.lcrs.3y)
  }
  if (5 %in% nyears) {
    r.lcrs.5y <- ifelse(data$smkstatus == 0, 1 - 0.998060**hr.lcrs,
                        ifelse(data$smkstatus == 1 | data$smkstatus == 2, 1 - 0.995305**hr.lcrs, NA)
    )
    r.lcrs <- cbind(r.lcrs, lcrs.5y = r.lcrs.5y)
  }
  if (6 %in% nyears) {
    r.lcrs.6y <- ifelse(data$smkstatus == 0, 1 - 0.997527**hr.lcrs,
                        ifelse(data$smkstatus == 1 | data$smkstatus == 2, 1 - 0.993970**hr.lcrs, NA)
    )
    r.lcrs <- cbind(r.lcrs, lcrs.6y = r.lcrs.6y)
  }
  if (10 %in% nyears) {
    r.lcrs.10y <- ifelse(data$smkstatus == 0, 1 - 0.994893**hr.lcrs,
                         ifelse(data$smkstatus == 1 | data$smkstatus == 2, 1 - 0.987736**hr.lcrs, NA)
    )
    r.lcrs <- cbind(r.lcrs, lcrs.10y = r.lcrs.10y)
  }
  # Final data frame with calculated risks for specified years
  r.lcrs <- data.frame(r.lcrs[, -1])
  colnames(r.lcrs) <- paste0("lcrs.", nyears, "y")
  r.lcrs
}

#' Model: OWL,
#' Reference：Pan Z, Zhang R, Shen S, et al. OWL: an optimized and independently validated machine learning prediction model for lung cancer screening based on the UK Biobank, PLCO, and NLST populations. EBioMedicine. 2023;88:104443. doi:10.1016/j.ebiom.2023.104443
#' Published years: 2023
#' Country: UK
#' Targeted population: General population
#' Prediction time range: 1-8 years
risk.owl <- function(data, nyears) {
  # Load two packages.
  if (!require(Matrix)) {
    install.packages("Matrix")
  } else {
    library(Matrix)
  }
  
  if (!require(xgboost)) {
    install.packages("xgboost")
  } else {
    library(xgboost)
  }
  
  owl.data <- data[, 1:15]
  colnames(owl.data) <- extra_data$vindependents
  owl.data$age <- data$age
  owl.data$gender <- ifelse(data$female == 0, 1, 0)
  owl.data$edu_level <- ifelse(data$edu == 0, 1, data$edu)
  owl.data$bmi <- data$bmi
  owl.data$smoke_status <- data$smkstatus
  owl.data$age_begin_smoke <- ifelse(owl.data$smoke_status == 0, 99, data$smoke.age.start)
  owl.data$smoke_last_yr <- ifelse(owl.data$smoke_status == 0, 0, as.numeric(data$smkyears) + 1)
  owl.data$smoke_quit_year <- ifelse(owl.data$smoke_status == 1, data$qtyears, 0)
  owl.data$ncig <- ifelse(owl.data$smoke_status == 0, 0, data$cpd)
  owl.data$packyr <- ifelse(owl.data$smoke_status == 0, 0, owl.data$ncig * owl.data$smoke_last_yr / 20)
  owl.data$diabetes <- data$diabetes
  owl.data$copd <- data$copd
  owl.data$emphysema <- data$emp
  owl.data$bronchitis <- data$bron
  owl.data$FDRLC <- data$fam.smoke.cancer
  
  vindependets <- extra_data$vindependents
  beta <- extra_data$beta
  baseline_hazard <- extra_data$baseline_hazard
  
  xgb.owl <- list(data = Matrix(data.matrix(owl.data[, vindependets]),
                                sparse = TRUE
  ))
  
  risk <- predict(owl,
                  newdata = xgb.DMatrix(data = xgb.owl$data)
  )
  
  # Calculate the OWL risk score
  riskscore <- beta * log(risk)
  
  # Calculate the OWL absolute risk within 8-years follow-up.
  surv_prob <- exp(exp(riskscore) %*% -matrix(baseline_hazard[, 1], nrow = 1))
  r.owl <- data.frame(1 - surv_prob[, findInterval(nyears, baseline_hazard[, 2])])
  
  colnames(r.owl) <- paste0("owl.", nyears, "y")
  
  r.owl
}

# data form--------------------------------
#' The data must be a data frame or matrix containing individuals' covariate values.
#'  Covariates should be in the following columns and numerical formats:
#'
#'  \itemize{
#'  \item column 1 - current age (numeric);
#'  \item column 2 - gender (1=Female, 0=Male);
#'  \item column 3 - smoke status (0=never smoker, 1=former smoker, 2=current smoker)
#'  \item column 4 - years smoked (numeric);
#'  \item column 5 - years quit (numeric);
#'  \item column 6 - cigarettes per day (numeric);
#'  \item column 7 - race (0=Non-hispanic white,
#'                   1=Non-hispanic Black/African American, 
#'                   2=Hispanic, 
#'                   3=Other Ethnicity);
#'  \item column 8 - lung disease (1=Emphysema, 0=No Emphysema);
#'  \item column 9 - number of first degree relatives with lung cancer (0,1,2);
#'  \item column 10 - BMI;
#'  \item column 11 - highest education level(0=elementary school or lower,
#'                    1=middle school, 
#'                    2=HS graduate, 
#'                    3=post hs, no college, 
#'                    4=associate degree/some college, 
#'                    5=bachelors degree,
#'                    6=graduate school);
#'  \item column 12 - asbestos exposure (1=Yes,0=No);
#'  \item column 13 - prior history of pneumonia (1=Yes,0=No);
#'  \item column 14 - prior history of cancer (1=Yes,0=No);
#'  \item column 15 - family history of lung cancer (0=none, 1=early onset, 2=late onset);
#'  \item column 16 - prior history of copd (1=Yes,0=No);
#'  \item column 17 - Dust exposure  (1=Yes,0=No);
#'  \item column 18 - 2 or more first degree relatives with cancer (binary indicator);
#'  \item column 19 - 1 or more first degree relatives with smoking cancer (binary indicator);
#'  \item column 20 - no hay fever (1=No Hay Fever,0=Yes Hay Fever);
#'  \item column 21 - asian ethnicity (1=Yes,0=No);
#'  \item column 22 - islander ethnicity (1=Yes,0=No);
#'  \item column 23 - American indian ethnicity (1=Yes,0=No);
#'  \item column 24 - environment tabacco smoke (1=Yes,0=No);
#'  \item column 25 - indoor smoke exposure in hours (numeric);
#'  \item column 26 - History of emphysema and/or chronic bronchitis (1=Yes,0=No);
#'  \item column 27 - physical activity (0=No;
#'                    1 = ≤4 times/week at <30 minutes/session; 
#'                    2 = 2-4 times/week at ≥30 minutes/session or ≥5 times/week at <30 minutes/session; 
#'                    3 = ≥5 times/week at ≥30 minutes/session);
#'  \item column 28 - exercising for less than a minimum of 30 minutes 3 times a week(0=No; 1=Yes);
#'  \item column 29 - prior history of diabetes (1=Yes,0=No);
#'  \item column 30 - fasting glucose levels, mg/dl (numeric);
#'  \item column 31 - cough daily during periods of the year (0=No; 1=Yes);
#'  \item column 32 - height, cm (numeric);
#'  \item column 33 - residential area (0=rural; 1=urban);
#'  \item column 34 - frequent coughing during the day or at night (lasting 3 months or more) in the past 12 months (0=No; 1=Yes);
#'  \item column 35 - smoking inhalation to the lungs (0=No; 1=Yes);
#'  }



# main function-----------------------------------

lcrisks <- function(data,
                    models = "all",
                    nyears.bach = 10,
                    nyears.lcrat = 5,
                    nyears.lcrs = 10,
                    nyears.owl = 8) {
  
  x<-data
  if (!length(models) | (length(models)==1 && models == "all")) {
    models <- c(
      "bach", "llpv3", "llpi", "plcom2012", "plcoall2014", 
      "pittsburgh", "lcrat", "koreanmen", "hunt", "lcrs", "owl"
    )
  }
  
  # Extract and process the columns from the input data
  age <- ifelse(is.na(x[, 1]) == 0 & x[, 1] > 0, x[, 1], NA)
  female <- ifelse(is.na(x[, 2]) == 0 & x[, 2] %in% c(0, 1), x[, 2], NA)
  smkstatus <- ifelse(is.na(x[, 3]) == 0 & x[, 3] %in% c(0, 1, 2), x[, 3], NA)
  smkyears <- ifelse(is.na(x[, 4]) == 0 & x[, 4] >= 0, x[, 4], NA)
  smkyears[!is.na(smkstatus) & smkstatus == 0] <- 0
  qtyears <- ifelse(is.na(x[, 5]) == 0 & x[, 5] >= 0, x[, 5], NA)
  qtyears[!is.na(smkstatus) & (smkstatus == 2)] <- 0
  qtyears[!is.na(smkstatus) & (smkstatus == 0)] <- NA
  cpd <- ifelse(is.na(x[, 6]) == 0 & x[, 6] >= 0, x[, 6], NA)
  cpd[!is.na(smkstatus) & smkstatus == 0] <- 0
  race <- as.factor(ifelse(is.na(x[, 7]) == 0 & x[, 7] %in% c(0, 1, 2, 3), x[, 7], NA))
  emp <- ifelse(is.na(x[, 8]) == 0 & x[, 8] %in% c(0, 1), x[, 8], NA)
  fam.lung.trend <- ifelse(is.na(x[, 9]) == 0 & x[, 9] %in% c(0, 1, 2), x[, 9], NA)
  bmi <- ifelse(is.na(x[, 10]) == 0 & x[, 10] > 0, x[, 10], NA)
  edu <- ifelse(is.na(x[, 11]) == 0 & x[, 11] %in% c(0, 1, 2, 3, 4, 5, 6), x[, 11], NA)
  pkyr.cat <- smkyears * cpd / 20
  age.stopped <- ifelse(is.na(age - qtyears) == 0 & age - qtyears >= 0, age - qtyears, NA)
  smoke.age.start <- ifelse(is.na(age - smkyears - qtyears) == 0 & age - smkyears - qtyears >= 0, age - smkyears - qtyears, NA)
  asb <- ifelse(is.na(x[, 12]) == 0 & x[, 12] %in% c(0, 1), x[, 12], NA)
  pneu <- ifelse(is.na(x[, 13]) == 0 & x[, 13] %in% c(0, 1), x[, 13], NA)
  prior.cancer <- ifelse(is.na(x[, 14]) == 0 & x[, 14] %in% c(0, 1), x[, 14], NA)
  fam.cancer.onset <- as.factor(ifelse(is.na(x[, 15]) == 0 & x[, 15] %in% c(0, 1, 2), x[, 15], NA))
  copd <- ifelse(is.na(x[, 16]) == 0 & x[, 16] %in% c(0, 1), x[, 16], NA)
  dust <- ifelse(is.na(x[, 17]) == 0 & x[, 17] %in% c(0, 1), x[, 17], NA)
  fam.cancer <- ifelse(is.na(x[, 18]) == 0 & x[, 18] %in% c(0, 1), x[, 18], NA)
  fam.smoke.cancer <- ifelse(is.na(x[, 19]) == 0 & x[, 19] %in% c(0, 1), x[, 19], NA)
  no.hayfever <- ifelse(is.na(x[, 20]) == 0 & x[, 20] %in% c(0, 1), x[, 20], NA)
  asian <- ifelse(is.na(x[, 21]) == 0 & x[, 21] %in% c(0, 1), x[, 21], NA)
  islander <- ifelse(is.na(x[, 22]) == 0 & x[, 22] %in% c(0, 1), x[, 22], NA)
  indian <- ifelse(is.na(x[, 23]) == 0 & x[, 23] %in% c(0, 1), x[, 23], NA)
  ets <- ifelse(is.na(x[, 24]) == 0 & x[, 24] %in% c(0, 1), x[, 24], NA)
  smkexp <- ifelse(is.na(x[, 25]) == 0 & x[, 25] >= 0, x[, 25], NA)
  bron <- ifelse(is.na(x[, 26]) == 0 & x[, 26] %in% c(0, 1), x[, 26], NA)
  activity <- ifelse(is.na(x[, 27]) == 0 & x[, 27] %in% c(0, 1, 2, 3), x[, 27], NA)
  low.activity <- ifelse(is.na(x[, 28]) == 0 & x[, 28] %in% c(0, 1), x[, 28], NA)
  diabetes <- ifelse(is.na(x[, 29]) == 0 & x[, 29] %in% c(0, 1), x[, 29], NA)
  glucose <- ifelse(is.na(x[, 30]) == 0 & x[, 30] > 0, x[, 30], NA)
  cough.daily <- ifelse(is.na(x[, 31]) == 0 & x[, 31] %in% c(0, 1), x[, 31], NA)
  height <- ifelse(is.na(x[, 32]) == 0 & x[, 32] > 0, x[, 32], NA)
  urban <- ifelse(is.na(x[, 33]) == 0 & x[, 33] %in% c(0, 1), x[, 33], NA)
  frequent.cough <- ifelse(is.na(x[, 34]) == 0 & x[, 34] %in% c(0, 1), x[, 34], NA)
  smk.into.lung <- ifelse(is.na(x[, 35]) == 0 & x[, 35] %in% c(0, 1), x[, 35], NA)
  
  # Combine the processed variate into a data frame
  covar <- data.frame(
    age = age,
    bmi = bmi,
    cpd = cpd,
    emp = emp,
    fam.lung.trend = fam.lung.trend,
    female = female,
    qtyears = qtyears,
    smkyears = smkyears,
    race = race,
    edu = edu,
    pkyr.cat = pkyr.cat,
    asb = asb,
    smoke.age.start = smoke.age.start,
    pneu = pneu,
    prior.cancer = prior.cancer,
    fam.cancer.onset = fam.cancer.onset,
    copd = copd,
    dust = dust,
    fam.cancer = fam.cancer,
    fam.smoke.cancer = fam.smoke.cancer,
    no.hayfever = no.hayfever,
    age.stopped = age.stopped,
    asian = asian,
    islander = islander,
    indian = indian,
    smkstatus = smkstatus,
    ets = ets,
    smkexp = smkexp,
    bron = bron,
    activity = activity,
    low.activity = low.activity,
    diabetes = diabetes,
    glucose = glucose,
    cough.daily = cough.daily,
    height = height,
    urban = urban,
    frequent.cough = frequent.cough,
    smk.into.lung = smk.into.lung
  )
  
  
  out <- data.frame(number = 1:nrow(covar))
  # Define a list of models
  model.list <- list(
    bach = list(fn = risk.bach, args = list(nyears = nyears.bach)),
    llpv3 = list(fn = risk.llpv3, args = list()),
    llpi = list(fn = risk.llpi, args = list()),
    plcom2012 = list(fn = risk.plcom2012, args = list()),
    plcoall2014 = list(fn = risk.plcoall2014, args = list()),
    pittsburgh = list(fn = risk.pittsburgh, args = list()),
    lcrat = list(fn = risk.lcrat, args = list(nyears = nyears.lcrat)),
    koreanmen = list(fn = risk.koreanmen, args = list()),
    hunt = list(fn = risk.hunt, args = list()),
    lcrs = list(fn = risk.lcrs, args = list(nyears = nyears.lcrs)),
    owl = list(fn = risk.owl, args = list(nyears = nyears.owl))
  )
  
  # Define model names used in the function
  model.names <- c(
    "bach", "llpv3", "llpi", "plcom2012", "plcoall2014", 
    "pittsburgh", "lcrat", "koreanmen", "hunt", "lcrs", "owl"
  )
  
  # Iterate through the models, call the appropriate function and merge the results
  for (model in model.names) {
    if (model %in% models) {
      model.fn <- model.list[[model]]$fn
      model.args <- model.list[[model]]$args
      model.result <- do.call(model.fn, c(list(covar), model.args))
      out <- cbind(out, model.result)
    }
  }
  out[, -1]
}

# example------------------------
age <- c(66, 58, 75, 72, 56, 50, 65)
bmi <- c(23, 28, 26, 27, 24, 23, 26)
smkstatus <- c(2, 2, 1, 1, 1, 2, 1)
cpd <- c(36, 36, 40, 24, 40, 10, 5)
emp <- c(0, 1, 1, 0, 1, 1, 0)
copd <- c(0, 1, 1, 0, 1, 1, 0)
fam.lung.trend <- c(0, 2, 0, 2, 0, 2, 0)
female <- c(0, 1, 0, 1, 0, 1, 0)
smkyears <- c(43, 37, 45, 42, 29, 15, 30)
qtyears <- c(0, 0, 9, 6, 6, 0, 5)
race <- c(0, 1, 2, 2, 3, 3, 3)
edu6 <- c(3, 5, 4, 5, 5, 4, 2)
asb <- c(0, 0, 0, 0, 0, 0, 0)
pneu <- c(0, 0, 0, 0, 0, 0, 0)
prior.cancer <- c(0, 0, 0, 0, 0, 0, 0)
fam.cancer.onset <- c(0, 1, 0, 2, 0, 0, 0)
dust <- c(0, 0, 0, 0, 0, 0, 0)
fam.cancer <- c(0, 1, 0, 1, 0, 0, 0)
fam.smoke.cancer <- c(0, 1, 0, 1, 0, 0, 0)
no.hayfever <- c(1, 1, 1, 1, 1, 1, 0)
asian <- c(0, 0, 0, 0, 1, 1, 1)
islander <- c(0, 0, 0, 0, 0, 0, 0)
indian <- c(0, 0, 0, 0, 0, 0, 0)
ets <- c(0, 0, 0, 0, 0, 0, 0)
smkexp <- c(2, 3, 5, 1, 0, 3, 2)
bron <- c(0, 0, 0, 0, 0, 0, 0)
activity <- c(1, 3, 2, 1, 2, 1, 1)
low.activity <- c(0, 1, 1, 0, 0, 0, 0)
diabetes <- c(1, 0, 0, 1, 1, 0, 1)
glucose <- c(130, 115, 101, 145, 126, 120, 135)
cough.daily <- c(0, 1, 0, 1, 0, 1, 1)
height <- c(160, 170, 183, 164, 172, 169, 180)
urban <- c(1, 1, 0, 1, 0, 0, 0)
frequent.cough <- c(0, 1, 0, 1, 0, 1, 1)
smk.into.lung <- c(0, 0, 1, 1, 0, 1, 1)

persons <- data.frame(
  age,
  female,
  smkstatus,
  smkyears,
  qtyears,
  cpd,
  race,
  emp,
  fam.lung.trend,
  bmi,
  edu6,
  asb,
  pneu,
  prior.cancer,
  fam.cancer.onset,
  copd,
  dust,
  fam.cancer,
  fam.smoke.cancer,
  no.hayfever,
  asian,
  islander,
  indian,
  ets,
  smkexp,
  bron,
  activity,
  low.activity,
  diabetes,
  glucose,
  cough.daily,
  height,
  urban,
  frequent.cough,
  smk.into.lung
)


lcrisks(persons,
        models = "all", 
        nyears.bach = c(5, 7, 9), 
        nyears.lcrat = 3,
        nyears.lcrs = 10, 
        nyears.owl = c(1:8)
)