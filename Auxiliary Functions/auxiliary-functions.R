################################################
# funcionts to calculate DALYS:

auxburden <- function(x,K = 0,C = .1658, beta = .04, r = .03, a){
  
  ## Auxiliary function to evaluate a DALY in one point, by Murray and Lopez
  
  ### Input:
  
  # K: binary constant: 1 if weighting by age applies, 0 if not
  # r: discount rate. WHO recommends 3%
  # a: age to which the disease will be assigned
  # x: current age, integration variable
  # C and beta are constants, see Murray and Lopez for details
  
  ### Output: 
  
  # Function evaluated to approximate the integral 
  
  K*C*x*exp(-beta*x)*exp(-r*(x-a))+(1-K)*exp(-r*(x-a))
  
}

burden <- function(N,DW,A,L,K= 0,r = .03,a){
 
  ## Calculate a  DALY during the course of the disease 
  
  ### Input:
  
  # N: number of people 
  # DW: disability weights
  # A: age at onset of the disease or age at death
  # L: duration of the disease or expected life years remaining at age of death
  # a: age to which the disease will be assigned
  
  N*DW*integrate(auxburden,lower = A,upper = A+L, K=K, r=r, a = a)$value
}

################################################
# Detection:
# function to determine the number of women detected in a determined program

detection <- function(percentage,periodicity,incidence,conapo.incidence,mean.sojourn.time,sensitivity,start = 2016,end = 2026,age.group.data.frame){
  
  ### Input:
  # percentage = coverage percentage by age groups
  # periodicity = periodicity of screening programs
  # indicencia = rate of incidence by age groups
  # conapo.incidence = population data set by group ages and new cases per year
  # mean.sojourn.time = mean sojourn time by age groups
  # sensitivity = sensitivity of the mammography by group ages
  # start  = starting year of the intervention 
  # end  = ending year of the intervention
  
  ### Output:
  # we = women examined by age group and year
  # det = women detected by age group and year
  
  ################################################
  
  # Years when the program takes place
  x <- seq(from = start, to = end, by = periodicity )
  
  # Define Lambda for sensitivity 
  lambda <- data.frame(lambda = 1/mean.sojourn.time$mean.sojourn.time, 
                       age.groups = age.group.data.frame)
  
  # Axuiliary function for adjusted sensitivity:
  f <- function(lambda,sensitivity,periodicity){
    sensitivity.adjusted <- (sensitivity*(1-exp(-lambda*periodicity)))/(lambda*periodicity*(1-(1-sensitivity)*exp((-lambda)*periodicity)))
  }
  
  # rate of detection for every 100,000 women
  for(j in 1:length(x)){
    if(j ==1 )
    {
      # prevalence
      rate.detection <- mean.sojourn.time$mean.sojourn.time*incidence$incidence*sensitivity$sensitivity 
    }
    else{
      # incidence
      sensitivity.mammography <- f(lambda$lambda,sensitivity$sensitivity,periodicity) 
      rate.detection <- cbind(rate.detection,periodicity*incidence$incidence*sensitivity.mammography)
    }
  }
  
  row.names(rate.detection) <- lambda$age.groups
  colnames(rate.detection) <- x
  rate.detection <- melt(rate.detection)
  colnames(rate.detection) <- c('age.group','period','rate.detection')
  
  # number of women examined:
  conapo.examined <- merge(conapo.incidence,percentage, by = 'age.group')
  conapo.examined$women.examined <- round(ifelse(conapo.examined$period %in% x, 
                                                       conapo.examined$total.group*conapo.examined$percentage,0))
  
  # number of women which are detected BrCa:
  
conapo.detected <- merge(conapo.examined,rate.detection,
                             by = intersect(names(rate.detection), names(conapo.examined)),
                             all = T)
  
conapo.detected$rate.detection[is.na(conapo.detected$rate.detection)] <- 0
  
conapo.detected$detected.women<-  round((conapo.detected$women.examined/100000)*conapo.detected$rate.detection)
  
  return(conapo.detected)
  
}

################################################
# function to determine new estimated BrCa cases based on populational data

population.incidence <- function(incidence,age.group,conapo,mean.age.group){
  
  
  # divide the data set in group ages and years that we're interested in
  conapo_2 <- subset(conapo,Age >= 25 & Age <= 74)
  conapo_3 <- conapo_2[,c(1,13:dim(conapo_2)[2])]
  
  # transform the data:
  conapo_4 <- conapo_3
  conapo_4$age.group <- cut(conapo_3$Age, age.group, right = FALSE)
  conapo_5 <- melt(conapo_4[,-1],by = colnames)
  conapo.acumulated <- ddply(conapo_5,c('age.group','variable'),summarize,total.group = sum(value))
  
  # merge incidence and population data set:
  conapo.incidence <- merge(conapo.acumulated,incidence, by = 'age.group')
  conapo.incidence$new.cases <- round((conapo.incidence$total.group/100000)*conapo.incidence$incidence) 
  
  # transform and remove variables
  conapo.incidence$period <- as.character(conapo.incidence$variable)
  conapo.incidence$period <- gsub('X','',conapo.incidence$period)
  conapo.incidence$period <- as.numeric(conapo.incidence$period)
  conapo.incidence <- conapo.incidence[,-2]
  
  # merge with the mean of each group age
  conapo.incidence <- merge(conapo.incidence,mean.age.group, by = 'age.group')
  
  return(conapo.incidence)
  
}

################################################
# cost effectiveness:
# function that quantifies the cost and benefit of each intervention

analisis.ce <- function(percentage,periodicity,incidence,conapo.incidence,mean.sojourn.time,
                        sensitivity,start,end,age.group.data.frame) {
  
  conapo.incidence <- detection(percentage,periodicity,incidence,conapo.incidence,mean.sojourn.time,
                                 sensitivity,start,end,age.group.data.frame)
  
  # Separate each case by stages:
  conapo.incidence$cases.undetected <- conapo.incidence$new.cases - conapo.incidence$detected.women
  
  conapo.incidence$Early <- round(conapo.incidence$cases.undetected*distribution[distribution$clinical.stage == 'Early',]$distribution) +
    conapo.incidence$detected.women
  
  conapo.incidence$LocallyAdvanced  <- round(conapo.incidence$cases.undetected*distribution[distribution$clinical.stage == 'Locally Advanced',]$distribution)
  
  conapo.incidence$Metastasic <- round(conapo.incidence$cases.undetected*distribution[distribution$clinical.stage == 'Metastasic',]$distribution)
  
  # cost of mammography
  
  conapo.incidence$cost.mammography <- conapo.incidence$women.examined*data$cost.mammography
  
  # cost of treatment
  conapo.incidence$cost.stage.treatment.early <- conapo.incidence$Early*costs.stage[costs.stage$clinical.stage == 'Early',]$costs
  
  conapo.incidence$cost.stage.treatment.Locally.Advanced <- conapo.incidence$LocallyAdvanced*costs.stage[costs.stage$clinical.stage == 'Locally Advanced',]$costs
  
  conapo.incidence$cost.stage.treatment.Metastasic <- conapo.incidence$Metastasic*costs.stage[costs.stage$clinical.stage == 'Metastasic',]$costs
  
  # DALYS
  conapo.incidence$DALYS <- conapo.incidence$Early*conapo.incidence$total.DALY.early+ 
    conapo.incidence$LocallyAdvanced*conapo.incidence$total.DALY.advanced +
    conapo.incidence$Metastasic*conapo.incidence$total.DALY.Metastasic
  
  return(conapo.incidence)
  
}

################################################
# function to export results

resume.ce <- function(percentage,periodicity,incidence,conapo.incidence,mean.sojourn.time,
                       sensitivity,start,end,age.group.data.frame) {
  
  
  data <- detection(percentage,periodicity,incidence,conapo.incidence,mean.sojourn.time,
                     sensitivity,start,end,age.group.data.frame)
  
  # Separate each case by stage:
  
  data$detected.women<- pmin(data$detected.women,data$new.cases)
  
  data$cases.undetected <- data$new.cases - data$detected.women
  
  data$Early <- round(data$cases.undetected*distribution[distribution$clinical.stage == 'Early',]$distribution) +
    data$detected.women
  
  data$LocallyAdvanced<- round(data$cases.undetected*distribution[distribution$clinical.stage == 'Locally Advanced',]$distribution)
  
  data$Metastasic <- round(data$cases.undetected*distribution[distribution$clinical.stage == 'Metastasic',]$distribution)
  
  # cost of mammography
  
  data$cost.total.mammography <- data$women.examined*data$cost.mammography
  
  # cost of treatment
  
  data$total.cost.early <- data$Early*data$cost.stage.Early 
  data$total.cost.locally.advanced <- data$LocallyAdvanced*data$cost.stage.locally.advanced
  data$total.cost.metastasic <- data$Metastasic*data$cost.stage.metastasic
  
  # DALYS
  
  data$DALYS <- data$Early*data$total.DALY.early+ 
    data$LocallyAdvanced*data$total.DALY.advanced +
    data$Metastasic*data$total.DALY.metastasic
  
  # Total dalys
  total.dalys <- sum(data$DALYS)
  
  # Total cost mammography
  total.cost.mammo <- sum(data$cost.total.mammography)
  
  # Total cost treatment
  total.treatment <- sum(data$total.cost.early ) + 
    sum(data$total.cost.locally.advanced ) + 
    sum(data$total.cost.metastasic)
  
  # Total cost 
  total.cost <- total.cost.mammo + total.treatment 
  
  resume <- data.frame(total.dalys = total.dalys, total.cost.masto = total.cost.mammo,
                        total.treatment = total.treatment, total.cost = total.cost)
  return(resume)
  
}











