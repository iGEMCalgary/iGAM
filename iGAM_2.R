######Genetic Algorithm for selection#########

set.seed(16)
install.packages("Peptides")
install.packages("stringi")
library(Peptides)
library(stringi)

##########################################
##   ~EDIT HERE~     
#Set gix to be the string of your proteins Amino acids in FASTA Format

gix <- "MLKLVVLISVYVCYVQARDYYLNPCLPNPCRYGGTCKSIGLFGYKCFCTNGYKGKNCQFNACTPNPCLNGGTCALIYGPPFYQCSCPYGYYGTKCEFKRHYYDRCGGCLNGGLCISDSYGKYVCRCKPGYYGKRCIDPYY"

###############
###Fitness Function
Eval <- function(f){ #f is a string of the FASTA we are looking to Change
  #Make it usable in a FASTA
  #Here is where you would put in any modifications present into your original sequence
  fu <- gix
  #GEt the Criteria Values

  ######################################

  #INSERT FITNESS FUNCTION HERE
  #MAKE IT RETURN UNDER THE VARIABLE Gcrit

  #Example test one
  Gcrit <- boman(fu)
  #####################################

  print(Gcrit)
  return(Gcrit)
}


################################
##     ~Everything below should be ok for changing 10 amino acids~

## If its not working and you would like a hand please contact andrew.symes1@ucalgary.ca

patynum <- c(1:20)
paty <- (c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","U","V","W"))






### population function generating 25 chromosomes

GenPop <- function(n){
  NewG <- data.frame(c(1:10))
  for(i in c(1:n)){
    Genecode <- sample(c(1:20),10)
    tempg <- c()
    for(ii in Genecode){
      tempg <- c(tempg,(patynum[ii]))
    }
    NewG <- data.frame(NewG,tempg)
  }
  NewG <- NewG[,2:n+1]
  return(NewG)
}

### Function for genetic Crossover

CrossOver <- function(Selpop,SelElites){ #selpop is a list of the selected population

  #loop through creating 20 different new ones along with the proElites
  crosspop <- data.frame(Selpop[,SelElites == 2]) # Population to go mutation step
  for(s in c(1:23)){
    crossnum1 <- as.numeric(sample(1:5,1)) # number that the crossover starts
    crossnum2 <- 10 - sample(1:5,1) # number that the crossover ends
    tempgen <- c()# Running generated gene
    gen1 <- sample(c(1:10),1) # gene 1 to select for start
    gen2 <- sample(c(1:10),1) # gene 2 to select to select for insertion

    tempgen <- Selpop[,gen1] #Select starting base gene
    tempgen <- as.numeric(replace(tempgen,c(crossnum1:crossnum2),Selpop[gen2,crossnum1:crossnum2])) # replace the part of the chromosome

    crosspop <- data.frame(crosspop,tempgen)
  }
  return(crosspop)
}

####Function for Evaluating the entire population

Evaluation <- function(pop){
  popEval <- c() # vector containing the population evaluation per chromosome
  Elites <- c() # Vector indicating the elites(highest scores of the population)

  #Loop through Entire Population to get the values
  for(p in c(1:length(names(pop)))){
    popEval <- c(popEval,Eval(pop[,p]))

  }
  top12 <- c(tail(sort(popEval),12))
  top2 <- c(tail(sort(popEval),2))

  numel <- 0
  numpro <- 0
  for(ii in c(1:length((popEval)))){
    if(popEval[ii] %in% top12){
      if(numel <= 11){
        Elites <- c(Elites,1)
        numel <- numel + 1
      }else{
        Elites <- c(Elites,0)
      }
    } else{
      Elites <- c(Elites,0)
    }

  }

  #Select Elites and the even better elites
  for(ii in c(1:length(popEval))){
    if(popEval[ii] %in% top2) {
      numpro <- numpro + 1
      if(numpro <= 2){
        Elites[ii] <- 2
          #replace(Elites,ii,2)
      }
      }

  }



  return(Elites)
}

###### Function for performing Mutation on the new cross population
Mutation <- function(crosspop,probb){#8 random
  #loop through creating 20 different new ones along with the proElites
  for(s in c(3:25)){
    crossnum1 <- as.numeric(sample(c(1:6),1)) # number that the crossover starts
    crossnum2 <- 10 - sample(c(1:6),1) # number that the crossover ends
    tempgen <- c()# Running generated gene
    gen1 <- as.numeric(sample(c(1:17),1)) # gene 1 to select for start
    gen2 <- as.numeric(sample(c(1:17),1)) # gene 2 to select to select for insertion
    TempMutt <- as.numeric(sample(patynum,10,replace = TRUE))

    tempgen <- crosspop[,gen1] #Select starting base gene
    tempgen <- (replace(tempgen,c(crossnum1:crossnum2),TempMutt[crossnum1:crossnum2])) # replace the part of the chromosome

    crosspop[s] <- (tempgen)
  }
  #introduce 8 random mutts each round
  #mutts <- GenPop(8)
  Mutpop <- data.frame(crosspop)#,mutts)

  #Assume 5 mutations per mutation success

  return(Mutpop)
}

########### Genetic Algorithm ################

#################Set up initial parameter #################
g <- 25 # size of the population
Newpop <- c()#New population generated through GenPop
EndCon <- FALSE #End Condition for the entire algorithm
CrossoverProb <- 0.2 #Probability for Crossover
MutationProb <- 0.8 # Probability for random mutation in the population
SelectionProb <- 12/25
Elitism <- 0.2 #Percent of top candidates that are passed on
genNum <- 0 #Number of generations that have passed
MaxGen <- 100 #Maximum number of Generations allowed
Selpop <- c() #Population selected on if they are elites
mutElites <- c(rep(1,5),rep(0,20)) # elites for the mutation population change
PlotFit <- c()
######################Initialization #############
pop <- (GenPop(g))

###Evaluation
Elites <- Evaluation(pop)

####While we are still going AKA while EndCon == FALSE
genNum = 0
while(EndCon == FALSE){

  ###Selection and Reproduction
  Selpop <- pop[,Elites ==1 | Elites ==2] # New population selected from evaluation
  SelElites <- Elites[Elites != 0] # New elites classifier for this selected populations

  ###Crossover
  crossPop <- CrossOver(Selpop, SelElites)

  ###Mutation
  pop <- Mutation(crossPop,MutationProb)
  PlotFit <- c(PlotFit,Eval(pop[,1]))
  #pop <- crossPop
  ###Evaluation
  Elites <- Evaluation(pop)
  ###End conditioning for the loop
  genNum <- genNum + 1
  if(genNum == MaxGen) EndCon <- TRUE
  print("##############################################################")
  print(genNum)
  print("###############################################################")
}
###Ends after breakpoints reached

plot(PlotFit,type = 'l',main = "Fitness Values Over Generations",ylab = "Fitness Value",xlab = "Generation #", xlim = c(0,500),lwd = 4)
legend(35,7,legend = c('Genetic Algorithm Score','Max Fitness Value'),col = c('black', 'red'),lty= 1:2)
#abline(v=101,lwd=2,col='red',lty =2)

