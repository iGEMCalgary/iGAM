Mutation <- function(crosspop,probb,SeqLength){#8 random
  #loop through creating 20 different new ones along with the proElites
  for(s in c(3:25)){
    crossnum1 <- as.numeric(sample(c(1:floor(n/2)),1)) # number that the mutationstarts
    crossnum2 <- 12 - sample(c(1:floor(n/2)),1) # number that the crossover ends
    tempgen <- c()# Running generated gene
    gen1 <- as.numeric(sample(c(1:length(crosspop)),1)) # gene 1 to select for start
    gen2 <- as.numeric(sample(c(1:length(crosspop)),1)) # gene 2 to select to select for insertion
    TempMutt <- as.numeric(sample(patynum,SeqLength,replace = TRUE))

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
