Evaluation <- function(pop,pronum,elinum,Eval){
  popEval <- c() # vector containing the population evaluation per chromosome
  Elites <- c() # Vector indicating the elites(highest scores of the population)

  #Loop through Entire Population to get the values
  for(p in c(1:length(names(pop)))){
    popEval <- c(popEval,Eval(pop[,p]))

  }
  top12 <- c(tail(sort(popEval),elinum))
  top2 <- c(tail(sort(popEval),pronum))

  numel <- 0
  numpro <- 0
  for(ii in c(1:length((popEval)))){
    if(popEval[ii] %in% top12){
      if(numel <= elinum - 1){
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
      if(numpro <= pronum){
        Elites[ii] <- 2
      }
    }

  }



  return(Elites)
}
