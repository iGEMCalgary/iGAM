##Function to generate n random chromosomes of length SeqLengh
GenPop <- function(n,SeqLength){
  #n(numeric value) number of chromosomes to develop
  #SeqLength is the length of each chromosome
  NewG <- data.frame(c(1:SeqLength))
  for(i in c(1:n)){
    Genecode <- sample(c(1:23),SeqLength)
    tempg <- c()
    for(ii in Genecode){
      tempg <- c(tempg,(patynum[ii]))
    }
    NewG <- data.frame(NewG,tempg)
  }
  NewG <- NewG[,2:n+1]
  return(NewG)
}
