GenLoop <- function(n,Seqlength,Muttprob,Eval,elinum,pronum,genNum){
  PlotFit <-c()
  EndCon = FALSE
  Genera
  while(EndCon == FALSE){
    pop <- GenPop(n,SeqLength)
    Elites <- Evaluation(pop,pronum,elinum,Eval)
    Selpop <- pop[,Elites ==1 | Elites ==2] # New population selected from evaluation
    SelElites <- Elites[Elites != 0] # New elites classifier for this selected populations
    crossPop <- CrossOver(Selpop, SelElites,n)
    pop <- Mutation(crossPop,MutationProb)
    PlotFit <- c(PlotFit,Eval(pop[,1]))
    print("##############################################################")
    print(genNum)
    print("###############################################################")
    Genera <- Genera + 1
    if(Genera == genNum) EndCon <- TRUE
  }
  plot(PlotFit,type = 'l',main = "Fitness Values Over Generations",ylab = "Fitness",xlab = "Generation #", xlim = c(0,200))

}
