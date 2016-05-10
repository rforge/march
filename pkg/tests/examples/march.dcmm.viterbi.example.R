# Estimation of a DCMM  
set.seed(111)
# Computation of a very simple DCMM model
model <- march.dcmm.construct(y=pewee,orderHC=2,orderVC=1,M=2,popSize=2,gen=1)
# Extraction of the best sequence of hidden states using the Viterbi algorithm.
bestSeq <- march.dcmm.viterbi(model,pewee)
print(bestSeq)
