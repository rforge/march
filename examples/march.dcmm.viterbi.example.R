# Estimation of a DCMM  
model <- march.dcmm.construct(y=pewee,orderHC=2,orderVC=3,M=4,popSize=2,gen=2)
# Extraction of the best sequence of hidden states using the Viterbi algorithm.
bestSeq <- march.dcmm.viterbi(model,pewee)
print(bestSeq)
