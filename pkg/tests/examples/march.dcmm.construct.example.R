# Construct a 2 hidden states DCMM for the pewee data
# with hidden and visible orders set to 1.
# The estimation procedure uses only the evolutionary algorithm.
march.dcmm.construct(y=pewee,orderHC=1,orderVC=1,M=2,popSize=2,gen=2)

# Construct a first-order, three hidden states HMM for the pewee data.
# The estimation procedure uses only the Bauw-Welch algorithm.
march.dcmm.construct(pewee,orderHC=1,orderVC=0,M=3,gen=1,popSize=1,iterBw=3,stopBw=0.001)

## Not run:

# Construct a DCMM using an order 2 hidden chain,
# a visible chain of order 3 and 3 hidden states.
# A first model is computed using both EA and Baum-Welch algorithms.
# Then the initial solution is improved through two rounds of Baum-Welch optimization.
models <- list()
models[[length(models)+1]] <- march.dcmm.construct(y=pewee,orderHC=2,
                                                   orderVC=3,M=3,popSize=2,gen=2)
models[[length(models)+1]] <- march.dcmm.construct(y=pewee,seedModel=models[[1]],
                                                   orderVC=3,iterBw=10,stopBw=0.001)
models[[length(models)+1]] <- march.dcmm.construct(y=pewee,seedModel=models[[2]],
                                                   orderVC=3,iterBw=10,stopBw=0.0001)
# Show performance indicators (ll, number of independent parameters,
# BIC and AIC) for all computed models.
r <- do.call(rbind,lapply(models,march.summary))
print(r)

# Construct a three hidden states, first-order HMM (hence OrderVC=0) for the sleep data.
# By setting gen=1 and popSize=1, the estimation procedure uses only the Baum-Welch algorithm.
HMM <- march.dcmm.construct(sleep,orderHC=1,orderVC=0,M=3,gen=1,popSize=1,iterBw=10,stopBw=0.0001)
print(HMM)
# Compute the best sequence of hidden states for each subject in the sleep data.
HS <- march.dcmm.viterbi(HMM,sleep)
# Display the hidden states for the first 10 subjects.
print(HS[1:10])

## End(Not run)
