import numpy as np
import scipy.special as special
import scipy.stats as stats
def computeBinomialPrior(alphas,P,N):
    spkRng = np.asmatrix(np.arange(N+1)).H
    if(np.isscalar(P)):
        uwts = np.multiply( np.power(P,spkRng),np.power(1-P,(N-spkRng)) )
        wts = stats.binom.pmf(spkRng,N,P)
        wts2 = np.multiply(wts,uwts)
    else:
        assert len(P)==N+1
        P = P/np.sum(P)
        wts = np.asmatrix(P).H
        uwts = np.log(wts) - special.gammaln(N+1) + special.gammaln(np.arange(N+1)+1) + special.gammaln(N-np.arange(N+1)+1)
        wts2 = np.exp(np.log(wts)+uwts)
        uwts2 = np.exp(uwts)
    uwts = np.asmatrix(uwts).T
    wts = np.asmatrix(wts).T
    wts2 = np.asmatrix(wts2).T
    Z = special.polygamma(3,alphas+1)
    prior = np.zeros(np.shape(alphas))
    for i in range(len(alphas)):
        prior[i] = Z[i] - np.asmatrix(wts2).H*special.polygamma(3,alphas[i]*uwts+1)
    return prior, wts, uwts

P = 0.5
N=10
alphas =  np.random.rand(N+1,N+1)
computeBinomialPrior(alphas,P,N)