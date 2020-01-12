import numpy as np
from scipy.special import factorial
from scipy.optimize import curve_fit

def gauss(x, *p):
    A,mu,sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#### for doing some convolution fitting ####
def gaussian(x, *p):
    A,mu,sig = p
    return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def poisson(*q):
    k,lamb=q
    return (lamb**k/factorial(k)) * np.exp(-lamb)

def simpleConv(x, A, mu, sig, lamb_d, b=1, maxK=20):
    mySumArr=np.float64(np.zeros_like(x))
    # fS="/home/frank/mess/damic/karthicFiles/fanoRed.tsv"
    for k in range(maxK):
        mySumArr+=np.float64(poisson(k,lamb_d)*gaussian(x,A,mu+b*k,sig))
        # mySumArr+=poisson(k,lamb_d)*gaussian(x,A,mu,sig+k)
        # n = np.subtract(yn, self.y_mean, out=yn, casting="unsafe")
        # mySumArr = np.add(mySumArr, poisson(k,lamb_d)*gaussian(x,A,mu,sig+k),\
        #                   out=mySumArr, casting="unsafe")
        # mySumArr += poisson(k,lamb_d)*gaussian(x,A,mu,sig+k).astype(mySumArr.dtype)
    return mySumArr

def simpleConv2(x, A, mu, sig, lamb_d, lamb, maxK=10, maxKK=10):
    mySumArr=np.float64(np.zeros_like(x))
    # fS="/home/frank/mess/damic/karthicFiles/fanoRed.tsv"
    # print(type(maxK))
    for k in range(maxK):
        for kk in range(maxKK):
            mySumArr+=np.float64(poisson(kk,lamb)*poisson(k,lamb_d)*gaussian(x,A,mu,sig+k))
    return mySumArr

#### Convolution fitting functs end ####

def getUppIntL(sigmaBVal,myImportantList):
    uppIntL=[]

    for e in myImportantList:
        uV,fC=getFreqCount(e.astype(int))
        myMaxIdx=np.argmax(fC)#locating the index
        myMaxV=uV[myMaxIdx]#guess for the mean
        myMaxC=fC[myMaxIdx]#guess for A
        mySigma=100#find better estimate
        popt,pcov = curve_fit(gauss,uV,fC,p0=[myMaxC, myMaxV, mySigma])
        A,mean,sigma=popt
        uppIntV=int(mean+sigmaBVal*sigma)
        uppIntL.append(uppIntV)

    return uppIntL


def getFreqCount(myArray):
    """Given a numpy array (can be n dimensional) it returns a unique
array with frequency counts"""
    myNewArr=myArray.flatten()#make it 1D array
    uniqVals,freqCounts=np.unique(myNewArr,return_counts=True)
    return uniqVals,freqCounts
