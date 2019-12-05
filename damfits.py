#!/usr/bin/python3
import sys
import os
from os.path import basename #For printing nicely argv[0]
import os.path
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import factorial

import matplotlib.pyplot as plt
from astropy.io import fits

#After the -p option are the extra options, just be careful with that
#-i b4 -p, it's useful outside too
accOpts=['-h','--help','--header',\
         '-r','--rectangle',\
         '--pValue','--gFit','--testing',\
         '--noPlot','--noLog',\
         '--xPlot','--xAve','--yAve',\
         '--dump','--upperB','--sOver',\
         '--save2Pdf','--range','-i','--pDist',\
         '--r2','--gFit','-p',\
         '-i', '-b','--side','--noLog',\
         '--color', '--conv1']

#A consistency dictionary
cDict={'--help':[], '--header':[],\
       '--rectangle': ['-i', '--xPlot',\
                       '--dump','--xAve','--yAve',\
                       '--pDist','--r2','--gFit',\
                       '--noLog','--noPlot',\
                       '--upperB','--sOver','--save2Pdf',\
                       '--conv1'],\
       '--xAve': ['-r','--rectangle','-i',\
                  '--dump','--upperB','--sOver',\
                  '--save2Pdf','--range'],\
       '--sOver': ['-r','--rectangle','-i',\
                   '--dump','--upperB','--xAve',\
                   # '--yAve','--save2Pdf','--range'],\
                   '--yAve','--save2Pdf','--range',\
                   '--pDist','--gFit', '--conv1',\
                   '--noPlot', '--noLog'],\
       '--pVal': ['-i','--save2Pdf'],\
       '-p': ['-i','-d','--side','--noLog','--color']}

#For the equivalent expressions
cDict['-r']=cDict['--rectangle']
cDict['--yAve']=cDict['--xAve']
#Making all consistent with themselves and fitsFiles
for cVal in cDict:
    cDict[cVal].append(cVal)
    cDict[cVal].append('fitsFiles')

def getXYAveTag(myOptDict):
    if '--xAve' in myOptDict:
        return '--xAve'
    if '--yAve' in myOptDict:
        return '--yAve'
    return None

def getTagsL(rVals,hdu_list):
    tagsL=[]
    xVals,yVals=rVals[:2],rVals[2:]
    # xLen=len(hdu_list[2].data[0])
    # yLen=len(hdu_list[2].data)
    imageDataShape=hdu_list[1].shape
    yLen,xLen=imageDataShape
    dSL=getDataSecList(hdu_list)
    dSX,dSY=dSL
    for xV in  xVals:
        if xV < 1 or xV > xLen:
            if 'xRangeError' not in tagsL:
                tagsL.append('xRangeError')
        if 1 <= xV <= dSX[0] or dSX[1] <= xV <= xLen:
            if 'xOverscan' not in tagsL:
                tagsL.append('xOverscan')
        if dSX[0] <= xV <= dSX[1]:
            if 'centerRegion' not in tagsL:
                tagsL.append('xCenterRegion')
        if xV <= xLen//2:
            if 'Left' not in tagsL:
                tagsL.append('Left')
        if xV > xLen//2:
            if 'Right' not in tagsL:
                tagsL.append('Right')

    for yV in yVals:
        if yV < 1 or yV > yLen:
            if 'yRangeError' not in tagsL:
                tagsL.append('yRangeError')
        if dSY[1] < yV :
            if 'yOverscan' not in tagsL:
                tagsL.append('yOverscan')
        if dSY[0] <= yV <= dSY[1]:
            if 'yCenterRegion' not in tagsL:
                tagsL.append('yCenterRegion')

    return tagsL

def isOk2Go(tagsL,xTraCheck=None):
    if 'xRangeError' in tagsL or 'yRangeError' in tagsL:
        return False
    if xTraCheck != None:
        if 'Left' in tagsL and 'Right' in tagsL and\
           xTraCheck == '--yAve':
            return False
        if 'yOverscan' in tagsL and xTraCheck == "--xAve":
            return False
    return True

def getDataSecList(hdu_list):
    #Note: I'm assumming that the values in hdu_list are always
    #formatted correctly
    dStr=hdu_list[0].header['DATASEC'][1:-1]
    xStr,yStr=dStr.split(',')
    # xRange=xStr.split(':') #[int(xVar) for xVar in xStr.split(':')]
    xRange=[int(xVar) for xVar in xStr.split(':')]

    yRange=[int(yVar) for yVar in yStr.split(':')]
    return [xRange,yRange]

def checkOptConsistency(myOptDict):
    """Checks if the options given from the command line can be used
together"""
    theOptL=[myOpt for myOpt in myOptDict]
    for myOptE in theOptL:
        if myOptE in cDict:
            accOptsList=cDict[myOptE]
            list2See=[x for x in theOptL if x not in accOptsList]
            if len(list2See) != 0:
                print("error: %s is inconsistent with" %(myOptE))
                for l2c in list2See:
                    print(l2c)
                return False
    return True

def createExtraOptionsDict(accOpts):
    extrOptDict={}
    for e in accOpts[accOpts.index('-p')+1:]:
        extrOptDict[e]=""

    extrOptDict['-i']="""Needs one integer as argument
 \t\tso it can select the corresponding image"""
    extrOptDict['-b']="""Needs one integer for defining the bining"""
    extrOptDict['--side']="""Needs one argument specifying the left or right
\t\tside of the image"""
    extrOptDict['--noLog']="""A simple flag for turning off the ylog scale"""
    extrOptDict['--color']="""Needs one argument defining the color"""
    return extrOptDict

extrOptDict=createExtraOptionsDict(accOpts)

def gauss(x, *p):
    A,mu,sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#### for doinf some convolution fitting ####
def gaussian(x, *p):
    A,mu,sig = p
    return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def poisson(*q):
    k,lamb=q
    return (lamb**k/factorial(k)) * np.exp(-lamb)

def simpleConv(x, A, mu, sig, lamb_d, maxK=20):
    mySumArr=np.float64(np.zeros_like(x))
    # fS="/home/frank/mess/damic/karthicFiles/fanoRed.tsv"
    for k in range(maxK):
        mySumArr+=np.float64(poisson(k,lamb_d)*gaussian(x,A,mu+k,sig))
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

def handleHeader(hdu_list,myOptDict,argv):
    if len(myOptDict['--header']) == 0:
        print("error: header needs an extention number.")
        return

    myHeaderStr=argv[myOptDict['--header'][0]]
    print("Header is %s" %(myHeaderStr))
    if not myHeaderStr.isdigit():
        print("error: header number neends to be a positive integer")
        return

    hIdx=int(myHeaderStr)
    checkBool,totNum=checkIfValidFitsIdx(hdu_list,hIdx)
    if not checkBool:
        print("error:  0 <= %s < %d should be satisfied" %(hIdx,totNum))
        return

    for e in hdu_list[hIdx].header:
        print("%s\t\t\t%s" %(e, hdu_list[hIdx].header[e]))

def getMyOptDict(myArgs):
    myOptDict={}
    myOptDict['fitsFiles']=[]

    tmpOpt=''
    for i in range(len(myArgs)):
        e=myArgs[i]
        if e[0] == '-': #only changes if new option is found
            myOptDict[e]=[]
            tmpOpt=e
            continue #Just skipping the option

        if e.endswith('.fits'):
            myOptDict['fitsFiles'].append(i)
            #Assuming that if a fits file is found then whatever
            #option was used has finished receiving it's arguments.
            tmpOpt = ''

        if tmpOpt != '':
            myOptDict[tmpOpt].append(i)

    return myOptDict

def handleMinusI(iNumStr,hdu_list):
    if not iNumStr.isdigit():
        print("error: %s should be a positive integer")
        return False
    iNum=int(iNumStr)
    checkBool,totNum=checkIfValidFitsImgNum(hdu_list,iNum)
    if not checkBool:
        print("error:  1 <= %s < %d should be satisfied" %(iNum,totNum))
        return False
    return True

def getRectangleList(myOptDict,argv):
    #Assuming the rectangle check has already been done
    if '-r' in myOptDict:
        myOptDict['--rectangle']=myOptDict['-r']
    #assuming that all the values are valid
    rectangleList=[int(argv[rVal]) for rVal in myOptDict['--rectangle'][0:4]]
    return rectangleList

def handleSubCases(hdu_list, myOptDict, argv):
    iNum=1 #Which image to plot
    bining=300 #default bining
    logB=True #log y
    side="right" #or left of the image
    color="blue" #for the histogram

    if "-b" in myOptDict:
        bining=argv[myOptDict["-b"][0]]
        if not bining.isdigit():
            print("error: bining needs a positive integer")
            return False
        bining=int(bining)
        if bining <= 0:
            print("error: bining needs a positive integer")
            return False

    if "--side" in myOptDict:
        side=argv[myOptDict["--side"][0]]
        if side == "l" or side == "left":
            side="left"
        elif side =="r" or side == "right":
            side="right"
        else:
            print("error: %s is an invalid side" %(side))
            print("valid sides are left or l and right or r")
            return False

    if "--noLog" in myOptDict:
        logB=False

    iNum=1
    if "-i" in myOptDict:
        #The image number
        iNumStr=argv[myOptDict["-i"][0]]
        if not handleMinusI(iNumStr,hdu_list):
            return False
        iNum=int(iNumStr)
        if '--color' in myOptDict:
            color=argv[myOptDict['--color'][0]]

    plotHistoData(hdu_list,iNum,bining,logB,side,color)
    return True

def printHelp(argv):
    print("%s [-h|--help]\n" %(basename(argv[0])))
    print("%s file0.fits [file1.fits ...] #displays fits file info\n" %(basename(argv[0])))
    print("%s --header number file0.fits [file1.fits ...] #displays fits header info\n" %(basename(argv[0])))
    print("%s (-r|--rectangle) xMin xMax yMin yMax [-i iNum] [--xPlot [--save2Pdf file.pdf]] file0.fits [file1.fits ...] #prints average pixel value in rectangle region (improve this...)\n" %(basename(argv[0])))
    print("%s (-r|--rectangle) xMin xMax yMin yMax [-i iNum] (--xAve|--yAve) [--upperB upperBound] [--dump | --save2Pdf file.pdf] [--sOver] file0.fits [file1.fits ...] #plots the average pixel values along axes, if dump is used then it prints the values\n" %(basename(argv[0])))
    print("%s (-r|--rectangle) xMin xMax yMin yMax [--r2 xMin2 xMax2 yMin2 yMax2] [-i iNum] --pDist [--gFit | --conv1] [--noPlot | --save2Pdf file.pdf]  [--sOver] [--upperB upperBound] file0.fits [file1.fits ...] #plots the pixel distribution values\n" %(basename(argv[0])))
    print("%s --pValue xVal yVal [-i iNum] file0.fits [file1.fits ...] #prints the pixel value\n" %(basename(argv[0])))
    # print("%s -p [extraOptions] file.fits #plots \n" %(basename(argv[0])))
    # print("extraOptions:\n")
    # for e in extrOptDict:
    #     print("\t%s:\t%s\n" %(e,extrOptDict[e]))

def getFreqCount(myArray):
    """Given a numpy array (can be n dimensional) it returns a unique
array with frequency counts"""
    myNewArr=myArray.flatten()#make it 1D array
    uniqVals,freqCounts=np.unique(myNewArr,return_counts=True)
    return uniqVals,freqCounts

def myXFill(imageStuff,side):
    halfXVal=len(imageStuff[0])//2 #// is integer division
    x=[]
    if side == 'right':
        for myRow in imageStuff:
            x+=list(myRow[:halfXVal+1]) #the right side
    else:
        for myRow in imageStuff:
            x+=list(myRow[halfXVal+1:]) #the left side

    return x

def plotHistoData(hdu_list,iNum=1,bining=300,\
                 logB=True,side="right",\
                 color="blue"):

    halfXVal=len(hdu_list[iNum].data[0])//2 #// is integer division
    imageStuff=hdu_list[iNum].data
    x=myXFill(imageStuff,side)
    n,bins,patches = plt.hist(x,bining,facecolor=color,alpha=0.5)
    if logB:
        plt.yscale('log', nonposy='clip')
    plt.show()

#deprecated
def openfits(b):
    hdu_list = fits.open(b)
    hdu_list.info()
    print(hdu_list[0])
    print(hdu_list[1])
    print(hdu_list[1].data[0])
    print(len(hdu_list[1].data))
    print(len(hdu_list[1].data[0]))

    i=3
    print(i)
    imageStuff=hdu_list[i].data
    print("half x value = ")
    halfXVal=len(hdu_list[i].data[0])//2
    print(halfXVal)
    # return
    x=[]
    for  myRow in imageStuff:
        # print("len(myRow) = %d " % (len(myRow)))
        # x+=list(myRow[halfXVal+1:]) #the left side
        x+=list(myRow[:halfXVal+1]) #the right side

    # plt.imshow(hdu_list[1], cmap='gray')
    # plt.colorbar()

    num_bins = 300
    n,bins,patches = plt.hist(x, num_bins, facecolor='blue', alpha=0.5)
    plt.yscale('log', nonposy='clip')
    plt.show()

def checkIfValidOpts(myOptDict, accOpts):
    for e in myOptDict:
        if e == 'fitsFiles':
            #just ommiting this one, it's not really an option
            continue
        if e not in accOpts:
            print("error: %s is not a valid option" %(e))
            return False
        if '--xAve' in myOptDict and '--yAve' in myOptDict:
            print("error: --xAve and --yAve cannot be used simultaneously")
            return False
    return True

def getNumOfFitsImgs(hdu_list):
    return len(hdu_list)

def checkIfValidFitsImgNum(hdu_list,iNum):
    totNum=getNumOfFitsImgs(hdu_list)
    if 0< iNum < totNum:
        return True,totNum
    return False,totNum

def checkIfValidFitsIdx(hdu_list,iNum):
    totNum=getNumOfFitsImgs(hdu_list)
    if 0<= iNum < totNum:
        return True,totNum
    return False,totNum

def checkIfIntArgs(simpleList,argv):
    """Checks if the elements in argv, with indeces in simpleList, are all
integers positive. Useful function for the rectangle option and others.

    """
    for e in simpleList:
        if not argv[e].isdigit():
            return False
        if int(argv[e]) == 0:
            return False

    return True

def oSRH(argv,hdu_list,myOptDict,fitsFileIdx,marginD=0):
    rList=getRectangleList(myOptDict,argv)
    # print("rList = ", rList)
    tagsL=getTagsL(rList,hdu_list)
    # print("tagsL = ", tagsL)
    oSRect=None
    xTraCheck=getXYAveTag(myOptDict)
    if '--sOver' in myOptDict and not isOk2Go(tagsL,xTraCheck):
        print('error: there are inconsistencies in the defined rectangle')
        sys.exit()

    dSL=getDataSecList(hdu_list)
    oSRect=getOverscanRectangle(rList, hdu_list, xTraCheck,marginD)
    xMin,xMax=oSRect[0]
    yMin,yMax=oSRect[1]
    if '--sOver' in myOptDict:
        myOptDict['--sOver']=oSRect
    iNum=1
    if '-i' in myOptDict:
        iNumStr=argv[myOptDict["-i"][0]]
        if not handleMinusI(iNumStr,hdu_list):
            return 990
        iNum=int(iNumStr)

    imageStuff=hdu_list[iNum].data
    croppedArr=imageStuff[yMin-1:yMax-1,xMin-1:xMax-1]
    return croppedArr

def rectangleHandling(argv,hdu_list,myOptDict,fitsFileIdx):
    """This function will return a number greater (or equal) than 600 if
there was an error. Else it will return a numpy array if... """
    if not checkRectangle(myOptDict,argv,fitsFileIdx):
        print("error: invalid rectangle")
        return False

    #Getting the integers from the command line
    # myRectVals=[int(argv[e]) for e in myRectList]
    myRectVals=getRectangleList(myOptDict,argv)
    xMin,xMax,yMin,yMax=myRectVals
    if not checkIfValidRect(myRectVals):
        print("error: invalid rectangle range")
        return 668
    iNum=1
    if '-i' in myOptDict:
        iNumStr=argv[myOptDict["-i"][0]]
        if not handleMinusI(iNumStr,hdu_list):
            return 990
        iNum=int(iNumStr)
    imageStuff=hdu_list[iNum].data
    uppB=np.amax(imageStuff)
    if '--upperB' in myOptDict:
        if len(myOptDict["--upperB"]) < 1:
            print("error: --upperB option needs an argument")
            return 788
        uppBStr=argv[myOptDict['--upperB'][0]]
        if not uppBStr.isdigit():
            print("error: --upperB needs a positive integer")
            return 789
        uppB=int(uppBStr)
    myOptDict['uppBInt']=uppB

    croppedArr=imageStuff[yMin-1:yMax-1,xMin-1:xMax-1]
    if '--xPlot' in myOptDict:
        xSum=croppedArr.sum(axis=0)
        xPos=np.arange(xMin,xMax)
        if '--noPlot' not in myOptDict:
            plt.xlabel("x pixel")
            plt.ylabel("Sum of y vals")
            plt.title("Cummulative vals of y vs x")
            plt.plot(xPos,xSum)
    elif '--xAve' in myOptDict or '--yAve' in myOptDict:
        return croppedArr
    elif '--pDist' in myOptDict:
        if '--r2' in myOptDict:
            myRectList2=myOptDict['--r2']
            if len(myRectList2) < 4:
                print("error: --r2 option needs 4 arguments")
                return 4004
            myRectList2=myRectList2[0:4]
            if not checkIfIntArgs(myRectList2,argv):
                print("error: rectangle arguments must be integers")
                return 4005
            myRectVals2=[int(argv[e]) for e in myRectList2]
            if not checkIfValidRect(myRectVals2):
                print("error: invalid --r2 range")
                return 4006
            xMin2,xMax2,yMin2,yMax2=myRectVals2
            croppedArr2=imageStuff[yMin2-1:yMax2-1,xMin2-1:xMax2-1]
            return [croppedArr,croppedArr2]
        else:
            return [croppedArr,[]]
    else:
        flatArr=croppedArr.flatten()
        myAverage=np.average(flatArr)
        print(myAverage)
    return 34

def getOverscanRectangle(rList,hdu_list,xTraCheck,marginD=0):
    """Gets the proper overscan rectangle region"""
    #The isOk2Go function has to be
    #previously run in order to
    #avoid inconsistencies.
    tagsL=getTagsL(rList,hdu_list)
    dSL=getDataSecList(hdu_list)

    xVRange,yVRange=dSL
    iDShape=hdu_list[1].shape
    yMax,xMax=iDShape

    # print('xTraCheck = ', xTraCheck)
    if xTraCheck == None:
        #assuming we got pDist and the operation we want is the same
        #as for --yAve
        xTraCheck = '--yAve'

    if xTraCheck == '--yAve':
        oSYRange=rList[2:4]
        if 'Right' in tagsL:
            oSXRange=[xVRange[1]+1+marginD,xMax-marginD]
        if 'Left' in tagsL:
            oSXRange=[1+marginD,xVRange[0]-1-marginD]

    if xTraCheck == '--xAve':
        oSXRange=rList[:2]
        oSYRange=[yVRange[1]+1+marginD,yMax-marginD]

    oSRange=[oSXRange,oSYRange]

    return oSRange

def getMaskArrAndSum(list4NumpyStuff, myOptDict):
    mySum=np.zeros(list4NumpyStuff[0].shape)
    myMaskArr=[]
    uppBInt=myOptDict['uppBInt']
    for i in range(len(list4NumpyStuff)):
        myMaskArr.append(list4NumpyStuff[i]<=uppBInt)
        try:
            mySum+=list4NumpyStuff[i]*myMaskArr[i]
        except:
            print("error: the pixel sizes don't appear to be the same in all files.")
            sys.exit()
    return myMaskArr,mySum

def simpleOSSub(myArray,myOSAverages,\
                myMaskArr, myAxis):
    """Returns the overscan subtracted array"""
    #Getting the average of the masks
    myMaskAve=sum(myMaskArr)/len(myMaskArr)
    if myAxis == 1: #'--yAve'
        #Doing 1d transpose
        subArr=myArray-myMaskAve*myOSAverages[:, np.newaxis]
    else:
        subArr=myArray-myMaskAve*myOSAverages

    return subArr

def getProcessedCroppedArrLists(croppedArrLists, overScanList,\
                                myOptDict,myAxis):
    #Just in case --upperB wasn't used
    myMaskArr=np.ones(croppedArrLists[0].shape)
    if '--upperB' in myOptDict:
        myMaskArr,croppedArrLists = getMaskArrAndSum(croppedArrLists,\
                                                     myOptDict)
    if '--sOver' in myOptDict:
        myAxis=1 #Doing it on the y part
        myOSAverages = getOverscanAverages(overScanList, myAxis)
        croppedArrLists = [simpleOSSub(croppedArrLists,\
                                       myOSAverages,\
                                       myMaskArr, myAxis)]
    return croppedArrLists

    
def getOverscanAverages(overScanList, myAxis):
    """Returns the average overscan values with respect to the columns or
rows"""
    oSSum=np.zeros(overScanList[0].shape)

    for oVal in overScanList:
        oSSum+=oVal

    oSCSum= oSSum.sum(axis=myAxis)
    numberOfOSEntries=len(overScanList)
    numberOfOSRows=len(overScanList[0])
    numberOfOSCols=len(overScanList[0][0])

    if myAxis==0:
        numOfOSEle=numberOfOSRows
    else:# myAxis==1
        numOfOSEle=numberOfOSCols

    #No number of entries
    #division has been done at
    #this point
    myOSAverages=np.array([float(oSCV)/(numOfOSEle)\
                           for oSCV in oSCSum])

    return myOSAverages


def getAverageList(list4NumpyStuff,myKey,myOptDict,overScanList=[]):
    myAxis=1 #'--yAve'
    if myKey=='--xAve':
        myAxis=0
    #Getting a zeros numpy array with the same shape
    mySum=np.zeros(list4NumpyStuff[0].shape)

    #Getting a summed array from the info in the other fits files.
    myMaskArr,mySum = getMaskArrAndSum(list4NumpyStuff, myOptDict)

    divEle=None
    for i in range(len(list4NumpyStuff)):
        if type(divEle).__module__ != np.__name__:
            divEle=myMaskArr[i].sum(axis=myAxis)
        else:
            divEle+=myMaskArr[i].sum(axis=myAxis)

    numberOfEntries=len(list4NumpyStuff)
    # numberOfRows=len(list4NumpyStuff[0])
    # numberOfCols=len(list4NumpyStuff[0][0])

    cumSum=mySum.sum(axis=myAxis)

    if '--sOver' in myOptDict:
        myOSAverages = getOverscanAverages(overScanList, myAxis)
        # print("mySum = ", mySum)
        subArr = simpleOSSub(mySum, myOSAverages, myMaskArr, myAxis)
        cumSubSum=subArr.sum(axis=myAxis)
        myAverage=[float(cumSS/divE) for divE,cumSS\
                   in zip(divEle,cumSubSum)]

        # myAverage=[float(cumS-(float(divE)/numberOfEntries)\
        #                  *myOSAVal)/(divE) for divE,cumS,myOSAVal\
        #            in zip(divEle,cumSum,myOSAverages)]
    else:
        myAverage=[float(cumS)/(divE)\
                   for divE,cumS in zip(divEle,cumSum)]

    return myAverage

def checkIfValidRect(myRectVals):
    """Checks if the ranges of the rectangle are valid (work in
progress)."""
    if any(t < 0 for t in myRectVals):
        return False
    xMin,xMax,yMin,yMax=myRectVals
    if xMin > xMax:
        return False
    if yMin > yMax:
        return False
    return True

def checkIfValidPixel(argv,hdu_list,myOptDict):
    """Still missing checking if pixel is inside the image"""
    if '--pValue' not in myOptDict:
        return False

    theLen=len(myOptDict['--pValue'])
    if theLen<2:
        print("error: --pValue needs 2 arguments")
        return False

    myRectList=myOptDict['--pValue'][0:2]

    if not checkIfIntArgs(myRectList,argv):
        print("error: pValue arguments must be positive integers")
        return False

    return True

def handlePValue(argv,hdu_list,myOptDict,fitsFileIdx):
    pBoolVal=checkIfValidPixel(argv,hdu_list,myOptDict)
    if not pBoolVal:
        return None

    iNum=1
    if '-i' in myOptDict:
        iNumStr=argv[myOptDict["-i"][0]]
        if not handleMinusI(iNumStr,hdu_list):
            return None
        iNum=int(iNumStr)

    myRectList=myOptDict['--pValue']
    myXP,myYP=int(argv[myRectList[0]]),\
        int(argv[myRectList[1]])

    imageStuff=hdu_list[iNum].data
    pValue=imageStuff[myYP-1,myXP-1]
    return pValue

def checkRectangle(myOptDict,argv,fitsFileIdx):
    if '--rectangle' in myOptDict:
        myRectList=myOptDict['--rectangle']
    else:
        myRectList=myOptDict['-r']
    rectNum=len(myRectList)
    #Maybe it's now unecessary (the fits check)
    if fitsFileIdx in myRectList:
        rectNum-=1
    if rectNum != 4: #are we ready to put an equal here?
        print("error: rectangle option needs exactly 4 arguments")
        return False
    myRectList=myRectList[0:4]
    if not checkIfIntArgs(myRectList,argv):
        print("error: rectangle arguments must be integers")
        return False
    return True

def checkRange(myOptDict,argv):
    if '--range' not in myOptDict:
        #This check was probably done before the function call
        return False

    rangeIdxList=myOptDict['--range']
    if len(rangeIdxList) < 2:
        print("error: --range option needs at least 2 arguments")
        return False

    if len(rangeIdxList) % 2 != 0:
        print("error: --range option needs even numbers of arguments")
        return False

    if not checkIfIntArgs(rangeIdxList,argv):
        print("error: --range option needs only positive integers")
        return False

    rangeArgL=[int(argv[idx]) for idx in rangeIdxList]

    for i in range(len(rangeArgL)//2):
        if rangeArgL[2*i] > rangeArgL[2*i+1]:
            print("error: --range iVal < fVal has to be satisfied")
            return False

    #Rectangle option has to be previously defined in order to use
    #also the range option. Along with either --xAve or --yAve
    rectangleList=getRectangleList(myOptDict,argv)
    if not ('--yAve' in myOptDict or '--xAve' in myOptDict):
        print("error: --range option needs either --yAve or --xAve")
        return False

    if '--xAve' in myOptDict:
        rMin,rMax=rectangleList[:2]
    else: #--yAve case
        rMin,rMax=rectangleList[2:]#rectangle list returns a 4 element
                                   #list so these are the last 2.

    #Now checking that the ranges are within the rMin and rMax range.
    for e in rangeArgL:
        if not rMin <= e <= rMax:
            print("error: all --range values have to be inside rectangle %d <= %d <= %d is False" %(rMin, e, rMax))
            return False
    return True

def main(argv):
    if len(argv) < 2:
        print("usage: "\
              +basename(argv[0])+\
              "[options] afile.fits# use -h for help")
        return 1
    myOptDict=getMyOptDict(argv)

    if not checkIfValidOpts(myOptDict, accOpts):
        return 4

    if '-h' in myOptDict or '--help' in myOptDict:
        printHelp(argv)
        return 5

    if '--range' in myOptDict:
        if not checkRange(myOptDict,argv):
            return 9384
        return 9990

    if '--save2Pdf' in myOptDict:
        #parsing should be done for this part, getting 1 argument etc.
        myPltFig=plt.figure()
        myPdfFile=argv[myOptDict['--save2Pdf'][0]]
        if not myPdfFile.endswith('.pdf'):
            print('error: pdfFilename has to end with ".pdf"!')
            return 999

    fitsFIdxs=myOptDict['fitsFiles']
    if len(fitsFIdxs) == 0:
        print("error: at least 1 fits file has to be provided")
        return 6

    iNum=1
    if '-i' in myOptDict:
        ccdNumLIdx=myOptDict['-i']
        if len(ccdNumLIdx) < 1:
            print("error: -i option cannot be left empty")
            return 660
        if not argv[ccdNumLIdx[0]].isdigit():
            print("error: -i option needs an integer")
            return 662
        iNum=int(argv[ccdNumLIdx[0]])
    if not checkOptConsistency(myOptDict):
        print("Check help for proper syntax")
        return 8

    list4NumpyStuff=[]#for getting the numpy stuff in case they exist.
    croppedArrLists=[]
    croppedArrLists2=[]
    overScanList=[]
    #Looping through all the fits file indeces given though the
    #command line
    for fitsIdx in fitsFIdxs:
        # myFitsF=argv[fitsFIdxs[-1]]
        myFitsF=argv[fitsIdx]
        if not os.path.isfile(myFitsF):
            print("error: "+myFitsF+" is not a file.")
            print("Maybe you are in the wrong directory.")
            return 2
        # openfits("d44_snolab_Int-200_Exp-100000_4283.fits")
        # fitsFileIdx=argv.index(argv[-1])
        fitsFileIdx=fitsFIdxs[-1]
        try:
            hdu_list = fits.open(myFitsF)
        except:
            print("error: "+myFitsF+" is not a valid fits file.")
            return 3

        #Handling the rectangle option
        if '-r' in myOptDict or '--rectangle' in myOptDict:
            #simply making main return the same thing as the function
            if '-r' in myOptDict:
                #making sure '--rectangle' exists (for simplicity)
                myOptDict['--rectangle']=myOptDict['-r']
            exitVal=rectangleHandling(argv,hdu_list,myOptDict,fitsFileIdx)
            if type(exitVal).__module__ == np.__name__:
                list4NumpyStuff.append(exitVal)
                #Assuming for now that the (x|y)Ave option was given
            elif type(exitVal) == type([]): #Horrible way
                #This --pDist was probably used
                newCroppedArr=exitVal[0]
                croppedArrLists.append(newCroppedArr)
                if '--r2' in myOptDict:
                    newSubRectArr=exitVal[1]
                    croppedArrLists2.append(newSubRectArr)
            elif exitVal >= 600:
                return 666
            # if fitsIdx == fitsFIdxs[-1]:
                #No error is shown if we are not plotting
                # if '--save2Pdf' not in myOptDict:
                #     plt.show()
                # else:
                #     print("saving to pdf...?")

            #Handling the overscan part (inneficient will improve)
            if '--sOver' not in myOptDict:
                continue
            # if  myOptDict['--sOver'] != []:
            #     continue
            marginD=0
            oSCropped=oSRH(argv,hdu_list,myOptDict,fitsFileIdx,marginD=6)
            overScanList.append(oSCropped)
            continue

        if '--pValue' in myOptDict:
            pValue=handlePValue(argv,hdu_list,myOptDict,fitsFileIdx)
            if pValue == None:
                return 1
            print(pValue)
            continue

        if '--testing' in myOptDict:
            print("Inside the testing part")
            print("Hello world!")
            print(len(hdu_list[2].data[0])) # xLen
            print(len(hdu_list[2].data)) # yLen
            imageDataShape=hdu_list[1].shape
            print("DATASEC")
            print(hdu_list[0].header['DATASEC'])
            print(hdu_list[0].header['TRIMSEC'])
            print(type(hdu_list[0].header['DATASEC']))
            print("imageSHape = ", imageDataShape)
            tagsL=getTagsL([3500,4600,5,148],hdu_list)
            print("isOk2Go = ", isOk2Go(tagsL))
            print("tagsL = ", tagsL)

            return hdu_list
            sys.exit()

        # openfits(myFitsF)
        if '-p' not in myOptDict:
            if '--header' not in myOptDict:
                hdu_list.info()
            else:
                handleHeader(hdu_list,myOptDict,argv)
                return 0
        else:
            handleSubCases(hdu_list, myOptDict, argv)

    xTraCheck=None
    #This is done outside the for (easier to manipulate here)
    if '--xPlot' in myOptDict:
        if "--save2Pdf" not in myOptDict:
            plt.show()
        else:
            plt.savefig(myPdfFile, bbox_inches='tight')
        return 0

    if '--xAve' in myOptDict or '--yAve' in myOptDict:
        myKey='--xAve'
        minIdx,maxIdx=myOptDict['--rectangle'][0:2]
        myXLabel='x pixel'
        myYLabel='average y value per pixel'
        dumpXLab='xPixel'
        dumpYLab='yAverage'
        xTraCheck='--xAve'
        if '--yAve' in myOptDict:
            myKey='--yAve'
            minIdx,maxIdx=myOptDict['--rectangle'][2:4]
            myXLabel='y pixel'
            myYLabel='average x value per pixel'
            dumpXLab='yPixel'
            dumpYLab='xAverage'
            xTraCheck='--yAve'

        minVal,maxVal=int(argv[minIdx]),int(argv[maxIdx])
        aPos=np.arange(minVal,maxVal)
        myAverageL=getAverageList(list4NumpyStuff, myKey, myOptDict,overScanList)
        if '--dump' in myOptDict:
            print("#%s\t%s" %(dumpXLab,dumpYLab))
            for xVal,yVal in zip(aPos,myAverageL):
                print("%d\t%0.3f" %(xVal,yVal))
        else:
            if '--noPlot' not in myOptDict:
                plt.xlabel(myXLabel)
                plt.ylabel(myYLabel)
                plt.plot(aPos,myAverageL, marker='^')
                if '--save2Pdf' not in myOptDict:
                    plt.show()
                else:
                    plt.savefig(myPdfFile, bbox_inches='tight')
        return 0

    if '--pDist' in myOptDict:
        myAxis=1 #--yAve case, rewrite this part later
        croppedArrLists = getProcessedCroppedArrLists(croppedArrLists,\
                                                      overScanList,\
                                                      myOptDict,\
                                                      myAxis)
        newFlatLists=[e.flatten() for e in croppedArrLists]

        myKey='--yAve'

        flatFlatArr=newFlatLists[0]
        for i in range(1,len(newFlatLists)):
            flatFlatArr=np.append(flatFlatArr,newFlatLists[i])
        #Important to make the flatFlatArr an integer array specially
        #if there were any floating point operations.
        uV,fC=getFreqCount(flatFlatArr.astype(int))
        # p0 my initial guess for the fitting coefficients (A, mu and sigma)
        if '--noPlot' not in myOptDict:
            plt.xlabel("pixel value")
            plt.ylabel("counts")

            plt.plot(uV,fC)
            if not '--noLog' in myOptDict:
                plt.yscale('log', nonposy='clip')
        if '--r2' in myOptDict:
            croppedArrLists2=getProcessedCroppedArrLists(croppedArrLists2,\
                                                      overScanList,\
                                                      myOptDict,\
                                                      myAxis)
            newFlatLists2=[e.flatten() for e in croppedArrLists2]
            flatFlatArr2=newFlatLists2[0]
            for i in range(1,len(newFlatLists2)):
                flatFlatArr2=np.append(flatFlatArr2,newFlatLists2[i])
            uV2,fC2=getFreqCount(flatFlatArr2)
            if '--noPlot' not in myOptDict:
                plt.plot(uV2,fC2)

        if '--gFit' in myOptDict:
            print("#rect\tA\tmean\tsigma")
            myMaxIdx=np.argmax(fC)#locating the index
            myMaxV=uV[myMaxIdx]#guess for the mean
            myMaxC=fC[myMaxIdx]#guess for A
            mySigma=100#find better estimate
            popt,pcov = curve_fit(gauss,uV,fC,p0=[myMaxC, myMaxV, mySigma])
            if '--noPlot' not in myOptDict:
                plt.plot(uV,gauss(uV,*popt),label='fit')
            A,mean,sigma=popt
            print("r1\t%0.2f\t%0.2f\t%0.2f" %(A,mean,sigma))
            if '--r2' in myOptDict:
                myMaxIdx2=np.argmax(fC2)
                myMaxV2=uV2[myMaxIdx2]
                myMaxC2=fC2[myMaxIdx2]
                popt2,pcov2 = curve_fit(gauss,uV2,fC2,p0=[myMaxC2, myMaxV2, 100.])
                if '--noPlot' not in myOptDict:
                    plt.plot(uV2,gauss(uV2,*popt2),label='fit2')
                A2,mean2,sigma2=popt2
                print("r2\t%0.2f\t%0.2f\t%0.2f" %(A2,mean2,sigma2))
                print("r1-r2\t%0.2f\t%0.2f\t%0.2f" %(A-A2,\
                                                     mean-mean2,\
                                                     sigma-sigma2))

        if '--conv1' in myOptDict:
            # return 10
            print("#rect\tA\tmean\tsigma\tlambda")
            myMaxIdx=np.argmax(fC)#locating the index
            myMaxV=uV[myMaxIdx]#guess for the mean
            myMaxC=fC[myMaxIdx]#guess for A
            mySigma=20#find better estimate
            lamb_d=2 #Find a better way to do this!
            lamb=4.3 #Find a better way to do this!
            popt,pcov = curve_fit(simpleConv,uV,fC,p0=[myMaxC, myMaxV, mySigma,lamb_d])
            if '--noPlot' not in myOptDict:
                plt.plot(uV,simpleConv(uV,*popt),label='fit')
            A,mean,sigma,lamb_d=popt
            print("r1\t%0.2f\t%0.2f\t%0.2f\t%0.2f" %(A,mean,sigma,lamb_d))

            # popt,pcov = curve_fit(simpleConv2,uV,fC,p0=[myMaxC, myMaxV, mySigma,lamb_d,lamb])
            # if '--noPlot' not in myOptDict:
            #     plt.plot(uV,simpleConv(uV,*popt),label='fit')
            # A,mean,sigma,lamb_d,lamb=popt
            # print("r1\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f" %(A,mean,sigma,lamb_d,lamb))

        if not '--noPlot' in myOptDict:
            if '--save2Pdf' not in myOptDict:
                plt.show()
            else:
                plt.savefig(myPdfFile, bbox_inches='tight')

        return 0

if __name__ == "__main__":
   main(sys.argv)
