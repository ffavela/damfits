import numpy as np
from myLibs.miscellaneous import *


def getMaskArrAndSum(list4NumpyStuff, myOptDict):
    mySum=np.zeros(list4NumpyStuff[0].shape)
    myMaskArr=[]
    # uppBInt=myOptDict['uppBInt']
    uppIntL=myOptDict['uppIntL']

    for i,uppBInt in zip(range(len(list4NumpyStuff)),uppIntL):
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

    if '--upperB' in myOptDict or 'uppIntL' in myOptDict:
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
        # print("Inside the -i handling 2")
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
