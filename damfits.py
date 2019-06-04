#!/usr/bin/python3
import sys
import os
from os.path import basename #For printing nicely argv[0]
import os.path
import numpy as np

import matplotlib.pyplot as plt
from astropy.io import fits

#After the -p option are the extra options, just be careful with that
#-i b4 -p, it's useful outside too
accOpts=['-h', '--help','--header',\
         '-r', '--rectangle',\
         '--xPlot','--xAve','--yAve',\
         '--dump','-i', '-p',\
         '-i', '-b', '--side', '--noLog',\
         '--color']

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

def handleHeader(hdu_list,myOptDict,argv):
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
            #leaving the tmpOpt conditional after this one for now

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
    print("%s (-r|--rectangle) xMin xMax yMin yMax [-i iNum] [--xPlot] file0.fits [file1.fits ...] #prints average pixel value in rectangle region (improve this...)\n" %(basename(argv[0])))
    print("%s (-r|--rectangle) xMin xMax yMin yMax [-i iNum] (--xAve|--yAve) [--dump]file0.fits [file1.fits ...] #plots the averages along axes, if dump is used then it prints the values\n" %(basename(argv[0])))
    print("%s -p [extraOptions] file.fits #plots \n" %(basename(argv[0])))
    print("extraOptions:\n")
    for e in extrOptDict:
        print("\t%s:\t%s\n" %(e,extrOptDict[e]))


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
integers. Useful function for the rectangle option.

    """
    for e in simpleList:
        if not argv[e].isdigit():
            return False
    return True

def rectangleHandling(argv,hdu_list,myOptDict,fitsFileIdx):
    if '--rectangle' in myOptDict:
        myRectList=myOptDict['--rectangle']
    else:
        myRectList=myOptDict['-r']
    rectNum=len(myRectList)
    if fitsFileIdx in myRectList:
        rectNum-=1
    if rectNum != 4:
        print("error: rectangle option needs exactly 4 arguments")
        return 666
    myRectList=myRectList[0:4]
    if not checkIfIntArgs(myRectList,argv):
        print("error: rectangle arguments must be integers")
        return 667
    #Getting the integers from the command line
    myRectVals=[int(argv[e]) for e in myRectList]
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
    croppedArr=imageStuff[yMin:yMax,xMin:xMax]
    if '--xPlot' in myOptDict:
        xSum=croppedArr.sum(axis=0)
        xPos=np.arange(xMin,xMax)
        plt.xlabel("x pixel")
        plt.ylabel("Sum of y vals")
        plt.title("Cummulative vals of y vs x")
        plt.plot(xPos,xSum)
    elif '--xAve' in myOptDict or '--yAve' in myOptDict:
        return croppedArr
    else:
        flatArr=croppedArr.flatten()
        myAverage=np.average(flatArr)
        print(myAverage)

    return 34

def getAverageList(list4NumpyStuff,myKey,myOptDict):
    #Getting a zeros numpy array with the same shape
    mySum=np.zeros(list4NumpyStuff[0].shape)
    #Getting a sumed array from the info in the other fits files.
    for i in range(len(list4NumpyStuff)):
        mySum+=list4NumpyStuff[i]
    myAxis=1 #'--yAve'
    if myKey=='--xAve':
        myAxis=0
    numberOfEntries=len(list4NumpyStuff)
    numberOfRows=len(list4NumpyStuff[0])
    numberOfCols=len(list4NumpyStuff[0][0])
    numOfEle=numberOfRows
    if myAxis==1:
        numOfEle=numberOfCols
    cumSum=mySum.sum(axis=myAxis)

    myAverage=[float(cumS)/(numberOfEntries*numOfEle)\
               for cumS in cumSum]
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

    fitsFIdxs=myOptDict['fitsFiles']
    if len(fitsFIdxs) == 0:
        print("error: at least 1 fits file has to be provided")
        return 6

    list4NumpyStuff=[]#for getting the numpy stuff in case they exist.

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
            elif exitVal >= 600:
                return 666
            if fitsIdx == fitsFIdxs[-1]:
                #No error is shown if we are not plotting
                plt.show()
            continue

        # openfits(myFitsF)
        if '-p' not in myOptDict:
            if '--header' not in myOptDict:
                hdu_list.info()
            else:
                handleHeader(hdu_list,myOptDict,argv)
        else:
            handleSubCases(hdu_list, myOptDict, argv)

    #This is done outside the for (easier to manipulate here)
    if '--xAve' in myOptDict or '--yAve' in myOptDict:
        myKey='--xAve'
        minIdx,maxIdx=myOptDict['--rectangle'][0:2]
        myXLabel='x pixel'
        myYLabel='average value over y'
        dumpXLab='xPixel'
        dumpYLab='yAverage'
        if '--yAve' in myOptDict:
            myKey='--yAve'
            minIdx,maxIdx=myOptDict['--rectangle'][2:4]
            myXLabel='y pixel'
            myYLabel='average value over x'
            dumpXLab='yPixel'
            dumpYLab='xAverage'

        minVal,maxVal=int(argv[minIdx]),int(argv[maxIdx])
        aPos=np.arange(minVal,maxVal)
        myAverageL=getAverageList(list4NumpyStuff, myKey, myOptDict)
        if '--dump' in myOptDict:
            print("#%s\t%s" %(dumpXLab,dumpYLab))
            for xVal,yVal in zip(aPos,myAverageL):
                print("%d\t%0.3f" %(xVal,yVal))
        else:
            plt.xlabel(myXLabel)
            plt.ylabel(myYLabel)
            plt.plot(aPos,myAverageL)
            plt.show()

if __name__ == "__main__":
   main(sys.argv)
