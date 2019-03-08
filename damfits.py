#!/usr/bin/python3

import sys
import os
from os.path import basename #For printing nicely argv[0]
import os.path

import matplotlib.pyplot as plt
from astropy.io import fits

accOpts=['-h', '--help','--header', '-p',\
         '-i', '-b', '--side', '--noLog',\
         '--color']

def createExtraOptionsDict(accOpts):
    extrOptDict={}
    for e in accOpts[4:]:
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

    tmpOpt=''
    for i in range(len(myArgs)):
        e=myArgs[i]
        if e[0] == '-':
            myOptDict[e]=[]
            tmpOpt=e
            continue #Just skipping the option

        if tmpOpt != '':
            myOptDict[tmpOpt].append(i)

    return myOptDict

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

    if "-i" in myOptDict:
        #The image number
        iNumStr=argv[myOptDict["-i"][0]]
        if not iNumStr.isdigit():
            print("error: %s should be a positive integer")
        iNum=int(iNumStr)
        checkBool,totNum=checkIfValidFitsImgNum(hdu_list,iNum)
        if not checkBool:
            print("error:  1 <= %s < %d should be satisfied" %(iNum,totNum))
            return False

        if '--color' in myOptDict:
            color=argv[myOptDict['--color'][0]]

    plotHistoData(hdu_list,iNum,bining,logB,side,color)
    return True

def printHelp(argv):
    print("%s [-h|--help]\n" %(basename(argv[0])))
    print("%s file.fits #displays fits file info\n" %(basename(argv[0])))
    print("%s --header number file.fits #displays fits header info\n" %(basename(argv[0])))
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
        if e not in accOpts:
            print("error: %s is not a valid option" %(e))
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

    # for e in myOptDict:
    #     print(e, myOptDict[e])

    myFitsF=argv[-1]
    if not os.path.isfile(myFitsF):
        print("error: "+myFitsF+" is not a file")
        return 2
    # openfits("d44_snolab_Int-200_Exp-100000_4283.fits")

    try:
        hdu_list = fits.open(myFitsF)
    except:
        print("error: "+myFitsF+" is not a valid fits file.")
        return 3

    # openfits(myFitsF)
    if '-p' not in myOptDict:
        if '--header' not in myOptDict:
            hdu_list.info()
        else:
            handleHeader(hdu_list,myOptDict,argv)
    else:
        handleSubCases(hdu_list, myOptDict, argv)

if __name__ == "__main__":
   main(sys.argv)
