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

#personal library
from myLibs.fitting import *
from myLibs.parsing import *
from myLibs.miscellaneous import *
from myLibs.processing import *


#Ignoring numpy errors (usually they com from the curve_fit function)
np.seterr(all='ignore') #Uncomment for debugging

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

    if '-c' in myOptDict:
        calConstIdx=myOptDict['-c']
        if len(calConstIdx) < 1:
            print("error: -c option cannot be left empty")
            return 661
        if not is_float(argv[calConstIdx[0]]):
            print("error: -c option needs a float argument")
            return 663
        calConst=float(argv[calConstIdx[0]])
        myOptDict['-cFloatVal']=calConst

    if '--sigmaB' in myOptDict:
        sigmaBIdx=myOptDict['--sigmaB']
        if len(sigmaBIdx) < 1:
            print("error: --sigmaB option cannot be left empty")
            return 666
        if not is_float(argv[sigmaBIdx[0]]):
            print("error: --sigmaB option needs a float argument")
            return 667
        sigmaBVal=float(argv[sigmaBIdx[0]])
        myOptDict['-sigmaBFloat']=sigmaBVal

    list4NumpyStuff=[]#for getting the numpy stuff in case they exist.
    croppedArrLists=[]
    croppedArrLists2=[]
    overScanList=[]
    #Looping through all the fits file indeces given though the
    #command line
    idxCounter=0
    for fitsIdx in fitsFIdxs:
        idxCounter+=1
        # myFitsF=argv[fitsFIdxs[-1]]
        # print("fitsIdx = ",fitsIdx)
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

    if 'uppBInt' in myOptDict:
        uppInt=myOptDict['uppBInt']
        #Just simplifying masking process when looping.
        uppIntL=[uppInt for dummy in fitsFIdxs]
        #If a mapping is done make sure it's though this list.
        myOptDict["uppIntL"]=uppIntL

    if '--sigmaB' in myOptDict:
        sigmaBVal=myOptDict["-sigmaBFloat"]
        if len(croppedArrLists) == 0:
            #Maybe should remove this part, I think it never enters...
            uppIntL=getUppIntL(sigmaBVal, list4NumpyStuff)
        else:
            #Maybe leave only this part?...
            uppIntL=getUppIntL(sigmaBVal, croppedArrLists)

        myOptDict["uppIntL"]=uppIntL

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
        # print(len(flatFlatArr.astype(int)))
        # print(flatFlatArr.astype(int))
        # print(len(croppedArrLists))
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
            lamb_d=2.1 #Find a better way to do this!
            lamb=4.3 #Find a better way to do this!
            b=4.2
            b=1 #The default value
            if '-cFloatVal' in myOptDict:
                b=myOptDict['-cFloatVal']
            #Doing the fitting but leaving b constant
            myConv = lambda x, A, mu, sig, lamb_d:\
                simpleConv(x, A, mu, sig, lamb_d, b)
            popt,pcov = curve_fit(myConv,\
                                  uV,fC,\
                                  p0=[myMaxC, myMaxV, mySigma,lamb_d])
            if '--noPlot' not in myOptDict:
                plt.plot(uV,myConv(uV,*popt),label='fit')
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
