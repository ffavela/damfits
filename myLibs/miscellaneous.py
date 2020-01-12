from myLibs.parsing import *

def is_float(n):
    try:
        float_n = float(n)
    except ValueError:
        return False
    else:
        return True

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

#This mightbe deprecated
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

def printHelp(argv):
    print("%s [-h|--help]\n" %(basename(argv[0])))
    print("%s file0.fits [file1.fits ...] #displays fits file info\n" %(basename(argv[0])))
    print("%s --header number file0.fits [file1.fits ...] #displays fits header info\n" %(basename(argv[0])))
    print("%s (-r|--rectangle) xMin xMax yMin yMax [-i iNum] [--xPlot [--save2Pdf file.pdf]] file0.fits [file1.fits ...] #prints average pixel value in rectangle region (improve this...)\n" %(basename(argv[0])))
    print("%s (-r|--rectangle) xMin xMax yMin yMax [-i iNum] (--xAve|--yAve) [--upperB upperBound | --sigmaB sigmaV] [--dump | --save2Pdf file.pdf] [--sOver] [-c calConst] file0.fits [file1.fits ...] #plots the average pixel values along axes, if dump is used then it prints the values\n" %(basename(argv[0])))
    print("%s (-r|--rectangle) xMin xMax yMin yMax [--r2 xMin2 xMax2 yMin2 yMax2] [-i iNum] --pDist [--gFit | --conv1] [--noPlot | --save2Pdf file.pdf]  [--sOver] [--upperB upperBound | --sigmaB sigmaV] [-c calConst] file0.fits [file1.fits ...] #plots the pixel distribution values\n" %(basename(argv[0])))
    print("%s --pValue xVal yVal [-i iNum] file0.fits [file1.fits ...] #prints the pixel value\n" %(basename(argv[0])))
    # print("%s -p [extraOptions] file.fits #plots \n" %(basename(argv[0])))
    # print("extraOptions:\n")
    # for e in extrOptDict:
    #     print("\t%s:\t%s\n" %(e,extrOptDict[e]))
