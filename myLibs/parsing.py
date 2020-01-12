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
         '--color', '--conv1', '-c',\
         '--sigmaB']

#A consistency dictionary
cDict={'--help':[], '--header':[],\
       '--rectangle': ['-i', '--xPlot',\
                       '--dump','--xAve','--yAve',\
                       '--pDist','--r2','--gFit',\
                       '--noLog','--noPlot',\
                       '--upperB','--sOver','--save2Pdf',\
                       '--conv1','-c','--sigmaB'],\
       '--xAve': ['-r','--rectangle','-i',\
                  '--dump','--upperB','--sOver',\
                  '--save2Pdf','--range','-c','--sigmaB'],\
       '--sOver': ['-r','--rectangle','-i',\
                   '--dump','--upperB','--xAve',\
                   # '--yAve','--save2Pdf','--range'],\
                   '--yAve','--save2Pdf','--range',\
                   '--pDist','--gFit', '--conv1',\
                   '--noPlot', '--noLog','-c','--sigmaB'],\
       '--pVal': ['-i','--save2Pdf'],\
       '--sigmaB': ['-i', '--xPlot','--rectangle',\
                       '--dump','--xAve','--yAve',\
                       '--pDist','--r2','--gFit',\
                       '--noLog','--noPlot',\
                       '--sOver','--save2Pdf',\
                       '--conv1','-c'],\
       '-p': ['-i','-d','--side','--noLog','--color']}

#For the equivalent expressions
cDict['-r']=cDict['--rectangle']
cDict['--yAve']=cDict['--xAve']
#Making all consistent with themselves and fitsFiles
for cVal in cDict:
    cDict[cVal].append(cVal)
    cDict[cVal].append('fitsFiles')

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

        if e.endswith('.fits') or e.endswith('.fits.fz'):
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
    #This function is deprecated
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
        # print("Inside the -i part handling")
        #The image number
        iNumStr=argv[myOptDict["-i"][0]]
        if not handleMinusI(iNumStr,hdu_list):
            return False
        iNum=int(iNumStr)
        if '--color' in myOptDict:
            color=argv[myOptDict['--color'][0]]

    plotHistoData(hdu_list,iNum,bining,logB,side,color)
    return True

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

def getXYAveTag(myOptDict):
    if '--xAve' in myOptDict:
        return '--xAve'
    if '--yAve' in myOptDict:
        return '--yAve'
    return None
