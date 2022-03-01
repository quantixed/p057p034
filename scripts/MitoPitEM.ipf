#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// IMOD models converted using model2point are the input here.
// This procedure is a fork of LiposomeEM.ipf

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	Submenu	"MPDV Analysis"
		"Load IMOD Models...", /Q, IMODModelAnalysis()
		"Start Over", /Q, CleanSlate()
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function IMODModelAnalysis()
	PreLoader()
End
// we stop beween these two. User input needed.
Function TheLoader()
	LoadIMODModels()
	ProcessAllModels()
	CollectAllMeasurements()
	MakeTheLayouts("MPDV",4,2, alphaSort = 1, saveIt = 0)
	MakeTheLayouts("p",5,4, rev = 1, saveIt = 0)
End

Function TheProcessor()
	ProcessAllModels()
	CollectAllMeasurements()
	MakeTheLayouts("MPDV",5,3, alphaSort = 1, saveIt = 0)
	MakeTheLayouts("p",5,3, rev = 1, saveIt = 0)
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function PreLoader()
	NewPath/O/Q/M="Please find disk folder" expDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	
	PathInfo/S expDiskFolder
	NewPath/O/Q ExpDiskFolder, S_path
	String FileList = IndexedFile(expDiskFolder,-1,".txt")
	FileList = ReplaceString(".txt",FileList,"")
	if(ItemsInList(FileList) == 0)
		DoAlert 0, "No point (txt) files found"
		return -1
	endif
	// all files are the same condition so we will just process them, unblinded
	Wave/T fileNameFWave = ListToTextWave(fileList,";")
	MoveWave fileNameFWave, root:fileNameFWave // save a copy in root
	TheLoader()
End

Function LoadIMODModels()
	
	Wave/T/Z fileNameFWave = root:fileNameFWave
	Variable nFiles = numpnts(fileNameFWave)
	// 461 pixels is 200 nm at 25kX so pixel size is 
	Make/O/N=(nFiles)/D pixelSizeFWave = 200 / 461
	
	NewDataFolder/O/S root:data
	String ThisFile
	
	Variable i
 	
	// now load files into correct data folders
	for (i = 0; i < nFiles; i += 1)
		ThisFile = fileNameFWave[i] + ".txt"
		LoadWave/A/J/D/O/Q/K=1/V={" "," $",0,0}/L={0,0,0,1,0}/P=expDiskFolder ThisFile
		MakeObjectContourWaves(i,pixelSizeFWave[i])
	endfor
	
	SetDataFolder root:
End

Function MakeObjectContourWaves(id,pxSize)
	Variable id // this is the row number from fileNameFWave - will be used as a UID
	Variable pxSize // nm per pixel
	
	Concatenate/O/KILL wavelist("wave*",";",""), matA
	WaveStats/Q/RMD=[][0] matA
	// Scale the coordinates to real values
	matA[][2,4] *= pxSize
	Variable nObjects = V_max + 1 // objects in IMOD is 1-based
	Variable nContours, contourVar
	String wName
	
	Variable i,j
	
	for (i = 0; i < nObjects; i += 1)
		MatrixOP/O filtObj = col(matA,0)
		filtObj[] = (filtObj[p] == i) ? matA[p][1] : NaN
		WaveTransform zapnans filtObj
		FindDuplicates/RN=uniqueContours filtObj
		nContours = numpnts(uniqueContours)
		// zero-indexed list of contours in this object
		for (j = 0; j < nContours; j += 1)
			contourVar = uniqueContours[j]
			// find the rows that correspond to each contour
			Duplicate/O/FREE matA,matB
			matB[][] = (matB[p][0] == i && matB[p][1] == contourVar) ? matB[p][q] : NaN
			MatrixOp/O/FREE xW = col(matB,2)
			MatrixOp/O/FREE yW = col(matB,3)
			// no need to take the z column here
			WaveTransform zapnans xW
			WaveTransform zapnans yW
			// Now make ObjectContour waves
			wName = "vs_" + num2str(id) + "_" + num2str(i) + "_" + num2str(contourVar)
			Concatenate/O/NP=1 {xW,yW}, $wName
			// do not close contour if it's not already closed 
			Wave cW = $wName // contour wave
			if(cW[DimSize(cW,0)-1][0] == cW[0][0] && cW[DimSize(cW,0)-1][1] == cW[0][1])
				DeletePoints DimSize(cW,0) - 1, 1, cW
			endif
			// delete if it's just 0-3 pixels?
			if(DimSize(cW,0) < 4)
				KillWaves cW
			endif
		endfor
	endfor
	KillWaves/Z filtObj,UniqueContours,MatA
End

// this function goes into the datafolder and runs some code on the contours in there
Function ProcessAllModels()
	SetDataFolder root:data:	// relies on earlier load
	String plotName = "MPDV_allVsPlot"
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	plotName = "MPDV_allRotVsPlot"
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	
	DetermineInnerAndOuterMembranes(0)
	DetermineInnerAndOuterMembranes(2)
	// plot summary of intermembrane distances
	plotName = "p_intermembrane"
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	Wave distList_0,distList_2
	Make/T/O/N=2 labelW={"MPDV","Mito"}
	AppendBoxPlot/W=$plotName distList_0 vs labelW
	AddWavesToBoxPlot/W=$plotName distList_2
	Label/W=$plotName left "Distance (nm)"
	SetAxis/A/N=1/E=1/W=$plotName left
	ModifyBoxPlot/W=$plotName trace=distList_0,whiskerMethod=4,markers={19,-1,19},markerSizes={2,2,2}
	ModifyBoxPlot/W=$plotName trace=distList_0,medianMarkerColor=(34952,34952,34952,32768),dataColor=(34952,34952,34952,32768)
	ModifyBoxPlot/W=$plotName trace=distList_0,outlierColor=(34952,34952,34952,32768),farOutlierColor=(34952,34952,34952,32768)
	// stats on intermembrane distance
	StatsTTest/T=0 distList_0,distList_2
	
	// now we need to retrieve the properties of the outer membrane of MPDV
	FindVesicleCentresAndPlotOut()
	
	SetDataFolder root:
End

Function DetermineInnerAndOuterMembranes(outerObj)
	Variable outerObj
	String oList = WaveList("Vs_*_" + num2str(outerObj) + "_*",";","")
	oList = RemoveFromList(WaveList("Vs_*_r",";",""),oList)
	Variable outer = ItemsInList(oList)
	Make/O/T/N=(outer,2) $("pairList_" + num2str(outerObj))
	Wave/T w = $("pairList_" + num2str(outerObj))
	String iName, oName, iList, ocList
	Variable inner
	String expr = "vs\\w([[:digit:]]+)\\w([[:digit:]]+)\\w([[:digit:]]+)"
	String img, obj, ctr
	
	Variable i,j
	
	for(i = 0; i < outer; i += 1)
		oName = StringFromList(i,oList)
		w[i][0] = oName
		SplitString/E=expr oName, img, obj, ctr
		iList = WaveList("vs_" + img + "_" + num2str(str2num(obj) + 1) + "_*",";","")
		ocList = WaveList("vs_" + img + "_" + obj + "_*",";","")
		inner = ItemsInList(iList)
		if(inner == 0)
			w[i][1] = ""
			continue
		elseif(inner == 1 && ItemsInList(ocList) == 1)
			w[i][1] = StringFromList(0, iList)
		else
			w[i][1] = MatchContour(oName,ocList,iList)
		endif
	endfor
	
	Make/O/N=(outer) $("distList_" + num2str(outerObj))
	Wave dW = $("distList_" + num2str(outerObj))
	
	// record distances between the two
	for(i = 0; i < outer; i += 1)
		oName = w[i][0]
		iName = w[i][1]
		Wave/Z oW = $oName
		Wave/Z iW = $iName
		if(!WaveExists(iW))
			dW[i] = NaN
			continue
		else
			dW[i] = FindMedianDistance(oW,iW)
		endif
	endfor
End

STATIC Function/S MatchContour(oName,ocList,iList)
	String oName, ocList, iList
	Variable ocNum = ItemsInList(ocList)
	Variable iNum = ItemsInList(iList)
	Make/O/N=(ocNum,iNum)/FREE distMat, ocMat, iMat
	ocMat[][] = p
	iMat[][] = q
	String matchStr = ""
	
	Variable i,j
	
	for(i = 0; i < ocNum; i += 1)
		Wave ocW = $(StringFromList(i,ocList))
		MatrixOp/O/FREE ocCentroid = sumcols(ocW) / numrows(ocW)
		
		for(j = 0; j < iNum; j += 1)
			Wave iW = $(StringFromList(j,iList))
			MatrixOp/O/FREE iCentroid = sumcols(iW) / numrows(iW)
			distMat[i][j] = abs(ocCentroid[0][0] - iCentroid[0][0]) + abs(ocCentroid[0][1] - iCentroid[0][1])
		endfor
	endfor
	Redimension/N=(ocNum * iNum) distMat, ocMat, iMat
	Sort distMat, distMat, ocMat, iMat
	Variable pairs = min(ocNum,iNum)
	DeletePoints pairs, (ocNum * iNum) - pairs, distMat, ocMat, iMat
	Variable ocIndex = WhichListItem(oName,ocList)
	if(ocIndex == -1)
		matchStr = ""
		return matchStr
	endif
	FindValue/V=(ocIndex) ocMat
	if(V_Value == -1)
		matchStr = ""
	else
		matchStr = StringFromList(iMat[V_Value],iList)
	endif
		
	return matchStr
End

STATIC Function FindMedianDistance(m0,m1)
	Wave m0,m1
	Variable nRowsA = dimsize(m0,0)
	Variable nRowsB = dimsize(m1,0)
	
	MatrixOp/O/FREE ax = col(m0,0)
	MatrixOp/O/FREE ay = col(m0,1)
	MatrixOp/O/FREE bx = col(m1,0)
	MatrixOp/O/FREE by = col(m1,1)
	MatrixOp/O/FREE matAx = colRepeat(ax,nRowsB)
	MatrixOp/O/FREE matAy = colRepeat(ay,nRowsB)
	MatrixOp/O/FREE matBx = rowRepeat(bx,nRowsA)
	MatrixOp/O/FREE matBy = rowRepeat(by,nRowsA)
	MatrixOp/O/FREE distanceX = matAx - matBx
	MatrixOp/O/FREE distanceY = matAy - matBy
	MatrixOp/O/FREE matDist = sqrt(distanceX * distanceX + distanceY * distanceY)
	// matDist has m0 points as rows and distances to every point in m1 as columns
	// find the minimum for each point of m0
	MatrixOp/O/FREE minDist = minRows(matDist)
	
	return Median(minDist)
End

Function FindVesicleCentresAndPlotOut()
	// finding the centre of each vesicle and placing coords into a wave called VsWave
	WAVE/Z/T pairList_0
	Variable nWaves = DimSize(pairList_0,0)
	Variable nCol = 2 // x y 
	Make/O/N=(nWaves,nCol) Img_VsCentre
	Make/O/N=(nWaves) Img_VsMinAxis,Img_VsMajAxis,Img_VsPerimeter,Img_VsArea,Img_VsMean2R
	String wName, iName, plotName
	
	Variable i,j
	
	for (i = 0; i < nWaves; i += 1)
		wName = pairList_0[i][0]
		iName = pairList_0[i][1]
		Wave w0 = $wName
		Wave/Z iW = $(pairList_0[i][1])
		for (j = 0; j < nCol; j += 1)
			WaveStats/Q/M=1/RMD=[][j] w0
			Img_VsCentre[i][j] = V_avg
		endfor
		plotName = "MPDV_allVsPlot"
		AppendToGraph/W=$plotName w0[][1] vs w0[][0]
		// offset to origin
		ModifyGraph/W=$plotName offset($wName)={-Img_VsCentre[i][0],-Img_VsCentre[i][1]}
		ModifyGraph/W=$plotName rgb($wName)=(100*257,49*257,142*257,6554)
		// add inner membrane if present
		if(WaveExists(iW))
			AppendToGraph/W=$plotName iW[][1] vs iW[][0]
			ModifyGraph/W=$plotName offset($iName)={-Img_VsCentre[i][0],-Img_VsCentre[i][1]}
			ModifyGraph/W=$plotName rgb($iName)=(34*257,119*257,185*257,6554)
		endif
		TidyGraph(plotName)
		// find eigenvectors and rotate vesicle coords (also offset them)
		Wave w1 = FindEV(w0)
		plotName = "MPDV_allRotVsPlot"
		AppendToGraph/W=$plotName w1[][1] vs w1[][0]
		ModifyGraph/W=$plotName rgb($NameOfWave(w1))=(100*257,49*257,142*257,6554)
		// need to rotate the inner wave (if there is one)
		WAVE/Z M_C
		Wave/Z iW = $(pairList_0[i][1])
		if(WaveExists(iW))
			Wave w2 = RotatePartner(iW,M_C)
			AppendToGraph/W=$plotName w2[][1] vs w2[][0]
			ModifyGraph/W=$plotName rgb($NameOfWave(w2))=(34*257,119*257,185*257,6554)
		endif
		TidyGraph(plotName)
		
		Img_VsMinAxis[i] = VesicleAxisLength(w1,1)
		Img_VsMajAxis[i] = VesicleAxisLength(w1,0)
		Img_VsPerimeter[i] = FindLengthOfXYCoords(w1)
		MatrixOp/O/FREE w1c0 = col(w1,0)
		MatrixOp/O/FREE w1c1 = col(w1,1)
		Img_VsArea[i] = PolygonArea(w1c0,w1c1)
		// use polar coords to get mean radius, convert to diam
		Make/O/N=(dimsize(w1,0)-1)/FREE rW
		rW[] = sqrt(w1[p][0]^2 + w1[p][1]^2)
		Img_VsMean2R[i] = 2 * mean(rW)
	endfor
	if(numpnts(Img_VsArea) > 1)
		MatrixOp/O/NTHR=0 Img_VsAspectRatio = Img_VsMinAxis / Img_VsMajAxis
		MatrixOp/O/NTHR=0 Img_VsCircularity = (4 * pi * Img_VsArea) / (Img_VsPerimeter * Img_VsPerimeter)
	elseif(numpnts(Img_VsArea) == 1)
		MatrixOp/O/NTHR=0 Img_VsAspectRatio = Img_VsMinAxis / Img_VsMajAxis
		MatrixOp/O/FREE/NTHR=0 tempMat = Img_VsPerimeter * Img_VsPerimeter
		MatrixOp/O/NTHR=0 Img_VsCircularity = 4 * pi * tempMat
	else
		Print ""
	endif
End

///	@param	m1	2D wave of xy coords
Function/WAVE FindEV(m1)
	Wave m1
	MatrixOp/O/FREE xCoord = col(m1,0)
	MatrixOp/O/FREE yCoord = col(m1,1)
	
	// translate to origin
	Variable offX = mean(xCoord)
	Variable offY = mean(yCoord)
	xCoord[] -= offX
	yCoord[] -= offY
	// do PCA. Rotated points are in M_R
	PCA/ALL/SEVC/SRMT/SCMT xCoord,yCoord
	WAVE M_R
	String mName = NameOfWave(m1) + "_r"
	Duplicate/O M_R, $mName
	Wave m2 = $mName
	
	Return m2
End

///	@param	m1	2D wave of xy coords
Function/WAVE RotatePartner(m1,columnMatrix)
	Wave m1,columnMatrix
	MatrixOp/O/FREE xCoord = col(m1,0)
	MatrixOp/O/FREE yCoord = col(m1,1)
	
	// translate to origin
	Variable offX = mean(xCoord)
	Variable offY = mean(yCoord)
	xCoord[] -= offX
	yCoord[] -= offY
	// PCA was already done. We have columnMatrix
	Duplicate/O/FREE columnMatrix, tempMat
	tempMat[0][1] *= -1
	tempMat[1][0] *= -1
	Concatenate/O/NP=1/FREE {xCoord,yCoord}, tempCoord
	
	String mName = NameOfWave(m1) + "_r"
	MatrixOp/O $mName = tempCoord x tempMat
	Wave m2 = $mName
	
	Return m2
End

///	@param	m1	2D wave of xy coords
///	@param	colNo	column number to use for search
STATIC Function VesicleAxisLength(m0,colNo)
	Wave m0
	Variable colNo
	Variable nRow = DimSize(m0,0)
	Duplicate/O/FREE m0, m1
	InsertPoints nRow,1, m1
	m1[nRow - 1][] = m1[0][q]
	
	MatrixOp/O/FREE m1c0 = col(m1,0)
	MatrixOp/O/FREE m1c1 = col(m1,1)
	Variable V_Value,len
	if(colNo == 1)
		FindLevel/Q/EDGE=1/P m1c0, 0
		len = abs(m1c1(V_LevelX))
		FindLevel/Q/EDGE=2/P m1c0, 0
		len += abs(m1c1(V_LevelX))
	else
		FindLevel/Q/EDGE=1/P m1c1, 0
		len = abs(m1c0(V_LevelX))
		FindLevel/Q/EDGE=2/P m1c1, 0
		len += abs(m1c0(V_LevelX))
	endif
	return len
End

///	@param	m1	2D wave of xy coords
STATIC Function FindLengthOfXYCoords(m0)
	Wave m0
	
	Variable nRow = DimSize(m0,0)
	Duplicate/O/FREE m0, m1
	InsertPoints nRow,1, m1
	m1[nRow - 1][] = m1[0][q]
	// make new 2D wave of xy coords
	Duplicate/O/FREE m1,tempDist
	// offset to zero
	tempDist[][0] -= m1[0][0]
	tempDist[][1] -= m1[0][1]
	// Differentiate, backward difference
	Differentiate/METH=2 tempDist
	// find norm, cumlative distance
	MatrixOp/O/FREE/NTHR=0 tempNorm = sqrt(sumRows(tempDist * tempDist))
	tempNorm[0] = 0 // first point is garbage
	// return the sum of distances
	return sum(tempNorm)
End


Function TidyGraph(plotName)
	String plotName
	SetAxis/W=$plotName left -100,100
	SetAxis/W=$plotName bottom -100,100
	ModifyGraph/W=$plotName width={Aspect,1}
	ModifyGraph/W=$plotName mirror=1
End


STATIC Function CollectAllMeasurements()
	SetDataFolder root:data:	// relies on earlier load
	String wList = "Img_VsArea;Img_VsAspectRatio;Img_VsCircularity;Img_VsMajAxis;Img_VsMinAxis;Img_VsMean2R;Img_VsPerimeter;"
	Variable nWaves = ItemsInList(wList)
	String wName, plotName
	Make/O/N=(1)/T bLabelWave = "MPDV"
	
	Variable i
		
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		plotName =  "p_" + wName
		Wave w = $wName
		KillWindow/Z $plotName
		Display/N=$plotName
		AppendBoxPlot/W=$plotName w vs bLabelWave
		ModifyBoxPlot/W=$plotName trace=$wName,markers={19,-1,19},markerSizes={2,2,2}
		ModifyBoxPlot/W=$plotName trace=$wName,whiskerMethod=4
		ModifyBoxPlot/W=$plotName trace=$wName,medianMarkerColor=(34952,34952,34952,32768),dataColor=(34952,34952,34952,32768)
		ModifyBoxPlot/W=$plotName trace=$wName,outlierColor=(34952,34952,34952,32768),farOutlierColor=(34952,34952,34952,32768)
		SetAxis/A/N=1/E=1/W=$plotName left
		ModifyGraph/W=$plotName toMode=-1
		ModifyGraph/W=$plotName margin(left)=40
	endfor
	
	// Label y-axes
	Label/W=p_Img_VsArea left "Area (nm\S2\M)"
	Label/W=p_Img_VsminAxis left "Minor axis (nm)"
	Label/W=p_Img_VsmajAxis left "Major axis (nm)"
	Label/W=p_Img_VsAspectRatio left "Aspect Ratio"
	Label/W=p_Img_VsCircularity left "Circularity"
	Label/W=p_Img_VsMean2R left "Diameter (nm)"
	Label/W=p_Img_VsPerimeter left "Perimeter (nm)"
	SetDataFolder root:
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////
Function CleanSlate()
	String fullList = WinList("*", ";","WIN:71")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
	
	KillDataFolder/Z root:data:
		
	// Kill waves in root
	KillWaves/A/Z
	// Look for data folders and kill them
	DFREF dfr = GetDataFolderDFR()
	allItems = CountObjectsDFR(dfr, 4)
	for(i = 0; i < allItems; i += 1)
		name = GetIndexedObjNameDFR(dfr, 4, i)
		KillDataFolder $name		
	endfor
End

STATIC Function KillTheseWaves(wList)
	String wList
	Variable nWaves = ItemsInList(wList)
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w0 = $wName
		KillWaves/Z w0
	endfor
End


STATIC Function MakeTheLayouts(prefix,nRow,nCol,[iter, filtVar, rev, alphaSort, saveIt, orient])
	String prefix
	Variable nRow, nCol
	Variable iter	// this is if we are doing multiple iterations of the same layout
	Variable filtVar // this is the object we want to filter for
	Variable rev // optional - reverse plot order
	Variable alphaSort // optional - do alphanumeric sort
	Variable saveIt
	Variable orient //optional 1 = landscape, 0 or default is portrait
	if(ParamIsDefault(filtVar) == 0)
		String filtStr = prefix + "_*_" + num2str(filtVar) + "_*"	// this is if we want to filter for this string from the prefix
	endif
	
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	String allList = WinList(prefix+"*",";","WIN:1") // edited this line from previous version
	String modList = allList
	Variable nWindows = ItemsInList(allList)
	String plotName
	
	Variable i
	
	if(ParamIsDefault(filtVar) == 0)
		modList = "" // reinitialise
		for(i = 0; i < nWindows; i += 1)
			plotName = StringFromList(i,allList)
			if(stringmatch(plotName,filtStr) == 1)
				modList += plotName + ";"
			endif
		endfor
	endif
	
	if(ParamIsDefault(alphaSort) == 0)
		if(alphaSort == 1)
			modList = SortList(modList)
		endif
	endif
	
	nWindows = ItemsInList(modList)
	Variable PlotsPerPage = nRow * nCol
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/nCol)) + ","+num2str(nCol)+")"
	
	Variable pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		if(ParamIsDefault(rev) == 0)
			if(rev == 1)
				plotName = StringFromList(nWindows - 1 - i,modList)
			else
				plotName = StringFromList(i,modList)
			endif
		else
			plotName = StringFromList(i,modList)
		endif
		AppendLayoutObject/W=$layoutName/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			if(ParamIsDefault(orient) == 0)
				if(orient == 1)
					LayoutPageAction size(-1)=(842,595), margins(-1)=(18, 18, 18, 18)
				endif
			else
				// default is for portrait
				LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			endif
			ModifyLayout/W=$layoutName units=0
			ModifyLayout/W=$layoutName frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=$layoutName appendpage
				pgNum += 1
				LayoutPageAction/W=$layoutName page=(pgNum)
			endif
		endif
	endfor
	
	String fileName
	// if anthing is passed here we save an iteration, otherwise usual name
	if(!ParamIsDefault(iter))
		fileName = layoutName + num2str(iter) + ".pdf"
	else
		fileName = layoutName + ".pdf"
	endif
	// if anthing is passed here we save the filtered version
	if(ParamIsDefault(filtVar) == 0)
		fileName = ReplaceString(".pdf",fileName, "_" + num2str(filtVar) + ".pdf")
	endif
	if(ParamIsDefault(saveIt) == 0)
		if(saveIt == 1)
			SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
		endif
	else
		// default is to save
		SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
	endif
End


// for axis scaling
///	@param	value	this is the input value that requires rounding up
///	@param	roundto	round to the nearest...
STATIC Function RoundFunction(value,roundTo)
	Variable value, roundTo
	
	value /= roundTo
	Variable newVal = ceil(value)
	newVal *= roundTo
	return newVal
End