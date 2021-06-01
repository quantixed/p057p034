#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// CSVs from Analyze Particles

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	"Spot Size...",  ParticleAnalysis()
	"Start Over", CleanSlate()
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function ParticleAnalysis()
	PreLoader()
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function PreLoader()
	// this workflow requires a condWave which has unique strings that are used to parse the csv files
	WAVE/Z/T condWave
	if(!WaveExists(condWave))
		DoAlert 0, "For this workflow to run, you need to make a textwave called condWave"
		return -1
	endif
	
	NewPath/O/Q/M="Please find folder with particle csv files" expDiskFolder1
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	String FileList1 = IndexedFile(expDiskFolder1,-1,".csv")
	if(ItemsInList(FileList1) == 0)
		DoAlert 0, "No csv files found"
		return -1
	endif
	Wave/T fileName1Wave = ListToTextWave(FileList1,";")
	MoveWave fileName1Wave, root:fileName1Wave
	
	LoadandProcessCSVs(fileName1Wave, "expDiskFolder1", condWave)
End

Function LoadAndProcessCSVs(fWaveToLoad, pathNameString, condWave)
	Wave/T fWaveToLoad
	String pathNameString
	Wave/T condWave 
	
	Variable nFiles = numpnts(fWaveToLoad)
	String thisFile, matName
	Variable nCond = numpnts(condWave)
	Variable cond
	
	Variable i, j
	
	// set up folder
	NewDataFolder/O/S root:data
 	
	for (i = 0; i < nFiles; i += 1)
		ThisFile = fWaveToLoad[i]
		for(j = 0; j < nCond; j += 1)
			if(stringmatch(ThisFile,"*"+condWave[j]+"*"))
				cond = j
				break
			else
				cond = nCond
			endif
		endfor
		matName = FindNextName("cond" + num2str(cond) + "_")
		LoadWave/A/D/J/O/Q/K=1/L={0,1,0,1,0}/P=$pathNameString/W ThisFile
		// this part was modified from mitochondrial aggregation so the circularity is calculated
		// although we don't really need it
		Wave/Z AreaW, Perim_
//		// cut below 0.2 um^2
//		Perim_[] = (AreaW[p] > 0.2) ? Perim_[p] : NaN
//		AreaW[] = (AreaW[p] > 0.2) ? AreaW[p] : NaN
		WaveTransform zapnans AreaW
		WaveTransform zapnans Perim_
		if(numpnts(AreaW) > 0)
			MatrixOP/O/FREE tempW = 4 * pi * (AreaW / (Perim_ * Perim_))
			Concatenate/O/NP=1 {AreaW,Perim_,tempW}, $matName
		endif
		
		KillAllExcept(WaveList("cond*",";",""))
	endfor
	
	// now move the data and summaries to root
	String wList, wName
	Variable nWaves
	
	for(i = 0; i <= nCond; i += 1)
		wList = WaveList("cond" + num2str(i) + "_*",";","")
		if(ItemsInList(wList) == 0)
			continue
		endif
		matName = "root:cond" + num2str(i) + "_all"
		Concatenate/O/NP=0 wList, $matName
		// now we have stored all the data in root
		// let's get the summaries
		nWaves = ItemsInList(wList)
		Make/O/N=(nWaves) $("root:medianArea" + num2str(i))
		Wave w0 = $("root:medianArea" + num2str(i))
		Make/O/N=(nWaves) $("root:spotsPerCell" + num2str(i))
		Wave w1 = $("root:spotsPerCell" + num2str(i))
		Make/O/N=(nWaves) $("root:totalSpotArea" + num2str(i))
		Wave w2 = $("root:totalSpotArea" + num2str(i))
	
		for (j = 0; j < nWaves; j += 1)
			wName = StringFromList(j, wList)
			Wave w = $wName
			MatrixOP/O/FREE tempW = col(w,0)
			w0[j] = median(tempW)
			w1[j] = DimSize(tempW,0)
			w2[j] = sum(tempW)
		endfor
	endfor
	
	SetDataFolder root:
	
	// make a box plot of median spot area per cell
	MakeTheBoxPlot("medianArea","Area (\u03BCm\S2\M)")
	// make a box plot of spot count per cell
	MakeTheBoxPlot("spotsPerCell","Spots per cell")
	// make a box plot of spot count per cell
	MakeTheBoxPlot("totalSpotArea","Total spot area")
	
	MakeTheLayouts("p",5,3, rev = 1, saveIt = 0)
End

STATIC Function MakeTheBoxPlot(rootName, leftLabel)
	String rootName, leftLabel
	String plotName = "p_" + rootName
	
	WAVE/Z/T condWave
	Variable nCond = numpnts(condWave)
	WAVE/Z/T labelWave
	if(!waveExists(labelWave))
		Duplicate/O/T condWave, labelWave
	endif
	
	KillWindow/Z $plotName
	Display/N=$plotName
	String wName = rootName + num2str(0)
	Wave w0 = $wName
	AppendBoxPlot/W=$plotName w0 vs labelWave
	
	Variable i
	
	for(i = 1; i < nCond; i += 1)
		// any mismatches won't get plotted
		Wave w1 = $(rootName + num2str(i))
		if(WaveExists(w1))
			AddWavesToBoxPlot/W=$plotName w1
		endif
	endfor
	Label/W=$plotName left leftLabel
	SetAxis/W=$plotName/A/N=1/E=1 left
	ModifyGraph/W=$plotName margin(left)=42
	ModifyBoxPlot/W=$plotName trace=$wName,whiskerMethod=4
	ModifyBoxPlot/W=$plotName trace=$wName,markers={19,-1,19},medianMarkerColor=(65535,0,0,32768)
	ModifyBoxPlot/W=$plotName trace=$wName,dataColor=(65535,0,0,32768),outlierColor=(65535,0,0,32768)
	ModifyBoxPlot/W=$plotName trace=$wName,farOutlierColor=(65535,0,0,32768),markerSizes={2,2,2}
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

STATIC Function KillAllExcept(exceptList)
	String exceptList
	String wList = WaveList("*",";","")
	wList = RemoveFromList(exceptList, wList)
	Variable nWaves = ItemsInList(wList)
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		KillWaves/Z $wName
	endfor
End

STATIC Function KillTheseWaves(wList)
	String wList
	Variable nWaves = ItemsInList(wList)
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		KillWaves/Z $wName
	endfor
End
// Colours are taken from Paul Tol SRON stylesheet
// Colours updated. Brighter palette for up to 6 colours, then palette of 12 for > 6
// Define colours
StrConstant SRON_1 = "0x4477aa;"
StrConstant SRON_2 = "0x4477aa;0xee6677;"
StrConstant SRON_3 = "0x4477aa;0xccbb44;0xee6677;"
StrConstant SRON_4 = "0x4477aa;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_5 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_6 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;0xaa3377;"
StrConstant SRON_7 = "0x332288;0x88ccee;0x44aa99;0x117733;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_8 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_9 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_10 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_11 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_12 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0xaa4466;0x882255;0xaa4499;"

/// @param hex		variable in hexadecimal
Function hexcolor_red(hex)
	Variable hex
	return byte_value(hex, 2) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_green(hex)
	Variable hex
	return byte_value(hex, 1) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_blue(hex)
	Variable hex
	return byte_value(hex, 0) * 2^8
End

/// @param data	variable in hexadecimal
/// @param byte	variable to determine R, G or B value
STATIC Function byte_value(data, byte)
	Variable data
	Variable byte
	return (data & (0xFF * (2^(8*byte)))) / (2^(8*byte))
End

/// @param	nRow	variable for number of conditions
/// @param	wName	string to name the resulting colorwave
/// @param	[alpha]	optional variable for alpha column, n.b. 16-bit integer
Function MakeColorWave(nRow, wName, [alpha])
	Variable nRow
	String wName
	Variable alpha
	
	// Pick colours from SRON palettes
	String pal
	if(nRow == 1)
		pal = SRON_1
	elseif(nRow == 2)
		pal = SRON_2
	elseif(nRow == 3)
		pal = SRON_3
	elseif(nRow == 4)
		pal = SRON_4
	elseif(nRow == 5)
		pal = SRON_5
	elseif(nRow == 6)
		pal = SRON_6
	elseif(nRow == 7)
		pal = SRON_7
	elseif(nRow == 8)
		pal = SRON_8
	elseif(nRow == 9)
		pal = SRON_9
	elseif(nRow == 10)
		pal = SRON_10
	elseif(nRow == 11)
		pal = SRON_11
	else
		pal = SRON_12
	endif
	
	Variable color
	String colorWaveFullName = "root:" + wName
	if(ParamisDefault(alpha) == 1)
		Make/O/N=(nRow,3) $colorWaveFullName
		WAVE w = $colorWaveFullName
	else
		Make/O/N=(nRow,4) $colorWaveFullName
		WAVE w = $colorWaveFullName
		w[][3] = alpha
	endif
	
	Variable i
	
	for(i = 0; i < nRow; i += 1)
		// specify colours
		color = str2num(StringFromList(mod(i, 12),pal))
		w[i][0] = hexcolor_red(color)
		w[i][1] = hexcolor_green(color)
		w[i][2] = hexcolor_blue(color)
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

STATIC Function TidyAndSave(prefix)
	String prefix
	String layoutName = "all"+prefix+"Layout"
	// go to first page
	LayoutPageAction/W=$layoutName page=(1)
	// build the key
	WAVE/Z/T labelWave = root:labelWave
	WAVE/Z colorWave = root:colorWave
	Variable cond = numpnts(labelWave)
	String boxString = ""
	
	Variable i
	
	for(i = 0; i < cond; i += 1)
		// add text colour for condition
		boxString += "\\K(" + num2str(colorWave[i][0]) + "," + num2str(colorWave[i][1]) + "," + num2str(colorWave[i][2])
		boxString += ")" + labelWave[i]
		if (i < cond - 1)
			boxString += "\r"
		endif
	endfor
	TextBox/W=$layoutName/C/N=text0/F=0/A=RB/X=5.00/Y=5.00 boxString
	String fileName = layoutName + ".pdf"
	SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
End

STATIC Function DecideOpacity(nTrace)
	Variable nTrace
	Variable alpha
	if(nTrace < 10)
		alpha = 1
	elseif(nTrace < 50)
		alpha = 0.5
	elseif(nTrace < 100)
		alpha = 0.3
	else
		alpha = 0.2
	endif
	alpha = round(65535 * alpha)
	return alpha
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

STATIC Function/S FindNextName(proposedName)
	String proposedName
	
	Variable start = 0
	String result
	
	do
		result = proposedName + num2str(start)
		start += 1
	while (CheckName(result,1) != 0)
	
	return result
End