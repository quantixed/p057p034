#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Waves Average>

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	Submenu	"MitoPits"
		"Load MitoPits...",  MitoPitAnalysis()
		"Start Over", CleanSlate()
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function MitoPitAnalysis()
	if(LoadCSVFiles() != -1)
		if(GraphOutAllWaves() == -1)
			return -1
		endif
		FindZeroDist()
		GraphOutAllAlignedWaves()
		MakeTheAveragesAndPlot()
		FitThePits()
		MakeTheLayouts("p_", 5, 3, rev = 1)
	endif
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function LoadCSVFiles()
	NewPath/O/Q/M="Please find disk folder" expDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		Return -1
	endif
	
	PathInfo/S expDiskFolder
	Make/O/N=1/T pathWave = {S_path}	
	
	String FileList, ThisFile
	Variable FileLoop
	
	FileList = IndexedFile(expDiskFolder,-1,".csv")
	Variable nFiles = ItemsInList(FileList)
	Make/O/N=(nFiles)/T fileNameWave
	MakeColorWave()
	
	NewDataFolder/O/S root:data
	
	for (FileLoop = 0; FileLoop < nFiles; FileLoop += 1)
		ThisFile = StringFromList(FileLoop, FileList)
		fileNameWave[fileLoop] = ReplaceString(".csv",ThisFile,"")
		LoadWave/A=lineData/D/J/K=1/M/L={0,1,0,1,0}/P=expDiskFolder/Q ThisFile
	endfor
	
	SetDataFolder root:
End

Function GraphOutAllWaves()
	WAVE/Z colorWave = root:colorWave
	SetDataFolder root:data
	String wList = WaveList("lineData*",";","")
	Variable nWaves = ItemsInList(wList)
	
	// how many channels? Channels is cols - 1
	Wave w = $(StringFromList(0,wList))
	Variable channels = DimSize(w,1) - 1
	
	// make the graph windows
	String plotList
	if(channels == 2)
		plotList = "p_linesGreen;p_linesRed;"
	elseif(channels == 3)
		plotList = "p_linesGreen;p_linesRed;p_linesFarRed;"
	else
		DoAlert 0, "Number of channels is not 2 or 3"
		return -1
	endif
	
	String plotName
	Variable i
	
	for(i = 0; i < ItemsInList(plotList); i += 1)
		plotName = StringFromList(i, plotList)
		KillWindow/Z $plotName
		Display/N=$plotName
	endfor
	
	String wName
	Variable j
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w = $wName
		
		for(j = 0; j < ItemsInList(plotList); j += 1)
			plotName = StringFromList(j, plotList)
			AppendToGraph/W=$plotName w[][j+1] vs w[][0]
		endfor
	endfor
	
	Variable alpha = DecideOpacity(nWaves)
	
	for(i = 0; i < ItemsInList(plotList); i += 1)
		plotName = StringFromList(i, plotList)
		ModifyGraph/W=$plotName rgb=(colorWave[i][0],colorWave[i][1],colorWave[i][2],alpha)
		Label/W=$plotName left "Intensity";Label/W=$plotName bottom "Distance (μm)";SetAxis/A/E=1/N=1/W=$plotName left
		TextBox/W=$plotName/C/N=text0/F=0/S=3/B=1/X=0.00/Y=0.00 "Raw data"
	endfor
	
	SetDataFolder root:
End

Function FindZeroDist()
	SetDataFolder root:data
	String wList = WaveList("lineData*",";","")
	// how many channels? Channels is cols - 1
	Wave w = $(StringFromList(0,wList))
	Variable channels = DimSize(w,1) - 1
	
	Variable nWaves = ItemsInList(wList)
	String wName
	Variable subTime
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w = $wName
		// make free waves of distance and green fluorescence for fitting
		Duplicate/O/RMD=[][0] w, xW
		Duplicate/O/RMD=[][1] w, yW
		// make 1D
		Redimension/N=-1 xW,yW
		// find the distance for Gauss peak (high resolution)
		subTime = TheFitter(yW,xW)
		if(numType(subTime) != 0 || subTime < 0)
			// if fit has failed, just use max value of yW
			WaveStats/Q yW
			subTime = xW[V_maxRowLoc]
		endif
		xW -= subTime
		// we need to make 1D waves for averaging
		SplitAndScaleWaves(w,xW[0],xW[1],channels)
		KillWaves/Z xW, yW
	endfor
	
	SetDataFolder root:
	KillVariables/A/Z
End

///	@param	w0	matrix wave
///	@param	aa	variable, first point of x wave
///	@param	bb	variable, second point of x wave - allows for start/delta calc
///	@param	cc	variable, how many channels
STATIC Function SplitAndScaleWaves(w0,aa,bb,cc)
	Wave w0
	Variable aa,bb,cc
	
	String wName = "s_" + NameOfWave(w0) + "_"
	String newName
	
	Variable i
	
	for(i = 1; i <= cc; i += 1) // loop to do cols 1 and 2 or 1, 2 and 3
		newName = wName + num2str(i)
		Duplicate/O/RMD=[][i] w0, $newName
		Wave w = $newName
		Redimension/N=-1 w
		SetScale/P x aa,(bb - aa),"", w
	endfor
End

Function GraphOutAllAlignedWaves()
	WAVE/Z colorWave = root:colorWave
	SetDataFolder root:data
	
	String wList = WaveList("lineData*",";","")
	// how many channels? Channels is cols - 1
	Wave w = $(StringFromList(0,wList))
	Variable channels = DimSize(w,1) - 1
	
	wList = WaveList("s_lineData*",";","")
	Variable nWaves = ItemsInList(wList)
	
	// make the graph windows
	String plotList
	if(channels == 2)
		plotList = "p_s_linesGreen;p_s_linesRed;"
	elseif(channels == 3)
		plotList = "p_s_linesGreen;p_s_linesRed;p_s_linesFarRed;"
	endif

	String plotName
	Variable i
	for(i = 0; i < channels; i += 1)
		plotName = StringFromList(i, plotList)
		KillWindow/Z $plotName
		Display/N=$plotName
	endfor
	
	String wName
	Variable col
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		// find if it is from col 1, 2 or 3
		col = str2num(wName[strlen(wName)-1])
		plotName = StringFromList(col - 1, plotList)
		AppendToGraph/W=$plotName $wName
	endfor
	
	Variable alpha = DecideOpacity(ceil(nWaves / 3))
	
	for(i = 0; i < channels; i += 1)
		plotName = StringFromList(i, plotList)
		ModifyGraph/W=$plotName rgb=(colorWave[i][0],colorWave[i][1],colorWave[i][2],alpha)
		Label/W=$plotName left "Intensity";Label/W=$plotName bottom "Distance (μm)";SetAxis/A/E=1/N=1/W=$plotName left
		TextBox/W=$plotName/C/N=text0/F=0/S=3/B=1/X=0.00/Y=0.00 "Ch1-aligned"
	endfor

	SetDataFolder root:
End

Function MakeTheAveragesAndPlot()
	WAVE/Z colorWave = root:colorWave
	SetDataFolder root:data:
	
	String aveList, aveName, errName
	String plotList = WinList("p_s_lines*",";","")
	if(ItemsInList(plotList) == 2)
		plotList = "p_s_linesGreen;p_s_linesRed;"
	else
		plotList = "p_s_linesGreen;p_s_linesRed;p_s_linesFarRed;"
	endif
	String plotName
	KillWindow/Z p_s_linesAll
	Display/N=p_s_linesAll
	
	Variable i
	
	for(i = 0; i < ItemsInList(plotList); i += 1)
		plotName = StringFromList(i, plotList)
		aveList = WaveList("*",";","WIN:"+plotName)
		aveName = ReplaceString("p_s_lines",plotName,"Mean_")
		errName = ReplaceString("p_s_lines",plotName,"SD_")
		fWaveAverage(aveList, "", 1, 1, aveName, errName)
		AppendToGraph/W=p_s_linesAll $aveName
		ErrorBars/W=p_s_linesAll $aveName SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($errName,$errName)
		ModifyGraph/W=p_s_linesAll rgb($aveName)=(colorWave[i][0],colorWave[i][1],colorWave[i][2])
	endfor
	ModifyGraph/W=p_s_linesAll lsize=2
	Label/W=p_s_linesAll left "Intensity";Label/W=p_s_linesAll bottom "Distance (μm)";SetAxis/A/N=1/W=p_s_linesAll left
	
	SetDataFolder root:
End

Function FitThePits()
	SetDataFolder root:data:
	
	String wList = WaveList("lineData*",";","")
	// how many channels? Channels is cols - 1
	Wave w = $(StringFromList(0,wList))
	Variable channels = DimSize(w,1) - 1
	
	wList = WaveList("s_line*",";","")
	wList = SortList(wList)
	Variable nPits = ItemsInList(wList)
	Make/O/N=(nPits)/T allFitNames = StringFromList(p, wList)
	// make a wave to hold the fit data
	// rows are pits of all channels, columns coefficients
	Make/O/N=(nPits,4) root:allFitData
	Wave allFitData = root:allFitData
	String wName
	
	Variable i
	
	for(i = 0; i < nPits; i += 1)
		wName = StringFromList(i, wList)
		Wave w = $wName
		Wave resultW = TheFitterForScaledData(w)
		allFitData[i][] = resultW[q]
	endfor
	
	// now we check the quality of the fits
	// this function checks for outliers in the fit
	OutlierDetector(channels)
	SetDataFolder root:
	// we will assume that wavelist takes s_line* in order 0_1,0_2,0_3,1_1 and so on
	// plot out the width and height of channels as box plots
	PlotWidthsAndAmpsOfChannels(channels)
	PlotIntensities()
	// now plot out the distances from channel 1 peak to each of the peaks in 2 and 3
	PlotPeakDistances(channels)
End

///	@param	cc	variable, how many channels
STATIC Function OutlierDetector(cc)
	Variable cc
	// this one uses 1.5 * IQR
		
	Wave w = root:allFitData
	
	Variable nRows = DimSize(w,0)
	Make/O/N=(nRows,4) root:qualW = NaN // 4 cols because 4 Gauss co-effs
	Wave qualW = root:qualW
	
	Variable i
	
	for(i = 0; i < 4; i += 1)
		Duplicate/O/RMD=[][i]/FREE w, zW
		Wave resultW = IQRCrunch(zW)
		qualW[][i] = resultW[p]
	endfor
	
	Wave qualW = root:qualW
	MatrixOP/O/FREE qSumW = sumrows(qualW)
	// if more than 1 co-efficient was an outlier we won't use the doublet/triplet
	qSumW[] = (qSumW[p] > 1) ? 0 : 1 // 0 is fail, 1 is OK. Doublet should = 2, Triplet should = 3
	Redimension/N=(cc, nRows / cc) qSumW
	MatrixOp/O tripletOKWave = sumcols(qSumW)^t
	tripletOKWave[] = (tripletOKWave[p] == cc) ? 1 : 0 // 1 is OK, 0 is fail
End

STATIC Function/WAVE IQRCrunch(w)
	Wave w
	Duplicate/O/FREE w, fW
	StatsQuantiles/Q/ALL w
	WAVE/Z W_StatsQuantiles
	fW[] = (w[p] < W_StatsQuantiles[%lowerInnerFence] || w[p] > W_StatsQuantiles[%upperInnerFence]) ? 1 : 0
	
	return fW
End

///	@param	cc	variable, how many channels
Function PlotWidthsAndAmpsOfChannels(cc)
	Variable cc
	
	Wave allFitData = root:allFitData
	Wave colorWave = root:colorWave
	
	Variable nPits = DimSize(allFitData,0) / cc
	Make/O/N=(cc)/T root:labelAllChWave = "Ch " + num2str(p + 1)
	Wave/T labelAllChWave = root:labelAllChWave
	String wName, plotName
	
	Variable i
	
	// get width data for each channel
	for(i = 0; i < cc; i += 1)
		wName = "ch" + num2str(i+1) + "Widths"
		Make/O/N=(nPits) $wName
		Wave w = $wName
		w[] = allFitData[i + (p * cc)][3] // the 4th column contains width
		w[] *= 1000 // convert to nm
		// print out some stats
		Print "Channel "+ num2str(i+1) + " width: median", StatsMedian(w)
		WaveStats/Q w
		Print "    Mean", V_avg, "standard deviation", V_sdev
	endfor
	
	// plot width data
	plotName = "p_widths"
	KillWindow/Z $plotName
	// safe to assume that we have these waves
	WAVE/Z ch1Widths, ch2Widths
	if(cc > 2)
		WAVE/Z ch3Widths
	endif
	Display/N=$plotName; AppendBoxPlot ch1Widths vs labelAllChWave
	AddWavesToBoxPlot/W=$plotName ch2Widths
	if(cc > 2)
		AddWavesToBoxPlot/W=$plotName ch3Widths
	endif
	ModifyBoxPlot/W=$plotName trace=ch1Widths,whiskerMethod=4
	for(i = 0; i < cc; i += 1)
		ModifyBoxPlot/W=$plotName trace=ch1Widths,markers[i]={19,-1,19},markerSizes[i]={2,2,2},medianMarkerColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch1Widths,dataColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767),outlierColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch1Widths,farOutlierColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767)
	endfor
	Label/W=$plotName left "Width (nm)"
	SetAxis/A/N=1/E=1/W=$plotName left
	
	// now get height data for each channel
	for(i = 0; i < cc; i += 1)
		wName = "ch" + num2str(i+1) + "Intensities"
		Make/O/N=(nPits) $wName
		Wave w = $wName
		w[] = allFitData[i + (p * cc)][1] - allFitData[i + (p * cc)][0] // 0 = y0, 1 = A
		// print out some stats
		Print "Channel "+ num2str(i+1) + " intensity: median", StatsMedian(w)
		WaveStats/Q w
		Print "    Mean", V_avg, "standard deviation", V_sdev
	endfor
	
	// plot width data
	plotName = "p_intensities"
	KillWindow/Z $plotName
	// safe to assume that we have these waves
	WAVE/Z ch1Intensities, ch2Intensities
	if(cc > 2)
		WAVE/Z ch3Intensities
	endif
	Display/N=$plotName; AppendBoxPlot ch1Intensities vs labelAllChWave
	AddWavesToBoxPlot/W=$plotName ch2Intensities
	if(cc > 2)
		AddWavesToBoxPlot/W=$plotName ch3Intensities
	endif
	ModifyBoxPlot/W=$plotName trace=ch1Intensities,whiskerMethod=4
	for(i = 0; i < 3; i += 1)
		ModifyBoxPlot/W=$plotName trace=ch1Intensities,markers[i]={19,-1,19},markerSizes[i]={2,2,2},medianMarkerColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch1Intensities,dataColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767),outlierColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch1Intensities,farOutlierColor[i]=(colorWave[i][0],colorWave[i][1],colorWave[i][2],32767)
	endfor
	Label/W=$plotName left "Intensity"
	SetAxis/A/N=1/E=1/W=$plotName left
End

Function PlotIntensities()
	SetDataFolder root:
	WAVE/Z ch1Intensities,ch2Intensities
	WAVE/Z colorWave
	// note that function requires (and assumes) exact correpsondence between rows of two waves
	// they originate from AllFitData
	Duplicate/O ch1Intensities, ch1Intensities_n
	Duplicate/O ch2Intensities, ch2Intensities_n
	ch1Intensities_n[] = ch1Intensities[p] / WaveMax(ch1Intensities)
	ch2Intensities_n[] = ch2Intensities[p] / WaveMax(ch2Intensities)
	String plotName = "p_intcomp"
	KillWindow/Z $plotName
	Display/N=$plotName ch1Intensities_n vs ch2Intensities_n
	ModifyGraph/W=$plotName mode=3,marker=19,rgb=(0,0,0,32767)
	SetAxis/W=$plotName left 0,1
	SetAxis/W=$plotName bottom 0,1
	Label/W=$plotName left "Ch1 Intensity"
	Label/W=$plotName bottom "Ch2 Intensity"
	ModifyGraph/W=$plotName width={Aspect,1}
	ModifyGraph/W=$plotName axRGB(bottom)=(colorWave[1][0],colorWave[1][1],colorWave[1][2]),tlblRGB(bottom)=(colorWave[1][0],colorWave[1][1],colorWave[1][2]),alblRGB(bottom)=(colorWave[1][0],colorWave[1][1],colorWave[1][2])
	ModifyGraph/W=$plotName axRGB(left)=(colorWave[0][0],colorWave[0][1],colorWave[0][2]),tlblRGB(left)=(colorWave[0][0],colorWave[0][1],colorWave[0][2]),alblRGB(left)=(colorWave[0][0],colorWave[0][1],colorWave[0][2])
End

///	@param	cc	variable, how many channels
Function PlotPeakDistances(cc)
	Variable cc
	
	WAVE allFitData = root:allFitData
	WAVE tripletOKWave = root:data:tripletOKWave
	
	Variable nPits = DimSize(allFitData,0) / cc
	Make/O/N=(nPits) ch2Peak
	ch2Peak[] = allFitData[(p * cc) + 1][2] // the 2nd column contains x0
	// now delete the triplets that were not OK
	ch2Peak[] = (tripletOKWave[p] == 1) ? ch2Peak[p] * 1000 : NaN
	WaveTransform zapnans ch2Peak
	
	if(cc == 3) // optional farRed channel
		Make/O/N=(nPits) ch3Peak
		ch3Peak[] = allFitData[(p * cc) + 2][2] // the 2nd column contains x0
		ch3Peak[] = (tripletOKWave[p] == 1) ? ch3Peak[p] * 1000 : NaN
		WaveTransform zapnans ch3Peak
	endif
	
	if(cc == 2)
		Sort/R ch2Peak, ch2Peak
	else
		Sort/R ch2Peak, ch2Peak, ch3Peak
	endif
	
	Make/O/N=(numpnts(ch2Peak)) chXPeak = p
	
	// plot out pit by pit representation of the distances
	String plotName = "p_peaks"
	KillWindow/Z $plotName
	Display/N=$plotName ch2Peak vs chXPeak
	if(cc == 3)
		AppendToGraph/W=$plotName ch3Peak vs chXPeak
	endif
	ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0
	ModifyGraph/W=$plotName swapXY=1
	SetAxis/A/R/W=$plotName left
	ModifyGraph/W=$plotName noLabel(left)=2
	ModifyGraph/W=$plotName zero(bottom)=4
	SetAxis/A/N=1 bottom
	Label/W=$plotName bottom "Distance from Ch1 (nm)"
	Wave colorWave = root:colorWave
	ModifyGraph/W=$plotName rgb(ch2Peak)=(colorWave[1][0],colorWave[1][1],colorWave[1][2],floor(65535 / 2))
	if(cc == 3)
		ModifyGraph/W=$plotName rgb(ch3Peak)=(colorWave[2][0],colorWave[2][1],colorWave[2][2],floor(65535 / 2))
	endif
	
	// plot out channel3 distances vs channel 2, if cc ==3
	if(cc == 3)
		plotName = "p_twoPeaks"
		KillWindow/Z $plotName
		Display/N=$plotName ch3Peak vs ch2Peak
		ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0,rgb=(0,0,0,32768)
		SetAxis/A/N=1/W=$plotName bottom
		SetAxis/A/N=1/W=$plotName left
		ModifyGraph/W=$plotName width={Aspect,1}
		ModifyGraph/W=$plotName zero=4
		Label/W=$plotName bottom "Ch2 distance from Ch1 (nm)"
		Label/W=$plotName left "Ch3 distance from Ch1 (nm)"

		// box plot of the distances
		plotName =  "p_distances"
		KillWindow/Z $plotName
		Display/N=$plotName
		Make/O/N=(2)/T labelWave={"Ch 2","Ch 3"}
		AppendBoxPlot/W=$plotName ch2peak vs labelWave
		AddWavesToBoxPlot/W=$plotName ch3peak
		ModifyBoxPlot/W=$plotName trace=ch2peak,markers={19,-1,19},markerSizes={2,2,2},whiskerMethod=4
		ModifyBoxPlot/W=$plotName trace=ch2peak,medianMarkerColor[0]=(colorWave[1][0],colorWave[1][1],colorWave[1][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch2peak,dataColor[0]=(colorWave[1][0],colorWave[1][1],colorWave[1][2],32767),outlierColor[0]=(colorWave[1][0],colorWave[1][1],colorWave[1][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch2peak,farOutlierColor[0]=(colorWave[1][0],colorWave[1][1],colorWave[1][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch2peak,medianMarkerColor[1]=(colorWave[2][0],colorWave[2][1],colorWave[2][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch2peak,dataColor[1]=(colorWave[2][0],colorWave[2][1],colorWave[2][2],32767),outlierColor[1]=(colorWave[2][0],colorWave[2][1],colorWave[2][2],32767)
		ModifyBoxPlot/W=$plotName trace=ch2peak,farOutlierColor[1]=(colorWave[2][0],colorWave[2][1],colorWave[2][2],32767)
	
		SetAxis/A/N=1/E=3/W=$plotName left
		Label/W=$plotName left "Distance to peak (nm)"
	endif
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

// This is a function for fitting
///	@param	yW	the data wave
///	@param	xW	the time wave (X)
STATIC Function TheFitter(yW,xW)
	Wave yW,xW
	
	String yName = NameOfWave(yW)
	String gName = yName + "_gauss"
	Make/O/N=(numpnts(yW))/D $gName=NaN
	// single Gauss fit to data
	CurveFit/Q gauss yW /X=xW /D=$gName
	// store coef
	WAVE/Z W_coef,W_sigma
	Variable x0 = W_coef[2]
	// tidy up
	Wave fitWave = $gName
	KillWaves/Z W_coef,W_sigma,fitWave
	
	return x0
End

// This is another function for fitting (single Y wave that is scaled)
// Output is a FREE copy of W_coef
///	@param	yW	the data wave
STATIC Function/WAVE TheFitterForScaledData(yW)
	Wave yW
	
	String yName = NameOfWave(yW)
	String gName = yName + "_gauss"
	Make/O/N=(numpnts(yW))/D $gName=NaN
	// single Gauss fit to data
	CurveFit/Q gauss yW /D=$gName
	// store coef
	WAVE/Z W_coef,W_sigma
	Duplicate/O/FREE W_coef, outputW
	// tidy up
	Wave fitWave = $gName
	KillWaves/Z W_coef,W_sigma,fitWave
	
	return outputW
End


// Colours are taken from Royle Lab stylesheet
// Define colours
StrConstant pal = "0x00a651;0xed1c24;0x64318e;"

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

Function MakeColorWave()
	
	Variable cond = ItemsInList(pal)
	Variable color
	Make/O/N=(cond,3) root:colorWave
	WAVE colorWave = root:colorWave
	Variable i
	
	for(i = 0; i < cond; i += 1)
		// specify colours
		color = str2num(StringFromList(i,pal))
		colorwave[i][0] = hexcolor_red(color)
		colorwave[i][1] = hexcolor_green(color)
		colorwave[i][2] = hexcolor_blue(color)
	endfor
End

STATIC Function MakeTheLayouts(prefix,nRow,nCol,[iter, filtVar, rev, saveIt, orient])
	String prefix
	Variable nRow, nCol
	Variable iter	// this is if we are doing multiple iterations of the same layout
	Variable filtVar // this is the object we want to filter for
	Variable rev // optional - reverse plot order
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