import shutil

def findColumn(filename,columnname):
	firstline = open(filename,"r").readline()
	i=0
	xvals = firstline.split(",")
	while xvals[i] != columnname:
	  i+=1
	  if(i > len(xvals)):
		i = -1
		break

	return i
	
def writeTikzArray(filename,arrayname,values,iAppend):
	linename = "\\def\\" + arrayname
	myfile = open(filename,"r")
	mylines = myfile.readlines()
	myfile.close()
	i = 0
	iline = 0;

	for myline in mylines:
		vals = myline.split("{")
#		print vals[0]
		if(vals[0] == linename):
			iline = i
			break
		i+=1
	
	
	
	templine = mylines[iline].split("{")
	tempvals = templine[0]
	
	if(isinstance(values[0],float)):
		valuesout = "{0:2.4g}".format(values[0])
	else:
		valuesout = str(values[0])
		
	for j in range(len(values)-1):
		tval=""
		if(isinstance(values[j+1],float)):
			tval = "{0:2.4g}".format(values[j+1])
		else:
			tval = str(values[j+1])
		valuesout += ","+str(tval)
	
	if(iAppend):
		lineout = linename+"{"+tempvals+","+valuesout+"}\n"
	else:
		lineout = linename+"{"+valuesout+"}\n"
	print iline
	
	mylines[iline] = lineout
	
	fileout = open(filename,"w")
	fileout.writelines(mylines)
	fileout.close()

def getValues(benchmarkName,valueName,typeout):
	
	icolumn = findColumn(benchmarkName,valueName)
	myfile = open(filename,"r").readlines()
	times = [None]*(len(myfile)-1)
	
	i = 0
	for myline in myfile:
		if(i > 0):
			vals = myline.split(",")
			times[i-1] = typeout(vals[icolumn])
		i+=1
	
	return times

def writeTimings(benchmarkName,outputName):
	
	PushTime=getValues(benchmarkName,"Push_time(ms)",float)
	OPiccard=getValues(benchmarkName,"Step_time(ms)",float)
	Tstep=getValues(benchmarkName,"Total_time(ms)",float)
	CommTime=getValues(benchmarkName,"Comm_time(ms)",float)
	LOSolve=getValues(benchmarkName,"LOSolve_time(ms)",float)
	PPiccard = getValues(benchmarkName,"Ppicard_time(ms)",float)
	AccelTime = getValues(benchmarkName,"Accel_time(ms)",float)
	TallyTime = getValues(benchmarkName,"Tally_time(ms)",float)
	CrossTime = getValues(benchmarkName,"Crossing_time(ms)",float)
	DtauTime = getValues(benchmarkName,"Dtau_est_time(ms)",float)
	PLoadTime = getValues(benchmarkName,"PLoadStore_time(ms)",float)
	POther = getValues(benchmarkName,"PPiccardOther_time(ms)",float)
	OutputID = getValues(benchmarkName,"OutputID",str)
	
	NNodes = getValues(benchmarkName,"num_nodes",int)
	Nsteps = getValues(benchmarkName,"nsubsteps_total",float)
	NCores = getValues(benchmarkName,"num_cores",int)
	Nptcls = getValues(benchmarkName,"nptcls_total",int)
	Nx = getValues(benchmarkName,"ncells",int)
	VecL = getValues(benchmarkName,"vector_length",int)
	
	# Convert ms to s
	for i in range(len(PushTime)):
		scale = 1.0e6/Nsteps[i]
		PushTime[i] *= scale
		OPiccard[i] *= scale
		Tstep[i] *= scale
		CommTime[i] *= scale
		LOSolve[i] *= scale
		PPiccard[i] *= scale
		AccelTime[i] *= scale
		TallyTime[i] *= scale
		CrossTime[i] *= scale
		DtauTime[i] *= scale
		PLoadTime[i] *= scale
		POther[i] *= scale
		
		Nx[i] *= Nx[i]
		NCores[i] *= NNodes[i]
		Nptcls[i] = int(Nptcls[i]*1.0/Nx[i])
		
		OutputID[i] = OutputID[i].split("/")[2]
		
	# Get the Run Parameters
	
	
	
	
	
	
	Arrays = [(Tstep,"A"),(OPiccard,"B"),(PushTime,"C"),(PPiccard,"D"),(CommTime,"E"),
	(LOSolve,"F"),(NNodes,"NNode"),(NCores,"NCore"),(Nx,"NX"),(Nsteps,"NStep"),
	(OutputID,"ID"),(Nptcls,"Nptcls"),
	(AccelTime,"Accel"),(TallyTime,"Current"),(CrossTime,"Crossing"),
	(DtauTime,"DtauEst"),(VecL,"VecL"),(PLoadTime,"PRead"),(POther,"POther")]
	
	
	
	shutil.copyfile("timing_breakdown.tex",outputName)
	
	for (myarray,arrayname) in Arrays:
		print myarray[0], arrayname
		writeTikzArray(outputName,"Array"+arrayname,myarray,False)
		
	writeTikzArray(outputName,"SourceFile",[benchmarkName],False)

###############################################

filename = "../benchmarks/benchmark21827.dat"
outputname = "timing_breakdown21827.tex"

writeTimings(filename,outputname)


