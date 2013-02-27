import shutil
import numpy as np
import pylab as pylab
import matplotlib.pyplot as plt

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

def getValues(benchmarkName,valueName,typeout,offset=0):
	
	icolumn = findColumn(benchmarkName,valueName)+offset
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
	
	
def extract_timings(benchmarkName):
	timing_names = ("Particle\n Read",
	"Time Step\n Estimator",
	"Cell\n Crossing",
	"Intrp\n Accel",
	"Current\n Tally",
	"Picard\n Math",
	"Comm",
	"LO-Solve")
	
	timing_ref_names = ["PLoadStore_time(ms)",
	"Dtau_est_time(ms)",
	"Crossing_time(ms)",
	"Accel_time(ms)",
	"Tally_time(ms)",
	"PPiccardOther_time(ms)",
	"Comm_time(ms)",
	"LOSolve_time(ms)"]
	
	timing_ref_names = timing_ref_names[::-1]
	timing_names = timing_names[::-1]
	
	
	times = [[None]]*(len(timing_names))
	errors = [[None]]*len(timing_names)
	

	
	Nsteps = getValues(benchmarkName,"nsubsteps_total",float)
	NNodes = getValues(benchmarkName,"num_nodes",int)
	
	for i in range(len(timing_names)):
		times[i] = getValues(benchmarkName,timing_ref_names[i],float)
		errors[i] = getValues(benchmarkName,timing_ref_names[i],float,1)
	
	for j in range(len(timing_names)):
		for i in range(len(times[j])):
			scale = 1.0e6*NNodes[i]/Nsteps[i]
			times[j][i] *= scale
			errors[j][i] *= scale
	
		print len(times[0])
	
	return [timing_names,times,errors] 
	
def make_rect(timing_names,times_in,errors_in,ax,irow,ibar,width=0.35):
	times = [None]*len(timing_names)
	errors = [None]*len(timing_names)
	colors = ['r','y','g','b','purple']
#	ecolors = ['pink','yellow','green','blue','purple']
	
	for i in range(len(timing_names)):
		times[i] = times_in[i][irow]
		errors[i] = errors_in[i][irow]
		
	ind = np.arange(len(timing_names))

	
	rects = ax.barh(ind+width*ibar,times,align='center',height=width,
					color=colors[ibar])
	
	return rects
	
	
  

###############################################




filename = "../benchmarks/benchmark21827.dat"
outputname = "timing_breakdown21827.tex"

width = 0.35


[timing_names,times,errors] = extract_timings(filename)
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(len(timing_names))
plt.yticks(ind,timing_names)
rects1 = make_rect(timing_names,times,errors,ax,0,0)


ax.axvline(0,color='k',lw=3)   # poor man's zero level

ax.set_xlabel('Time per particle-subcycle (ns)')
ax.set_title('2D3V Initial Performance Analysis')

x1,x2,y1,y2 = plt.axis()

plt.subplots_adjust(left=0.18, right=0.82)
plt.axis((x1,x2,y1-0.5,y2+0.5))
pylab.savefig('2D3V_Perf.pdf',transparent=True)
#plt.show()

##################################################

LegendNames = ('Opt=000','Opt=100','Opt=110','Opt=111','Opt=101')
width=0.15
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ind = np.arange(len(timing_names))
plt.yticks(ind+1.5*width,timing_names)

rects = [None]*5

for i in range(5):
	rects[i] = make_rect(timing_names,times,errors,ax2,i,i,0.15)


plt.legend((rects[0][0],rects[1][0],rects[2][0],rects[3][0],rects[4][0]),LegendNames)
ax2.axvline(0,color='k',lw=3)   # poor man's zero level

ax2.set_xlabel('Time per particle-subcycle (ns)')
ax2.set_title('2D3V Optimized Performance Analysis')

x1,x2,y1,y2 = plt.axis()

plt.subplots_adjust(left=0.18, right=0.82)
plt.axis((x1,x2,y1-0.2,y2+0.2))
pylab.savefig('2D3V_OptPerf.pdf',transparent=True)
#plt.show()

###############################################




filename = "../benchmarks/benchmark22201m.dat"


width = 0.35


[timing_names,times,errors] = extract_timings(filename)
fig = plt.figure()
ax3 = fig.add_subplot(111)
ind = np.arange(len(timing_names))
plt.yticks(ind,timing_names)
rects1 = make_rect(timing_names,times,errors,ax3,1,0)


ax3.axvline(0,color='k',lw=3)   # poor man's zero level

ax3.set_xlabel('Time per particle-subcycle (ns)')
ax3.set_title('Two stream instability on 32 nodes 16 Cores per Node')

x1,x2,y1,y2 = plt.axis()

plt.subplots_adjust(left=0.18, right=0.82)
plt.axis((x1,x2,y1-0.5,y2+0.5))
pylab.savefig('1D1V_Perf.pdf',transparent=True)


