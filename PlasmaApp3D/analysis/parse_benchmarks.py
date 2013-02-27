def getMPIScaling(filename):
	myfile = open(filename,"r").readlines()

	i = 0
	xval = [None]*len(myfile)
	yval = [None]*len(myfile)
	for myline in myfile:
		if(i > 0):
			vals = myline.split(",")
			xval[i-1] = float(vals[8])
			yval[i-1] = float(vals[33])
			print xval[i-1], yval[i-1]
		i+=1

	scale = yval[0]
	for j in range(len(xval)-1):
		yval[j] = scale/(yval[j]*xval[j])
		print xval[j], yval[j]

	values = [xval[:(len(xval)-1)],yval]
	return values

def getOMPScaling(filename):
	myfile = open(filename,"r").readlines()

	i = 0
	tdict = {}
	for myline in myfile:
		if(i > 0):
			vals = myline.split(",")
			temp = int(vals[7])
			tdict[temp] = float(vals[17])
			print tdict[temp], temp
		i+=1

	xval = tdict.keys()
	xval.sort()
	yval = tdict.values()
	k = 0
	for x in xval:
		yval[k] = tdict[x]
		k+=1

	scale = yval[0]
	for j in range(len(xval)):
		print xval[j], yval[j]
		yval[j] = scale/(yval[j]*xval[j])
		
#	xval.remove(len(xval)-1)
	values = [xval,yval]
	return values

def getSIMDScaling(filename):
	myfile = open(filename,"r").readlines()

	i = 0
	tdict = {}
	for myline in myfile:
		if(i > 0):
			vals = myline.split(",")
			temp = int(vals[6])
			tdict[temp] = float(vals[17])
			print tdict[temp], temp
		i+=1

	xval = tdict.keys()
	xval.sort()
	yval = tdict.values()
	k = 0
	for x in xval:
		yval[k] = tdict[x]
		k+=1

	scale = yval[0]
	for j in range(len(xval)):
		print xval[j], yval[j]
		yval[j] = scale/(yval[j])
		
#	xval.remove(len(xval)-1)
	values = [xval,yval]
	return values

def writeVals(fileout,xval,yval):
	fout = open(fileout,"w")

	for i in range(len(xval)):
		print xval[i], yval[i]
		value = (xval[i]," ",yval[i])
		s = '{} {}\n'.format(xval[i],yval[i])
		fout.write(s)
	fout.close()
###############################################


filename = "./benchmarks_d/benchmark1.dat"
values = getMPIScaling(filename)

fileout = "mpiscaling_d.dat"

xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);

filename = "./benchmarks_d/benchmark2.dat"
values = getMPIScaling(filename)

fileout = "mpiscaling2_d.dat"

xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);

filename = "./benchmarks_m/benchmark1.dat"
values = getMPIScaling(filename)

fileout = "mpiscaling_m.dat"


xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);

filename = "./benchmarks_m/benchmark2.dat"
values = getMPIScaling(filename)

fileout = "mpiscaling2_m.dat"


xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);
######################################
print 'Getting SIMD Scalings'
filename = "./benchmarks_d/benchmark5.dat"
values = getSIMDScaling(filename)

fileout = "simdscaling_d.dat"


xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);

filename = "./benchmarks_m/benchmark3.dat"
values = getSIMDScaling(filename)

fileout = "simdscaling_m.dat"


xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);

######################################
print 'Getting OMP Scalings'
filename = "./benchmarks_d/benchmark4.dat"
values = getOMPScaling(filename)

fileout = "ompscaling_d.dat"


xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);

filename = "./benchmarks_m/benchmark4.dat"
values = getOMPScaling(filename)

fileout = "ompscaling_m.dat"


xval = values[0]
yval = values[1]

writeVals(fileout,xval,yval);
