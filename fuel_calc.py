'''
Created on 14 Jun 2013

@author: Robert
'''

import math
import matplotlib.pyplot as plt
#import sys
from MASS_ZBO_SA import calc_mass

'-------------------------Function Definitions---------------------------------------------'
'Function Generates a list of delta-Vs between a minim and maximum value'
def gendeltav(mini,maxi,step):
	'''Generates a list of deltav's based up a minimum, maximum and step size entered.
	Returns as an ordered list (lowest to highest)

	Syntax: gendeltav(minimum,maximum,stepsize)

	minimum = minimum deltav to be calculated
	maximum = maximum deltav to be calculated
	stepsize = difference between steps of deltav 

	Returns: List
	'''

	#Removes a Delta V of 0
	if mini == 0:
		deltav = []
	else:
		deltav = [mini]


	#Generates List of Delta V's from the minimum in increments to the maximum value
	current = mini
	while current < maxi:
		current = current + step
		current = round(current, 3)
		deltav.append(current)
	return deltav

'Function Calculates the mass of a fuel tank based upon its geometry'
def tank_geometryc(volume,diameter,inpudict):
	'''Calculates the mass and dimensions of a fuel tank from an inputted volume and diameter

	Syntax: tank_geometryc(volume,diameter)

	volume = volume of hydrogen to be calculated
	diameter = internal diameter of the tank

	'''

	'Initialise Values'
	#Density of propellant (e.g. hydrogen)
	hydrogendensity = 70.85
	
	#Volume of tank, fuel mass (from density)
	Volt = volume
	fuelmass = hydrogendensity*Volt
	
	#Density of propellant (e.g. Aluminium-Lithium 2195)
	density = 2685
	
	#Stress Factor from Graph
	bigK = 0.79
	
	#Welding Efficiency
	ew = 0.6
	
	#Internal Radius of the Tank
	a = diameter/2
	
	#Ultimate Tensile Strength
	Fu = 447000000
	
	#Yield Strength
	Fy = 312000000
	
	#Maximum Operating Pressure
	pt = 230000

	#Max Working Stress 1
	Sw1 = Fy/1.33
	
	#Max Working Stress 2
	Sw2 = Fu/1.65
	
	#Pick the lowest Max Working Stress
	if Sw1 < Sw2:
		Sw = Sw1
	else:
		Sw = Sw2

	'Calculate for Ellipsoidal End Caps'
	#Chosen Ellipsoidal Radius
	k = 1.395 
	
	#Calculated Ellipsoidal Height
	b = a/k 
	
	#Tank End Crown Radius
	R = k*a 
	
	#Volume of Ellipsoidal Ends
	Volends = (2*math.pi*math.pow(a,2)*b)/3
	
	#Thickness at Knuckle
	tk = (bigK*pt*a)/(Sw*ew) 
	
	#Thickness at Crown
	tcr = (pt*R)/(2*Sw*ew) 
	
	#Equivalent Wall thickness (two methods for comparison)
	te1 = (tk + tcr)/2 
	te2 = (pt*a*(bigK+k/2))/(2*Sw*ew)
	
	#Eccentricity (two methods for comparison)
	e1 = math.pow((math.pow(a,2)-math.pow(b,2)),0.5)/a
	e2 = math.pow((1-(1/math.pow(k,2))),0.5)

	#Surface Area of Ellipsoidal & Insulation Mass
	SAe = math.pow(a,2) + ((math.pi*math.pow(b,2)*math.log((1+e1)/(1-e1))))/(2*e1)
	Wei = SAe*(0.78+0.015*inpudict['mlilayers'])

	#Design Factor 
	desfactor = 2*k + (1/math.pow((math.pow(k,2)-1),0.5))*math.log((k+math.pow(math.pow(k,2)-1,0.5))/(k-math.pow(math.pow(k,2)-1,0.5)))
	
	#Mass of each ellipsoidal end
	We = (math.pi*math.pow(a,2)*te1*desfactor*density)/(2*k) + Wei #Mass of Ellipsoidal Ends
	
	'Calculate for Cylindrical Body Tank'
	#Volume of Cylinder
	Volcyl = Volt - 2*Volends        

	#Length of Cylinder
	lengthc = Volcyl/(math.pi*math.pow(a,2))

	#Thickness of Cylinder
	tc = (pt*a)/(Sw*ew) 

	#Surface Area & Insulation mass
	SAc = 2*math.pi*a*lengthc
	Wci = SAc*(0.78+0.015*inpudict['mlilayers'])

	#Cylindrical Section Mass
	Wc = 2*math.pi*a*lengthc*tc*density + Wci

	'''
	"-----------------------------Summary-----------------------------"
	print('\n------------------------Summary-------------------------|')
	print('--------------------------------------------------------|')
	print('Fuel:\t\tHydrogen')
	print('Fuel Mass:     ', fuelmass, 'Kg')
	print('Fuel Volume:   ', round(Volt,2), 'M^3')
	print('--------------------------------------------------------|')
	print('Ellipsoidal Ends: ')
	print('\t\tAverage thickness:\t', round(te1*1000,2),'\t\tmm')
	print('\t\tRadius of End:\t\t', round(a,2),'\t\tm')
	print('\t\tHeight of End:\t\t', round(b,2),'\t\tm')
	print('\t\tMass of End:\t\t', round(We,2),'\tKg')
	print('--------------------------------------------------------|')
	print('Cylindrical Section:')
	print('\t\tAverage thickness:\t', round(tc*1000,2),'\t\tmm')
	print('\t\tRadius of Cylinder:\t', round(a,2),'\t\tm')
	print('\t\tLength of Cylinder:\t', round(lengthc,2),'\t\tm')
	print('\t\tMass of Cylinder:\t', round(Wc,2),'\tKg')
	print('--------------------------------------------------------|')
	print('Totals: (1x Cylinder + 2x Ends)')
	print('\t\tTotal Length:\t\t',round(lengthc+2*b,2),'\t\tm')
	print('\t\tTotal Mass\t\t', round(2*We + Wc,2),'\tKg')
    
	print('\n--------------------------------------------------------|')
	print('------------------------End-----------------------------|')
	print('--------------------------------------------------------|\n')
	'''

	#Total Mass
	Wt = 2*We + Wc

	#Total Surface Area
	SA = SAc + 2*SAe
	
	#Total Length
	totlength = lengthc+2*b

	#Calculate total thickness of material (1.243 is average thickness of mlilayer according to DRA)
	if tc > te1:
		totalthickness = 25.4 + tc + inpudict['mlilayers']*1.243
	else:
		totalthickness = 25.4 + te + inpudict['mlilayers']*1.243
	
	extradius = (a + (totalthickness/1000))

	#Return Important Values
	return Wt,We,Wc,a,extradius,totalthickness,te1,tc,b,lengthc,totlength,SA

'Function used by fuelmassc to iterate through fuel/tank mass'
def sconvergemass(inpudict,diameter):
	#Density of Liquid Hydrogen
	dlh2 = 70.85
	
	#Define Initial Tankmass
	tankmass = 0	

	#Get ISP and structuremass
	isp = inpudict['Isp']
	structuremass = inpudict['structuremass']
	deltav = inpudict['deltav']
	payloadmass = inpudict['payloadmass']

	#Find Total Initial Mass
	systemmass = payloadmass + structuremass + tankmass
	
	#Define Iteration Counter
	inter = 0

	print('Delta-V: ', deltav, 'Payloadmass: ', payloadmass, 'Inner Diameter: ', diameter)
	
	#Calculate First Guess Fuelmasss and Volume
	#*1.03 includes ullage,1.15 includes 15% error margin (and for propellant cooldown)
	fuelmass = ((systemmass*(math.exp((deltav*1000)/(isp*9.81))-1))*1.03)*1.15
	volume = fuelmass/dlh2
	
	#Calculate fuel Geometry from first guess
	(masstank,ellipmass,cylmass,intradius,extradius,totalthickness,ellipthick,cylthick,ellipheight,cyllength,totlength,totalsa) = tank_geometryc(volume,diameter,inpudict)
	
	#Calculate ZBO & Power System Masses
	(Marray,Mzbo,loss) = calc_mass(totalsa,inpudict['mlilayers'],inpudict['ZBO'],inpudict['days'],inpudict['spow'])
	masstank = masstank + Marray + Mzbo

	lasttankmass = tankmass
	tankmass = masstank
	diff = tankmass - lasttankmass
	print('Converging',end="")

	#Convergent Loop with limiter of 100 iterations
	while ((tankmass-lasttankmass) > 1) and (inter < 100) :
		inter = inter + 1
		systemmass = payloadmass + structuremass + tankmass	
		fuelmass = ((systemmass*(math.exp((deltav*1000)/(isp*9.81))-1))*1.03)*1.15
		volume = fuelmass/dlh2
		(masstank,ellipmass,cylmass,intradius,extradius,totalthickness,ellipthick,cylthick,ellipheight,cyllength,totlength,totalsa) = tank_geometryc(volume,diameter,inpudict)
		
		#Calculate ZBO & Power System Masses
		(Marray,Mzbo,loss) = calc_mass(totalsa,inpudict['mlilayers'],inpudict['ZBO'],inpudict['days'],inpudict['spow'])
		masstank = masstank + Marray + Mzbo

		lasttankmass = tankmass	
		tankmass = masstank
		newsystemass = payloadmass + structuremass + tankmass + fuelmass
		diff = tankmass-lasttankmass
		print('.',end="")
		
	if inter == 100:
		print('\nFailure to converge, data ignored\n')
		return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	else: 
		print('Converged.\n')
		newsystemass = payloadmass + structuremass + tankmass + fuelmass
		tankmass = tankmass - (Mzbo + Marray)
		return newsystemass, Mzbo, Marray, fuelmass, tankmass, ellipmass, cylmass, intradius, extradius, totalthickness, ellipthick, cylthick, ellipheight, cyllength, totlength, totalsa  

'Returns data'
def fuelmassc(diameterlist,inpudict):
	'''Calculates the fuel masses based upon the ideal rocket equation. Uses 
	input of a list of deltav's and a list of payloadmass of equal length. These
	can be generated using gendeltav and genpayloadmass routines

	Syntax: fuelmass(deltav,payloadmass,diameter)

	deltav = list of deltav
	payloadmass = list of payloadmasses
	diameter = list of diameters to be inspected
	isp = isp of vehicle
	'''
	
	print('Calculating mass of fuels for given payload mass and delta v...')
	
	datalist = []
	rlist = importrocketdata()
	
	for i in diameterlist:
		(newsystemass, Mzbo, Marray, fuelmass, 
			tankmass, ellipmass, cylmass, intradius, extradius, 
			totalthickness, ellipthick, cylthick, ellipheight, 
			cyllength, totlength, totalsa) = sconvergemass(inpudict,i)  
		
		if (cyllength > 0) and (fuelmass is not 0):
			clist = []
			
			for a in range(len(rlist)):
				if (newsystemass-inpudict['structuremass']) < rlist[a][3]:
					clist.append(rlist[a])
			
			clist = sorted(clist, key=lambda x: x[4])
			
			datalist.append([inpudict['deltav'], inpudict['payloadmass'], newsystemass, inpudict['structuremass'], 
				Mzbo, Marray, fuelmass, tankmass, ellipmass, cylmass, 
				intradius, extradius, totalthickness, ellipthick, cylthick, 
				ellipheight, cyllength, totlength, totalsa, clist])
	
	datalist = sorted(datalist, key=lambda x: x[2])

	return datalist

'Function saving data to file'
def savefueldata(testfuel,tformat=0,outputfile='output.dat'):
	'''Function to save output data from fuelmass in a .dat file
	
	Syntax: savefueldata(testfuel,optional=tformat,optional=outputfile)
	
	testfuel = list returned from fuelmass function
	tformat = toggle. 1 = formatted, 0 = unformatted. Default = 0.
	outputfile = string with filename to be outputted, default output.dat

	'''

	#Open File and produce initial format
	f = open(outputfile,'a')
	f.write('======================================================================================\n')
	f.write('                                          NEW\n')
	f.write('======================================================================================\n')

	#Select whether formatted data required
	if tformat == 1:
		for a in range(len(testfuel)):
			string1 = ''.join('Delta-V: ' + str(testfuel[a][0]) + 'km/s\n')
			string2 = ''.join('Mass Breakdown:\n')
			string3 = ''.join('\t|Total Mass: ' + str(round(testfuel[a][2],2)) + 'kg\n')
			string4 = ''.join('\t\t|Payload Mass: ' + str(round(testfuel[a][1],2)) + ' kg\n')
			string5 = ''.join('\t\t|NTR Mass: ' + str(round(testfuel[a][3],2)) + ' kg\n')
			string6 = ''.join('\t\t|Fuel Tank Mass: ' + str(round((testfuel[a][1]+testfuel[a][4]+testfuel[a][5]),2)) + ' kg\n')
			string7 = ''.join('\t\t\t|Solar Array Mass: ' + str(round(testfuel[a][5])) + ' kg\n')
			string8 = ''.join('\t\t\t|ZBO System Mass: ' + str(round(testfuel[a][4])) + ' kg\n')
			string9 = ''.join('\t\t\t|Tank Structure Mass: ' + str(round(testfuel[a][7])) + ' kg\n')
			stringa = ''.join('\t\t|Fuel Mass: ' + str(round(testfuel[a][6],2)) + ' kg\n')
			launchmass = testfuel[a][2] - testfuel[a][3]
			stringb = ''.join('\t|Mission Launch Mass: ' + str(round(launchmass,2)) + ' kg\n')
			stringc = ''.join('Tank Parameters:\n')
			stringd = ''.join('\t|Total Length: ' + str(round(testfuel[a][17],2)) + ' m\n')
			stringe = ''.join('\t\t|Cylinder Length: ' + str(round(testfuel[a][16],2)) + ' m\n')
			stringf = ''.join('\t\t|End Height: ' + str(round(testfuel[a][15],2)) + ' m\n')
			stringg = ''.join('\t|Radius:\n')
			stringh = ''.join('\t\t|Inner: ' + str(round(testfuel[a][10],2)) + 'm\n')
			stringi = ''.join('\t\t|Outer: ' + str(round(testfuel[a][11],2)) + 'm\n')
			f.write(string1)
			f.write(string2)
			f.write(string3)
			f.write(string4)
			f.write(string5)
			f.write(string6)
			f.write(string7)
			f.write(string8)
			f.write(string9)
			f.write(stringa)
			f.write(stringb)
			f.write(stringc)
			f.write(stringd)
			f.write(stringe)
			f.write(stringf)
			f.write(stringg)
			f.write(stringh)
			f.write(stringi)
			f.write('Launch Vehicles:\n')
			for x in range(len(testfuel[a][19])):
				string = ''.join('\t' + str(testfuel[a][19][x][0]) + ' ' + str(testfuel[a][19][x][1]) + ':\t' + str(testfuel[a][19][x][3]) + '\t' + str(testfuel[a][19][x][4]) + '\n')
				f.write(string)
			f.write('\n')

			f.write('--------------------------------------------------------------------------------------\n')
	else:
		for a in range(len(testfuel)):
			string = ', '.join(str(x) for x in testfuel[a])
			f.write(string)
			f.write('\n')
	
	f.close()

'Imports rocket data from file in local directory (name: rocketdatafile)'
def importrocketdata():
	f = open('rocketdatafile','r')

	rlist = f.readlines()

	srlist = []

	for x in range(len(rlist)):
		rlist[x] = rlist[x].rstrip('\n')
		srlist.append(rlist[x].split(','))

	for i in range(len(srlist)):
		srlist[i][3] = float(srlist[i][3])
		srlist[i][4] = float(srlist[i][4])

	return srlist


'-------------------------Main-------------------------------------------------------------'		


#sys.stdout = open('file.log', 'w')

#Main Input Deck
inpudict = {'deltav':2.1,
			'payloadmass':6400,
			'days':240,
			'Isp':900,
			'structuremass':4000,
			'mlilayers':60,
			'ZBO':0,
			'spow':1000,
			'tloss':0}

#Input Deck for Iterative Generation Fuctions
iterdict = {'mindiam': 		1,
			'maxdiam':		10,
			'stepdiam': 	0.5}

diameterlist = gendeltav(iterdict['mindiam'],iterdict['maxdiam'],iterdict['stepdiam'])


'Perform Operations'
testfuel = fuelmassc(diameterlist,inpudict)
savefueldata(testfuel,1)