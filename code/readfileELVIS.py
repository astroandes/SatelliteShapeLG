import numpy as np
import scipy as sp
import scipy.stats
import pdb
import matplotlib.pyplot as plt		# paquete para graficar
import matplotlib.cm as cm


#PairName ="Romulus&Remus"
#PairNumber = 0
#PairName ="Scylla&Charybdis"
#PairNumber = 1
#PairName ="Kek&Kauket"
#PairNumber = 2
#PairName ="Hamilton&Burr"
#PairNumber = 3
#PairName ="Orion&Taurus"
#PairNumber = 4
#PairName ="Siegfried&Roy"
#PairNumber = 5
#PairName ="Sonny&Cher"
#PairNumber = 6
#PairName ="Lincoln&Douglas"
#PairNumber = 7
#PairName ="Thelma_and_Louise"
#PairNumber = 8
#PairName ="Zeus&Hera"
#PairNumber = 9
#PairName ="Hall&Oates"
#PairNumber = 10
PairName ="Serena&Venus"
PairNumber = 11 


##plotting routines
def HaloPositionsplot(HaloX_H1, HaloY_H1, HaloZ_H1, HaloX_H2, HaloY_H2, HaloZ_H2):
	plt.figure()
	plt.plot(HaloX_H1,HaloY_H1,'.', markersize=5,color="blue")
	plt.plot(HaloX_H2,HaloY_H2,'.', markersize=5,color="red")
	plt.xlabel('X (Mpc)', fontsize=15) 
	plt.ylabel('Y (Mpc)', fontsize=15)
	plt.savefig('HalosXY'+PairName+'.pdf')
	
	plt.figure()
	plt.plot(HaloX_H1,HaloZ_H1,'.', markersize=5,color="blue")
	plt.plot(HaloX_H2,HaloZ_H2,'.', markersize=5,color="red")
	plt.xlabel('X (Mpc)', fontsize=15) 
	plt.ylabel('Z (Mpc)', fontsize=15)
	plt.savefig('HalosXZ'+PairName+'.pdf')
	
	plt.figure()
	plt.plot(HaloY_H1,HaloZ_H1,'.', markersize=5,color="blue")
	plt.plot(HaloY_H2,HaloZ_H2,'.', markersize=5,color="red")
	plt.xlabel('Y (Mpc)', fontsize=15) 
	plt.ylabel('Z (Mpc)', fontsize=15)
	plt.savefig('HalosYZ'+PairName+'.pdf')

def Halo_K_L_plot(posHaloH1, L_haloH1, posHaloH2, L_haloH2):
		
##	posM31Sat =[0.06752411772, 
##	0.08597745082, 
##	0.27778807473,
##	0.22235086371,
##	0.12998620354,
##	0.10268162264,
##	0.16318735258,
##	0.31951258637,
##	0.06706897363, 
##	0.11605411258,
##	0.12902180114,
##	0.13578819570,
##	0.28098256183,
##	0.12831966729,
##	0.16943994399,
##	0.09065788310, 
##	0.39249727963,
##	0.11791725543]
##	
##	posM31SatPlane =[0.0675241177239,
##	0.0859774508276,
##	0.102681622647 ,
##	0.163187352585 ,
##	0.31951258637  ,
##	0.0670689736352,
##	0.0906578831035,
##	0.117917255434]
##		
##	LM31Sat = [24.0308601868,
##	24.2624703503,
##	35.8501368509,
##	15.4310423895,
##	2.75518707717,
##	19.121089961 ,
##	31.8689008181,
##	15.6958941041,
##	24.4130565045,
##	23.5620534898,
##	28.6241226452,
##	5.02868270358,
##	58.306195574 ,
##	0.18339024743,
##	24.9534853556,
##	27.4851732919,
##	7.34478011433,
##	29.6213335464]
##		
##	LM31SatPlane = [24.0308601868,
##	24.2624703503,
##	19.121089961 ,
##	31.8689008181,
##	15.6958941041,
##	24.4130565045,
##	27.4851732919,
##	29.6213335464]

	posM31SatPlaneTAN = [ 0.11]

	#LM31SatPlaneTAN =[44.8696208889]
	LM31SatPlaneTAN =[10]

#http://arxiv.org/pdf/1503.07176v1.pdf
#Fornax    140   -31.8  1.7    196.0  29.0   1.00
#LeoI      261   167.9  2.8    101.0  34.4   1.00
#LMC       49    93.2   3.7    346.0  8.5    0.99
#SMC       60    6.8    2.4    259.0  17.0   0.99
#Sculptor  87    79.0   6.0    198.0  50.0   1.00
#Draco     92    -98.5  2.6    210.0  25.0   1.00


	posMWSat = [140,261,49,60,87,92]
	vrMWSat =[-31.8,167.9,93.2 ,6.8,79.0,-98.5]
	vtanMWSat =[196.0,101.0,346.0,259.0,198.0,210.0]	
	vtotMWSat = []
	LMWSat = []
	
	for i in xrange(6):
		posMWSat[i] = 0.001*posMWSat[i]
		LMWSat_i = posMWSat[i]*vtanMWSat[i]
		LMWSat.append(LMWSat_i)

	plt.figure()
	plt.plot(K_haloH1, L_haloH1,'.', markersize=5,color="blue")
	plt.plot(K_haloH2, L_haloH2,'.', markersize=5,color="red")
	#plt.plot(posM31SatPlane, LM31SatPlane ,'*', markersize=8,color="black")
	#plt.plot(posM31SatPlaneTAN, LM31SatPlaneTAN ,'*', markersize=12,color="grey")
	#plt.plot(posM31Sat, LM31Sat ,'.', markersize=8,color="black")
	#plt.plot(posMWSat, LMWSat ,'*', markersize=8,color="green")
	plt.xlabel('Kinetic energy per unit mass (km^2/s^2)', fontsize=15) 
	plt.ylabel('Angular momentum per unit mass (Mpc.km/s)', fontsize=15)
	plt.savefig('Halo_K_L'+PairName+'.pdf')

def Halo_Pos_L_plot(posHaloH1, L_haloH1, posHaloH2, L_haloH2):
		
	posM31Sat =[0.06752411772, 
	0.08597745082, 
	0.27778807473,
	0.22235086371,
	0.12998620354,
	0.10268162264,
	0.16318735258,
	0.31951258637,
	0.06706897363, 
	0.11605411258,
	0.12902180114,
	0.13578819570,
	0.28098256183,
	0.12831966729,
	0.16943994399,
	0.09065788310, 
	0.39249727963,
	0.11791725543]
	
	posM31SatPlane =[0.0675241177239,
	0.0859774508276,
	0.102681622647 ,
	0.163187352585 ,
	0.31951258637  ,
	0.0670689736352,
	0.0906578831035,
	0.117917255434]
		
	LM31Sat = [24.0308601868,
	24.2624703503,
	35.8501368509,
	15.4310423895,
	2.75518707717,
	19.121089961 ,
	31.8689008181,
	15.6958941041,
	24.4130565045,
	23.5620534898,
	28.6241226452,
	5.02868270358,
	58.306195574 ,
	0.18339024743,
	24.9534853556,
	27.4851732919,
	7.34478011433,
	29.6213335464]
		
	LM31SatPlane = [24.0308601868,
	24.2624703503,
	19.121089961 ,
	31.8689008181,
	15.6958941041,
	24.4130565045,
	27.4851732919,
	29.6213335464]

	posM31SatPlaneTAN = [ 0.11]

	#LM31SatPlaneTAN =[44.8696208889]
	LM31SatPlaneTAN =[10]

#http://arxiv.org/pdf/1503.07176v1.pdf
#Fornax    140   -31.8  1.7    196.0  29.0   1.00
#LeoI      261   167.9  2.8    101.0  34.4   1.00
#LMC       49    93.2   3.7    346.0  8.5    0.99
#SMC       60    6.8    2.4    259.0  17.0   0.99
#Sculptor  87    79.0   6.0    198.0  50.0   1.00
#Draco     92    -98.5  2.6    210.0  25.0   1.00


	posMWSat = [140,261,49,60,87,92]
	vrMWSat =[-31.8,167.9,93.2 ,6.8,79.0,-98.5]
	vtanMWSat =[196.0,101.0,346.0,259.0,198.0,210.0]	
	vtotMWSat = []
	LMWSat = []
	
	for i in xrange(6):
		posMWSat[i] = 0.001*posMWSat[i]
		LMWSat_i = posMWSat[i]*vtanMWSat[i]
		LMWSat.append(LMWSat_i)

	plt.figure()
	plt.plot(posHaloH1, L_haloH1,'.', markersize=5,color="blue")
	plt.plot(posHaloH2, L_haloH2,'.', markersize=5,color="red")
	plt.plot(posM31SatPlane, LM31SatPlane ,'*', markersize=8,color="black")
	plt.plot(posM31SatPlaneTAN, LM31SatPlaneTAN ,'*', markersize=12,color="grey")
	plt.plot(posM31Sat, LM31Sat ,'.', markersize=8,color="black")
	plt.plot(posMWSat, LMWSat ,'*', markersize=8,color="green")
	plt.xlabel('Distance to parent halo (Mpc)', fontsize=15) 
	plt.ylabel('Angular momentum per unit mass (Mpc.km/s)', fontsize=15)
	plt.savefig('Halo_pos_L'+PairName+'.pdf')

def cosangle_prob_plot(cosAngleH1,cosAngleH2, eigenvector, SubHalo1, SubHalo2):
	plt.figure()
	x = np.linspace(0, 20, 1000)  # 100 evenly-spaced values from 0 to 50
	y = x

	#cosAngleRandom = np.random.normal(size=300) 
	countsH1, startH1, dxH1, _ = scipy.stats.cumfreq(cosAngleH1, numbins=50)
	countsH2, startH2, dxH2, _ = scipy.stats.cumfreq(cosAngleH2, numbins=50)
	#countsR, startR, dxR, _ = scipy.stats.cumfreq(cosAngleRandom, numbins=15)
	xH1 = np.arange(countsH1.size) * dxH1 + startH1
	xH2 = np.arange(countsH2.size) * dxH2 + startH2
	#xR = np.arange(countsR.size) * dxR + startR
	countsH1 = countsH1/SubHalo1
	countsH2 = countsH2/SubHalo2
	#countsR = countsR/300
	
	plt.plot(xH1, countsH1, linewidth=3, color = "red")#'ro')
	plt.plot(xH2, countsH2, linewidth=3, color = "blue")#'bo')
	plt.plot(x, y, color="black")
	#plt.plot(xR, countsR, 'go')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('Cos theta '+eigenvector)
	plt.ylabel('Cumulative Frequency')
	plt.savefig('freq_pos_dot_'+eigenvector+'_'+PairName+'.pdf')
	
	#plt.show()

def Angle_HaloPos_eigenvector_plot(cosAngleH1,cosAngleH2, eigenvector, SubHalo1, SubHalo2):
	#cosAngleH1 = np.asarray(cosAngleH1)/SubHalo1
	#cosAngleH2 = np.asarray(cosAngleH2)/SubHalo2
	plt.figure()
	fig = plt.figure(figsize=(6, 6))
	ax = fig.add_subplot(1,1,1) # one row, one column, first plot
	#plt.hist(cosAngleH1,  bins =25, color ="blue", alpha =0.7 )
	#plt.hist(cosAngleH2,  bins =25, color ="red", alpha =0.7 )
	ax.hist((cosAngleH1, cosAngleH2),  bins =15, color =("red", "blue"))
	plt.xlabel('cos angle pos '+eigenvector, fontsize=15) 
	plt.ylabel('N Halos', fontsize=15)
	plt.savefig('Hist_pos_dot_'+eigenvector+'_'+PairName+'.pdf')
	
def FileOpenRead():
	count = 0
	f = open('../ELVIS_Halo_Catalogs/PairedCatalogs/'+PairName+".txt", 'r')

	header1 = f.readline()
	header2 = f.readline()
	header3 = f.readline()
	halo=[]
	data=[]
	for line in f:
		data = line.strip().split()
		cols = []
		for col in data:
			cols.append(float(col))
		halo.append(cols)
		count= count +1
	print "Total number of Halos",count
	halo = np.asarray(halo)
	return halo, count

def FileOpenReadEigenvectors():
	count3 = 0
	f = open('../Jaime/Eigenvectorsylambdas.tex', 'r')
	header1 = f.readline()
	LambdaEigenvec=[]
	data=[]
	for line in f:
		data = line.strip().split()
		cols = []
		for col in data:
			cols.append(float(col))
		LambdaEigenvec.append(cols)
		count3= count3 +1
	#LambdaEigenvec = np.asarray(LambdaEigenvec)
	return LambdaEigenvec, count3

halo, count = FileOpenRead()
LambdaEigenvec, count3 = FileOpenReadEigenvectors()

count2 = 0
count4 = 0
count5 = 0
cosAngleH1e1=[]
cosAngleH2e1=[]
cosAngleH1e2=[]
cosAngleH2e2=[]
cosAngleH1e3=[]
cosAngleH2e3=[]
cosAngleH1randomVec=[]
cosAngleH2randomVec=[]
L_haloH1=[]
L_haloH2=[]
K_haloH1=[]
K_haloH2=[]
posHaloH1=[]
posHaloH2=[]

HaloX_H1=[] 
HaloY_H1=[] 
HaloZ_H1=[] 
HaloX_H2=[] 
HaloY_H2=[] 
HaloZ_H2=[] 


HaloVX_H1=[] 
HaloVY_H1=[] 
HaloVZ_H1=[] 
HaloVX_H2=[] 
HaloVY_H2=[] 
HaloVZ_H2=[] 

for j in xrange(int(count)):
	if j > 1:
		if halo[j][10] > 250000000:
			count2 = count2 + 1
			distHalo1 = np.sqrt((halo[j][1]-halo[0][1])*(halo[j][1]-halo[0][1])+(halo[j][2]-halo[0][2])*(halo[j][2]-halo[0][2])+(halo[j][3]-halo[0][3])*(halo[j][3]-halo[0][3]))  
			distHalo2 = np.sqrt((halo[j][1]-halo[1][1])*(halo[j][1]-halo[1][1])+(halo[j][2]-halo[1][2])*(halo[j][2]-halo[1][2])+(halo[j][3]-halo[1][3])*(halo[j][3]-halo[1][3]))  
			if halo[j][17] == halo[0][0]: 
				posH1_X=halo[j][1]-halo[0][1]
				posH1_Y=halo[j][2]-halo[0][2]
				posH1_Z=halo[j][3]-halo[0][3]
				HaloX_H1.append(posH1_X)
				HaloY_H1.append(posH1_Y)
				HaloZ_H1.append(posH1_Z)
				e3= np.array([ LambdaEigenvec[PairNumber][12],LambdaEigenvec[PairNumber][13],LambdaEigenvec[PairNumber][14]]) 
				posHalo = np.array([ posH1_X, posH1_Y, posH1_Z])
				posHalo_modulus = np.sqrt((posHalo*posHalo).sum())

				e3_modulus = np.sqrt((e3*e3).sum())
				dotH1e3 = np.dot(posHalo, e3)
				cos_angleH1e3 = np.absolute (dotH1e3 / posHalo_modulus / e3_modulus) # cosine of angle between poshalo and e3
				#angleH1e3 = np.arccos(cos_angleH1e3)
				#angleH1e3 = angleH1e3* 180 / np.pi
				cosAngleH1e3.append(cos_angleH1e3)

				e1= np.array([ LambdaEigenvec[PairNumber][6],LambdaEigenvec[PairNumber][7],LambdaEigenvec[PairNumber][8]]) 
				e1_modulus = np.sqrt((e1*e1).sum())
				dotH1e1 = np.dot(posHalo, e1)
				cos_angleH1e1 = np.absolute (dotH1e1 / posHalo_modulus / e1_modulus) # cosine of angle between poshalo and e1
				#angleH1e1 = np.arccos(cos_angleH1e1)
				#angleH1e1 = angleH1e1* 180 / np.pi
				cosAngleH1e1.append(cos_angleH1e1)

				e2= np.array([ LambdaEigenvec[PairNumber][9],LambdaEigenvec[PairNumber][10],LambdaEigenvec[PairNumber][11]]) 
				e2_modulus = np.sqrt((e2*e2).sum())
				dotH1e2 = np.dot(posHalo, e2)
				cos_angleH1e2 = np.absolute (dotH1e2 / posHalo_modulus / e2_modulus) # cosine of angle between poshalo and e2
				#angleH1e2 = np.arccos(cos_angleH1e2)
				#angleH1e2 = angleH1e2* 180 / np.pi
				cosAngleH1e2.append(cos_angleH1e2)

				randomVec= np.array([ 0, 0, 1]) 
				randomVec_modulus = np.sqrt((randomVec*randomVec).sum())
				dotH1randomVec = np.dot(posHalo, randomVec)
				cos_angleH1randomVec = np.absolute (dotH1randomVec / posHalo_modulus / randomVec_modulus) # cosine of angle between poshalo and e2
				#angleH1randomVec = np.arccos(cos_angleH1randomVec)
				#angleH1randomVec= angleH1randomVec* 180 / np.pi
				cosAngleH1randomVec.append(cos_angleH1randomVec)
				#print AngleH1randomVec
				count4 = count4+1
				print "4", count4

				#looking for structure in phase space
				velH1_X=halo[j][4]-halo[0][4]
				velH1_Y=halo[j][5]-halo[0][5]
				velH1_Z=halo[j][6]-halo[0][6]
				HaloVX_H1.append(velH1_X)
				HaloVY_H1.append(velH1_Y)
				HaloVZ_H1.append(velH1_Z)
				velHalo = np.array([ velH1_X, velH1_Y, velH1_Z])
				# Angular momentum per unit mass
				L_halo = np.cross(posHalo, velHalo)	
				L_halo_modulus = np.sqrt((L_halo*L_halo).sum())
				HaloVY_H1.append(velH1_Y)
				posHaloH1.append(posHalo_modulus)
				L_haloH1.append(L_halo_modulus)
				# kinetic energy per unit mass
				VelVel = ((velHalo*velHalo).sum())	
				K_halo = 0.5*VelVel
				K_haloH1.append(K_halo)

			if halo[j][17] == halo[1][0]: 
				posH2_X=halo[j][1]-halo[1][1]
				posH2_Y=halo[j][2]-halo[1][2]
				posH2_Z=halo[j][3]-halo[1][3]
				HaloX_H2.append(posH2_X)
				HaloY_H2.append(posH2_Y)
				HaloZ_H2.append(posH2_Z)
				posHalo = np.array([ posH2_X, posH2_Y, posH2_Z])
				posHalo_modulus = np.sqrt((posHalo*posHalo).sum())

				e3= np.array([ LambdaEigenvec[PairNumber][21],LambdaEigenvec[PairNumber][22],LambdaEigenvec[PairNumber][23]]) 
				e3_modulus = np.sqrt((e3*e3).sum())
				dotH2e3 = np.dot(posHalo, e3)
				cos_angleH2e3 = np.absolute (dotH2e3 / posHalo_modulus / e3_modulus) # cosine of angle between poshalo and e3
				#angleH2e3 = np.arccos(cos_angleH2e3)
				#angleH2e3 = angleH2e3* 180 / np.pi
				cosAngleH2e3.append(cos_angleH2e3)

				e1= np.array([ LambdaEigenvec[PairNumber][15],LambdaEigenvec[PairNumber][16],LambdaEigenvec[PairNumber][17]])
				e1_modulus = np.sqrt((e1*e1).sum())
				dotH2e1 = np.dot(posHalo, e1)
				cos_angleH2e1 = np.absolute (dotH2e1 / posHalo_modulus / e1_modulus) # cosine of angle between poshalo and e3
				#angleH2e1 = np.arccos(cos_angleH2e1)
				#angleH2e1 = angleH2e1* 180 / np.pi
				cosAngleH2e1.append(cos_angleH2e1)

				e2= np.array([ LambdaEigenvec[PairNumber][18],LambdaEigenvec[PairNumber][19],LambdaEigenvec[PairNumber][20]])
				e2_modulus = np.sqrt((e2*e2).sum())
				dotH2e2 = np.dot(posHalo, e2)
				cos_angleH2e2 = np.absolute (dotH2e2 / posHalo_modulus / e2_modulus) # cosine of angle between poshalo and e2
				#angleH2e2 = np.arccos(cos_angleH2e2)
				#angleH2e2 = angleH2e2* 180 / np.pi
				cosAngleH2e2.append(cos_angleH2e2)

				randomVec= np.array([ 0, 0, 1]) 
				randomVec_modulus = np.sqrt((randomVec*randomVec).sum())
				dotH2randomVec = np.dot(posHalo, randomVec)
				cos_angleH2randomVec = np.absolute (dotH2randomVec / posHalo_modulus / randomVec_modulus) # cosine of angle between poshalo and e2
				#angleH2randomVec = np.arccos(cos_angleH2randomVec)
				#angleH2randomVec= angleH2randomVec* 180 / np.pi
				cosAngleH2randomVec.append(cos_angleH2randomVec)

				count5 = count5+1
				print "5", count5

				#looking for structure in phase space
				velH2_X=halo[j][4]-halo[1][4]
				velH2_Y=halo[j][5]-halo[1][5]
				velH2_Z=halo[j][6]-halo[1][6]
				HaloVX_H2.append(velH2_X)
				HaloVY_H2.append(velH2_Y)
				HaloVZ_H2.append(velH2_Z)
				velHalo = np.array([ velH2_X, velH2_Y, velH2_Z])
				# Angular momentum per unit mass
				L_halo = np.cross(posHalo, velHalo)	
				L_halo_modulus = np.sqrt((L_halo*L_halo).sum())
				posHaloH2.append(posHalo_modulus)
				L_haloH2.append(L_halo_modulus)
				# kinetic energy per unit mass
				VelVel = ((velHalo*velHalo).sum())	
				K_halo = 0.5*VelVel
				K_haloH2.append(K_halo)

SubHalo1 = count4
SubHalo2 = count5
				
print "Number of Halos after mass cut", count2	
#print "Number of Halos after mass and distance cuts", count4	
print "."	
#print AngleH1randomVec
print "."	
#print AngleH2randomVec
print "."	
#print AngleH2e2


HaloPositionsplot(HaloX_H1, HaloY_H1, HaloZ_H1, HaloX_H2, HaloY_H2, HaloZ_H2)
Angle_HaloPos_eigenvector_plot(cosAngleH1e3,cosAngleH2e3, "e3", SubHalo1, SubHalo2)
Angle_HaloPos_eigenvector_plot(cosAngleH1e1,cosAngleH2e1, "e1", SubHalo1, SubHalo2)
Angle_HaloPos_eigenvector_plot(cosAngleH1e2,cosAngleH2e2, "e2", SubHalo1, SubHalo2)
Angle_HaloPos_eigenvector_plot(cosAngleH1randomVec,cosAngleH2randomVec, "random", SubHalo1, SubHalo2)
Halo_Pos_L_plot(posHaloH1, L_haloH1, posHaloH2, L_haloH2)
Halo_K_L_plot(posHaloH1, L_haloH1, posHaloH2, L_haloH2)
cosangle_prob_plot(cosAngleH1e3,cosAngleH2e3, "e3" , SubHalo1, SubHalo2)
cosangle_prob_plot(cosAngleH1e2,cosAngleH2e2, "e2" , SubHalo1, SubHalo2)
cosangle_prob_plot(cosAngleH1e1,cosAngleH2e1, "e1" , SubHalo1, SubHalo2)
cosangle_prob_plot(cosAngleH1randomVec,cosAngleH2randomVec, "random" , SubHalo1, SubHalo2)

#pdb.set_trace()




