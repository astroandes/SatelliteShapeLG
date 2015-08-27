import numpy as np
import scipy as sp
import scipy.stats
import pdb
import matplotlib.pyplot as plt		# paquete para graficar
import matplotlib.cm as cm


GroupNum ="36"
PairNumber = 36 


##plotting routines
def Positionsplot(x_H1, y_H1, z_H1, x_satH1, y_satH1, z_satH1, x_H2, y_H2, z_H2, x_satH2, y_satH2, z_satH2, x_DMsatH1, y_DMsatH1, z_DMsatH1, x_DMsatH2, y_DMsatH2, z_DMsatH2):
	plt.figure()
	plt.plot(x_H1, y_H1,'.', markersize=15,color="blue")
	plt.plot(x_DMsatH1, y_DMsatH1,'.', markersize=5,color="blue")
	plt.plot(x_satH1, y_satH1,'*', markersize=5,color="blue")
	plt.plot(x_H2, y_H2,'.', markersize=10,color="red")
	plt.plot(x_DMsatH2, y_DMsatH2,'.', markersize=5,color="red")
	plt.plot(x_satH2, y_satH2,'*', markersize=5,color="red")
	plt.xlabel('X (Kpc)', fontsize=15) 
	plt.ylabel('Y (Kpc)', fontsize=15)
	plt.savefig('HalosXY'+GroupNum+'.pdf')
	
	plt.figure()
	plt.plot(x_H1, z_H1,'.', markersize=15,color="blue")
	plt.plot(x_DMsatH1, z_DMsatH1,'.', markersize=5,color="blue")
	plt.plot(x_satH1, z_satH1,'*', markersize=5,color="blue")
	plt.plot(x_H2, z_H2,'.', markersize=10,color="red")
	plt.plot(x_DMsatH2, z_DMsatH2,'.', markersize=5,color="red")
	plt.plot(x_satH2, z_satH2,'*', markersize=5,color="red")
	plt.xlabel('X (Kpc)', fontsize=15) 
	plt.ylabel('Z (Kpc)', fontsize=15)
	plt.savefig('HalosXZ'+GroupNum+'.pdf')
	
	plt.figure()
	plt.plot(y_H1, z_H1,'.', markersize=15,color="blue")
	plt.plot(y_DMsatH1, z_DMsatH1,'.', markersize=5,color="blue")
	plt.plot(y_satH1, z_satH1,'*', markersize=5,color="blue")
	plt.plot(y_H2, z_H2,'.', markersize=10,color="red")
	plt.plot(y_DMsatH2, z_DMsatH2,'.', markersize=5,color="red")
	plt.plot(y_satH2, z_satH2,'*', markersize=5,color="red")
	plt.xlabel('Y (Kpc)', fontsize=15) 
	plt.ylabel('Z (Kpc)', fontsize=15)
	plt.savefig('HalosYZ'+GroupNum+'.pdf')
	

def Inertiaplots(All_AxisRatioH1, All_AxisRatioH1DM,All_AxisRatioH2, All_AxisRatioH2DM):
	
	plt.figure()
	x = np.linspace(0, 20, 1000)  # 100 evenly-spaced values from 0 to 50
	y = x

	plt.plot(All_AxisRatioH1DM, All_AxisRatioH1,'.', markersize=5,color="blue")
	plt.plot(All_AxisRatioH2DM, All_AxisRatioH2,'.', markersize=5,color="red")
	plt.plot(x, y, color="black")
	plt.ylabel('Axis ratio', fontsize=15)
	plt.ylim([0,1])
	plt.xlim([0,1])
	plt.xlabel('Axis ratio DM', fontsize=15)
	plt.ylabel('Axis ratio Luminous', fontsize=15)
	plt.savefig('AxisRatio_All.pdf')


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
	plt.savefig('Halo_pos_L'+GroupNum+'.pdf')


# Juan Nicolas me mostro lo siguiente (mucho mas eficiente que como yo leia los datos!)
def loading_snapshot(snap_name):
	data = np.loadtxt(snap_name)
        x = data[:,1]
        y = data[:,2]
        z = data[:,3]
        vx = data[:,4]
        vy = data[:,5]
        vz = data[:,6]
        vmax = data[:,7]
        Mag  = data[:,8]
        return x, y, z, vx, vy, vz, vmax, Mag

def MainHalos(x, y, z, vx, vy, vz, vmax, Mag):
        max_vmax = np.amax(vmax)
        indexH = np.where(vmax >= max_vmax)
        index = np.where(vmax < max_vmax)
	print "indexH", indexH #[0]
        index = index[0]
	return x[index], y[index], z[index], vx[index], vy[index], vz[index], vmax[index], Mag[index], x[indexH], y[indexH], z[indexH], vx[indexH], vy[indexH], vz[indexH], vmax[indexH], Mag[indexH]  
	
def AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, dist_a, dist_b, DistTreshold, MagTresholdUPPER, MagTresholdLOWER):
        index = np.where((dist_a < dist_b)*(dist_a < DistTreshold)*(Mag<MagTresholdUPPER)*(Mag>MagTresholdLOWER)) 
        index = index[0]
        return x[index], y[index], z[index], vx[index], vy[index], vz[index], vmax[index], Mag[index] 

def inertiaTensor(x, y, z):
	I=[]
	for index in range(9):
        	I.append(0)

	Ixx= np.sum(y*y+z*z) 
	Ixy= np.sum(-y*x)
	Ixz= np.sum(-x*z)
	Iyx= np.sum(-y*x)
	Iyy= np.sum(x*x+z*z)
	Iyz= np.sum(-y*z)
	Izx= np.sum(-z*x)
	Izy= np.sum(-z*y)
	Izz= np.sum(x*x+y*y)

	I[0] = Ixx
        I[1] = Ixy
        I[2] = Ixz
        I[3] = Iyx
        I[4] = Iyy
        I[5] = Iyz
        I[6] = Izx
        I[7] = Izy
        I[8] = Izz

	tensor = np.array([(I[0:3]), (I[3:6]), (I[6:9])])
	vals, vects = np.linalg.eig(tensor)  # they come out unsorted, so the command below is needed
	eig_ord = np.argsort(vals)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.
	ord_vals = vals[eig_ord]
	ord_vects = vects[:, eig_ord].T

	TriaxParam = (ord_vals[2]*ord_vals[2]-ord_vals[1]*ord_vals[1])/(ord_vals[2]*ord_vals[2]-ord_vals[0]*ord_vals[0])
	AxisRatio = ord_vals[0]/ord_vals[2]
	print AxisRatio
	return ord_vals, ord_vects, TriaxParam, AxisRatio


#MAIN

All_TriaxParamH1  =[] 
All_AxisRatioH1   =[] 
All_TriaxParamH1DM=[]  
All_AxisRatioH1DM =[]  
All_TriaxParamH2  =[] 
All_AxisRatioH2   =[] 
All_TriaxParamH2DM=[]   
All_AxisRatioH2DM =[]  

for i in range (53):
	#print "Grupo =", i
	snap_name ="../data/Illustris_group_"+str(i)+".dat"
	x, y, z, vx, vy, vz, vmax, Mag = loading_snapshot(snap_name)
	#Finds the two main Halos
	print "MAIN HALO 1"
	x_aux, y_aux, z_aux, vx_aux, vy_aux, vz_aux, vmax_aux, Mag_aux, x_H1, y_H1, z_H1, vx_H1, vy_H1, vz_H1, vmax_H1, Mag_H1 = MainHalos(x, y, z, vx, vy, vz, vmax, Mag)
	print "MAIN HALO 2"
	x_aux, y_aux, z_aux, vx_aux, vy_aux, vz_aux, vmax_aux, Mag_aux, x_H2, y_H2, z_H2, vx_H2, vy_H2, vz_H2, vmax_H2, Mag_H2 = MainHalos(x_aux, y_aux, z_aux, vx_aux, vy_aux, vz_aux, vmax_aux, Mag_aux)
	#distance from all the subhalos to the Main halos	
	distH1 = np.sqrt((x-x_H1)*(x-x_H1) + (y-y_H1)*(y-y_H1) + (z-z_H1)*(z-z_H1))
	distH2 = np.sqrt((x-x_H2)*(x-x_H2) + (y-y_H2)*(y-y_H2) + (z-z_H2)*(z-z_H2))
	
	#if distH1 < distH2,  distH1 < distTreshold, luminous/not luminous, :
	x_satH1, y_satH1, z_satH1, vx_satH1, vy_satH1, vz_satH1, vmax_satH1, Mag_satH1 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH1, distH2, 200. , -6., -1000000.)
	x_DMsatH1, y_DMsatH1, z_DMsatH1, vx_DMsatH1, vy_DMsatH1, vz_DMsatH1, vmax_DMsatH1, Mag_DMsatH1 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH1, distH2, 500. , 10e45, 10000.)
	x_satH1, y_satH1, z_satH1 = x_satH1 - x_H1, y_satH1 - y_H1, z_satH1 - z_H1
	x_DMsatH1, y_DMsatH1, z_DMsatH1 = x_DMsatH1 - x_H1, y_DMsatH1 - y_H1, z_DMsatH1 - z_H1
	
	#if distH2 < distH1,  distH1 < distTreshold, luminous/not luminous, :
	x_satH2, y_satH2, z_satH2, vx_satH2, vy_satH2, vz_satH2, vmax_satH2, Mag_satH2 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH2, distH1, 200. , -6., -1000000.)
	x_DMsatH2, y_DMsatH2, z_DMsatH2, vx_DMsatH2, vy_DMsatH2, vz_DMsatH2, vmax_DMsatH2, Mag_DMsatH2 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH2, distH1, 500. , 10e45, 10000.)
	x_satH2, y_satH2, z_satH2 = x_satH2 - x_H2, y_satH2 - y_H2, z_satH2 - z_H2
	x_DMsatH2, y_DMsatH2, z_DMsatH2 = x_DMsatH2 - x_H2, y_DMsatH2 - y_H2, z_DMsatH2 - z_H2
	x_H1, y_H1, z_H1, x_H2, y_H2, z_H2 = 0.,0.,0.,0.,0.,0.
	#print x_satH1	
	#print x_satH2	
	#Positionsplot(x_H1, y_H1, z_H1, x_satH1, y_satH1, z_satH1, x_H2, y_H2, z_H2, x_satH2, y_satH2, z_satH2, x_DMsatH1, y_DMsatH1, z_DMsatH1, x_DMsatH2, y_DMsatH2, z_DMsatH2)
	
	eiValH1, eiVecH1 , TriaxParamH1, AxisRatioH1 = inertiaTensor(x_satH1, y_satH1, z_satH1)
	eiValH1DM, eiVecH1DM , TriaxParamH1DM, AxisRatioH1DM = inertiaTensor(x_DMsatH1, y_DMsatH1, z_DMsatH1)
	eiValH2, eiVecH2 , TriaxParamH2, AxisRatioH2 = inertiaTensor(x_satH2, y_satH2, z_satH2)
	eiValH2DM, eiVecH2DM , TriaxParamH2DM, AxisRatioH2DM = inertiaTensor(x_DMsatH2, y_DMsatH2, z_DMsatH2)
	
	#print "H1  ",TriaxParamH1, AxisRatioH1 
	#print "H1DM",TriaxParamH1DM, AxisRatioH1DM 
	#print "H2  ",TriaxParamH2, AxisRatioH2 
	#print "H2DM",TriaxParamH2DM, AxisRatioH2DM 
	

	All_TriaxParamH1.append(TriaxParamH1) 
	All_AxisRatioH1.append(AxisRatioH1)
	#All_eigenvalcH1.append(eiValH1[0])
	#All_eigenvalaH1.append(eiValH1[2])
	All_TriaxParamH1DM.append(TriaxParamH1DM)
	All_AxisRatioH1DM.append(AxisRatioH1DM)
	#All_eigenvalcH1DM.append(eiValH1DM[0])
	#All_eigenvalaH1DM.append(eiValH1DM[2])
	All_TriaxParamH2.append(TriaxParamH2)
	All_AxisRatioH2.append(AxisRatioH2 )  
	All_TriaxParamH2DM.append(TriaxParamH2DM)
	All_AxisRatioH2DM.append(AxisRatioH2DM ) 



#print "_AxisRatioH1", All_AxisRatioH1
#print "_AxisRatioH1DM", All_AxisRatioH1DM
Inertiaplots(All_AxisRatioH1, All_AxisRatioH1DM, All_AxisRatioH2, All_AxisRatioH2DM)



###	#looking for structure in phase space
###	velH1_X=halo[j][4]-halo[H1_line_Num][4]
###	velH1_Y=halo[j][5]-halo[H1_line_Num][5]
###	velH1_Z=halo[j][6]-halo[H1_line_Num][6]
###	HaloVX_H1.append(velH1_X)
###	HaloVY_H1.append(velH1_Y)
###	HaloVZ_H1.append(velH1_Z)
###	velHalo = np.array([ velH1_X, velH1_Y, velH1_Z])
###	# Angular momentum per unit mass
###	L_halo = np.cross(posHalo, velHalo)	
###	L_halo_modulus = np.sqrt((L_halo*L_halo).sum())
###	HaloVY_H1.append(velH1_Y)
###	posHaloH1.append(posHalo_modulus)
###	L_haloH1.append(L_halo_modulus)
###	# kinetic energy per unit mass
###	VelVel = ((velHalo*velHalo).sum())	
###	K_halo = 0.5*VelVel
###	K_haloH1.append(K_halo)



##################################      
##################################      			if halo[j][17] == halo[1][0]: 
##################################      				posH2_X=halo[j][1]-halo[1][1]
##################################      				posH2_Y=halo[j][2]-halo[1][2]
##################################      				posH2_Z=halo[j][3]-halo[1][3]
##################################      				HaloX_H2.append(posH2_X)
##################################      				HaloY_H2.append(posH2_Y)
##################################      				HaloZ_H2.append(posH2_Z)
##################################      				posHalo = np.array([ posH2_X, posH2_Y, posH2_Z])
##################################      				posHalo_modulus = np.sqrt((posHalo*posHalo).sum())
##################################      




