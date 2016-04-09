import numpy as np
import scipy as sp
import scipy.stats
import pdb
import matplotlib.pyplot as plt		# paquete para graficar
import matplotlib.cm as cm
from scipy.interpolate import interp1d

def loading_snapshot(snap_name,i):
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

#finds the 2 main halos 
def MainHalos(x, y, z, vx, vy, vz, vmax, Mag):
        max_vmax = np.amax(vmax)
        max_vmax = np.amax(vmax)
        indexH = np.where(vmax >= max_vmax)
        index = np.where(vmax < max_vmax)
        index = index[0]
	return x[index], y[index], z[index], vx[index], vy[index], vz[index], vmax[index], Mag[index], x[indexH], y[indexH], z[indexH], vx[indexH], vy[indexH], vz[indexH], vmax[indexH], Mag[indexH]  

# determines which subhalos "belong" to which main halo
def AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, dist_a, dist_b, DistTreshold, MagTresholdUPPER, MagTresholdLOWER):
        index = np.where((dist_a < dist_b)*(dist_a < DistTreshold)*(Mag<MagTresholdUPPER)*(Mag>MagTresholdLOWER)) 
        index = index[0]
	return x[index], y[index], z[index], vx[index], vy[index], vz[index], vmax[index], Mag[index] 

#selects subhalos above a certain vmax
def vmaxCut(x, y, z, vx, vy, vz, vmax, Mag, dist_a, dist_b, DistTreshold, MagTresholdUPPER, MagTresholdLOWER, vmaxTresholdLOWER):
        index = np.where((dist_a < dist_b)*(dist_a < DistTreshold)*(Mag<MagTresholdUPPER)*(Mag>MagTresholdLOWER)*(vmax>vmaxTresholdLOWER)) 
        index = index[0]
	return x[index], y[index], z[index], vx[index], vy[index], vz[index], vmax[index], Mag[index] 

def inertiaTensor(x, y, z):
	I=[]
	for index in range(9):
        	I.append(0)
	I[0] = np.sum(y*y+z*z) 
	I[1] = np.sum(-y*x)    
	I[2] = np.sum(-x*z)    
	I[3] = np.sum(-y*x)    
	I[4] = np.sum(x*x+z*z) 
	I[5] = np.sum(-y*z)    
	I[6] = np.sum(-z*x)    
	I[7] = np.sum(-z*y)    
	I[8] = np.sum(x*x+y*y) 

	tensor = np.array([(I[0:3]), (I[3:6]), (I[6:9])])
	vals, vects = np.linalg.eig(tensor)  # they come out unsorted, so the command below is needed
	eig_ord = np.argsort(vals)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.
	ord_vals = vals[eig_ord]
	ord_vects = vects[:, eig_ord].T

	TriaxParam = (ord_vals[2]*ord_vals[2]-ord_vals[1]*ord_vals[1])/(ord_vals[2]*ord_vals[2]-ord_vals[0]*ord_vals[0])
	AxisRatio = ord_vals[0]/ord_vals[2]
	return ord_vals, ord_vects, TriaxParam, AxisRatio


def Alignment_H1H2_haloShape(eiVecH1, eiVecH2, eiVecH1DM, eiVecH2DM, HaloVec, countH1, countH2, countH1H2):
		#All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2, All_cos_angleH1DM_H1H2, All_cos_angleH2DM_H1H2, All_cos_angleH1DMH2DM = Alignment_H1H2_haloShape(eiVecH1, eiVecH2, eiVecH1DM, eiVecH2DM, HaloVec)
	vecH1H2 = HaloVec 
	print "vecH1H2", vecH1H2
	print "eiVecH1", eiVecH1 
	#vecH1H2_modulus = np.sqrt((HaloVec*HaloVec).sum())  
	vecH1H2_modulus = np.sqrt(vecH1H2[0]*vecH1H2[0]+vecH1H2[1]*vecH1H2[1]+vecH1H2[2]*vecH1H2[2])  
	SemimajorAxisH1=eiVecH1[0] 
	SemimajorAxisH2=eiVecH2[0]
	SemimajorAxisH1DM=eiVecH1DM[0] 
	SemimajorAxisH2DM=eiVecH2DM[0]
	SemimajorAxisH1_modulus = np.sqrt((SemimajorAxisH1*SemimajorAxisH1).sum())
	SemimajorAxisH2_modulus = np.sqrt((SemimajorAxisH2*SemimajorAxisH2).sum())
	SemimajorAxisH1DM_modulus = np.sqrt((SemimajorAxisH1DM*SemimajorAxisH1DM).sum())
	SemimajorAxisH2DM_modulus = np.sqrt((SemimajorAxisH2DM*SemimajorAxisH2DM).sum())

	###e3_modulus = np.sqrt((e3*e3).sum())
	dotH1_H1H2 = np.dot(SemimajorAxisH1, vecH1H2)
	cos_angleH1_H1H2 = np.absolute (dotH1_H1H2 / vecH1H2_modulus / SemimajorAxisH1_modulus) # cosine of angle between poshalo and e3
	###cosAngleH1_H1H2.append(cos_angleH1e3)
	dotH2_H1H2 = np.dot(SemimajorAxisH2, vecH1H2)
	cos_angleH2_H1H2 = np.absolute (dotH2_H1H2 / vecH1H2_modulus / SemimajorAxisH2_modulus) # cosine of angle between poshalo and e3
	###cosAngleH1H2.append(cos_angleH1e3)
	dotH1H2 = np.dot(SemimajorAxisH1, SemimajorAxisH2)
	cos_angleH1H2 = np.absolute (dotH1H2 / SemimajorAxisH1_modulus / SemimajorAxisH2_modulus) # cosine of angle between poshalo and e3

	###e3_modulus = np.sqrt((e3*e3).sum())
	dotH1DM_H1H2 = np.dot(SemimajorAxisH1DM, vecH1H2)
	cos_angleH1DM_H1H2 = np.absolute (dotH1DM_H1H2 / vecH1H2_modulus / SemimajorAxisH1DM_modulus) # cosine of angle between poshalo and e3
	###cosAngleH1_H1H2.append(cos_angleH1e3)
	dotH2DM_H1H2 = np.dot(SemimajorAxisH2DM, vecH1H2)
	cos_angleH2DM_H1H2 = np.absolute (dotH2DM_H1H2 / vecH1H2_modulus / SemimajorAxisH2DM_modulus) # cosine of angle between poshalo and e3
	###cosAngleH1H2.append(cos_angleH1e3)
	dotH1H2DM = np.dot(SemimajorAxisH1DM, SemimajorAxisH2DM)
	cos_angleH1DMH2DM = np.absolute (dotH1H2DM / SemimajorAxisH1DM_modulus / SemimajorAxisH2DM_modulus) # cosine of angle between poshalo and e3

	if num_x_satH1 > SatNumCutoff_alignment :
		All_cos_angleH1_H1H2.append(cos_angleH1_H1H2)
		All_cos_angleH1DM_H1H2.append(cos_angleH1DM_H1H2)
		countH1 = countH1 + 1
	if num_x_satH2 > SatNumCutoff_alignment :
		All_cos_angleH2_H1H2.append(cos_angleH2_H1H2)
		All_cos_angleH2DM_H1H2.append(cos_angleH2DM_H1H2)
		countH2 = countH2 + 1
	if num_x_satH1 > SatNumCutoff_alignment :
		if num_x_satH2 > SatNumCutoff_alignment :
			All_cos_angleH1H2.append(cos_angleH1H2)
			All_cos_angleH1DMH2DM.append(cos_angleH1DMH2DM)
		countH1H2 = countH1H2 + 1
	return All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2, All_cos_angleH1DM_H1H2, All_cos_angleH2DM_H1H2, All_cos_angleH1DMH2DM, countH1, countH2, countH1H2 	

def Alignment_DM_luminous(eiVecH, eiVecHDM, countH, countHDM):
		#All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2, All_cos_angleH1DM_H1H2, All_cos_angleH2DM_H1H2, All_cos_angleH1DMH2DM = Alignment_H1H2_haloShape(eiVecH1, eiVecH2, eiVecH1DM, eiVecH2DM, HaloVec)
	#vecH1H2 = HaloVec 
	#print "vecH1H2", vecH1H2
	#print "eiVecH1", eiVecH1 
	#vecH1H2_modulus = np.sqrt((HaloVec*HaloVec).sum())  
	#vecH1H2_modulus = np.sqrt(vecH1H2[0]*vecH1H2[0]+vecH1H2[1]*vecH1H2[1]+vecH1H2[2]*vecH1H2[2])  
	SemimajorAxisH1=eiVecH[0] 
	SemimajorAxisHDM=eiVecHDM[0]
	#SemimajorAxisH1DM=eiVecH1DM[0] 
	#SemimajorAxisH2DM=eiVecH2DM[0]
	SemimajorAxisH_modulus = np.sqrt((SemimajorAxisH*SemimajorAxisH).sum())
	#SemimajorAxisH2_modulus = np.sqrt((SemimajorAxisH2*SemimajorAxisH2).sum())
	SemimajorAxisHDM_modulus = np.sqrt((SemimajorAxisHDM*SemimajorAxisHDM).sum())
	#SemimajorAxisH2DM_modulus = np.sqrt((SemimajorAxisH2DM*SemimajorAxisH2DM).sum())

	###cosAngleH1H2.append(cos_angleH1e3)
	dotHHDM = np.dot(SemimajorAxisH, SemimajorAxisHDM)
	cos_angleHHDM = np.absolute (dotHHDM / SemimajorAxisH_modulus / SemimajorAxisHDM_modulus) # cosine of angle between poshalo and e3

	if num_x_satH1 > SatNumCutoff_alignment :
		All_cos_angleH1_H1H2.append(cos_angleH1_H1H2)
		All_cos_angleH1DM_H1H2.append(cos_angleH1DM_H1H2)
		countH1 = countH1 + 1
	if num_x_satH2 > SatNumCutoff_alignment :
		All_cos_angleH2_H1H2.append(cos_angleH2_H1H2)
		All_cos_angleH2DM_H1H2.append(cos_angleH2DM_H1H2)
		countH2 = countH2 + 1
	if num_x_satH1 > SatNumCutoff_alignment :
		if num_x_satH2 > SatNumCutoff_alignment :
			All_cos_angleH1H2.append(cos_angleH1H2)
			All_cos_angleH1DMH2DM.append(cos_angleH1DMH2DM)
		countH1H2 = countH1H2 + 1
	return All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2, All_cos_angleH1DM_H1H2, All_cos_angleH2DM_H1H2, All_cos_angleH1DMH2DM, countH1, countH2, countH1H2 	

def Angle_HaloPos_eigenvector_plot(All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2):
	plt.figure()
	fig = plt.figure(figsize=(6, 6))
	ax = fig.add_subplot(1,1,1) # one row, one column, first plot
	#plt.hist(cosAngleH1,  bins =25, color ="blue", alpha =0.7 )
	#plt.hist(cosAngleH2,  bins =25, color ="red", alpha =0.7 )
	ax.hist((All_cos_angleH1_H1H2, All_cos_angleH2_H1H2),  bins =5, color =("red", "blue"))
	plt.xlabel('cos angle pos ', fontsize=15) 
	plt.ylabel('N Halos', fontsize=15)
	plt.savefig('Hist_H1H2_shape.pdf')

	plt.figure()
	fig = plt.figure(figsize=(6, 6))
	#ax = fig.add_subplot(1,1,1) # one row, one column, first plot
	#plt.hist(cosAngleH1,  bins =25, color ="blue", alpha =0.7 )
	#plt.hist(cosAngleH2,  bins =25, color ="red", alpha =0.7 )
	plt.hist(All_cos_angleH1H2,  bins =10, color ="black")
	plt.xlabel('cos angle pos ', fontsize=15) 
	plt.ylabel('N Halos', fontsize=15)
	plt.savefig('Hist_shape.pdf')

def cosangle_prob_plot(p_All_cos_angleH1_H1H2, p_All_cos_angleH2_H1H2, p_All_cos_angleH1H2, SubHalo1, SubHalo2, name, numSat):
	plt.figure()
	
	x = np.linspace(0, 20, 1000)  # 100 evenly-spaced values from 0 to 50
	y = x

	a1,b1,c1,d1 = get_randoms(n_points=countH1)
	a2,b2,c2,d2 = get_randoms(n_points=countH2)
	a12,b12,c12,d12 = get_randoms(n_points=countH2)

	countsH1, startH1, dxH1, _ = scipy.stats.cumfreq(p_All_cos_angleH1_H1H2, numbins=50)
	countsH2, startH2, dxH2, _ = scipy.stats.cumfreq(p_All_cos_angleH2_H1H2, numbins=50)
	countsH1H2, startH1H2, dxH1H2, _ = scipy.stats.cumfreq(p_All_cos_angleH1H2, numbins=50)
	xH1 = np.arange(countsH1.size) * dxH1 + startH1
	xH2 = np.arange(countsH2.size) * dxH2 + startH2
	xH1H2 = np.arange(countsH1H2.size) * dxH1H2 + startH1H2
	countsH1 = countsH1/len(p_All_cos_angleH1_H1H2)
	countsH2 = countsH2/len(p_All_cos_angleH2_H1H2) 
	countsH1H2 = countsH1H2/len(p_All_cos_angleH1H2) 
	
	plt.plot(xH1, countsH1, linewidth=3, color = "red")
	plt.plot(xH2, countsH2, linewidth=3, color = "blue")
	plt.plot(xH1H2, countsH1H2, linewidth=3, color = "green")
	plt.plot(a1, b1, linewidth=3, color = "red")
	plt.plot(a1, c1, linewidth=3, color = "red")
	plt.plot(a1, d1, linewidth=3, color = "red")
	plt.plot(a2, b2, linewidth=3, color = "blue")
	plt.plot(a2, c2, linewidth=3, color = "blue")
	plt.plot(a2, d2, linewidth=3, color = "blue")
	plt.plot(a12, b12, linewidth=3, color = "green")
	plt.plot(a12, c12, linewidth=3, color = "green")
	plt.plot(a12, d12, linewidth=3, color = "green")
	plt.plot(x, y, color="black")
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('Cos theta ')
	plt.ylabel('Cumulative Frequency')
	plt.savefig('freq_pos_dot_'+str(name)+'_'+str(numSat)+'.pdf')

def get_randoms(n_points=10):
	n_random = 5000
	n_points_x = 80
	y_axis = np.linspace(0.0,1.0,n_points)
	x_axis_dense = np.linspace(0.0,1.0,n_points_x)
	y_axis_dense_all = np.zeros([n_random,n_points_x])
	for i in range(n_random):
	    ran_range = np.random.random(n_points)
	    ran_range = np.sort(ran_range)
	    ran_range[0] = 0.0
	    ran_range[-1] = 1.0
	    f = interp1d(ran_range, y_axis)
	    y_axis_dense = f(x_axis_dense)
	    y_axis_dense_all[i,:] = y_axis_dense[:]
	
	y_axis_median = np.ones(n_points_x)
	y_axis_low = np.ones(n_points_x)
	y_axis_high = np.ones(n_points_x)
	for i in range(n_points_x):
	    y_axis_median[i] = np.median(y_axis_dense_all[:,i])
	    y_axis_low[i] = np.percentile(y_axis_dense_all[:,i],5)
	    y_axis_high[i] = np.percentile(y_axis_dense_all[:,i],95)
	
	return x_axis_dense, y_axis_median, y_axis_low, y_axis_high


def Inertiaplots_VS_random(All_AxisRatioH1, All_AxisRatioH1All ,All_AxisRatioH2, All_AxisRatioH2All, All_AxisRatioH1random, All_AxisRatioH2random, All_AxisRatioH1random_mean,  All_AxisRatioH1random_mean_err, All_AxisRatioH2random_mean, All_AxisRatioH2random_mean_err):

	MW_AxisRatio = MW_data_AxisRatio()

	plt.figure()
	x = np.linspace(0, 20, 1000)  # 100 evenly-spaced values from 0 to 50
	y = x
	
	#x = MW_AxisRatio
	#point size depending on the ratio (number of sats)/(number of dark subhalos)
	#plt.scatter(All_AxisRatioH1random_mean, All_AxisRatioH1, marker='o', c='r', s=SatRatio_H1)
	#plt.scatter(All_AxisRatioH2random_mean, All_AxisRatioH2, marker='o', c='b', s=SatRatio_H2)

	#point size depending on the number of sats
	plt.scatter(All_AxisRatioH1random_mean, All_AxisRatioH1, marker='o', c='b', s=All_SatNumH1)
	plt.scatter(All_AxisRatioH2random_mean, All_AxisRatioH2, marker='o', c='r', s=All_SatNumH2)
	#plt.scatter(MW_AxisRatio,MW_AxisRatio, marker='o', c='grey', s=10*All_SatNumH2)
	plt.axhline(y=MW_AxisRatio, xmin=0., xmax=1., linewidth=2, color = "grey")

	plt.errorbar(All_AxisRatioH1random_mean,All_AxisRatioH1, xerr=All_AxisRatioH1random_mean_err, yerr = 0. , color="blue", fmt='.')
	plt.errorbar(All_AxisRatioH2random_mean,All_AxisRatioH2, xerr=All_AxisRatioH2random_mean_err, yerr = 0. , color="red", fmt='.')

	plt.plot(x, y, color="black")
	plt.ylabel('Axis ratio', fontsize=15)
	plt.ylim([0,1])
	plt.xlim([0,1])
	plt.xlabel('Mean axis ratio of subsets', fontsize=15)
	plt.ylabel('Axis ratio Luminous', fontsize=15)
	plt.savefig('AxisRatio_RandomMean_VS_Luminous.pdf')

	plt.figure()
	x = np.linspace(0, 20, 1000)  # 100 evenly-spaced values from 0 to 50
	y = x
	plt.plot(All_AxisRatioH1All, All_AxisRatioH1random_mean,'.', markersize=15,color="blue")
	plt.plot(All_AxisRatioH2All, All_AxisRatioH2random_mean,'.', markersize=15,color="red")
	plt.plot(x, y, color="black")
	plt.ylabel('Axis ratio', fontsize=15)
	plt.ylim([0,1])
	plt.xlim([0,1])
	plt.xlabel('Axis ratio all subhalos', fontsize=15)
	plt.ylabel('Axis ratio Luminous', fontsize=15)
	plt.savefig('AxisRatio_All_VS_Luminous_VS_RandomMean.pdf')

def plots_Mass_VS_Mag(vmax_satH1, Mag_satH1, vmax_DMsatH1, Mag_DMsatH1, vmax_satH2, Mag_satH2, vmax_DMsatH2, Mag_DMsatH2,i, vmax_DMsatH1_vmax, Mag_DMsatH1_vmax ):
	plt.figure()
	#print "vmax, Mag", vmax, Mag
	plt.plot(vmax_satH1, Mag_satH1, marker='*', c='b')
	plt.plot(vmax_satH2, Mag_satH2, marker='o', c='b')
	plt.plot(vmax_DMsatH1, np.random.uniform(low=2.001, high=2.003, size=len(Mag_DMsatH1)), marker='*', c='green')
	plt.plot(vmax_DMsatH1_vmax, np.random.uniform(low=1.5001, high=1.5001, size=len(Mag_DMsatH1_vmax)), marker='*', c='red')
	plt.plot(vmax_DMsatH2, np.random.uniform(low=1.00, high=1.001, size=len(Mag_DMsatH2)), marker='o', c='green')
	plt.plot(vmax_DMsatH2_vmax, np.random.uniform(low=0.500, high=0.5001, size=len(Mag_DMsatH2_vmax)), marker='o', c='red')
	plt.xlabel('Mass', fontsize=15)
	plt.ylabel('Magnitude', fontsize=15)
	#plt.xlim([0,1])
	plt.ylim([-25,2.5])
	plt.savefig('Magnitude_vs_Mass_'+str(i)+'.pdf')


def MW_data_AxisRatio():
	x_MW_sats = np.array([17.1, -0.6, 16.5, -4.3, -22.2, -5.2, -36.7, -25.0, -41.3, -77.3, -123.6])
	y_MW_sats = np.array([2.5, -41.8, -38.5, 62.2, 52.0, -9.8, -56.9, -95.9, -51.10, -58.3, -119.3])
	z_MW_sats = np.array([-6.4, -27.5, -44.7, -43.2, 53.5, -85.3, 57.8, -39.8, -134.91, 215.2, 191.7])
	eiValMW, eiVecMW, TriaxParamMW, AxisRatioMW = inertiaTensor(x_MW_sats, y_MW_sats, z_MW_sats)

	print AxisRatioMW
	return AxisRatioMW



#MAIN


Npairs = 53 #53 #dm 53
count=0
N=100
#*ratioAxisRatioH1=[]
#*ratioAxisRatioH2=[]
All_AxisRatioH1   =[] 
All_AxisRatioH1DM =[] 
All_AxisRatioH2   =[] 
All_AxisRatioH2DM =[] 

All_AxisRatioH1All =[] 
All_AxisRatioH2All =[] 

All_AxisRatioH1random =[] 
All_AxisRatioH2random =[] 

All_AxisRatioH1random_mean =[] 
All_AxisRatioH1random_mean_err =[] 
All_AxisRatioH2random_mean =[]
All_AxisRatioH2random_mean_err =[]

All_SatNumH1 =[] 
All_SatNumH2 =[]

All_cos_angleH1_H1H2 =[]
All_cos_angleH2_H1H2 =[]
All_cos_angleH1H2 =[]
All_cos_angleH1DM_H1H2 =[]
All_cos_angleH2DM_H1H2 =[]
All_cos_angleH1DMH2DM =[]

countHaloPairs =0
#number of satelites cutoff for the distribution flatness comparison
SatNumCutoff = 3 
#number of satelites cutoff for the distribution and alignment between H1 and H2 luminous or dark satellites 
SatNumCutoff_alignment = 3   
countH1 = 0
countH2 = 0
countH1H2 = 0



for i in range (Npairs):
	AxisRatioH1_r=np.zeros([N])
	AxisRatioH2_r=np.zeros([N])
	snap_name ="../data/dm_selected/Illustris_group_"+str(i)+".dat"
	x, y, z, vx, vy, vz, vmax, Mag = loading_snapshot(snap_name,i)
	#print i , x

	#Finds the two main Halos
	#print "MAIN HALO 1"
	x_aux, y_aux, z_aux, vx_aux, vy_aux, vz_aux, vmax_aux, Mag_aux, x_H1, y_H1, z_H1, vx_H1, vy_H1, vz_H1, vmax_H1, Mag_H1 = MainHalos(x, y, z, vx, vy, vz, vmax, Mag)
	#print "MAIN HALO 2"
	x_aux, y_aux, z_aux, vx_aux, vy_aux, vz_aux, vmax_aux, Mag_aux, x_H2, y_H2, z_H2, vx_H2, vy_H2, vz_H2, vmax_H2, Mag_H2 = MainHalos(x_aux, y_aux, z_aux, vx_aux, vy_aux, vz_aux, vmax_aux, Mag_aux)
	#print i, x_H1, x_H2

	#distance from all the subhalos to the Main halos	
	distH1 = np.sqrt((x-x_H1)*(x-x_H1) + (y-y_H1)*(y-y_H1) + (z-z_H1)*(z-z_H1))
	distH2 = np.sqrt((x-x_H2)*(x-x_H2) + (y-y_H2)*(y-y_H2) + (z-z_H2)*(z-z_H2))
	#distance between the 2 main halos
	distH1H2 = np.sqrt((x_H1-x_H2)*(x_H1-x_H2) + (y_H1-y_H2)*(y_H1-y_H2) + (z_H1-z_H2)*(z_H1-z_H2))

	if distH1H2 >=600:
		Radius = 300. #guess of the virial radius. NOTE: could be improved using Vmax
	else:
		Radius = 0.4*distH1H2
		print "i=", i, "Radius =", Radius

	x_AllH1, y_AllH1, z_AllH1, vx_AllH1, vy_AllH1, vz_AllH1, vmax_AllH1, Mag_AllH1 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH1, distH2, DistTreshold = Radius , MagTresholdUPPER = 10e45, MagTresholdLOWER = -1000000.)
	x_satH1, y_satH1, z_satH1, vx_satH1, vy_satH1, vz_satH1, vmax_satH1, Mag_satH1 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH1, distH2, DistTreshold = Radius , MagTresholdUPPER =-6., MagTresholdLOWER = -1000000.)
	x_DMsatH1, y_DMsatH1, z_DMsatH1, vx_DMsatH1, vy_DMsatH1, vz_DMsatH1, vmax_DMsatH1, Mag_DMsatH1 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH1, distH2, DistTreshold = Radius, MagTresholdUPPER = 10e45, MagTresholdLOWER = 10000.)
	x_DMsatH1_vmax, y_DMsatH1_vmax, z_DMsatH1_vmax, vx_DMsatH1_vmax, vy_DMsatH1_vmax, vz_DMsatH1_vmax, vmax_DMsatH1_vmax, Mag_DMsatH1_vmax = vmaxCut(x, y, z, vx, vy, vz, vmax, Mag, distH1, distH2, DistTreshold = Radius, MagTresholdUPPER = 10e45, MagTresholdLOWER = 10000., vmaxTresholdLOWER = np.min(vmax_satH1))
	
	#if distH2 < distH1,  distH1 < distTreshold, luminous/not luminous, :
	x_AllH2, y_AllH2, z_AllH2, vx_AllH2, vy_AllH2, vz_AllH2, vmax_allH2, Mag_AllH2 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH2, distH1, DistTreshold = Radius , MagTresholdUPPER =10e45, MagTresholdLOWER = -1000000.)
	x_satH2, y_satH2, z_satH2, vx_satH2, vy_satH2, vz_satH2, vmax_satH2, Mag_satH2 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH2, distH1, DistTreshold = Radius , MagTresholdUPPER =-6., MagTresholdLOWER = -1000000.)
	x_DMsatH2, y_DMsatH2, z_DMsatH2, vx_DMsatH2, vy_DMsatH2, vz_DMsatH2, vmax_DMsatH2, Mag_DMsatH2 = AllCriteria(x, y, z, vx, vy, vz, vmax, Mag, distH2, distH1, DistTreshold = Radius, MagTresholdUPPER = 10e45, MagTresholdLOWER = 10000.)
	x_DMsatH2_vmax, y_DMsatH2_vmax, z_DMsatH2_vmax, vx_DMsatH2_vmax, vy_DMsatH2_vmax, vz_DMsatH2_vmax, vmax_DMsatH2_vmax, Mag_DMsatH2_vmax = vmaxCut(x, y, z, vx, vy, vz, vmax, Mag, distH2, distH1, DistTreshold = Radius, MagTresholdUPPER = 10e45, MagTresholdLOWER = 10000., vmaxTresholdLOWER = np.min(vmax_satH2))


	#Centers all the coordinates in the center of the Parent Halo (H1 or H2)
	x_AllH1, y_AllH1, z_AllH1 = x_AllH1 - x_H1, y_AllH1 - y_H1, z_AllH1 - z_H1
	x_satH1, y_satH1, z_satH1 = x_satH1 - x_H1, y_satH1 - y_H1, z_satH1 - z_H1
	x_DMsatH1, y_DMsatH1, z_DMsatH1 = x_DMsatH1 - x_H1, y_DMsatH1 - y_H1, z_DMsatH1 - z_H1
	x_AllH2, y_AllH2, z_AllH2 = x_AllH2 - x_H2, y_AllH2 - y_H2, z_AllH2 - z_H2
	x_satH2, y_satH2, z_satH2 = x_satH2 - x_H2, y_satH2 - y_H2, z_satH2 - z_H2
	x_DMsatH2, y_DMsatH2, z_DMsatH2 = x_DMsatH2 - x_H2, y_DMsatH2 - y_H2, z_DMsatH2 - z_H2
	HaloVec = (x_H1[0]- x_H2[0], y_H1[0]- y_H2[0], z_H1[0] - z_H2[0])
	#print "x, y ,z halos ", x_H1, x_H2, y_H1, y_H2, z_H1, z_H2
	#*print "HaloVec", HaloVec, HaloVec[0], HaloVec[1], HaloVec[2]
	x_H1, y_H1, z_H1, x_H2, y_H2, z_H2 = 0.,0.,0.,0.,0.,0.
	#GroupNum =np.str(i)
	#PairNumber = i
	#print i, x_AllH1, x_AllH2, x_satH2

	num_x_H1   =len(x_AllH1)
	#print "i, num_x_H1  ", i, num_x_H1
	num_x_DMsatH1   =len(x_DMsatH1)
	#print "i ,num_x_DMsatH1", i, num_x_DMsatH1
	num_x_satH1   =len(x_satH1)
	#print "i, num_x_satH1" , i, num_x_satH1
	num_x_H2   =len(x_AllH2)
	num_x_DMsatH2   =len(x_DMsatH2)
	num_x_satH2   =len(x_satH2)

	#plots to see how many dark sub-halos are more massive than the least massive luminous sub-halo
	plots_Mass_VS_Mag(vmax_satH1, Mag_satH1, vmax_DMsatH1, Mag_DMsatH1, vmax_satH2, Mag_satH2, vmax_DMsatH2, Mag_DMsatH2, i, vmax_DMsatH1_vmax, Mag_DMsatH1_vmax)

	#counts how many of the halos have at least 16 luminous satellites
	#if num_x_satH1>=16:
	#	count=count+1
		#print "i=", i , "num_x_H1=", num_x_H1, "num_x_satH1", num_x_satH1 , "count" , count

	#For each halo: from all the subhalos, I create N random subset of num_x_satH1+1 satellites, where num_x_satH1+1 is the number of luminous subhalos.
	for j in range (N):
	        index_randH1=np.random.randint(0, num_x_H1, num_x_satH1+1) #genero num_x_satH1 numeros aleatorios
	        ran_x_satH1 = x_AllH1[index_randH1]
	        ran_y_satH1 = y_AllH1[index_randH1]
	        ran_z_satH1 = z_AllH1[index_randH1]

		#calculates the axis ratio of each subset
		eiValH1_r1, eiVecH1_r1 , TriaxParamH1_r1, AxisRatioH1_r1 = inertiaTensor(ran_x_satH1, ran_y_satH1, ran_z_satH1)
		AxisRatioH1_r[j]=AxisRatioH1_r1
		#print " eiVecH1_r1[0]", eiVecH1_r1#[0][0]
		#eiVecH1_r[:j]=eiVecH1_r1[0]

	#calculates the mean axis ratio of all the subsets
	AxisRatioH1_r_mean = scipy.stats.tmean(AxisRatioH1_r, limits=None, inclusive=(True, True))
	AxisRatioH1_r_mean_err = np.var(AxisRatioH1_r, axis=None, dtype=None, out=None, ddof=0)	

	#Calculates the axis ratio of the simulated subhalos for the total number of subhalos and for the luminous ones.	
	eiValH1, eiVecH1 , TriaxParamH1, AxisRatioH1 = inertiaTensor(x_satH1, y_satH1, z_satH1)
	#print "eiVecH1", eiVecH1
	eiValH1DM, eiVecH1DM , TriaxParamH1DM, AxisRatioH1DM = inertiaTensor(x_DMsatH1, y_DMsatH1, z_DMsatH1)
	eiValH1All, eiVecH1All , TriaxParamH1All, AxisRatioH1All = inertiaTensor(x_AllH1, y_AllH1, z_AllH1)

	#*H1_ratioAxisRatioH1=AxisRatioH1/AxisRatioH1_r_mean
	#*ratioAxisRatioH1.append(H1_ratioAxisRatioH1)

	#selects halos with a minimum number of luminous satellites:
	if num_x_satH1>=SatNumCutoff:
		All_AxisRatioH1.append(AxisRatioH1)
		All_AxisRatioH1DM.append(AxisRatioH1DM)

		All_AxisRatioH1All.append(AxisRatioH1All)

		All_AxisRatioH1random.append(AxisRatioH1_r)

		All_AxisRatioH1random_mean.append(AxisRatioH1_r_mean)
		All_AxisRatioH1random_mean_err.append(AxisRatioH1_r_mean_err)
		All_SatNumH1.append(4.*num_x_satH1)

	#counts how many of the halos have at least 16 luminous satellites
	#if num_x_satH2>=16:
	#	count=count+1
		#print "i=", i , "num_x_H2=", num_x_H2, "num_x_satH2", num_x_satH2 , "count" , count

	#For each halo: from all the subhalos, I create N random subset of num_x_satH2+1 satellites, where num_x_satH2+1 is the number of luminous subhalos.
	for j in range (N):
	        index_randH2=np.random.randint(0, num_x_H2, num_x_satH2+1) #genero num_x_satH2 numeros aleatorios
	        ran_x_satH2 = x_AllH2[index_randH2]
	        ran_y_satH2 = y_AllH2[index_randH2]
	        ran_z_satH2 = z_AllH2[index_randH2]

		#calculates the axis ratio of each subset
		eiValH2_r2, eiVecH2_r2 , TriaxParamH2_r2, AxisRatioH2_r2 = inertiaTensor(ran_x_satH2, ran_y_satH2, ran_z_satH2)
		AxisRatioH2_r[j]=AxisRatioH2_r2
		####****AxisRatioH2_r[j]=AxisRatioH2_r2
		#print "EI vecr2" , eiVecH2_r2 

	#calculates the mean axis ratio of all the subsets
	AxisRatioH2_r_mean = scipy.stats.tmean(AxisRatioH2_r, limits=None, inclusive=(True, True))
	AxisRatioH2_r_mean_err = np.var(AxisRatioH2_r, axis=None, dtype=None, out=None, ddof=0)	

	#Calculates the axis ratio of the simulated subhalos for the total number of subhalos and for the luminous ones.	
	eiValH2, eiVecH2 , TriaxParamH2, AxisRatioH2 = inertiaTensor(x_satH2, y_satH2, z_satH2)
	#print "eivecH2",  eiVecH2
	#print "eivecH2[0][0]",  eiVecH2[0][0]
	#print "eivecH2[0][1]",  eiVecH2[0][1]
	#print "eivecH2[1][0]",  eiVecH2[1][0]
	eiValH2DM, eiVecH2DM , TriaxParamH2D0M, AxisRatioH2DM = inertiaTensor(x_DMsatH2, y_DMsatH2, z_DMsatH2)
	eiValH2All, eiVecH2All , TriaxParamH2All, AxisRatioH2All = inertiaTensor(x_AllH2, y_AllH2, z_AllH2)

	#*H2_ratioAxisRatioH2=AxisRatioH2/AxisRatioH2_r_mean
	#*ratioAxisRatioH2.append(H2_ratioAxisRatioH2)

	#selects halos with a minimum number of luminous satellites:
	if num_x_satH2>=SatNumCutoff:
		All_AxisRatioH2.append(AxisRatioH2 )  
		All_AxisRatioH2DM.append(AxisRatioH2DM ) 

		All_AxisRatioH2All.append(AxisRatioH2All ) 

		All_AxisRatioH2random.append(AxisRatioH2_r) 

		All_AxisRatioH2random_mean.append(AxisRatioH2_r_mean) 
		All_AxisRatioH2random_mean_err.append(AxisRatioH2_r_mean_err) 
		All_SatNumH2.append(4.*num_x_satH2) #el 4 es solo para que en el plot se vean mejor
	#print "eiVecH1", eiVecH1
	#if num_x_satH1>= SatNumCutoff_alignment & num_x_satH2>= SatNumCutoff_alignment: 	
	All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2, All_cos_angleH1DM_H1H2, All_cos_angleH2DM_H1H2, All_cos_angleH1DMH2DM, countH1, countH2, countH1H2 = Alignment_H1H2_haloShape(eiVecH1 , eiVecH2, eiVecH1DM, eiVecH2DM, HaloVec, countH1, countH2, countH1H2)
		#All_cos_angleH1DM_H1H2, All_cos_angleH2DM_H1H2, All_cos_angleH1DMH2DM = Alignment_H1H2_haloShape(eiVecH1DM, eiVecH2DM, HaloVec)
		#countHaloPairs = countHaloPairs  +1 
print "xxxxxxxxxxxxxxxx i, countHaloPairs , SatNumCutoff xxxxxxxxxxxxxx", i, countHaloPairs, SatNumCutoff 

SubHalo1, SubHalo2 =  num_x_satH1,  num_x_satH2
cosangle_prob_plot(All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2, SubHalo1, SubHalo2, name = 'luminous', numSat = str(SatNumCutoff_alignment))
cosangle_prob_plot(All_cos_angleH1DM_H1H2, All_cos_angleH2DM_H1H2, All_cos_angleH1DMH2DM, SubHalo1, SubHalo2, name = 'DM', numSat = str(SatNumCutoff_alignment))

#*Angle_HaloPos_eigenvector_plot(All_cos_angleH1_H1H2, All_cos_angleH2_H1H2, All_cos_angleH1H2)

Inertiaplots_VS_random(All_AxisRatioH1, All_AxisRatioH1All, All_AxisRatioH2, All_AxisRatioH2All, All_AxisRatioH1random, All_AxisRatioH2random, All_AxisRatioH1random_mean,  All_AxisRatioH1random_mean_err, All_AxisRatioH2random_mean, All_AxisRatioH2random_mean_err)
