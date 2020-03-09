import pandas as pd 
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

#Q1
OIS=pd.read_excel('IR Data.xlsx', sheet_name='OIS', header=0, index_col=0).Rate
#calculate f0_6m, f0_1y, Do0_6m, Do0_1y
f01=lambda x: (1/(1+x/360)**180)*OIS[0]*0.5-(1/(1+x/360)**180)*((1+x/360)**180-1)
f0_6m=brentq(f01,0,1)
Do0_6m=1/(1+f0_6m/360)**180

f02=lambda x: Do0_6m*(1/(1+x/360)**180)*OIS[1]- Do0_6m*(1/(1+x/360)**180)*((1+f0_6m/360)**180*(1+x/360)**180-1)
f6m_1y=brentq(f02,0,1)
Do0_1y=Do0_6m*(1/(1+f6m_1y/360)**180)

#get rid of 6m OIS, because the loop strat from 2y OIS, but need 1y OIS
OIS=OIS[1:]

#turn index into number
OIS.index=[int(i.strip('y')) for i in OIS.index]

def OISbootstrap(OIS):
#initialize OIS discount factor, f rate and floatingleg using 1y OIS
	dis_OIS=[Do0_1y]
	F=[f6m_1y]
	Floatleg=[Do0_1y*((1+f0_6m/360)**180*(1+f6m_1y/360)**180-1)]

	def fixed(x, OIS, i):
		dis_fixed=sum(dis_OIS)
		for j in range(1, OIS.index[i]-OIS.index[i-1]+1):
			dis_fixed+=dis_OIS[-1]*(1/(1+x/360)**360)**j
		fixedleg=(dis_fixed)*OIS.iloc[i]
		return fixedleg

	def float(x, OIS, i):
		floatleg=Floatleg[-1]
		for j in range(1, OIS.index[i]-OIS.index[i-1]+1):
			floatleg+=dis_OIS[-1]*(1/(1+x/360)**360)**j*((1+x/360)**360-1)
		return floatleg

	for i in range(1,len(OIS)):
		f03= lambda x: fixed(x, OIS, i)-float(x,OIS,i)
		f=brentq(f03,0,1)
		F.append(f)
		floatleg=Floatleg[-1]
		for j in range(1, OIS.index[i]-OIS.index[i-1]+1):
			dis=dis_OIS[OIS.index[i-1]-1]*(1/(1+f/360)**360)**j
			floatleg+=dis_OIS[OIS.index[i-1]-1]*(1/(1+f/360)**360)**j*((1+f/360)**360-1)
			dis_OIS.append(dis)
			Floatleg.append(floatleg)
	return [F, dis_OIS]

dis_OIS=OISbootstrap(OIS)[1]
#Add on Do0_0, Do0_6m 
dis_OIS=[1]+[Do0_6m]+dis_OIS

F=OISbootstrap(OIS)[0]
#Add on f0_6m
F=[f0_6m]+F

#turn OIS discount factor and Fed fund rate into DataFrame
T1=[0]+[0.5]+list(range(1,31))
OISD=pd.DataFrame(dis_OIS, index=T1, columns=['OIS discount factor'])
FF=pd.DataFrame(F, index=[0.5]+list(OIS.index), columns=['F rate'])

#list f rate for all the years, in this way, the f rate curve will display flat in each period
FF_all_halfyear=[]
for i in range(2,len(FF.index)):
	f= int((FF.index[i]-FF.index[i-1]))*2*[FF.loc[FF.index[i],'F rate']]
	FF_all_halfyear.append(f)
FF_all_halfyear=[j for sub in FF_all_halfyear for j in sub]
FF_all_halfyear=[f0_6m]+[f6m_1y]+FF_all_halfyear
FF_all_halfyear=pd.DataFrame(FF_all_halfyear, index=np.arange(0.5,30.5,0.5),columns=['F rate'])

if __name__ == "__main__":
	OISD.to_excel('OISD.xlsx')
	FF.to_excel('FF.xlsx')
	print(OISD)
	print(FF)
	#plot OIS discount curve
	OISD.plot()
	plt.xlabel('T')
	plt.ylabel('OIS discount factor')
	plt.title('OIS discount curve')
	plt.show()

	FF_all_halfyear.plot()
	plt.xlabel('T')
	plt.ylabel('Forward Fed fund o/n rate')
	plt.title('Forward Fed fund o/n rate')
	plt.show()


#Q2
IRS=pd.read_excel('IR Data.xlsx', sheet_name='IRS', header=0, index_col=0).Rate
#calculate L0_6m, L6m_1y, D0_6m, D0_1y
L0_6m=IRS[0]
f04=lambda x: L0_6m-2*((1-x)/x)
D0_6m=brentq(f04, 0.1, 1)

f05=lambda x: (dis_OIS[1]+dis_OIS[2])*IRS[1]/2- (dis_OIS[1]*L0_6m/2+dis_OIS[2]*x/2)
L6m_1y=brentq(f05, 0.00001,0.5)
f06=lambda x: L6m_1y-2*((D0_6m-x)/x)
D0_1y=brentq(f06, 0.1,1)

#get rid of L0_6m, because the loop strat from 2y IRS, but need 1y IRS
IRS=IRS[1:]
#turn index into number
IRS.index=[int(i.strip('y')) for i in IRS.index]

#calculate half year discount factors, loop starts from Do0_1y
#firstly list out all the f rate corresponding to each year from 1y-30y, there should be 29 f rate
FF_all_year=[]
for i in range(2,len(FF.index)):
	f= int((FF.index[i]-FF.index[i-1]))*[FF['F rate'][FF.index[i]]]
	FF_all_year.append(f)
FF_all_year=[j for sub in FF_all_year for j in sub]
FF_all_year=pd.DataFrame(FF_all_year, index=range(1,30),columns=['OIS discount factor'])
#Turn f rate into mutiplier coresponding to each year
mutiplier= 1/(1+FF_all_year/360)**180
#OISD2 is from Do(0,1y) to Do(0,29y), 29 discount factors
OISD2=OISD.iloc[2:31,:]
#OISD_half is all the half year discount factors from Do(0,1.5y) to Do(0,29.5y)
OISD_half=OISD2*mutiplier
#put OISD2 and OISD_half to get all OIS discount factors
OIS_all=[]
for i in range(29):
	OIS_all.append(np.array(OISD2)[i,0])
	OIS_all.append(np.array(OISD_half)[i,0])
#add on Do(0,6m) and Do(0,30)
OIS_all= [Do0_6m]+OIS_all+[dis_OIS[-1]]
OIS_all=pd.DataFrame(OIS_all, index=np.arange(0.5,30.5,0.5),columns=['OIS discount factors'])

def IRSbootstrap(IRS):
#initialize Libor discount factor, Libor rate and floatingleg using 1y IRS
	dis_IRS=[D0_1y]
	Libor=[L6m_1y]
	Floatleg=[Do0_6m*L0_6m/2+Do0_1y*L6m_1y/2]

	def fixed(IRS, i):
		fixedleg=sum(OIS_all.iloc[:IRS.index[i]*2,0])*IRS.iloc[i]/2
		return fixedleg

	def float(x, IRS, i):
		floatleg=Floatleg[-1]
		n=2*(IRS.index[i]-IRS.index[i-1])
		for j in range(1, n+1):
			floatleg+=OIS_all.iloc[IRS.index[i-1]*2+j-1,0]*((dis_IRS[-1]-x)/((n-j)*dis_IRS[-1]+j*x))
		return floatleg

	for i in range(1,len(IRS)):
		f07= lambda x: fixed(IRS, i)-float(x,IRS,i)
		dis=brentq(f07,0.1,1)
		floatleg=Floatleg[-1]
		n=2*(IRS.index[i]-IRS.index[i-1])
		for j in range(1, n+1):
			l=2*(dis_IRS[-1]-dis)/((n-j)*dis_IRS[-1]+j*dis)
			Libor.append(l)
			floatleg+=OIS_all.iloc[IRS.index[i-1]*2+j-1,0]*((dis_IRS[-1]-dis)/((n-j)*dis_IRS[-1]+j*dis))
		Floatleg.append(floatleg)
		dis_IRS.append(dis)	
	return [Libor, dis_IRS]

Libor=IRSbootstrap(IRS)[0]

#add on L0_6m
Libor=[L0_6m]+Libor

#add on D0_0 and D0_6m 
dis_IRS=IRSbootstrap(IRS)[1]
dis_IRS=[1]+[D0_6m]+dis_IRS

#interpolate IRS discount factors
T2=[0]+[0.5]+list(IRS.index)
f08=interpolate.interp1d(T2,dis_IRS)
dis_IRS=f08(T1)

#turn IRS discount factor and Libor rate into DataFrame
IRSD=pd.DataFrame(dis_IRS, index=T1, columns=['Libor discount factor'])
Libor=pd.DataFrame(Libor, index=np.arange(0.5,30.5,0.5), columns=['Forward Libor rate'])
if __name__=='__main__':
	IRSD.to_excel('IRSD.xlsx')
	Libor.to_excel('Libor.xlsx')
	print(IRSD)
	print(Libor)
	#plot IRS discount factors
	IRSD.plot()
	plt.xlabel('T')
	plt.ylabel('Libor discount factor')
	plt.title('Libor discount curve')
	plt.show()

	Libor.plot()
	plt.xlabel('T')
	plt.ylabel('Forward Libor rate')
	plt.title('Forward Libor rate')
	plt.show()


#Q3
#input swaps
Fswap=np.array([[1,1],[1,2],[1,3],[1,5],
	[1,10],[5,1],[5,2],[5,3],[5,5],[5,10],
	[10,1],[10,2],[10,3],[10,5],[10,10]])

#list out all the Libor rate corresponding to each half year from 0-30y, there should be 60 libor rate

def calculate_Fswap(Fswap, OISD, Libor ):
	def fixed(i, OISD, Libor,S):
		fixedleg=sum(OISD.iloc[int(i[0]*2): int(i[0]*2+i[1]*2),0])*S
		return fixedleg
	def float(i, OISD, Libor):
		floatleg=0
		for j in range(1,int(i[1]*2+1)):
			floatleg+=OISD.iloc[int(i[0]*2+j-1),0]*Libor.iloc[int(i[0]*2+j-1),0]
		return floatleg
	Fswap_rate=[]
	for i in Fswap:
		f09= lambda S: fixed(i, OISD, Libor,S)-float(i, OISD, Libor)
		rate=brentq(f09,0,1)
		Fswap_rate.append(rate)
	return Fswap_rate

#turn Forward swap rate into DataFrame
Fswap_rate=calculate_Fswap(Fswap, OIS_all, Libor)
Fswap_rate=pd.DataFrame(np.array(Fswap_rate).reshape(3,5),index=['1Y','5Y','10Y'],columns=['1Y','2Y','3Y','5Y','10Y'])
if __name__ == '__main__':
	Fswap_rate.to_excel('Fswap_rate.xlsx')
	print(Fswap_rate)







































	








