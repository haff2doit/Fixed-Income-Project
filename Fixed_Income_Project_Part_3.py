import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
import matplotlib.pylab as plt
from scipy import interpolate, integrate
from scipy.optimize import least_squares
from enum import Enum
from Fixed_Income_Project_Part_1 import OIS_all, Fswap_rate, Libor,calculate_Fswap, FF_all_halfyear, IRSD
from Fixed_Income_Project_Part_2 import SABR, Black76,Type
from scipy.misc import derivative

Alpha=pd.read_excel('alpha.xlsx',header=0,index_col=0)
Rho=pd.read_excel('rho.xlsx',header=0,index_col=0)
Nu=pd.read_excel('nu.xlsx',header=0,index_col=0)

def IRR(Swap_rate,m,N):
	irr=0
	for i in range(1,int(N*m+1)):
		irr+=1/(1+Swap_rate/m)**i 
	return irr

def CMSrate(Swap_rate, IRR, Black76, SABR, alpha, rho, nu, beta, N,m, T,g):
	def integrant(K):
		sabrvol=SABR(Swap_rate,K, T,alpha,beta,rho,nu)
		def h(g,x):
			return g(x)/IRR(x,m,N)
		ddhk=derivative(lambda x: h(g,x),K,dx=0.0001,n=2)
		BlackReceiver=Black76(Swap_rate, K,1, sabrvol, T,Type.Receiver)
		BlackPayer=Black76(Swap_rate, K,1, sabrvol, T,Type.Payer)
		return [ddhk*BlackReceiver,ddhk*BlackPayer]
	F=Swap_rate
	CMS=g(F)+IRR(Swap_rate,m,N)*(integrate.quad(lambda K:integrant(K)[0],0,F)[0]+integrate.quad(lambda K:integrant(K)[1],F,1)[0])
	return CMS

def Liborquarter(IRSD):
	fl=interpolate.interp1d(IRSD.index, IRSD.iloc[:,0])
	TL_new=list(np.arange(0.25,30.25,0.25))
	IRSD_quarter=fl(TL_new)
	Libor_quater=[]
	for i in range(len(IRSD_quarter)-2):
		L=2*(IRSD_quarter[i]-IRSD_quarter[i+2])/IRSD_quarter[i+2]
		Libor_quater.append(L)
	return Libor_quater

if __name__ == '__main__':	
	T1=list(np.arange(0.5,5.5,0.5))
	T2=list(np.arange(0.25,10.25,0.25))
	xalpha=[1,5,10]
	yalpha=[Alpha.loc['1Y','10Y'],Alpha.loc['5Y','10Y'],Alpha.loc['10Y','10Y']]
	falpha=interpolate.interp1d(xalpha, yalpha, fill_value='extrapolate')
	alpha1=falpha(T1)
	alpha2=falpha(T2)

	xrho=[1,5,10]
	yrho=[Rho.loc['1Y','10Y'],Rho.loc['5Y','10Y'],Rho.loc['10Y','10Y']]
	frho=interpolate.interp1d(xrho, yrho,fill_value='extrapolate')
	rho1=frho(T1)
	rho2=frho(T2)

	xnu=[1,5,10]
	ynu=[Nu.loc['1Y','10Y'],Nu.loc['5Y','10Y'],Nu.loc['10Y','10Y']]
	fnu=interpolate.interp1d(xnu, ynu,fill_value='extrapolate')
	nu1=fnu(T1)
	nu2=fnu(T2)

	#calculate quarterly OIS discount factor
	mutiplier= np.array((1+FF_all_halfyear/360)**90).reshape(-1,1)
	OIS_quarter=OIS_all*mutiplier
	OIS_all_quarter=[]
	for i in range(60):
		OIS_all_quarter.append(OIS_quarter.iloc[i,0])
		OIS_all_quarter.append(OIS_all.iloc[i,0])
	OIS_all_quarter=pd.DataFrame(OIS_all_quarter, index=np.arange(0.25,30.25,0.25),columns=['OIS discount factors'])


	Swap1=list(zip(T1,[10]*10))
	Swap_rate1=calculate_Fswap(Swap1, OIS_all, Libor)

	Libor_quarter=Liborquarter(IRSD)
	Libor_quarter=pd.DataFrame(Libor_quarter, index=np.arange(0.25,29.75,0.25),columns=['Libor quarter'])

	Swap2=list(zip(T2,[2]*40))
	Swap_rate2=calculate_Fswap(Swap2, OIS_all, Libor)
	N1=[10]*len(T1)
	m1=2
	beta=0.9
	n1=5
	pf1=2
	g=lambda x: x
	PV1=0
	for i in range(len(Swap_rate1)):
		CMS=CMSrate(Swap_rate1[i], IRR, Black76, SABR, alpha1[i], rho1[i], nu1[i], beta, N1[i], m1, T1[i],g)
		PV1+=OIS_all.iloc[i,0]*(1/pf1)*CMS
	print('PV1',PV1)

	N2=[2]*len(T2)
	m2=2
	n2=10
	pf2=4
	PV2=0
	for i in range(len(Swap_rate2)):
		CMS=CMSrate(Swap_rate2[i], IRR, Black76, SABR, alpha2[i], rho2[i], nu2[i], beta, N2[i], m2, T2[i],g)
		PV2+=OIS_all_quarter.iloc[i,0]*(1/pf2)*CMS
	print('PV2',PV2)

	Fswap=[[1,1],[1,2],[1,3],[1,5],
		[1,10],[5,1],[5,2],[5,3],[5,5],[5,10],
		[10,1],[10,2],[10,3],[10,5],[10,10]]
	Swap_rate3=Fswap_rate.values.tolist()
	Swap_rate3=[j for sub in Swap_rate3 for j in sub]
	alpha3=Alpha.values.tolist()
	alpha3=[j for sub in alpha3 for j in sub]
	rho3=Rho.values.tolist()
	rho3=[j for sub in rho3 for j in sub]
	nu3=Nu.values.tolist()
	nu3=[j for sub in nu3 for j in sub]

	N3=[1,2,3,5,10]*3
	m3=2
	T3=[i[0] for i in Fswap]

	CMS_rate=[]
	for i in range(len(Swap_rate3)):
		CMS=CMSrate(Swap_rate3[i], IRR, Black76, SABR, alpha3[i], rho3[i], nu3[i], beta, N3[i], m3, T3[i],g)
		CMS_rate.append(CMS)
	expiry=['1Y','5Y','10Y']
	tenor=['1Y','2Y','3Y','5Y','10Y']
	CMS_rate=pd.DataFrame(np.array(CMS_rate).reshape(3,5),index=expiry,columns=tenor)
	print('CMS rate\n',CMS_rate)
	print('Fswap rate\n',Fswap_rate)
	convexity=CMS_rate-Fswap_rate
	print('convexity\n',convexity)

	print(IRSD)



















