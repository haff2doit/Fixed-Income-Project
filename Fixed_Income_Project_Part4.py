import pandas as pd
import numpy as np
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
import matplotlib.pylab as plt
from scipy import interpolate, integrate
from scipy.optimize import least_squares
from enum import Enum
from Fixed_Income_Project_Part_1 import OIS_all, Fswap_rate, Libor,calculate_Fswap, FF_all_halfyear
from Fixed_Income_Project_Part_2 import SABR, Black76,Type
from Fixed_Income_Project_Part_3 import IRR
from scipy.misc import derivative


def Valuation(Swap_rate, IRR, Black76, SABR, alpha, rho, nu, beta, N,m, T,g, OISD, plus):
	def integrant(K):
		sabrvol=SABR(Swap_rate,K, T,alpha,beta,rho,nu)
		def h(g,x):
			return g(x)/IRR(x,m,N)
		ddhk=derivative(lambda x: h(g,x),K,dx=0.0001,n=2)
		BlackReceiver=Black76(Swap_rate, K,1, sabrvol, T,Type.Receiver)
		BlackPayer=Black76(Swap_rate, K,1, sabrvol, T,Type.Payer)
		return [ddhk*BlackReceiver,ddhk*BlackPayer]
	F=Swap_rate
	if plus==False:
		V=OISD.iloc[T*2-1,0]*(g(F)+IRR(Swap_rate,m,N)*(integrate.quad(lambda K:integrant(K)[0],0,F)[0]+integrate.quad(lambda K:integrant(K)[1],F,np.inf)[0]))
	if plus==True:
		sabrvol=SABR(Swap_rate,.04**(.5), T,alpha,beta,rho,nu)
		def h(g,x):
			return g(x)/IRR(x,m,N)
		dhk=derivative(lambda x: h(g,x),0.0016,dx=0.0001,n=1)
		BlackPayer=Black76(Swap_rate, 0.0016,1, sabrvol, T,Type.Payer)
		V=OISD.iloc[T*2-1,0]*IRR(Swap_rate,m,N)*(dhk*BlackPayer+integrate.quad(lambda K:integrant(K)[1],0.0016,np.inf)[0])
	return V

Swap_rate=Fswap_rate.loc['5Y','10Y']
Alpha=pd.read_excel('alpha.xlsx',header=0,index_col=0)
Rho=pd.read_excel('rho.xlsx',header=0,index_col=0)
Nu=pd.read_excel('nu.xlsx',header=0,index_col=0)
print(Swap_rate)

alpha=Alpha.loc['5Y','10Y']
rho=Rho.loc['5Y','10Y']
nu=Nu.loc['5Y','10Y']
beta=0.9

N=10
m=2
T=5
p=4
q=2
g1=lambda x: x**(1/p)-0.04**(1/q)
g2=lambda x: x**(1/p)-0.04**(1/q)

PV1=Valuation(Swap_rate, IRR, Black76, SABR, alpha, rho, nu, beta, N,m, T,g1,OIS_all,plus=False)
print(PV1)
PV2=Valuation(Swap_rate, IRR, Black76, SABR, alpha, rho, nu, beta, N,m, T,g2,OIS_all,plus=True)
print(PV2)





