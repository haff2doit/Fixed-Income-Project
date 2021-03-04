import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
import matplotlib.pylab as plt
from scipy import interpolate
from scipy.optimize import least_squares
from enum import Enum
from Fixed_Income_Project_Part_1 import OIS_all, Fswap_rate, Libor,calculate_Fswap

class Type(Enum):
    Payer = 0
    Receiver = 1

def Black76(F,K,P,sigma,T,Type):
    d1= (np.log(F/K)+(1/2)*sigma**2*T) / (sigma*np.sqrt(T))
    d2=d1 - sigma*np.sqrt(T)
    if Type==Type.Payer:
        return P*(F*norm.cdf(d1) - K*norm.cdf(d2))
    if Type==Type.Receiver:
        return P*(K*norm.cdf(-d2) - F*norm.cdf(-d1))
    else:
        raise Exception('payoffType not supported')

def impliedVolatility(F,K,P,price,T,Type):
    impliedVol = brentq(lambda x: price -Black76(F,K,P,x,T,Type), 1e-6, 1)
    return impliedVol

def Displaced_diffusion(F,K,P,sigma,T,beta,Type):
    Fd=F/beta
    sigmad=sigma*beta
    Kd=K+((1-beta)/beta)*F
    d1 = (np.log(Fd/Kd)+(sigmad**2/2)*T) / (sigmad*np.sqrt(T))
    d2 = d1 - sigmad*np.sqrt(T)
    if Type==Type.Payer:
        return P*(Fd*norm.cdf(d1) - Kd*norm.cdf(d2))
    if Type==Type.Receiver:
        return P*(Kd*norm.cdf(-d2) - Fd*norm.cdf(-d1))
    else:
        raise Exception('payoffType not supported')

def SABR(F, K, T, alpha, beta, rho, nu):
    X = K
    # if K is at-the-money-forward
    if abs(F - K) < 1e-12:
        numer1 = (((1 - beta)**2)/24)*alpha*alpha/(F**(2 - 2*beta))
        numer2 = 0.25*rho*beta*nu*alpha/(F**(1 - beta))
        numer3 = ((2 - 3*rho*rho)/24)*nu*nu
        VolAtm = alpha*(1 + (numer1 + numer2 + numer3)*T)/(F**(1-beta))
        sabrsigma = VolAtm
    else:
        z = (nu/alpha)*((F*X)**(0.5*(1-beta)))*np.log(F/X)
        zhi = np.log((((1 - 2*rho*z + z*z)**0.5) + z - rho)/(1 - rho))
        numer1 = (((1 - beta)**2)/24)*((alpha*alpha)/((F*X)**(1 - beta)))
        numer2 = 0.25*rho*beta*nu*alpha/((F*X)**((1 - beta)/2))
        numer3 = ((2 - 3*rho*rho)/24)*nu*nu
        numer = alpha*(1 + (numer1 + numer2 + numer3)*T)*z
        denom1 = ((1 - beta)**2/24)*(np.log(F/X))**2
        denom2 = (((1 - beta)**4)/1920)*((np.log(F/X))**4)
        denom = ((F*X)**((1 - beta)/2))*(1 + denom1 + denom2)*zhi
        sabrsigma = numer/denom

    return sabrsigma

def sabrcalibration(x, strikes, vols, F, T):
    err = 0.0
    for i, vol in enumerate(vols):
        if i==5:
            err += 10*(vol - SABR(F, strikes[i], T,
                               x[0], 0.9, x[1], x[2]))**2
        else:
            err += (vol - SABR(F, strikes[i], T,
                               x[0], 0.9, x[1], x[2]))**2

    return err

def DDcalibration(F,ATMK,ReceiverK,PayerK,vols,x,T,P):
    K_Receiver= ATMK+ReceiverK
    K_Receiver=list(K_Receiver)+[ATMK]
    implied_Receiver=[]
    for i in K_Receiver:
        DD_Receiver=Displaced_diffusion(F,i,P,x[0],T,x[1],Type.Receiver)
        implied_Rec=impliedVolatility(F,i,P,DD_Receiver,T,Type.Receiver)
        implied_Receiver.append(implied_Rec)

    K_Payer= ATMK+PayerK
    K_Payer=list(K_Payer)
    implied_Payer=[]
    for i in K_Payer:
        DD_Payer=Displaced_diffusion(F,i,P,x[0],T,x[1],Type.Payer)
        implied_Pay=impliedVolatility(F, i,P,DD_Payer,T,Type.Payer)
        implied_Payer.append(implied_Pay)

    impliedvol_DD=implied_Receiver+implied_Payer

    err=0.0
    for i, vol in enumerate(vols):
        if i==5:
            err += 10*(vol - impliedvol_DD[i])**2
        else:
            err += (vol - impliedvol_DD[i])**2
    return err

if __name__ == '__main__':
    marketvol=pd.read_excel('IR Data.xlsx', sheet_name='Swaption', header=0, index_col=[0,1], skiprows=[0,1])
    print(marketvol)

    strikes=list(marketvol.columns)
    strikes.remove('ATM')
    allK=[int(i.strip('bps'))/10000 for i in strikes]
    ReceiverK=np.array([i for i in allK if i<0])
    PayerK=np.array([i for i in allK if i>0])
    print(OIS_all)
    print('Fswap_rate\n', Fswap_rate)

    Beta=[]
    Vol=[]
    Alpha=[]
    Rho=[]
    Nu=[]

    initialGuess1=[0.3,0.6]
    initialGuess2=[0.3, -0.5, 0.5]

    for i in range(len(marketvol.index)):
    #calibrate DD model
        vols=marketvol.iloc[i,:].values/100
        S=Fswap_rate.loc[marketvol.index[i][0],marketvol.index[i][1]]
        ATMK=S
        T=int(marketvol.index[i][0].strip('Y'))
        n=int(marketvol.index[i][0].strip('Y'))*2
        N=n+int(marketvol.index[i][1].strip('Y'))*2
        P=sum(OIS_all.iloc[n:N,0])/2
        res1 = least_squares(lambda x: DDcalibration(S,ATMK,ReceiverK,PayerK,vols,x,T,P),initialGuess1,bounds=([0.01,0.01],[0.6,1]))
        vol=res1.x[0]
        beta1=res1.x[1]
        Vol.append(vol)
        Beta.append(beta1)
    #calibrate SABR model
        strikes_all=np.array(list(ReceiverK)+[0]+list(PayerK))+ATMK
        res2= least_squares(lambda x: sabrcalibration(x, strikes_all, vols, S, T), initialGuess2,max_nfev=3000)
        alpha=res2.x[0]
        rho=res2.x[1]
        nu=res2.x[2]
        Alpha.append(alpha)
        Rho.append(rho)
        Nu.append(nu)

    expiry=['1Y','5Y','10Y']
    tenor=['1Y','2Y','3Y','5Y','10Y']
    Vol=pd.DataFrame(np.array(Vol).reshape(3,5),index=expiry,columns=tenor)
    Beta=pd.DataFrame(np.array(Beta).reshape(3,5),index=expiry,columns=tenor)
    Alpha=pd.DataFrame(np.array(Alpha).reshape(3,5),index=expiry,columns=tenor)
    Rho=pd.DataFrame(np.array(Rho).reshape(3,5),index=expiry,columns=tenor)
    Nu=pd.DataFrame(np.array(Nu).reshape(3,5),index=expiry,columns=tenor)
    
    Vol.to_excel('sigma.xlsx')
    Beta.to_excel('beta.xlsx')
    Alpha.to_excel('alpha.xlsx')
    Rho.to_excel('rho.xlsx')
    Nu.to_excel('nu.xlsx')

    print('sigma\n',Vol)
    print('beta\n',Beta)
    print('alpha\n',Alpha)
    print('rho\n',Rho)
    print('nu\n',Nu)


    # Price Swaption
    Swap=[[2,10],[8,10]]
    K=list(np.arange(0.01,0.09,0.01))
    Swap_rate=calculate_Fswap(Swap, OIS_all, Libor)
    PnN=[]
    T=[2,8]
    for i in Swap:
        n=i[0]*2
        N=n+i[1]*2
        p=sum(OIS_all.iloc[n:N,0])/2
        PnN.append(p)

    #Use DD model to price
    xsigma=[1,5,10]
    ysigma=[Vol.loc['1Y','10Y'],Vol.loc['5Y','10Y'],Vol.loc['10Y','10Y']]
    fsigma=interpolate.interp1d(xsigma, ysigma)
    sigma1,sigma2=fsigma(T)
    print('sigma1',sigma1)
    print('sigma2',sigma2)

    xbeta=[1,5,10]
    ybeta=[Beta.loc['1Y','10Y'],Beta.loc['5Y','10Y'],Beta.loc['10Y','10Y']]
    fbeta=interpolate.interp1d(xbeta, ybeta)
    beta1,beta2=fbeta(T)
    print('beta1',beta1)
    print('beta2',beta2)

    DD1=[]
    DD2=[]
    for i in K:
            dd1=Displaced_diffusion(Swap_rate[0],i,PnN[0],sigma1,T[0],beta1,Type.Payer)
            DD1.append(dd1)
            dd2=Displaced_diffusion(Swap_rate[1],i,PnN[1],sigma2,T[1],beta2,Type.Receiver)
            DD2.append(dd2)
    DD1=pd.DataFrame(np.array(DD1).reshape(1,len(K)),index=['payer 2y10y'],columns=K)
    DD2=pd.DataFrame(np.array(DD2).reshape(1,len(K)),index=['receiver 8y10y'],columns=K)
    DD=DD1.append(DD2)
    print('price under DD model\n',DD)

    #use SABR model to price
    xalpha=[1,5,10]
    yalpha=[Alpha.loc['1Y','10Y'],Alpha.loc['5Y','10Y'],Alpha.loc['10Y','10Y']]
    falpha=interpolate.interp1d(xalpha, yalpha)
    alpha1,alpha2=falpha(T)
    print('alpha1',alpha1)
    print('alpha2',alpha2)

    xrho=[1,5,10]
    yrho=[Rho.loc['1Y','10Y'],Rho.loc['5Y','10Y'],Rho.loc['10Y','10Y']]
    frho=interpolate.interp1d(xrho, yrho)
    rho1, rho2=frho(T)
    print('rho1',rho1)
    print('rho2',rho2)

    xnu=[1,5,10]
    ynu=[Nu.loc['1Y','10Y'],Nu.loc['5Y','10Y'],Nu.loc['10Y','10Y']]
    fnu=interpolate.interp1d(xnu, ynu)
    nu1,nu2=fnu(T)
    print('nu1',nu1)
    print('nu2',nu2)

    SABR1=[]
    SABR2=[]
    for i in K:
            sabrsigma1=SABR(Swap_rate[0], i, T[0], alpha1, 0.9, rho1, nu1)
            sabr1=Black76(Swap_rate[0],i,PnN[0],sabrsigma1,T[0],Type.Payer)
            SABR1.append(sabr1)
            sabrsigma2=SABR(Swap_rate[1], i, T[1], alpha2, 0.9, rho2, nu2)
            sabr2=Black76(Swap_rate[1],i,PnN[1],sabrsigma2,T[1],Type.Receiver)
            SABR2.append(sabr2)
    SABR1=pd.DataFrame(np.array(SABR1).reshape(1,len(K)),index=['payer 2y10y'],columns=K)
    SABR2=pd.DataFrame(np.array(SABR2).reshape(1,len(K)),index=['receiver 8y10y'],columns=K)
    SABR3=SABR1.append(SABR2)
    print('price under SABR model\n',SABR3)







