import numpy as np
import math

def lambda_r(lambada,r_i,Radius):
    return lambada*r_i/Radius

def phi1(lambda_r):          # initial phi
    return 2./3*np.arctan(lambda_r)

def sigma1_i(B,c_i,r_i):
    return B*c_i*r_i

def a1_i(phi,sigma_i,Cl_design):
    return 1./(1+(2*np.sin(phi))**2/(sigma_i*Cl_design*np.cos(phi)))
''' a' '''
def A1_i(a_i):
    return (1-3*a_i)/(4*a_i-1)

''' iterating elements '''
def Phi1_i(a_i,A_i,lambda_r):        # iterating phi
    return np.arctan((1-a_i)/((1+a_i)*lambda_r))

def tiploss1_i(Phi_i,r_i,B,R):
    return (2./np.pi)*np.arccos(np.exp(-1*((B*0.5)*(1-r_i/R))/(r_i*np.sin(Phi_i)/R)))

def alpha1_i(Phi_i,theta_i):
    return Phi_i-theta_i

def Cl1(alpha_i):
    deg = alpha_i*(180/3.14)
    z=np.float((5.863*math.pow(10,-7))*(deg**4) - (1.081*math.pow(10,-4))*(deg**3) -(0.0009739*(deg**2)) + (0.1144*(deg)) + (0.4865))
    return z

def Cd1(alpha_i):
    deg = alpha_i*(180/3.14)
    return np.float((7.692*math.pow(10,-7))*(deg**4) - (1.832*math.pow(10,-6))*(deg**3) -(3.828*math.pow(10,-5)*(deg**2)) + (4.201*  math.pow(10,-4))*(deg) + (0.009237))

def Ct1(sigma_i,Cl,Cd,Phi_i,a_i):
    return (sigma_i*((1-a_i)**2)*(Cl*np.cos(Phi_i)+Cd*np.sin(Phi_i)))/(np.sin(Phi_i))**2

def a_i_new(Ct,Phi_i,Cl,sigma_i,tiploss_i,length):      # returns both a_i and a'_i
    a_i=0
    A_i=0
    for i in range(length):
        if Ct<0.96:
            a_i=1/(1+(4*tiploss_i*(np.sin(Phi_i))**2)/(sigma_i*Cl*np.cos(Phi_i)))
        else:
            a_i=(1/tiploss_i)*(0.143+np.sqrt(0.0203-0.6427*(0.889-Ct)))
        A_i=1/((4*tiploss_i*np.cos(Phi_i))/(sigma_i*Cl)-1)
    return a_i,A_i

length=19     #number of divisions
Radius=1.87
B=3
Cl_design=1
theta=np.array([ 0.7950,0.6713,0.5647,0.4751,0.4014,0.3410,0.2913,0.2502,0.2159,0.1869,0.1622,0.1410,0.1226,0.1065,0.0923,0.0797,0.0685,0.0585,0.0494])
c_i=np.array([0.022,0.034,0.039,0.04,0.038,0.036,0.034,0.031,0.029,0.027,0.025,0.023,0.022,0.02,0.019,0.017,0.015,0.012,0.009])
r_i=np.array([0.089,0.178,0.267,0.356,0.445,0.534,0.623,0.712,0.801,0.891,0.98,1.069,1.158,1.247,1.336,1.425,1.514,1.603,1.692])
lambada=np.linspace(2,9,10)

for x in lambada:
    lambada_i=lambda_r(x,r_i,Radius)            # your r_i goes here
    phi=phi1(lambada_i)
    sigma_i=sigma1_i(B,c_i,r_i)
    #print(sigma_i.shape)
    a_i=a1_i(phi,sigma_i,Cl_design)
    A_i=A1_i(a_i)
    #print(x)
    Cp=0

    for j in range(length):
        a_old=a_i[j]
        A_old=A_i[j]
        zzzz=0
        while True:
            Phi_i=Phi1_i(a_old,A_old,lambada_i[j])
            #print(Phi_i, type(Phi_i))
            tiploss_i=tiploss1_i(Phi_i,r_i[j],B,Radius)
            alpha_i=alpha1_i(Phi_i,theta[j])
            Cl=Cl1(alpha_i)
            Cd=Cd1(alpha_i)
            Ct=Ct1(sigma_i[j],Cl,Cd,Phi_i,a_old)
            a_new,A_new=a_i_new(Ct,Phi_i,Cl,sigma_i[j],tiploss_i,length)
            #print(a_old, A_old)
            if np.abs(a_old-a_new)<0.01 or zzzz>1000:
                #print(Phi_i, alpha_i)
                a_i[j]=a_new
                A_i[j]=A_new
                Cp+=8/(x*length)*(tiploss_i*(np.sin(Phi_i))**2)*(np.cos(Phi_i)-lambada_i[j]*np.sin(Phi_i))*(np.sin(Phi_i)+lambada_i[j]*np.cos(Phi_i))*(1-(Cd/Cl)*(1/np.tan(Phi_i)))*lambada_i[j]**2

                break
            else:
                a_old=a_new
                A_old=A_new
            #print(type(a_old), A_old)
            #print()
            zzzz+=1
    #print("\n")
    print(x,Cp)





'''     power coeffienent       '''

for x in lambada:
    lambada_i=lambda_r(x,r_i,Radius)            # your r_i goes here
    phi=phi1(lambada_i)
    sigma_i=sigma1_i(B,c_i,r_i)
    #print(sigma_i.shape)
    a_i=a1_i(phi,sigma_i,Cl_design)
    A_i=A1_i(a_i)
