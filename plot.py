import numpy as np
import matplotlib.pyplot as plt
nm=10**-9
mev=1.6*10**-22
res=0.01
phis=np.arange(-np.pi/2,np.pi/2,res)

kf=2*np.pi/(50*nm)

D=100*nm
E=83*mev
V0=285*mev
vf=10**6
h_bar=(6.626*10**-34/(2*np.pi))

V=[0,V0,0]
s=np.sign(E-np.array(V))
x=[0,D]

def phi(i,ang):
	return np.arccos(kx(i,ang)/k(i))
		
def k(i):
	return np.abs(E-V[i])/(h_bar*vf)

def kx(i,ang):
	return np.sqrt(k(i)**2-(k(0)*np.sin(ang))**2)

def M(i,j,ang):
	N=[[np.exp(1j*kx(i,ang)*x[j]),np.exp(-1j*kx(i,ang)*x[j])],[s[i]*np.exp(1j*kx(i,ang)*x[j]+1j*phi(i,ang)),-1*s[i]*np.exp(-1j*kx(i,ang)*x[j]-1j*phi(i,ang))]]
	return np.array(N)

T=[]

for ang in phis:
	P=[[1,0],[0,1]]
	for m in range(0,len(x)):
		S=np.dot(np.linalg.inv(M(m,m,ang)),M(m+1,m,ang))
		P=np.dot(P,S)
	T.append(np.abs(P[0,0])**-2)

plt.polar(phis,np.array(T))
plt.savefig('res/'+str(int(V0/mev))+'.png')
plt.show()
