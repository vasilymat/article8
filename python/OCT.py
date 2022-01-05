import numpy as np
from numpy import pi

def FTS(T):        
    return np.fft.fftshift( np.fft.fft(np.fft.fftshift(T)) )
    
def iFTS(T):
    return np.fft.ifftshift( np.fft.ifft(np.fft.ifftshift(T)) )    
    

class OCT2D:
    
    def __init__(self,Ascan_smpl=128,dz=2.,Xsmpl=512,dx=.53,
                 lam=.85,foc=200,SpecW=.07,Arad=40):        
        self.aspx = Ascan_smpl #number of pixels in A-scan
        self.lam  = lam
        self.Xsmpl = Xsmpl
        self.foc  = foc
        self.SpecW = SpecW #Spectrum wight
        self.dx = dx
        self.dz = dz
        self.Arad = Arad
        
        L = dx*Xsmpl
        self.L  = L
        k = 2*pi/lam
        self.k  = k
        
        self.x  = np.arange(-L/2, L/2, dx)
        self.kx = np.arange(-2*pi/(2*dx), 2*pi/(2*dx), 2*pi/L)
        self.kz = np.sqrt(k**2-self.kx**2) #Is it necessary?
        
        sig = 1/np.sqrt(2)*Arad
        self.UL = np.exp(-self.x**2/2/sig**2).astype(complex)
        Lens = k*self.x**2/(2*foc)    
        self.UL *= np.exp(-1j*Lens)*(np.abs(self.x) < Arad)
        
    def get_params(self):
        print('Ascan semples = ',self.aspx)
        print('Wavelength = ',self.lam)
        print('Semples number of X = ',self.Xsmpl)
        print('Focal lenght = ', self.foc)
        print('Spectrum wight = ', self.SpecW)
        print('dx size = ', self.dx)
        print('dz size = ', self.dz)
        print('Aperture radius = ', self.Arad)
    
    def calculate_field(self,K):
        self.K = K
        specWd2 = self.SpecW/2
        dz = self.dz
        dx = self.dx
        ld2 = np.round(self.Xsmpl/2)
        
        SWa = self.SpecW/self.aspx
        Spec = np.arange(self.lam-specWd2,self.lam+specWd2-SWa,SWa)
        Spec = 2*np.pi/Spec
        Cr = K.nonzero()
        
        A = np.zeros(self.aspx,dtype=complex) 
        sig = 1/np.sqrt(2)*self.Arad 
        
        self.B  = np.zeros((self.aspx,self.Xsmpl),dtype=complex)
        Bt = np.zeros((self.aspx,self.Xsmpl),dtype=complex)
        
        kx2 = self.kx**2
        g = np.sqrt(1/self.UL.shape[0])*FTS(self.UL)
        
        for izs,ixs,Kt in zip(Cr[0],Cr[1],K[Cr]):
            xs = (ixs-ld2)*dx
            zs = izs*dz
            for j,x0 in enumerate(self.x):    
                Sh  = np.exp(1j*self.kx*(xs-x0))

                for i,k in enumerate(Spec):
                    kz = np.sqrt(k**2-kx2)
                    h  = np.exp(1j*kz*zs)
                    Us = ( g*h*Sh ).sum()
                    A[i] = Kt*Us**2
                Bt[:,j] = iFTS(A)
            self.B += Bt
        return self.B