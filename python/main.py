import sys
sys.path.append('C:/MatkivskiyV/!Актуальная наука/Статья 7/python/')

import numpy as np
import matplotlib.pyplot as plt
import fun_oct as fo
import fields_synthesis as fs

case = 'abless'

tmpl1 = np.arange(-1, 1, 1/128)
tmpl2 = np.arange(-1, 1, 1/256)
[x,y] = np.meshgrid(tmpl1,tmpl2)
r = np.sqrt(x**2+y**2)

if case == 'abless':
#    filename = 'C:/MatkivskiyV/!Актуальная наука/статья 8/data/USAF_2.dat'
#    B = A[::2,97:101,:]
    filename = 'C:/MatkivskiyV/!Актуальная наука/статья 8/data/USAF_4.dat'    
    A = fo.import_SU_dat(filename)
    B = A[::2,50,:]
elif case == 'ab':
    filename = 'C:/MatkivskiyV/!Актуальная наука/статья 8/data/USAF_1_ab.dat'
    A = fo.import_SU_dat(filename)
    B = A[::2,190:195,:]
elif case == 'abtor':
    filename = 'C:/MatkivskiyV/!Актуальная наука/статья 8/data/USAF_ab_tor.dat'
    A = fo.import_SU_dat(filename)
    B = A[1::2,225:244,:]


#C = np.sum(B,axis=1)
C = B
fC = np.fft.fft2(C)
plt.figure(3); plt.imshow(np.abs( fC )); plt.colorbar()
fC = np.fft.fftshift(fC,axes=0)

if case == 'abless':
    plt.figure(1); plt.imshow(np.abs( fC*fs.circ(r,0.5) ), 
               aspect=0.5, vmax=4*10**6)
    plt.colorbar()
    
    C2ab = np.fft.ifft2( fC*fs.circ(r,0.5) )
    plt.figure(2); plt.imshow(np.abs( C2ab ),aspect=0.5,vmin=800,vmax=7000)
    plt.colorbar()
elif case == 'ab':
    plt.figure(3); plt.imshow(np.abs( fC*fs.circ(r,0.5) ), 
               aspect=0.5,vmax=10**6)
    plt.colorbar()
    
    C2ab = np.fft.ifft2( fC*fs.circ(r,0.5) )
    plt.figure(4); plt.imshow(np.abs( C2ab ),aspect=0.5,vmin=300,vmax=2000)
    plt.colorbar()
elif case == 'abtor':
    plt.figure(1); plt.imshow(np.abs( fC*fs.circ(r,0.5) ), aspect=0.5)
    plt.colorbar()
    
    C2ab = np.fft.ifft2( fC*fs.circ(r,0.5) )
    plt.figure(2); plt.imshow(np.abs( C2ab ),aspect=0.5, vmax=1000)
    plt.colorbar()
    
