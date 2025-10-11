import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


def read(path):
    df=pd.read_csv(path,header=None,sep='\s+')
    ener = pd.to_numeric(df[0],errors='coerce')
    kenrl = pd.to_numeric(df[1],errors='coerce')
    return ener,kenrl


# Thomson cross section
sigma_t = 6.6524587321e-25

# dir 
path = '../debug'

test_enregy = [1,100,500]
test_temper = [1e7,5e7,5e8]

lim = [(0,2),(0,200),(80,800)]

fig,ax = plt.subplots(3,3,figsize=(10,10))

for j,energy in enumerate(test_enregy):
    for i,temper in enumerate(test_temper):

        # exact kernel
        file_name_exact = f'k_{energy:.3f}keV_{np.log10(temper):.3f}kernel_exact5000.dat'

        # gaus kernel
        # file_name_gaus = f'k_{energy:.3f}keV_{np.log10(temper):.3f}gaus.dat'

        # read data
        # ener,gaus = read(os.path.join(path,file_name_gaus))
        ener,exact = read(os.path.join(path,file_name_exact))

        breakpoint()
        # plot
        # ax[i,j].plot(ener/1e3,gaus/sigma_t*1e5,linestyle='--',color='red',label='gaus')
        ax[i,j].plot(ener/1e3,exact/sigma_t*1e5,color='orange',label='exact')

        ax[i,j].set_xlim(lim[j])
    
        if i==0 and j==0:
            ax[i,j].legend()
        
        if i==2:
            ax[i,j].set_xlabel('Energy [keV]')
        
        if j==0:
            ax[i,j].set_ylabel(r'R[E$_i$,E$_j$] $\sigma_T$ [kev$^{-1}$]')

        if i==0:
            ax[i,j].set_title(f'E={energy:3.0f} keV')
        if j==0:
            ax[i,j].text(0.05,0.95,f'T={temper:.0e}'.replace('e+0', r'$\times 10^').replace('e+', r'$\times 10^').replace('e-', r'$\times 10^{-') + '}$ K',transform=ax[i,j].transAxes,ha='left',va='top')



plt.tight_layout()
plt.savefig('test_kernel.png',dpi=300)
        









