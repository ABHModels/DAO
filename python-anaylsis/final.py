import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import os

def read(path):
    df=pd.read_csv(path,header=None,sep='\s+')
    ener = pd.to_numeric(df[0],errors='coerce')
    top = pd.to_numeric(df[1],errors='coerce')
    bot = pd.to_numeric(df[2],errors='coerce')
    inte = pd.to_numeric(df[3],errors='coerce')
    return ener,top,bot,inte
temps = ['2']
file = ['spec_2406b573fe_kernel_test.dat','spec_2406b573fe_kernel_test.dat','spec_c1209a7858_kernel_test.dat']
file = ['spec_2406b573fe_kernel_test.dat']
evtoerg = 1.60217662e-12
colors = mpl.cm.viridis(np.linspace(0,1,3))
for i,temp in enumerate(temps):

    ener,top,bot,inte = read(f'../out/{file[i]}')

    if i==0: 
        plt.loglog(ener/1e3,ener*top*evtoerg,color='k',linestyle='--',label=f'top')
        plt.loglog(ener/1e3,ener*bot*evtoerg,color='k',linestyle='-.',label=f'bottom')
    # plt.loglog(ener/1e3,ginte*ener*evtoerg,linestyle='--',color=colors[i],label=f'{temp}keV')
    plt.loglog(ener/1e3,inte*ener*evtoerg,color=colors[i],label=f'{temp}keV')

print(min(ener),max(ener))
plt.xlabel('Energy [KeV]')
plt.ylabel(r'EI$^+$ [erg$^2$/s/cm$^2$/erg]')
plt.xlim(1e-4,1000)
plt.ylim(1e4,1e18)
tick_positions = [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]
tick_labels = [r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',r'$1$',r'$10$',r'$10^2$',r'$10^3$']
plt.xticks(ticks=tick_positions, labels=tick_labels)
plt.legend()
plt.tight_layout()
plt.savefig('fin.png',dpi=300)