import os 
import pandas as pd
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def normflux(ener, flux):
    rflux = flux[(ener>10) & (ener<100)]
    rener = ener[(ener>10) & (ener<100)]
    normfactor = np.trapz(rflux, rener)
    return flux/normfactor

def read(path):
    data = pd.read_csv(path, sep='\s+', header=None).iloc[-1000:]
    ener = pd.to_numeric(data[0],errors='coerce')
    ener = ener/1e3
    top = pd.to_numeric(data[1],errors='coerce')
    flux = pd.to_numeric(data[3],errors='coerce')
    # return ener, normflux(ener, flux*ener)
    return ener,ener*flux

os.chdir('/data/home/yimin2/reflection')

colors = mpl.cm.plasma(np.linspace(0.2,0.8,8))
fig,ax = plt.subplots(1,2)
path = 'metadata'
metadata_files = os.listdir(path)
i = 0 
j= 0 
for metadata_file in metadata_files:
    with open(os.path.join(path, metadata_file), 'r') as f:
        metadata = json.load(f)
        dp = metadata['parameters']['op_spec']
        if metadata['parameters']['diy_tag']=='dxi':
            if metadata['parameters']['kernltype']=='gaus':
                i = i+1 
                ener,eemo = read(dp)
                ax[0].loglog(ener,eemo,color=colors[i])
            if metadata['parameters']['kernltype']=='exact':
                j = j+1
                ener,eemo = read(dp)
                ax[1].loglog(ener,eemo,color=colors[j])
    ax[0].set_ylabel(r'EF$_E$ (erg$^2$/cm$^2$/s/erg)')


    ax[0].set_title('Gaussian Kernel')
    ax[1].set_title('Exact Kernel')

    ax[0].set_xlabel('Energy (keV)')
    ax[1].set_xlabel('Energy (keV)')
    ax[1].set_yticks([])

plt.tight_layout()
plt.savefig('/data/home/yimin2/reflection/python-anaylsis/dxi.png',dpi=300)


 