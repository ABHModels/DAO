import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
top = '../testbb.dat'
# bot = '../illumination/top.txt'

def read(path):
    df=pd.read_csv(path,header=None,sep='\s+')
    ener = pd.to_numeric(df[0],errors='coerce')
    flux = pd.to_numeric(df[1],errors='coerce')
    return ener,flux

ener,topb = read(top)
# ener,botb = read(bot)

plt.figure()
plt.loglog(ener/1e3,topb*ener,label='top')
# plt.loglog(ener,botb*ener,label='bottom')
plt.legend()

plt.tight_layout()
plt.savefig('illumination.png',dpi=300)