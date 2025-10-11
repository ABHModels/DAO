import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

refsknfile = '../skntemp.dat'
pursknfile = '/data/home/yimin2/PureScat/skntest.dat'

def read(file):
    df = pd.read_csv(file,sep='\s+',header=None)
    ener = pd.to_numeric(df[0],errors='coerce')
    skn = pd.to_numeric(df[1],errors='coerce')
    return ener,skn

refener,refskn = read(refsknfile)
pursener,purskn = read(pursknfile)

plt.plot(refener,refskn/refskn,label='ref')
plt.xscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('tskn.png')

