import numpy as np
import pandas as pd


gamma = [1.2]

Kte = 60 


for g in gamma:
    path = f'spec_source_incident_spectra/nthcomp_g{g}_t{Kte}_spec.txt'

    try:
        data = pd.read_csv(path,skiprows=3,header=None,sep='\s+')
    except:
        print(f'Fail at {path}')
        continue
    ener = data[0]
    spec = data[2]

    header = "5000\n0\n1"
    datatowrite = np.column_stack((ener,spec))
    filename = f'nthcomp_g{g}_t{Kte}.txt'
    np.savetxt(filename,datatowrite,header=header,comments='',fmt='%.8e')