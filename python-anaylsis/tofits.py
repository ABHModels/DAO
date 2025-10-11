import numpy as np
import pandas as pd
import xspec
from datetime import datetime, time
import rdata
import os 
import json
import subprocess as sp

def extract_hash_with_split(file_path: str) -> str:
    filename = os.path.basename(file_path)
    parts = filename.split('_')
    if len(parts) >= 3:
        return parts[1]
    else:
        return None

def XSPEC_PLOT(path):
    xspec.Plot.device = '/null'
    xspec.Plot.add = True
    xspec.Plot.addCommand(f'wd {path}')
    xspec.Plot('eemodel') 
    xspec.Plot.addCommand(f'exit')
    xspec.Plot.commands = ()
    xspec.AllModels.clear()

os.chdir('/data/home/yimin2/reflection/') 
METADATA_DIR = 'metadata'
if not os.path.isdir(METADATA_DIR):
    raise FileNotFoundError(f"元数据目录 '{METADATA_DIR}' 不存在。")
datas = os.listdir(METADATA_DIR)


plot_data_list = []
allfits_name = []
ener_SET= True
xspec.AllModels.lmod("relxill", '/data/models/relxill_v2.3') 
xspec.AllModels.setEnergies("0.001 1000 5000 log")

for data_filename in datas:
    with open(os.path.join(METADATA_DIR, data_filename), 'r') as f:
        try:
            data = json.load(f)
            params = data['parameters']
            if params.get('diy_tag') == 'dinci_e5000_nmu30':

                timeinfo = datetime.fromisoformat(data['run_info']['start_time_iso'])
                target_timestamp_str = '2025-08-11T01:05'
                target_dt = datetime.fromisoformat(target_timestamp_str)
                if timeinfo>target_dt:
                    ener,top,bot,inte,flux = rdata.read_spec5000(params['op_spec'])
                    # From ev to keV

                    if ener_SET:                        
                        mask  = ener>0.01
                        in_mask = ener<0.01

                        ener_hi = np.array(ener[mask])
                        num_flux = len(ener_hi)
                        ener_lo = np.zeros(num_flux)
                        ener_lo[0] = 0.009981
                        ener_lo[1:] = ener_hi[:len(ener_hi)-1]
                        delta_energy = ener_hi-ener_lo
                        new_ener = (ener_lo+ener_hi)/2.
                        ener_SET = False

                        enerlo_reshaped = ener_lo.reshape(-1, 1)
                        enerhi_reshaped = ener_hi.reshape(-1, 1)

                    # FLUX in unit of ev/cm2/s/ev
                    flux = np.array(flux/ener)
                    flux = flux[mask]
                    rescale_inte = flux*delta_energy
                    rescale_inte = rescale_inte.reshape(-1,1)
                    # Set name of fits mother files with hash
                    hashvalue = extract_hash_with_split(params['op_spec'])
                    mfile_name = os.path.join("reldao",f"{hashvalue}.dat")

                    # Save rescaled data to files 
                    combined_data = np.hstack((enerlo_reshaped, enerhi_reshaped,rescale_inte))
                    np.savetxt(mfile_name, combined_data, fmt='%.6e', delimiter='   ')
                    print(f"文件 '{mfile_name}' 已成功保存！")

                    # To fits file 
                    fitsname = os.path.join("reldao",f"{hashvalue}.fits")
                    command = f"ftflx2tab {mfile_name} dao {fitsname}"

                    try:
                        sp.check_call(command, shell=True)
                        print(f"fits file is created: {fitsname}")
                    except sp.CalledProcessError as e:
                        print(f"Error creating fits file: {e}")
                        continue

                    # load xspec model and to plot 
                    cust = f"relconv*atable{{{fitsname}}}"
                    model = xspec.Model(cust)
                    

                    # reltivistic model data path
                    relfile = os.path.join("reldao",f"rel_{hashvalue}_{params['diy_tag']}.dat")
                    XSPEC_PLOT(relfile)
                    print(f"reltivistic model [EF_E] has been save to {relfile}")
                    xspec.AllModels.clear()
                    xspec.AllData.clear()

        except (json.JSONDecodeError, KeyError, TypeError) as e:
            print(f"警告: 跳过元数据文件 {data_filename}，因为解析错误: {e}")
            continue