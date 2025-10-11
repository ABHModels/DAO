# read data 
import pandas as pd
ERGSEV = 1.60217662e-12 
def read_spec(path):
    try:
        data = pd.read_csv(path, sep=r'\s+', header=None, comment='#').iloc[-1000:]
        ener = pd.to_numeric(data[0], errors='coerce')
        incitop = pd.to_numeric(data[1], errors='coerce') # [erg /cm2/s/erg/ster]
        incibot = pd.to_numeric(data[2], errors='coerce') # isotropic
        inte = pd.to_numeric(data[3], errors='coerce')  # intensity at incident angle
        flux = pd.to_numeric(data[4], errors='coerce')
        
        # drop nan 
        valid_indices = ener.notna() & flux.notna()
        ener = ener[valid_indices]
        flux = flux[valid_indices]
        
        # Convert to keV
        ener_kev = ener / 1e3

        # To EEMODEL 
        top = incitop * ener_kev * ERGSEV
        bot = incibot * ener_kev * ERGSEV
        intenisty = inte * ener_kev * ERGSEV
        emergent_flux = flux * ener_kev * ERGSEV
        
        return ener_kev, top,bot,intenisty,emergent_flux
    except FileNotFoundError:
        print(f"警告: 找不到数据文件 {path}")
        return None, None

def read_spec2(path):
    try:
        data = pd.read_csv(path, sep=r'\s+', header=None, comment='#').iloc[1000:2000]
        ener = pd.to_numeric(data[0], errors='coerce')
        incitop = pd.to_numeric(data[1], errors='coerce') # [erg /cm2/s/erg/ster]
        incibot = pd.to_numeric(data[2], errors='coerce') # isotropic
        inte = pd.to_numeric(data[3], errors='coerce')  # intensity at incident angle
        flux = pd.to_numeric(data[4], errors='coerce')
        
        # drop nan 
        valid_indices = ener.notna() & flux.notna()
        ener = ener[valid_indices]
        flux = flux[valid_indices]
        
        # Convert to keV
        ener_kev = ener / 1e3

        # To EEMODEL 
        top = incitop * ener_kev * ERGSEV
        bot = incibot * ener_kev * ERGSEV
        intenisty = inte * ener_kev * ERGSEV
        emergent_flux = flux * ener_kev * ERGSEV
        
        return ener_kev, top,bot,intenisty,emergent_flux
    except FileNotFoundError:
        print(f"警告: 找不到数据文件 {path}")
        return None, None


# --- 数据读取函数 ---
def read_temp(path):
    try:
        data = pd.read_csv(path, sep=r'\s+', header=None, comment='#').iloc[-200:]
        ener = pd.to_numeric(data[0], errors='coerce')
        temp = pd.to_numeric(data[1], errors='coerce')
        
        return ener, temp
    except FileNotFoundError:
        print(f"警告: 找不到数据文件 {path}")
        return None, None
def read_spec5000(path):
    try:
        data = pd.read_csv(path, sep=r'\s+', header=None, comment='#').iloc[-5000:]
        ener = pd.to_numeric(data[0], errors='coerce')
        incitop = pd.to_numeric(data[1], errors='coerce') # [erg /cm2/s/erg/ster]
        incibot = pd.to_numeric(data[2], errors='coerce') # isotropic
        inte = pd.to_numeric(data[3], errors='coerce')  # intensity at incident angle
        flux = pd.to_numeric(data[4], errors='coerce')
        
        # drop nan 
        valid_indices = ener.notna() & flux.notna()
        ener = ener[valid_indices]
        flux = flux[valid_indices]
        
        # Convert to keV
        ener_kev = ener / 1e3

        # To EEMODEL 
        top = incitop * ener_kev * ERGSEV
        bot = incibot * ener_kev * ERGSEV
        intenisty = inte * ener_kev * ERGSEV
        emergent_flux = flux * ener_kev * ERGSEV
        
        return ener_kev, top,bot,intenisty,emergent_flux
    except FileNotFoundError:
        print(f"警告: 找不到数据文件 {path}")
        return None, None
def read_spec333(path):
    try:
        data = pd.read_csv(path, sep=r'\s+', header=None, comment='#').iloc[:1000]
        ener = pd.to_numeric(data[0], errors='coerce')
        incitop = pd.to_numeric(data[1], errors='coerce') # [erg /cm2/s/erg/ster]
        incibot = pd.to_numeric(data[2], errors='coerce') # isotropic
        inte = pd.to_numeric(data[3], errors='coerce')  # intensity at incident angle
        flux = pd.to_numeric(data[4], errors='coerce')
        
        # drop nan 
        valid_indices = ener.notna() & flux.notna()
        ener = ener[valid_indices]
        flux = flux[valid_indices]
        
        # Convert to keV
        ener_kev = ener / 1e3

        # To EEMODEL 
        top = incitop * ener_kev * ERGSEV
        bot = incibot * ener_kev * ERGSEV
        intenisty = inte * ener_kev * ERGSEV
        emergent_flux = flux * ener_kev * ERGSEV
        
        return ener_kev, top,bot,intenisty,emergent_flux
    except FileNotFoundError:
        print(f"警告: 找不到数据文件 {path}")
        return None, None