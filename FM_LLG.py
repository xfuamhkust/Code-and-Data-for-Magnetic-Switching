'''
Author: Fu Xizhi
Date: 2025/11/10
Modifications:
    1. Polarizations of FL-SOT and DL-SOT are separate. 2025/01/04
    2. Time period of FL-SOT and DL-SOT are separate. 2025/01/05
    3. Modified output format for MATLAB compatibility. 2025/11/10
'''

# import sys
import time
import numpy as np
gamma = 1.76085963023e11
muB   = 9.27400999e-24

def llg_heun(Para):
    '''
    Parameters
    '''
    alpha = Para["alpha"]; m     = Para["m"]
    K     = Para["K"];     k     = Para["k"];
    pF    = Para["pF"];    pD    = Para["pD"]; 
    tF    = Para["tF"];    tD    = Para["tD"];
    m0    = Para["m0"];    dt    = Para["dt"];
    ttF   = Para["ttF"];   ttD   = Para["ttD"]; 
    l   = alpha
    muS = m * muB
    k   = np.array(k)/np.linalg.norm(k)
    pF  = np.array(pF)/np.linalg.norm(pF)
    pD  = np.array(pD)/np.linalg.norm(pD)
    m0  = np.array(m0)/np.linalg.norm(m0)
    ttF1, ttF2, ttF3 = ttF[0], ttF[0]+ttF[1], ttF[0]+ttF[1]+ttF[2]
    ttD1, ttD2, ttD3 = ttD[0], ttD[0]+ttD[1], ttD[0]+ttD[1]+ttD[2]
    t3 = np.max([ttF3,ttD3])
    gamma1 = gamma/(1+l**2)
    k_x, k_y, k_z = k; pF_x, pF_y, pF_z = pF; pD_x, pD_y, pD_z = pD
    
    '''
    Effective field
    '''
    def Hxyz(m_x, m_y, m_z, t):
        if t>ttF1 and t<=ttF2:
            IfFLSOT = 1
        else:
            IfFLSOT = 0
        if t>ttD1 and t<=ttD2:
            IfDLSOT = 1
        else:
            IfDLSOT = 0
        mk   = m_x* k_x + m_y* k_y + m_z* k_z
        H_x = -1/muS*(-2*K*mk*k_x) + IfFLSOT*tF*pF_x + IfDLSOT*tD*(m_y*pD_z - m_z*pD_y)
        H_y = -1/muS*(-2*K*mk*k_y) + IfFLSOT*tF*pF_y + IfDLSOT*tD*(m_z*pD_x - m_x*pD_z)
        H_z = -1/muS*(-2*K*mk*k_z) + IfFLSOT*tF*pF_z + IfDLSOT*tD*(m_x*pD_y - m_y*pD_x)
        return H_x, H_y, H_z
    
    def dm(m_x, m_y, m_z, t):
        H_x, H_y, H_z = Hxyz(m_x, m_y, m_z, t)
        mH_x = m_y*H_z - m_z*H_y
        mH_y = m_z*H_x - m_x*H_z
        mH_z = m_x*H_y - m_y*H_x
        mmH_x = m_y*mH_z - m_z*mH_y
        mmH_y = m_z*mH_x - m_x*mH_z
        mmH_z = m_x*mH_y - m_y*mH_x
        dm_x = -gamma1*(mH_x + l*mmH_x)
        dm_y = -gamma1*(mH_y + l*mmH_y)
        dm_z = -gamma1*(mH_z + l*mmH_z)
        return dm_x, dm_y, dm_z
    
    '''
    Heun method
    '''
    time_out = np.arange(t3+1)
    m_out = np.zeros((t3+1,3))
    m_out[0] = m0
    m0_x, m0_y, m0_z = m0
    for t in range(1,t3+1):
        dm0_x, dm0_y, dm0_z = dm(m0_x, m0_y, m0_z, t)
        m_x, m_y, m_z = m0_x + dm0_x*dt, m0_y + dm0_y*dt, m0_z + dm0_z*dt
        m_norm = (m_x*m_x+m_y*m_y+m_z*m_z)**-.5
        m_x, m_y, m_z = m_x*m_norm, m_y*m_norm, m_z*m_norm
        dm_x,  dm_y,  dm_z = dm(m_x,  m_y,  m_z, t)
        m_x, m_y, m_z = m0_x + (dm0_x+dm_x)*dt*.5, m0_y + (dm0_y+dm_y)*dt*.5, m0_z + (dm0_z+dm_z)*dt*.5
        m_norm = (m_x*m_x+m_y*m_y+m_z*m_z)**-.5
        m_x, m_y, m_z = m_x*m_norm, m_y*m_norm, m_z*m_norm
        m_out[t] = [m_x, m_y, m_z]
        m0_x, m0_y, m0_z = m_x, m_y, m_z
    
    return time_out,m_out

if __name__ == "__main__":
    
    '''
    Parameters
    '''
    # Damping constant
    alpha = 0.1
    # Atomic spin moment (uB)
    m = 1.72
    
    # Magnetic anisotropic energy (J)
    K = 1.0e-23
    # Magnetic anisotropic direction
    k = [0,0,1]
    # Polarizations of FL-SOT & DL-SOT
    pF = [0,0.02,1]
    pD = [0,0.02,1]
    # FL-SOT & DL-SOT (relative to Exchange coupling energy)
    tF = -1.0
    tD = 0.0
    muS = m * muB
    tF *= 10e-22/muS
    tD *= 10e-22/muS
    
    # Initial spin direction
    m0 = [ 0,0,1]
    
    # Time step (s)
    dt = 1.0E-15
    # Time steps
    ttF = [0,8000,0]
    ttD = [0,8000,0]
    
    Para = {"alpha": alpha, "m":         m,
            "K":         K, "k":         k,
            "pF":       pF, "pD":       pD,
            "tF":       tF, "tD":       tD,
            "m0":       m0, "dt":       dt, 
            "ttF":     ttF, "ttD":     ttD,
            }
    t0 = time.time()
    time_out,m_out = llg_heun(Para)
    t1 = time.time()
    print(f"Simulation time: {t1-t0:.2f} seconds")
    
    # Save without headers for MATLAB compatibility
    np.savetxt('magnetization_data.txt', m_out)
    np.savetxt('time_data.txt', time_out)
    
    print("Data saved as magnetization_data.txt and time_data.txt (no headers)")