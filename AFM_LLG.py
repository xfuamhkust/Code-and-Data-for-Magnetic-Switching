'''
Author: Fu Xizhi
Date: 2025/11/10
Modifications:
    1. Polarizations of FL-SOT and DL-SOT are separate. 2025/01/04
    2. Time period of FL-SOT and DL-SOT are separate. 2025/01/05
    3. Modified output format for MATLAB compatibility. 2025/11/10
    4. Modified Hamiltonian as hN * p·(N×d). 2025/11/10
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
    J     = Para["J"];     d     = Para["d"];     K     = Para["K"];     k     = Para["k"];
    pF    = Para["pF"];    pD    = Para["pD"]; 
    tF1   = Para["tF1"];   tF2   = Para["tF2"];   tD1   = Para["tD1"];   tD2   = Para["tD2"];   
    m10   = Para["m10"];   m20   = Para["m20"];
    dt    = Para["dt"];    ttF   = Para["ttF"];   ttD   = Para["ttD"];   hN    = Para["hN"];
    l   = alpha
    muS = m * muB
    k   = np.array(k)/np.linalg.norm(k)
    pF  = np.array(pF)/np.linalg.norm(pF)
    pD  = np.array(pD)/np.linalg.norm(pD)
    m10 = np.array(m10)/np.linalg.norm(m10)
    m20 = np.array(m20)/np.linalg.norm(m20)
    ttF1, ttF2, ttF3 = ttF[0], ttF[0]+ttF[1], ttF[0]+ttF[1]+ttF[2]
    ttD1, ttD2, ttD3 = ttD[0], ttD[0]+ttD[1], ttD[0]+ttD[1]+ttD[2]
    t3 = np.max([ttF3,ttD3])
    gamma1 = gamma/(1+l**2)
    d_x, d_y, d_z = d; k_x, k_y, k_z = k; pF_x, pF_y, pF_z = pF; pD_x, pD_y, pD_z = pD
    
    '''
    Effective field
    '''
    def Hxyz(m1_x, m1_y, m1_z, m2_x, m2_y, m2_z, t):
        if t>ttF1 and t<=ttF2:
            IfFLSOT = 1
        else:
            IfFLSOT = 0
        if t>ttD1 and t<=ttD2:
            IfDLSOT = 1
        else:
            IfDLSOT = 0
        m1k   = m1_x* k_x + m1_y* k_y + m1_z* k_z
        m2k   = m2_x* k_x + m2_y* k_y + m2_z* k_z
        H1_x = -1/muS*(-J*m2_x + (m2_y*d_z - m2_z*d_y) - 2*K*m1k*k_x) + IfFLSOT*tF1*pF_x + IfDLSOT*tD1*(m1_y*pD_z - m1_z*pD_y)
        H1_y = -1/muS*(-J*m2_y + (m2_z*d_x - m2_x*d_z) - 2*K*m1k*k_y) + IfFLSOT*tF1*pF_y + IfDLSOT*tD1*(m1_z*pD_x - m1_x*pD_z)
        H1_z = -1/muS*(-J*m2_z + (m2_x*d_y - m2_y*d_x) - 2*K*m1k*k_z) + IfFLSOT*tF1*pF_z + IfDLSOT*tD1*(m1_x*pD_y - m1_y*pD_x)
        H2_x = -1/muS*(-J*m1_x - (m1_y*d_z - m1_z*d_y) - 2*K*m2k*k_x) + IfFLSOT*tF2*pF_x + IfDLSOT*tD2*(m2_y*pD_z - m2_z*pD_y)
        H2_y = -1/muS*(-J*m1_y - (m1_z*d_x - m1_x*d_z) - 2*K*m2k*k_y) + IfFLSOT*tF2*pF_y + IfDLSOT*tD2*(m2_z*pD_x - m2_x*pD_z)
        H2_z = -1/muS*(-J*m1_z - (m1_x*d_y - m1_y*d_x) - 2*K*m2k*k_z) + IfFLSOT*tF2*pF_z + IfDLSOT*tD2*(m2_x*pD_y - m2_y*pD_x)
        # Modification 4: hN * p·(N×d)
        # print(abs(hN/(2*muS)*(d_x*pD_y - d_y*pD_x))/abs(H1_z))
        H1_x += -IfDLSOT*hN/(2*muS)*(d_y*pD_z - d_z*pD_y)
        H1_y += -IfDLSOT*hN/(2*muS)*(d_z*pD_x - d_x*pD_z)
        H1_z += -IfDLSOT*hN/(2*muS)*(d_x*pD_y - d_y*pD_x)
        H2_x += +IfDLSOT*hN/(2*muS)*(d_y*pD_z - d_z*pD_y)
        H2_y += +IfDLSOT*hN/(2*muS)*(d_z*pD_x - d_x*pD_z)
        H2_z += +IfDLSOT*hN/(2*muS)*(d_x*pD_y - d_y*pD_x)
        return H1_x, H1_y, H1_z, H2_x, H2_y, H2_z
    
    def dm(m1_x, m1_y, m1_z, m2_x, m2_y, m2_z, t):
        H1_x, H1_y, H1_z, H2_x, H2_y, H2_z = Hxyz(m1_x, m1_y, m1_z, m2_x, m2_y, m2_z, t)
        mH1_x = m1_y*H1_z - m1_z*H1_y
        mH1_y = m1_z*H1_x - m1_x*H1_z
        mH1_z = m1_x*H1_y - m1_y*H1_x
        mmH1_x = m1_y*mH1_z - m1_z*mH1_y
        mmH1_y = m1_z*mH1_x - m1_x*mH1_z
        mmH1_z = m1_x*mH1_y - m1_y*mH1_x
        dm1_x = -gamma1*(mH1_x + l*mmH1_x)
        dm1_y = -gamma1*(mH1_y + l*mmH1_y)
        dm1_z = -gamma1*(mH1_z + l*mmH1_z)
        mH2_x = m2_y*H2_z - m2_z*H2_y
        mH2_y = m2_z*H2_x - m2_x*H2_z
        mH2_z = m2_x*H2_y - m2_y*H2_x
        mmH2_x = m2_y*mH2_z - m2_z*mH2_y
        mmH2_y = m2_z*mH2_x - m2_x*mH2_z
        mmH2_z = m2_x*mH2_y - m2_y*mH2_x
        dm2_x = -gamma1*(mH2_x + l*mmH2_x)
        dm2_y = -gamma1*(mH2_y + l*mmH2_y)
        dm2_z = -gamma1*(mH2_z + l*mmH2_z)
        return dm1_x, dm1_y, dm1_z, dm2_x, dm2_y, dm2_z
    
    '''
    Heun method
    '''
    time_out = np.arange(t3+1)
    m1m2_out = np.zeros((t3+1,6))
    m1m2_out[0,0:3] = m10; m1m2_out[0,3:6] = m20
    m10_x, m10_y, m10_z = m10; m20_x, m20_y, m20_z = m20
    for t in range(1,t3+1):
        dm10_x, dm10_y, dm10_z, dm20_x, dm20_y, dm20_z = dm(m10_x, m10_y, m10_z, m20_x, m20_y, m20_z, t)
        m1_x, m1_y, m1_z = m10_x + dm10_x*dt, m10_y + dm10_y*dt, m10_z + dm10_z*dt
        m2_x, m2_y, m2_z = m20_x + dm20_x*dt, m20_y + dm20_y*dt, m20_z + dm20_z*dt
        m1_norm = (m1_x*m1_x+m1_y*m1_y+m1_z*m1_z)**-.5
        m2_norm = (m2_x*m2_x+m2_y*m2_y+m2_z*m2_z)**-.5
        m1_x, m1_y, m1_z = m1_x*m1_norm, m1_y*m1_norm, m1_z*m1_norm
        m2_x, m2_y, m2_z = m2_x*m2_norm, m2_y*m2_norm, m2_z*m2_norm
        dm1_x,  dm1_y,  dm1_z,  dm2_x,  dm2_y,  dm2_z  = dm(m1_x,  m1_y,  m1_z,  m2_x,  m2_y,  m2_z, t)
        m1_x, m1_y, m1_z = m10_x + (dm10_x+dm1_x)*dt*.5, m10_y + (dm10_y+dm1_y)*dt*.5, m10_z + (dm10_z+dm1_z)*dt*.5
        m2_x, m2_y, m2_z = m20_x + (dm20_x+dm2_x)*dt*.5, m20_y + (dm20_y+dm2_y)*dt*.5, m20_z + (dm20_z+dm2_z)*dt*.5
        m1_norm = (m1_x*m1_x+m1_y*m1_y+m1_z*m1_z)**-.5
        m2_norm = (m2_x*m2_x+m2_y*m2_y+m2_z*m2_z)**-.5
        m1_x, m1_y, m1_z = m1_x*m1_norm, m1_y*m1_norm, m1_z*m1_norm
        m2_x, m2_y, m2_z = m2_x*m2_norm, m2_y*m2_norm, m2_z*m2_norm
        m1m2_out[t,0:3] = [m1_x, m1_y, m1_z];  m1m2_out[t,3:6] = [m2_x, m2_y, m2_z]
        m10_x, m10_y, m10_z, m20_x, m20_y, m20_z = m1_x, m1_y, m1_z, m2_x, m2_y, m2_z
    
    return time_out,m1m2_out

if __name__ == "__main__":
    
    '''
    Parameters
    '''
    # Damping constant
    alpha = 0.1
    # Atomic spin moment (uB)
    m = 1.72
    
    # Exchange coupling energy (J)
    J = -10e-22
    # DMI (J)
    d = [1.0e-22,0,0]
    # Magnetic anisotropic energy (J)
    K = 1.0e-23
    # Magnetic anisotropic direction
    k = [0,0,1]
    # Polarizations of FL-SOT & DL-SOT
    pF = [0,0,1]
    pD = [0,1,0]
    # FL-SOT & DL-SOT (relative to Exchange coupling energy)
    tF1 = 0.0
    tF2 = 0.0
    tD1 = 0.01
    tD2 = 0.01
    muS = m * muB
    tF1 *= abs(J)/muS
    tF2 *= abs(J)/muS
    tD1 *= abs(J)/muS
    tD2 *= abs(J)/muS
    # Modification 4: hN
    hN = 0.5*tD1
    
    # Initial spin direction
    m10 = [0,0.0493241086, 0.998782825]
    m20 = [0,0.0493241086,-0.998782825]
    
    # Time step (s)
    dt = 1.0E-15
    # Time steps
    ttF = [0,   0,6000]
    ttD = [0,3000,3000]
    
    Para = {"alpha": alpha, "m":         m,
            "J":         J, "d":         d, "K":         K, "k":         k,
            "pF":       pF, "pD":       pD,
            "tF1":     tF1, "tF2":     tF2, "tD1":     tD1, "tD2":     tD2,
            "m10":     m10, "m20":     m20, 
            "dt":       dt, "ttF":     ttF, "ttD":     ttD, "hN":       hN,
            }
    t0 = time.time()
    time_out,m1m2_out = llg_heun(Para)
    n_out = (m1m2_out[:,0:3] - m1m2_out[:,3:6])/2
    t1 = time.time()
    print(f"Simulation time: {t1-t0:.2f} seconds")
    
    # Save without headers for MATLAB compatibility
    np.savetxt('magnetization_data.txt', m1m2_out)
    np.savetxt('neelvector_data.txt', n_out)
    np.savetxt('time_data.txt', time_out)
    
    print("Data saved as neelvector_data.txt and time_data.txt (no headers)")