import numpy as np
MMx = 1.0  
MMa = 1.0  
MMo = 1.0  
MMs = 1.0  
MMCO2 = 1.0  
Y_EG = 1.0  
Y_XG = 1.0  
Y_OG = 1.0  
Yo_EG = 1.0  
Yf_EG = 1.0  
Y_AG = 1.0  
Y_EA = 1.0  
Y_XA = 1.0  
Y_OA = 1.0  
Yo_EA = 1.0  
qsemol = 10.0  
qO2p = 10.0  
qr1max = 10.0  
qAcp = 10.0  

modes = np.zeros(5)

stoichiometric_matrix = np.array([
    [1, 0, 0, 0, 0],  # Mode A: Maintenance
    [1, 1, 0, 0, 0],  # Mode B: Oxidative growth on glucose
    [0, 1, 0, 0, 0],  # Mode C: Fermentative growth on glucose
    [0, 0, 1, 0, 0],  # Mode D: Overflow
    [0, 0, 0, 1, 0]   # Mode E: Oxidative growth on acetate
])

qmode = np.zeros(5)
qmode[0] = min(qsemol, qO2p / 6)
qsemol -= qmode[0]
qO2p -= qmode[0]
qmode[1] = min((Y_EG + Yo_EG) / Yo_EG * qr1max, qsemol,
               (Y_EG + Yo_EG) / (Y_EG * Y_OG) * qO2p)
qsr = qsemol - qmode[1]
qO2r = qO2p - (Y_EG * Y_OG) / (Y_EG + Yo_EG) * qmode[1]
vpr = (Y_EG + Yo_EG) / Yo_EG * qr1max - qmode[1]
qmode[2] = min((Y_EG + Yf_EG) / Yf_EG * vpr, qsr)
qsr -= qmode[2]
qmode[3] = min(qsr, qO2r)
qmode[4] = min((Y_EA + Yo_EA) / Yo_EA * vpr, qAcp,
               (Y_EA + Yo_EA) / (Yo_EA * Y_OA) * qO2r)
qO2r -= (Yo_EA * Y_OA) / (Y_EA + Yo_EA) * qmode[4]
mu = (Y_XG * Yo_EG / (Y_EG + Yo_EG) * qmode[1] +
      Y_XG * Yf_EG / (Y_EG + Yf_EG) * qmode[2] +
      Y_XA * Yo_EA / (Y_EA + Yo_EA) * qmode[4]) * MMx

qO2e = -(Y_EG * Y_OG / (Y_EG + Yo_EG) * qmode[1] +
          Y_EA * Y_OA / (Y_EA + Yo_EA) * qmode[4] +
          qmode[3]) * MMo

qAce = (Y_AG * Y_EG / (Y_EG + Yf_EG) * qmode[2] +
        Y_AG * qmode[3] - qmode[4]) * MMa

qsegr = -qsemol * MMs
qCO2e = -(qO2e / MMo) * MMCO2