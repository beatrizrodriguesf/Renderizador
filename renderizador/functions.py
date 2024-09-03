import numpy as np

def quater_rotation(values):
    angulo = values[3]/2
    u = np.sqrt(values[0]**2 + values[1]**2 + values[2]**2)
    qr = np.cos(angulo)
    qx = np.sin(angulo)*(values[0]/u)
    qy = np.sin(angulo)*(values[1]/u)
    qz = np.sin(angulo)*(values[2]/u)

    matriz = np.array([[(1-2*(qy**2 + qz**2)), (2*(qx*qy - qz*qr)), (2*(qx*qz + qy*qr)), 0],
              [(2*(qx*qy + qz*qr)), (1-2*(qx**2 + qz**2)), (2*(qy*qz - qx*qr)), 0],
              [(2*(qx*qz - qy*qr)), (2*(qy*qz + qx*qr)), (1-2*(qx**2 + qy**2)), 0],
              [0, 0, 0, 1]])
    
    return matriz