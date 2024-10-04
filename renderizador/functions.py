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

def points_in_triangle(vertices):
    # Função que recebe os vértices de um triângulo e retorna uma lista de tuplas [x,y] dos pontos que estão
    # dentro do triângulo
        pontos = []
        p0 = [vertices[0], vertices[1]]
        p1 = [vertices[2], vertices[3]]
        p2 = [vertices[4], vertices[5]]

        n0 = [(p0[1] - p2[1]), -(p0[0] - p2[0])]
        n1 = [(p1[1] - p0[1]), -(p1[0] - p0[0])]
        n2 = [(p2[1] - p1[1]), -(p2[0] - p1[0])]

        xmin = int(min([p0[0], p1[0], p2[0]]))
        xmax = int(max([p0[0], p1[0], p2[0]]))
        ymin = int(min([p0[1], p1[1], p2[1]]))
        ymax = int(max([p0[1], p1[1], p2[1]]))

        for x in range(xmin, xmax+1):
            y_lista = []
            if (p0[0] != p2[0]):
                y = int((p2[1]-p0[1])*(x-p0[0])/(p2[0]-p0[0]) + p0[1])
                if y >= ymin and y <= ymax:
                    y_lista.append(y)
            if (p0[0] != p1[0]):
                y = int((p1[1]-p0[1])*(x-p0[0])/(p1[0]-p0[0]) + p0[1])
                if y >= ymin and y <= ymax:
                    y_lista.append(y)
            if (p1[0] != p2[0]):
                y = int((p2[1]-p1[1])*(x-p1[0])/(p2[0]-p1[0]) + p1[1])
                if y >= ymin and y <= ymax:
                    y_lista.append(y)
            if len(y_lista) == 0:
                yi = ymin
                yf = ymin
            else:
                yi = min(y_lista)
                yf = max(y_lista)

            for y in range(yi, yf+1):
                r0 = (x+0.5-p2[0])*n0[0] + (y+0.5-p2[1])*n0[1]
                r1 = (x+0.5-p0[0])*n1[0] + (y+0.5-p0[1])*n1[1]
                r2 = (x+0.5-p1[0])*n2[0] + (y+0.5-p1[1])*n2[1]
                if (r0 >= 0 and r1 >= 0 and r2 >= 0):
                    pontos.append([x,y])
        return pontos

def triangle_projection(point, matrizes, width, height):
        z_points = []
        triangles = []
        for i in range(0, len(point)-2, 3):
            coords = np.array([[point[i]], [point[i+1]], [point[i+2]], [1]])
            # Triângulo posicionado e rotacionado pelo transform
            coordsTransformed = np.matmul(matrizes['transform'][-1], coords)

            # Triângulo no ponto de vista da camera projetado no NDC e normalizado
            coordsCamera = np.matmul(matrizes['viewCamera'], coordsTransformed)
            z_points.append(coordsCamera[2][0])

            coordsNDC = np.matmul(matrizes['NDC'], coordsCamera)
            coordsNormalized = coordsNDC/coordsNDC[3]

            # Triângulo mapeado para coordenadas da tela
            mTela = np.array([[width/4, 0, 0, width/4],
                               [0, -height/4, 0, height/4],
                               [0, 0, 1, 0],
                               [0, 0, 0, 1]])
            coordsTela = np.matmul(mTela, coordsNormalized)
            triangles.append(coordsTela[0][0])
            triangles.append(coordsTela[1][0])
        return triangles, z_points