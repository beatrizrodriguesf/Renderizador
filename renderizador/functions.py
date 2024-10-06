import numpy as np
import math

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
        ymin = min([p0[1], p1[1], p2[1]])
        ymax = max([p0[1], p1[1], p2[1]])

        for x in range(xmin, xmax+1):
            y_lista = []
            if (p0[0] != p2[0]):
                y = (p2[1]-p0[1])*(x-p0[0])/(p2[0]-p0[0]) + p0[1]
                if y >= ymin and y <= ymax:
                    y_lista.append(y)
            if (p0[0] != p1[0]):
                y = (p1[1]-p0[1])*(x-p0[0])/(p1[0]-p0[0]) + p0[1]
                if y >= ymin and y <= ymax:
                    y_lista.append(y)
            if (p1[0] != p2[0]):
                y = (p2[1]-p1[1])*(x-p1[0])/(p2[0]-p1[0]) + p1[1]
                if y >= ymin and y <= ymax:
                    y_lista.append(y)
            if len(y_lista) == 0:
                yi = int(ymin)
                yf = int(ymin)
            else:
                yi = int(min(y_lista))
                yf = int(max(y_lista))

            for y in range(yi, yf+2):
                r0 = (x+0.5-p2[0])*n0[0] + (y+0.5-p2[1])*n0[1]
                r1 = (x+0.5-p0[0])*n1[0] + (y+0.5-p0[1])*n1[1]
                r2 = (x+0.5-p1[0])*n2[0] + (y+0.5-p1[1])*n2[1]
                if (r0 >= 0 and r1 >= 0 and r2 >= 0):
                    pontos.append([x,y])
        return pontos

def triangle_projection(point, matrizes, width, height):
        z_points = []
        z_points_NDC = []
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
            z_points_NDC.append(coordsTela[2][0])
        return triangles, z_points, z_points_NDC

def nivel_image(image):
    new_image = []
    shape = [len(image), len(image[0])]
    for i in range(0,shape[0]-1, 2):
        linha = []
        for j in range(0,shape[1]-1,2):
            color1 = image[i][j]
            color2 = image[i+1][j]
            color3 = image[i][j+1]
            color4 = image[i+1][j+1]
            R = (int(color1[0]) + int(color2[0]) + int(color3[0]) + int(color4[0]))/4
            G = (int(color1[1]) + int(color2[1]) + int(color3[1]) + int(color4[1]))/4
            B = (int(color1[2]) + int(color2[2]) + int(color3[2]) + int(color4[2]))/4
            linha.append([R,G,B])
        new_image.append(linha)
    return new_image

def calcula_tex(point, vertices, z_points, tex, image_shape):
    p0 = [vertices[0], vertices[1]]
    p1 = [vertices[2], vertices[3]]
    p2 = [vertices[4], vertices[5]]

    tex0 = tex[0]
    tex1 = tex[1]
    tex2 = tex[2]

    x = point[0]
    y = point[1]

    a00 = (-(x+0.5-p1[0])*(p2[1]-p1[1]) + (y+0.5-p1[1])*(p2[0]-p1[0]))/(-(p0[0]-p1[0])*(p2[1]-p1[1]) + (p0[1]-p1[1])*(p2[0]-p1[0]))
    b00 = (-(x+0.5-p2[0])*(p0[1]-p2[1]) + (y+0.5-p2[1])*(p0[0]-p2[0]))/(-(p1[0]-p2[0])*(p0[1]-p2[1]) + (p1[1]-p2[1])*(p0[0]-p2[0]))
    c00 = 1-a00-b00
    z00 = 1/((a00/z_points[0])+(b00/z_points[1])+(c00/z_points[2]))
    u00 = z00*(tex0[0]*(a00/z_points[0]) + tex1[0]*(b00/z_points[1]) + tex2[0]*(c00/z_points[2]))
    v00 = z00*(tex0[1]*(a00/z_points[0]) + tex1[1]*(b00/z_points[1]) + tex2[1]*(c00/z_points[2]))

    a10 = (-(x+1+0.5-p1[0])*(p2[1]-p1[1]) + (y+0.5-p1[1])*(p2[0]-p1[0]))/(-(p0[0]-p1[0])*(p2[1]-p1[1]) + (p0[1]-p1[1])*(p2[0]-p1[0]))
    b10 = (-(x+1+0.5-p2[0])*(p0[1]-p2[1]) + (y+0.5-p2[1])*(p0[0]-p2[0]))/(-(p1[0]-p2[0])*(p0[1]-p2[1]) + (p1[1]-p2[1])*(p0[0]-p2[0]))
    c10 = 1-a10-b10
    z10 = 1/((a10/z_points[0])+(b10/z_points[1])+(c10/z_points[2]))
    u10 = z10*(tex0[0]*(a10/z_points[0]) + tex1[0]*(b10/z_points[1]) + tex2[0]*(c10/z_points[2]))
    v10 = z10*(tex0[1]*(a10/z_points[0]) + tex1[1]*(b10/z_points[1]) + tex2[1]*(c10/z_points[2]))

    a01 = (-(x+0.5-p1[0])*(p2[1]-p1[1]) + (y+1+0.5-p1[1])*(p2[0]-p1[0]))/(-(p0[0]-p1[0])*(p2[1]-p1[1]) + (p0[1]-p1[1])*(p2[0]-p1[0]))
    b01 = (-(x+0.5-p2[0])*(p0[1]-p2[1]) + (y+1+0.5-p2[1])*(p0[0]-p2[0]))/(-(p1[0]-p2[0])*(p0[1]-p2[1]) + (p1[1]-p2[1])*(p0[0]-p2[0]))
    c01 = 1-a01-b01
    z01 = 1/((a01/z_points[0])+(b01/z_points[1])+(c01/z_points[2]))
    u01 = z01*(tex0[0]*(a01/z_points[0]) + tex1[0]*(b01/z_points[1]) + tex2[0]*(c01/z_points[2]))
    v01 = z01*(tex0[1]*(a01/z_points[0]) + tex1[1]*(b01/z_points[1]) + tex2[1]*(c01/z_points[2]))

    dudx = (u10 - u00)*image_shape[0]
    dudy = (u01 - u00)*image_shape[0]
    dvdx = (v10 - v00)*image_shape[1]
    dvdy = (v01 - v00)*image_shape[1]

    L = max(np.sqrt(dudx**2 + dvdx**2), np.sqrt(dudy**2 + dvdy**2))
    D = int(math.log2(L))

    return u00, v00, D