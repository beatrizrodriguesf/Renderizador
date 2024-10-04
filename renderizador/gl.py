#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: <SEU NOME AQUI>
Disciplina: Computação Gráfica
Data: <DATA DE INÍCIO DA IMPLEMENTAÇÃO>
"""

import time         # Para operações com tempo
import gpu          # Simula os recursos de uma GPU
import math         # Funções matemáticas
import numpy as np  # Biblioteca do Numpy
from functions import *

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 8000   # largura da tela
    height = 6000  # altura da tela
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante
    
    matrizes = {'transform': [np.identity(4)], 'viewCamera': np.identity(4), 'NDC': np.identity(4)}

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width*2
        GL.height = height*2
        GL.near = near
        GL.far = far

    @staticmethod
    def polypoint2D(point, colors):
        """Função usada para renderizar Polypoint2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polypoint2D
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
        # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
        # você pode assumir inicialmente o desenho dos pontos com a cor emissiva (emissiveColor).

        for i in range(0, len(point), 2):
            color = []
            for value in colors['emissiveColor']:
                color.append(value*255)
            pos_x = int(point[i])
            pos_y = int(point[i+1])
            if (pos_x < GL.width and pos_y < GL.height and pos_x >= 0 and pos_y >= 0):
                gpu.GPU.draw_pixel([pos_x, pos_y], gpu.GPU.RGB8, color)
        
    @staticmethod
    def polyline2D(lineSegments, colors):
        """Função usada para renderizar Polyline2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polyline2D
        # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
        # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
        # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
        # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
        # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
        # vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        color = []
        for value in colors['emissiveColor']:
            color.append(int(value*255))

        for i in range(0, len(lineSegments)-2, 2):

            x1 = lineSegments[i]
            y1 = lineSegments[i+1]
            x2 = lineSegments[i+2]
            y2 = lineSegments[i+3]

            if x1 > x2:
                xmax = int(x1)
                xmin = int(x2)
            else:
                xmax = int(x2)
                xmin = int(x1)
            
            if (y1 > y2):
                ymax = int(y1)
                ymin = int(y2)
            else:
                ymax = int(y2)
                ymin = int(y1)

            if x1 != x2:
                m = (y2-y1)/(x2-x1)
                if abs(m) <= 1:
                    for x in range(xmin, xmax):
                        y = int(m*(x+0.5-x1) + y1)
                        if (x < GL.width and y < GL.height and x >= 0 and y >= 0):
                            gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, color)
                else:
                    m = 1/m
                    for y in range(ymin, ymax):
                        x = int(m*(y+0.5-y1) + x1)
                        if (x < GL.width and y < GL.height and x >= 0 and y >= 0):
                            gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)
            else:
                x = int(x1)
                for y in range(ymin, ymax):
                    if (x < GL.width and y < GL.height and x >= 0 and y >= 0):
                        gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, color)

    @staticmethod
    def circle2D(radius, colors):
        """Função usada para renderizar Circle2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Circle2D
        # Nessa função você receberá um valor de raio e deverá desenhar o contorno de
        # um círculo.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Circle2D
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

        color = []
        for value in colors['emissiveColor']:
            color.append(int(value*255))

        x_anterior = radius
        for y in range(0, round(radius)+1):
            x = round(np.sqrt(radius**2 - y**2))
            gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, color)

            if (x_anterior - x) > 1:
                for value in range(x+1, x_anterior):
                    gpu.GPU.draw_pixel([value,y-1], gpu.GPU.RGB8, color)
            x_anterior = x

    @staticmethod
    def triangleSet2D(vertices, colors):
        """Função usada para renderizar TriangleSet2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#TriangleSet2D
        # Nessa função você receberá os vertices de um triângulo no parâmetro vertices,
        # esses pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o
        # valor da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto.
        # Já point[2] é a coordenada x do segundo ponto e assim por diante. Assuma que a
        # quantidade de pontos é sempre multiplo de 3, ou seja, 6 valores ou 12 valores, etc.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        for i in range(len(vertices)):
            vertices[i] = vertices[i]*2
        
        color = []
        for value in colors['emissiveColor']:
            color.append(int(value*255))
        
        for i in range(0, len(vertices)-5, 6):
            pontos = points_in_triangle(vertices[i:i+6])
            for ponto in pontos:
                x = ponto[0]
                y = ponto[1]
                if (x < GL.width and y < GL.height and x >= 0 and y >= 0):
                    if colors['transparency'] != 0:
                        colorPoint = [0,0,0]
                        colorOld = gpu.GPU.read_pixel([x,y], gpu.GPU.RGB8)
                        for i in range(len(color)):
                            colorPoint[i] = color[i]*(1-colors['transparency']) + colorOld[i]*colors['transparency']
                        gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, colorPoint)
                    else:
                        gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, color)
    @staticmethod
    def triangleSet(point, colors):
        """Função usada para renderizar TriangleSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleSet
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
        # assim por diante.
        # No TriangleSet os triângulos são informados individualmente, assim os três
        # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
        # triângulo, e assim por diante.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, você pode assumir
        # inicialmente, para o TriangleSet, o desenho das linhas com a cor emissiva
        # (emissiveColor), conforme implementar novos materias você deverá suportar outros
        # tipos de cores.

        triangles = []
        for i in range(0, len(point)-2, 3):
            coords = np.array([[point[i]], [point[i+1]], [point[i+2]], [1]])
            # Triângulo posicionado e rotacionado pelo transform
            coordsTransformed = np.matmul(GL.matrizes['transform'][-1], coords)

            # Triângulo no ponto de vista da camera projetado no NDC e normalizado
            coordsCamera = np.matmul(GL.matrizes['viewCamera'], coordsTransformed)

            coordsNDC = np.matmul(GL.matrizes['NDC'], coordsCamera)
            coordsNormalized = coordsNDC/coordsNDC[3]

            # Triângulo mapeado para coordenadas da tela
            mTela = np.array([[GL.width/4, 0, 0, GL.width/4],
                               [0, -GL.height/4, 0, GL.height/4],
                               [0, 0, 1, 0],
                               [0, 0, 0, 1]])
            coordsTela = np.matmul(mTela, coordsNormalized)
            triangles.append(coordsTela[0][0])
            triangles.append(coordsTela[1][0])
        GL.triangleSet2D(triangles, colors)

    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        # inverte translation
        invTranslation = np.array([[1,0,0,-position[0]],[0,1,0,-position[1]],[0,0,1,-position[2]],[0,0,0,1]])

        # inverte Rotation
        mRotation = quater_rotation(orientation)
        invRotation = np.linalg.inv(mRotation)

        # matriz de Visualização
        mView = np.matmul(invRotation, invTranslation)
        GL.matrizes['viewCamera'] = mView

        # Aplicando fieldOfView
        fovy = 2*math.atan(math.tan(fieldOfView/2)*(GL.height/math.sqrt(GL.height**2 + GL.width**2)))
        top = GL.near*math.tan(fovy)
        right = top*(GL.width/GL.height)
        P = np.array([[GL.near/right, 0, 0, 0],
             [0, GL.near/top, 0, 0],
             [0, 0, -(GL.far + GL.near)/(GL.far - GL.near), -2*GL.far*GL.near/(GL.far - GL.near)],
             [0, 0, -1, 0]])
        
        GL.matrizes['NDC'] = P

    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo para depois potencialmente usar em outras chamadas. 
        # Quando começar a usar Transforms dentre de outros Transforms, mais a frente no curso
        # Você precisará usar alguma estrutura de dados pilha para organizar as matrizes.

        # matriz de translação
        mTranslation = np.array([[1,0,0,translation[0]],[0,1,0,translation[1]],[0,0,1,translation[2]],[0,0,0,1]])

        # matriz de Escala
        mScale = np.array([[scale[0],0,0,0], [0,scale[1],0,0], [0,0,scale[2],0], [0,0,0,1]])

        # matriz de Rotação
        mRotation = quater_rotation(rotation)
        mTransform = np.matmul(mTranslation, np.matmul(mRotation, mScale))
        
        GL.matrizes['transform'].append(np.matmul(GL.matrizes['transform'][-1], mTransform))

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        GL.matrizes['transform'].pop()

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleStripSet
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        inicio = 0
        for value in stripCount:
            for t in range(value-2):
                i = inicio + t*3
                if t % 2 == 0:
                    list = point[i:i+9]
                else:
                    p1 = point[i:i+3]
                    p2 = point[i+3:i+6]
                    p3 = point[i+6:i+9]
                    list = p1 + p3 + p2
                GL.triangleSet(list, colors)
            inicio += value*3

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#IndexedTriangleStripSet
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        i = 0
        while i < len(index):
            j = 0
            list = []
            while index[i+2] != -1:
                p1 = point[(index[i]*3):(index[i]*3 + 3)]
                p2 = point[(index[i+1]*3):(index[i+1]*3 + 3)]
                p3 = point[(index[i+2]*3):(index[i+2]*3 + 3)]
                if j % 2 == 0:
                    list = p1 + p2 + p3
                else:
                    list = p1 + p3 + p2
                GL.triangleSet(list, colors)
                j += 1
                i += 1
            i += 3

    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#IndexedFaceSet
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão não possui uma ordem oficial, mas em geral se o primeiro ponto com os dois
        # seguintes e depois este mesmo primeiro ponto com o terçeiro e quarto ponto. Por exemplo: numa
        # sequencia 0, 1, 2, 3, 4, -1 o primeiro triângulo será com os vértices 0, 1 e 2, depois serão
        # os vértices 0, 2 e 3, e depois 0, 3 e 4, e assim por diante, até chegar no final da lista.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        i = 0
        while i < len(coordIndex):
            inicio = coordIndex[i]

            while coordIndex[i+2] != -1:
                p1 = coord[(inicio*3):(inicio*3 + 3)]
                p2 = coord[(coordIndex[i+1]*3):(coordIndex[i+1]*3 + 3)]
                p3 = coord[(coordIndex[i+2]*3):(coordIndex[i+2]*3 + 3)]
                list = p1 + p2 + p3

                vertices, z_points, z_points_NDC = triangle_projection(list, GL.matrizes, GL.width, GL.height)

                for j in range(len(vertices)):
                    vertices[j] = vertices[j]*2
                
                points = points_in_triangle(vertices)
                p0 = [vertices[0], vertices[1]]
                p1 = [vertices[2], vertices[3]]
                p2 = [vertices[4], vertices[5]]
                
                for k in range(len(z_points)):
                    z_points[k] = z_points[k]*2

                if colorPerVertex and color:
                    color0 = color[(inicio*3):(inicio*3) + 3]
                    color1 = color[(coordIndex[i+1]*3):(coordIndex[i+1]*3 + 3)]
                    color2 = color[(coordIndex[i+2]*3):(coordIndex[i+2]*3 + 3)]

                    for point in points:
                        x = point[0]
                        y = point[1]
                        a = (-(x+0.5-p1[0])*(p2[1]-p1[1]) + (y+0.5-p1[1])*(p2[0]-p1[0]))/(-(p0[0]-p1[0])*(p2[1]-p1[1]) + (p0[1]-p1[1])*(p2[0]-p1[0]))
                        b = (-(x+0.5-p2[0])*(p0[1]-p2[1]) + (y+0.5-p2[1])*(p0[0]-p2[0]))/(-(p1[0]-p2[0])*(p0[1]-p2[1]) + (p1[1]-p2[1])*(p0[0]-p2[0]))
                        c = 1-a-b
                        z = 1/((a/z_points[0])+(b/z_points[1])+(c/z_points[2]))
                        colorR = z*(color0[0]*(a/z_points[0]) + color1[0]*(b/z_points[1]) + color2[0]*(c/z_points[2]))
                        colorG = z*(color0[1]*(a/z_points[0]) + color1[1]*(b/z_points[1]) + color2[1]*(c/z_points[2]))
                        colorB = z*(color0[2]*(a/z_points[0]) + color1[2]*(b/z_points[1]) + color2[2]*(c/z_points[2]))
                        colorPoint = [int(colorR*255), int(colorG*255), int(colorB*255)]

                        if (x < GL.width and y < GL.height and x >= 0 and y >= 0):
                            gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, colorPoint)

                if texCoord and current_texture:
                    tex0 = texCoord[(inicio*2):(inicio*2) + 2]
                    tex1 = texCoord[(coordIndex[i+1]*2):(coordIndex[i+1]*2 + 2)]
                    tex2 = texCoord[(coordIndex[i+2]*2):(coordIndex[i+2]*2 + 2)]

                    image = gpu.GPU.load_texture(current_texture[0])

                    x = points[0][0]
                    y = points[0][1]

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

                    dudx = (u10 - u00)*image.shape[0]
                    dudy = (u01 - u00)*image.shape[0]
                    dvdx = (v10 - v00)*image.shape[1]
                    dvdy = (v01 - v00)*image.shape[1]

                    L = max(np.sqrt(dudx**2 + dvdx**2), np.sqrt(dudy**2 + dvdy**2))
                    D = int(math.log2(L))
                    
                    for nivel in range(D):
                        new_image = nivel_image(image)
                        image = new_image
                    
                    shape = [len(image), len(image[0])]
                    
                    for point in points:
                        x = point[0]
                        y = point[1]
                        a = (-(x+0.5-p1[0])*(p2[1]-p1[1]) + (y+0.5-p1[1])*(p2[0]-p1[0]))/(-(p0[0]-p1[0])*(p2[1]-p1[1]) + (p0[1]-p1[1])*(p2[0]-p1[0]))
                        b = (-(x+0.5-p2[0])*(p0[1]-p2[1]) + (y+0.5-p2[1])*(p0[0]-p2[0]))/(-(p1[0]-p2[0])*(p0[1]-p2[1]) + (p1[1]-p2[1])*(p0[0]-p2[0]))
                        c = 1-a-b
                        z = 1/((a/z_points[0])+(b/z_points[1])+(c/z_points[2]))
                        u = z*(tex0[0]*(a/z_points[0]) + tex1[0]*(b/z_points[1]) + tex2[0]*(c/z_points[2]))
                        v = z*(tex0[1]*(a/z_points[0]) + tex1[1]*(b/z_points[1]) + tex2[1]*(c/z_points[2]))

                        colorPoint = image[int(u*shape[0])][int(shape[1]-v*shape[1])][0:3]
                        if (x < GL.width and y < GL.height and x >= 0 and y >= 0):
                            gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, colorPoint)
                else:
                    for point in points:
                        x = point[0]
                        y = point[1]
                        a = (-(x+0.5-p1[0])*(p2[1]-p1[1]) + (y+0.5-p1[1])*(p2[0]-p1[0]))/(-(p0[0]-p1[0])*(p2[1]-p1[1]) + (p0[1]-p1[1])*(p2[0]-p1[0]))
                        b = (-(x+0.5-p2[0])*(p0[1]-p2[1]) + (y+0.5-p2[1])*(p0[0]-p2[0]))/(-(p1[0]-p2[0])*(p0[1]-p2[1]) + (p1[1]-p2[1])*(p0[0]-p2[0]))
                        c = 1-a-b
                        z = 1/((a/z_points_NDC[0])+(b/z_points_NDC[1])+(c/z_points_NDC[2]))

                        colorPoint = []
                        for value in colors['emissiveColor']:
                            colorPoint.append(int(value*255))

                        if (x < GL.width and y < GL.height and x >= 0 and y >= 0):
                            gpu.GPU.read_framebuffer = 2
                            gpu.GPU.draw_framebuffer = 2
                            if gpu.GPU.read_pixel([x,y], gpu.GPU.DEPTH_COMPONENT32F) > z:
                                gpu.GPU.draw_pixel([x,y], gpu.GPU.DEPTH_COMPONENT32F, [z])
                                gpu.GPU.draw_framebuffer = 0
                                gpu.GPU.draw_pixel([x,y], gpu.GPU.RGB8, colorPoint)
                            gpu.GPU.draw_framebuffer = 0
                            gpu.GPU.read_framebuffer = 0
                i += 1
            i += 3

    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Box
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Box : size = {0}".format(size)) # imprime no terminal pontos
        print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixel([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Sphere
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Sphere : radius = {0}".format(radius)) # imprime no terminal o raio da esfera
        print("Sphere : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cone(bottomRadius, height, colors):
        """Função usada para renderizar Cones."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cone
        # A função cone é usada para desenhar cones na cena. O cone é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento bottomRadius especifica o
        # raio da base do cone e o argumento height especifica a altura do cone.
        # O cone é alinhado com o eixo Y local. O cone é fechado por padrão na base.
        # Para desenha esse cone você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cone : bottomRadius = {0}".format(bottomRadius)) # imprime no terminal o raio da base do cone
        print("Cone : height = {0}".format(height)) # imprime no terminal a altura do cone
        print("Cone : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cylinder(radius, height, colors):
        """Função usada para renderizar Cilindros."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cylinder
        # A função cylinder é usada para desenhar cilindros na cena. O cilindro é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da base do cilindro e o argumento height especifica a altura do cilindro.
        # O cilindro é alinhado com o eixo Y local. O cilindro é fechado por padrão em ambas as extremidades.
        # Para desenha esse cilindro você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cylinder : radius = {0}".format(radius)) # imprime no terminal o raio do cilindro
        print("Cylinder : height = {0}".format(height)) # imprime no terminal a altura do cilindro
        print("Cylinder : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/navigation.html#NavigationInfo
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#DirectionalLight
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#PointLight
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/environmentalEffects.html#Fog
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/time.html#TimeSensor
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nó TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#SplinePositionInterpolator
        # Interpola não linearmente entre uma lista de vetores 3D. O campo keyValue possui
        # uma lista com os valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantos vetores 3D quanto os
        # quadros-chave no key. O campo closed especifica se o interpolador deve tratar a malha
        # como fechada, com uma transições da última chave para a primeira chave. Se os keyValues
        # na primeira e na última chave não forem idênticos, o campo closed será ignorado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0.0, 0.0, 0.0]
        
        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#OrientationInterpolator
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
