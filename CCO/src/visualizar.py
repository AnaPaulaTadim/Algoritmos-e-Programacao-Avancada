#!/usr/bin/env python3

import sys
import os
import pyvista as pv
import numpy as np

# ============================================
# VERIFICA ARGUMENTOS DO TERMINAL
# ============================================
if len(sys.argv) < 2:
    print("Uso:")
    print("python3 visualizar.py arvore.vtk")
    sys.exit(1)

arquivo = sys.argv[1]

# ============================================
# LÊ O RAIO
# ============================================
R = 10.0

if os.path.exists("raio.txt"):
    with open("raio.txt", "r") as f:
        R = float(f.read())

if not os.path.exists(arquivo):
    print("Arquivo não encontrado:", arquivo)
    sys.exit(1)

# ============================================
# CARREGA A ÁRVORE
# ============================================
mesh = pv.read(arquivo)

# ============================================
# IDENTIFICA RAIZ E DEMAIS NÓS
# ============================================
pontos = mesh.points
num_pontos = mesh.n_points

dist_raiz = np.sqrt(
    pontos[:, 0]**2 +
    (pontos[:, 1] - R)**2
)

id_raiz = np.argmin(dist_raiz)

coordenadas_raiz = []
coordenadas_nos = []

for pid in range(num_pontos):

    if pid == id_raiz:
        coordenadas_raiz.append(pontos[pid])
    else:
        coordenadas_nos.append(pontos[pid])

# ============================================
# CRIA JANELA
# ============================================
plotter = pv.Plotter(
    window_size=[900, 900]
)

plotter.set_background("white")

# ============================================
# SEGMENTOS
# ============================================
cores_segmentos = [
    "red",
    "blue",
    "black"
]

for i in range(mesh.n_cells):

    segmento = mesh.extract_cells(i)

    cor = cores_segmentos[
        i % len(cores_segmentos)
    ]

    plotter.add_mesh(
        segmento,
        color=cor,
        line_width=5,
        render_lines_as_tubes=True
    )

# ============================================
# NÓ DE ENTRADA (RAIZ)
# ============================================
if len(coordenadas_raiz) > 0:

    plotter.add_points(
        np.array(coordenadas_raiz),
        color="orange",
        point_size=18,
        render_points_as_spheres=True
    )

# ============================================
# TODOS OS DEMAIS NÓS
# ============================================
if len(coordenadas_nos) > 0:

    plotter.add_points(
        np.array(coordenadas_nos),
        color="green",
        point_size=18,
        render_points_as_spheres=True
    )

# ============================================
# DOMÍNIO CIRCULAR
# ============================================
circulo = pv.Circle(
    radius=R,
    resolution=300
)

plotter.add_mesh(
    circulo,
    color="black",
    style="wireframe",
    line_width=3
)

# ============================================
# LEGENDA
# ============================================
plotter.add_legend(
    [
        ["Raiz", "orange"],
        ["Nós", "green"]
    ],
    bcolor="white",
    border=True,
    size=(0.11, 0.11),
    loc="upper right"
)

# ============================================
# TÍTULO
# ============================================
plotter.add_text(
    "MiniCCO - Estrutura Arterial Otimizada",
    position="upper_left",
    font_size=10,
    color="black"
)

# ============================================
# CONFIGURAÇÃO DA CÂMERA
# ============================================
plotter.enable_anti_aliasing()

plotter.camera.parallel_projection = True

dim_total = 2.0 * R

margem = dim_total * 0.05

plotter.camera_position = [
    (0, 0, 50),
    (0, 0, 0),
    (0, 1, 0)
]

plotter.camera.parallel_scale = (
    dim_total / 2.0
) + margem

plotter.view_xy()

# ============================================
# EXIBE A JANELA
# ============================================
plotter.show()









#==============================================================================================
'''#!/usr/bin/env python3

import sys
import os
import pyvista as pv

# ============================================
# VERIFICA ARGUMENTOS
# ============================================
if len(sys.argv) < 2:
    print("Uso:")
    print("python visualizar.py arvore.vtk")
    sys.exit(1)

arquivo = sys.argv[1]

# Valor padrão do raio
R = 10.0

if os.path.exists("raio.txt"):
    with open("raio.txt", "r") as f:
        R = float(f.read())

# ============================================
# VERIFICA EXISTÊNCIA DO ARQUIVO VTK
# ============================================
if not os.path.exists(arquivo):
    print("Arquivo não encontrado:", arquivo)
    sys.exit(1)

# ============================================
# LÊ O ARQUIVO VTK
# ============================================
mesh = pv.read(arquivo)

# ============================================
# CRIA VISUALIZAÇÃO
# ============================================
plotter = pv.Plotter(window_size=[900, 900])
plotter.set_background("white")

# ============================================
# CORES ALTERNADAS POR SEGMENTO
# ============================================
cores = ["red", "blue", "black"]

for i in range(mesh.n_cells):
    segmento = mesh.extract_cells(i)
    cor = cores[i % len(cores)]
    plotter.add_mesh(
        segmento,
        color=cor,
        line_width=5,
        render_lines_as_tubes=True
    )

# ============================================
# TÍTULO
# ============================================
plotter.add_text(
    "MiniCCO - Árvore Arterial Centrada",
    position="upper_left",
    font_size=9,
    color="black"
)

# ============================================
# DESENHA DOMÍNIO CIRCULAR (CENTRADO NA ORIGEM 0,0)
# ============================================
circulo = pv.Circle(
    radius=R,
    resolution=200
)

plotter.add_mesh(
    circulo,
    color="gray",
    style="wireframe",
    line_width=1
)

plotter.enable_anti_aliasing()

# ============================================
# ENQUADRAMENTO DA CÂMERA (CÍRCULO COMPLETO)
# ============================================
centro_x = 0
centro_y = 0  # Centrado na origem do domínio circular

# Dimensão total baseada no diâmetro do círculo com 5% de margem
dim_total = 2.0 * R
margem = dim_total * 0.05

plotter.camera.parallel_projection = True
plotter.camera_position = [
    (centro_x, centro_y, 50),
    (centro_x, centro_y, 0),
    (0, 1, 0)
]

plotter.camera.parallel_scale = (dim_total / 2.0) + margem
plotter.view_xy()

# ============================================
# EXIBE A JANELA
# ============================================
plotter.show()'''