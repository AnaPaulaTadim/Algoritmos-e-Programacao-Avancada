#!/usr/bin/env python3

import sys
import os
import pyvista as pv
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

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
# IDENTIFICA APENAS A RAIZ
# ============================================
pontos = mesh.points
num_pontos = mesh.n_points

dist_raiz = np.sqrt(
    pontos[:, 0]**2 +
    (pontos[:, 1] - R)**2
)

id_raiz = np.argmin(dist_raiz)

coordenadas_raiz = []

for pid in range(num_pontos):
    if pid == id_raiz:
        coordenadas_raiz.append(pontos[pid])

# ============================================
# CRIA JANELA
# ============================================
plotter = pv.Plotter(
    window_size=[900, 900]
)

plotter.set_background("white")

# ============================================
# SEGMENTOS (tubos proporcionais ao raio, coloridos pelo FLUXO)
# ============================================

# Le os raios de cada segmento exportados pelo C++ em CELL_DATA (Parte L)
if "raio" in mesh.cell_data:
    raios = mesh.cell_data["raio"]
else:
    # Fallback caso o VTK ainda nao tenha o campo raio
    raios = np.full(mesh.n_cells, 0.05 * R)

# --------------------------------------------------------------------
# FATOR DE ESCALA VISUAL
# --------------------------------------------------------------------
raio_maximo_fisico = float(np.max(raios)) if raios.size > 0 else 0.0
raio_maximo_desejado_na_tela = 0.025 * R  # tronco ocupa ~2.5% do raio do dominio

if raio_maximo_fisico > 0.0:
    fator_escala_visual = raio_maximo_desejado_na_tela / raio_maximo_fisico
else:
    fator_escala_visual = 1.0

# --------------------------------------------------------------------
# COR POR FLUXO (Alterado para o espectro de Vermelhos/Vinho)
# --------------------------------------------------------------------
raio_min = float(np.min(raios)) if raios.size > 0 else 0.0
raio_max = float(np.max(raios)) if raios.size > 0 else 1.0

if raio_max > raio_min:
    normalizador = mcolors.Normalize(vmin=raio_min, vmax=raio_max)
else:
    normalizador = mcolors.Normalize(vmin=0.0, vmax=1.0)

# 'Reds' gera o gradiente perfeito do vermelho claro (fino) ao vermelho vinho/escuro (grosso)
mapa_de_cor = cm.get_cmap("Reds")

for i in range(mesh.n_cells):
    segmento_bruto = mesh.extract_cells(i)
    p0 = segmento_bruto.points[0]
    p1 = segmento_bruto.points[-1]

    raio_segmento = float(raios[i]) * fator_escala_visual

    # Protege contra raio invalido (0, negativo ou NaN)
    if not np.isfinite(raio_segmento) or raio_segmento <= 0.0:
        raio_segmento = 0.005 * R

    cor_segmento = mapa_de_cor(normalizador(raios[i]))[:3]  # (R,G,B) 0..1

    linha_reta = pv.Line(p0, p1)
    linha = linha_reta.tube(radius=raio_segmento)

    plotter.add_mesh(
        linha,
        color=cor_segmento
    )

# ============================================
# NÓ DE ENTRADA (RAIZ)
# ============================================
if len(coordenadas_raiz) > 0:
    plotter.add_points(
        np.array(coordenadas_raiz),
        color="#1c4fd6",  # Azul escuro marcante para destacar a entrada na raiz
        point_size=18,
        render_points_as_spheres=True
    )

# ============================================
# DOMÍNIO CIRCULAR (MANTIDO)
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
# LEGENDA (ATUALIZADA SEM OS NÓS VERDES)
# ============================================
plotter.add_legend(
    [
        ["Raiz", "#1c4fd6"],
        ["Maior fluxo (tronco)", "maroon"],   # Vermelho vinho / escuro
        ["Menor fluxo (terminais)", "red"]     # Vermelho mais claro
    ],
    bcolor="white",
    border=True,
    size=(0.16, 0.11),
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