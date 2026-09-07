#!/usr/bin/env python3
"""
Visualizador 3D para Árvores Vasculares (CCO-TRICCO) utilizando PyVista.
Permite visualizar a árvore em formato de tubos 3D e colorir por:
  - raio (padrão)
  - fluxo
  - profundidade
"""

import sys
import os
from collections import Counter
import numpy as np
import pyvista as pv
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# ============================================================
# 1. VALIDAÇÃO DE ARGUMENTOS DA LINHA DE COMANDO
# ============================================================
if len(sys.argv) < 2:
    print("Uso:")
    print("  python3 visualizador.py arvore.vtk [campo_cor]")
    print("  campo_cor = raio (padrão) | fluxo | profundidade")
    sys.exit(1)

arquivo_vtk = sys.argv[1]
campo_cor = sys.argv[2].lower() if len(sys.argv) >= 3 else "raio"

if campo_cor not in ("raio", "fluxo", "profundidade"):
    print(f"[Erro] Campo de cor inválido: '{campo_cor}'. Use: raio, fluxo ou profundidade.")
    sys.exit(1)

if not os.path.exists(arquivo_vtk):
    print(f"[Erro] Arquivo VTK não encontrado: {arquivo_vtk}")
    sys.exit(1)

# ============================================================
# 2. CARREGAMENTO DO RAIO DO DOMÍNIO E MALHA VTK
# ============================================================
R = 10.0
if os.path.exists("raio.txt"):
    with open("raio.txt", "r") as f:
        try:
            R = float(f.read().strip())
        except ValueError:
            print("[Aviso] Erro ao ler raio.txt. Usando R = 10.0 como padrão.")

mesh = pv.read(arquivo_vtk)
pontos = mesh.points
num_pontos = mesh.n_points

# ============================================================
# 3. IDENTIFICAÇÃO TOPOLÓGICA (RAIZ, TERMINAIS E MULTIFURCAÇÕES)
# ============================================================
# A. Raiz: Ponto mais próximo da borda superior (0, R)
dist_raiz = np.sqrt(pontos[:, 0]**2 + (pontos[:, 1] - R)**2)
id_raiz = np.argmin(dist_raiz)
coordenadas_raiz = [pontos[id_raiz]]

# B. Posições Proximais (a) e Distais (b) dos segmentos
pontos_a = pontos[0::2]
pontos_b = pontos[1::2]

# C. Terminais (Folhas): Pontos 'b' que nunca aparecem como início 'a' de outro vaso
chaves_a = set((round(p[0], 6), round(p[1], 6)) for p in pontos_a)
coordenadas_terminais = [
    p for p in pontos_b
    if (round(p[0], 6), round(p[1], 6)) not in chaves_a
]

# D. Multifurcações (TRICCO): Pontos 'a' com 3 ou mais conexões de saída
contagem_a = Counter((round(p[0], 6), round(p[1], 6)) for p in pontos_a)
chaves_multifurcacao = {chave for chave, qtd in contagem_a.items() if qtd >= 3}

vistos = set()
coordenadas_multifurcacao = []
for p in pontos_a:
    chave = (round(p[0], 6), round(p[1], 6))
    if chave in chaves_multifurcacao and chave not in vistos:
        vistos.add(chave)
        coordenadas_multifurcacao.append(p)

# ============================================================
# 4. PROCESSAMENTO FÍSICO E MAPEAMENTO DE CORES
# ============================================================
# A. Extração dos Raios Físicos (Para Espessura dos Tubos)
if "raio" in mesh.cell_data:
    raios = np.asarray(mesh.cell_data["raio"], dtype=float)
else:
    print("[Aviso] Campo 'raio' não encontrado em cell_data. Usando raio genérico.")
    raios = np.full(mesh.n_cells, 0.05 * R)

# B. Extração dos Valores para Coloração
if campo_cor in mesh.cell_data:
    valores_cor = np.asarray(mesh.cell_data[campo_cor], dtype=float)
else:
    print(f"[Aviso] Campo '{campo_cor}' não encontrado no VTK. Alternando para 'raio'.")
    campo_cor = "raio"
    valores_cor = raios

# C. Escala Visual dos Tubos
raio_max_fisico = float(np.max(raios)) if raios.size > 0 else 1.0
raio_max_tela = 0.025 * R  # O tronco ocupará cerca de 2.5% do diâmetro do domínio
fator_escala = (raio_max_tela / raio_max_fisico) if raio_max_fisico > 0 else 1.0

# D. Normalização de Cores (Gradiente Reds)
cor_min = float(np.min(valores_cor)) if valores_cor.size > 0 else 0.0
cor_max = float(np.max(valores_cor)) if valores_cor.size > 0 else 1.0

normalizador = (
    mcolors.Normalize(vmin=cor_min, vmax=cor_max)
    if cor_max > cor_min
    else mcolors.Normalize(vmin=0.0, vmax=1.0)
)
mapa_de_cor = cm.get_cmap("Reds")

# ============================================================
# 5. RENDERIZAÇÃO DA CENA 3D (PYVISTA)
# ============================================================
plotter = pv.Plotter(window_size=[900, 900])
plotter.set_background("white")

# A. Adiciona Segmentos como Tubos 3D
for i in range(mesh.n_cells):
    segmento = mesh.extract_cells(i)
    p0, p1 = segmento.points[0], segmento.points[-1]

    # Ajusta e valida raio do tubo
    r_tubo = float(raios[i]) * fator_escala
    if not np.isfinite(r_tubo) or r_tubo <= 0.0:
        r_tubo = 0.005 * R

    rgb_cor = mapa_de_cor(normalizador(valores_cor[i]))[:3]
    tubo = pv.Line(p0, p1).tube(radius=r_tubo)

    plotter.add_mesh(tubo, color=rgb_cor)

# B. Destaque dos Nós de Interesse (Raiz, Terminais e Multifurcações)
if coordenadas_raiz:
    plotter.add_points(
        np.array(coordenadas_raiz),
        color="#1c4fd6",  # Azul marcante
        point_size=18,
        render_points_as_spheres=True
    )

if coordenadas_terminais:
    plotter.add_points(
        np.array(coordenadas_terminais),
        color="#1c9e4f",  # Verde
        point_size=10,
        render_points_as_spheres=True
    )

if coordenadas_multifurcacao:
    plotter.add_points(
        np.array(coordenadas_multifurcacao),
        color="blue",     # Azul brilhante para trifurcações (TRICCO)
        point_size=16,
        render_points_as_spheres=True
    )

# C. Domínio Circular Externo
circulo = pv.Circle(radius=R, resolution=300)
plotter.add_mesh(circulo, color="black", style="wireframe", line_width=3)

# ============================================================
# 6. LEGENDA, TÍTULO E CONFIGURAÇÃO DA CÂMERA
# ============================================================
rotulos_por_campo = {
    "raio":         ("Maior raio (tronco)", "Menor raio (terminais)"),
    "fluxo":        ("Maior fluxo (tronco)", "Menor fluxo (terminais)"),
    "profundidade": ("Mais profundo (folhas)", "Menos profundo (raiz)")
}
rotulo_alto, rotulo_baixo = rotulos_por_campo[campo_cor]

plotter.add_legend(
    [
        ["Raiz", "#1c4fd6"],
        ["Terminais (Folhas)", "#1c9e4f"],
        ["Multifurcação (3+ filhos)", "blue"],
        [rotulo_alto, "maroon"],
        [rotulo_baixo, "red"]
    ],
    bcolor="white",
    border=True,
    size=(0.24, 0.18),
    loc="upper right"
)

plotter.add_text(
    f"CCO-TRICCO - Estrutura Vascular (Coloração: {campo_cor.capitalize()})",
    position="upper_left",
    font_size=10,
    color="black"
)

# Configurações de Câmera 2D ortográfica paralela
plotter.enable_anti_aliasing()
plotter.camera.parallel_projection = True
plotter.camera_position = [(0, 0, 50), (0, 0, 0), (0, 1, 0)]
plotter.camera.parallel_scale = R * 1.05
plotter.view_xy()

# Exibe a janela interativa
plotter.show()