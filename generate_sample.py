"""
generate_sample.py — gera dados de EXEMPLO (SINTÉTICOS) para testar o EBSD Analyzer.

IMPORTANTE: os arquivos gerados NÃO são dados experimentais reais. São dados
sintéticos, criados por computador apenas para demonstrar o funcionamento do app.
Não devem ser usados como resultado de pesquisa.

Gera dois arquivos dentro de sample_data/:
  1) exemplo_graos.csv  -> dados por GRÃO (CSV), pronto para uso imediato
  2) exemplo_pixels.ctf -> dados por PIXEL (formato CTF Oxford/HKL), para
                           demonstrar a segmentação automática de grãos

Uso:
    python generate_sample.py
"""

import os
import numpy as np
import pandas as pd

SEED = 42
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_data")


def gerar_csv_por_grao(n=800):
    """Dados por grão: distribuição log-normal de tamanho + texturas e KAM."""
    rng = np.random.default_rng(SEED)

    # Diâmetro de grão (ECD) log-normal, média ~25 µm
    diametro = rng.lognormal(mean=np.log(22.0), sigma=0.45, size=n)

    # Desorientação média por grão (graus): mistura LAGB/HAGB
    lagb = rng.uniform(2, 15, size=n)
    hagb = rng.uniform(15, 62, size=n)
    is_hagb = rng.random(n) < 0.70
    misori = np.where(is_hagb, hagb, lagb)

    # KAM (graus)
    kam = np.abs(rng.normal(0.6, 0.3, size=n)).clip(0.01, 3.0)

    # Ângulos de Euler (graus) com leve preferência por fibra-gamma (Phi ~ 55)
    phi1 = rng.uniform(0, 360, size=n)
    Phi = rng.normal(55, 12, size=n).clip(0, 90)
    phi2 = rng.uniform(0, 90, size=n)

    # Fases: 95% Ferrita (BCC), 5% Austenita
    fase = np.where(rng.random(n) < 0.95, "Ferrite", "Austenite")

    # Qualidade de indexação
    ci = rng.uniform(0.2, 1.0, size=n).round(3)
    iq = rng.uniform(80, 220, size=n).round(0)

    # Injeta alguns outliers artificiais de tamanho de grão
    n_out = 20
    idx_out = rng.choice(n, size=n_out, replace=False)
    diametro[idx_out] = rng.uniform(110, 180, size=n_out)

    df = pd.DataFrame({
        "Grain Diameter": diametro.round(3),
        "Misorientation Angle": misori.round(3),
        "KAM": kam.round(4),
        "phi1": phi1.round(2),
        "Phi": Phi.round(2),
        "phi2": phi2.round(2),
        "Phase": fase,
        "CI": ci,
        "IQ": iq.astype(int),
    })
    return df


def gerar_ctf_por_pixel(ncols=60, nrows=50, step=1.0):
    """
    Gera um mapa de pixels em grade regular com ~6 grãos, no formato CTF
    (Oxford/HKL). Cada grão tem orientação de Euler própria; assim o app
    consegue segmentar os grãos automaticamente.
    """
    rng = np.random.default_rng(SEED)

    # Define centros de grãos e atribui cada pixel ao centro mais próximo (Voronoi)
    n_graos = 6
    cx = rng.uniform(0, ncols, size=n_graos)
    cy = rng.uniform(0, nrows, size=n_graos)
    # Orientação de Euler por grão
    g_phi1 = rng.uniform(0, 360, size=n_graos)
    g_Phi = rng.uniform(0, 90, size=n_graos)
    g_phi2 = rng.uniform(0, 90, size=n_graos)

    rows = []
    for j in range(nrows):
        for i in range(ncols):
            x = i * step
            y = j * step
            d = (cx - i) ** 2 + (cy - j) ** 2
            g = int(np.argmin(d))
            # pequeno ruído intragrão
            e1 = (g_phi1[g] + rng.normal(0, 0.8)) % 360
            e2 = float(np.clip(g_Phi[g] + rng.normal(0, 0.8), 0, 90))
            e3 = (g_phi2[g] + rng.normal(0, 0.8)) % 90
            bands = 8
            error = 0
            mad = round(abs(rng.normal(0.4, 0.1)), 3)
            bc = int(np.clip(rng.normal(150, 25), 0, 255))
            bs = int(np.clip(rng.normal(120, 25), 0, 255))
            phase = 1  # 1 = Ferrite
            rows.append((phase, x, y, bands, error,
                         round(e1, 3), round(e2, 3), round(e3, 3), mad, bc, bs))
    return rows, step


def escrever_ctf(path, rows, step):
    """Escreve um arquivo CTF mínimo, mas válido (cabeçalho + dados tab-sep)."""
    header = [
        "Channel Text File",
        "Prj\tEBSD Analyzer synthetic example",
        "Author\tgenerate_sample.py (DADOS SINTETICOS)",
        "JobMode\tGrid",
        f"XCells\t{int(max(r[1] for r in rows) / step) + 1}",
        f"YCells\t{int(max(r[2] for r in rows) / step) + 1}",
        f"XStep\t{step}",
        f"YStep\t{step}",
        "AcqE1\t0",
        "AcqE2\t0",
        "AcqE3\t0",
        "Euler angles refer to Sample Coordinate system (CS0)!",
        "Phases\t1",
        # lattice a;b;c;alpha;beta;gamma;SG;name
        "2.866;2.866;2.866;90;90;90;11;Ferrite;;",
        "Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    df = gerar_csv_por_grao()
    csv_path = os.path.join(OUT_DIR, "exemplo_graos.csv")
    df.to_csv(csv_path, index=False)
    print(f"[ok] CSV por grao gerado: {csv_path}  ({len(df)} graos)")

    rows, step = gerar_ctf_por_pixel()
    ctf_path = os.path.join(OUT_DIR, "exemplo_pixels.ctf")
    escrever_ctf(ctf_path, rows, step)
    print(f"[ok] CTF por pixel gerado: {ctf_path}  ({len(rows)} pixels)")

    print("\nATENCAO: arquivos SINTETICOS (gerados por computador) apenas para teste.")


if __name__ == "__main__":
    main()
