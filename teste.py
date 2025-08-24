import vpython as vp
import random

# ----------------- CONFIGURAÇÕES GERAIS -----------------
NUM_ATOMOS = 6
RAIO_ATOMO = 0.9
RAIO_CASCA = 1.25

CAIXA_SIZE = 22

DT = 0.1
NUM_PASSOS = 1000

R1 = 6.0
K = 0.5
K_REPULSAO = 0.6  

P_COLAGEM = 0.5
K_GLOBAL = 1

DAMPING_GLOBAL = 0.99
VEL_CUTOFF = 0   # corte explícito de velocidade

P_REBELDE = 0.05
REBELDE_MULT = 0.2  

# Lista de elementos e massas associadas
ELEMENTOS = [
    "S", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"
]

MASSAS = [
    32, 48, 51, 52, 55, 56, 59, 59, 63, 65,
    89, 91, 93, 96, 98, 101, 103, 106, 108, 112,
    175, 178, 181, 184, 186, 190, 192, 195, 197, 200
]

# ----------------- FUNÇÕES -----------------
def colisao_casca_caixa(a, massa):
    for axis in ['x','y','z']:
        limite = CAIXA_SIZE/2 - RAIO_CASCA
        if getattr(a.pos, axis) > limite:
            setattr(a.pos, axis, limite)
            setattr(a.vel, axis, -getattr(a.vel, axis))
        elif getattr(a.pos, axis) < -limite:
            setattr(a.pos, axis, -limite)
            setattr(a.vel, axis, -getattr(a.vel, axis))

def colisao_atomica(a1, a2, massa):
    r_vec = a2.pos - a1.pos
    dist = vp.mag(r_vec)

    if dist < 2*RAIO_ATOMO:  # colisão núcleo-núcleo
        n = vp.norm(r_vec)
        v_rel = a1.vel - a2.vel
        a1.vel -= (2*massa/(massa+massa)) * vp.dot(v_rel, n) * n
        a2.vel += (2*massa/(massa+massa)) * vp.dot(v_rel, n) * n
        overlap = 2*RAIO_ATOMO - dist
        a1.pos -= n * overlap/2
        a2.pos += n * overlap/2
        a1.cluster = [a1]
        a2.cluster = [a2]
        return

    # colagem probabilística
    min_col = RAIO_ATOMO + 0.05
    max_col = RAIO_CASCA
    if min_col < dist < max_col:
        if not any(x in a2.cluster for x in a1.cluster) and random.random() < P_COLAGEM:
            novo_cluster = list(set(a1.cluster + a2.cluster))
            for x in novo_cluster:
                x.cluster = novo_cluster

# ----------------- LOOP PRINCIPAL PARA TODOS OS ELEMENTOS -----------------
for SYS, MASSA in zip(ELEMENTOS, MASSAS):

    fator = 0.01
    RAIO_EFETIVO = RAIO_ATOMO - fator

    scene = vp.canvas(title=f"MD {SYS}", width=800, height=600, center=vp.vector(0,0,0))
    caixa = vp.box(pos=vp.vector(0,0,0), length=CAIXA_SIZE, height=CAIXA_SIZE, width=CAIXA_SIZE,
                   opacity=0.2, color=vp.color.white)
    esfera_central = vp.sphere(pos=vp.vector(0,0,0), radius=R1, color=vp.color.cyan, opacity=0.2)

    # Criar átomos
    atomos = []
    for i in range(NUM_ATOMOS):
        while True:
            pos = vp.vector(random.uniform(-CAIXA_SIZE/2 + RAIO_CASCA, CAIXA_SIZE/2 - RAIO_CASCA),
                            random.uniform(-CAIXA_SIZE/2 + RAIO_CASCA, CAIXA_SIZE/2 - RAIO_CASCA),
                            random.uniform(-CAIXA_SIZE/2 + RAIO_CASCA, CAIXA_SIZE/2 - RAIO_CASCA))
            if all(vp.mag(pos - a.pos) >= 2*RAIO_CASCA for a in atomos):
                break
        vel = vp.vector(random.uniform(-1,1), random.uniform(-1,1), random.uniform(-1,1))
        atomo = vp.sphere(pos=pos, radius=RAIO_EFETIVO, color=vp.color.red)
        casca = vp.sphere(pos=pos, radius=RAIO_CASCA, color=vp.color.yellow, opacity=0.3)
        atomo.vel = vel
        atomo.casca = casca
        atomo.cluster = [atomo]
        atomos.append(atomo)

    # Loop de simulação
    for step in range(NUM_PASSOS):
        vp.rate(100)

        factor_off = 1.0 if step < 0.9*NUM_PASSOS else max(0.0, (NUM_PASSOS-step)/(0.1*NUM_PASSOS))

        for i in range(NUM_ATOMOS):
            for j in range(i+1, NUM_ATOMOS):
                colisao_atomica(atomos[i], atomos[j], MASSA)

        for a in atomos:
            r_vec = a.pos - esfera_central.pos
            dist_centro = vp.mag(r_vec)
            if dist_centro < R1:
                F_repulsa = K_REPULSAO * r_vec * factor_off
                a.vel += (F_repulsa / MASSA) * DT   # massa influencia

            F_global = K_GLOBAL * (vp.vector(0,0,0) - a.pos)
            a.vel += (F_global / MASSA) * DT       # massa influencia

            a.vel *= DAMPING_GLOBAL                # amortecimento global

            if random.random() < P_REBELDE:
                direcao = vp.vector(random.uniform(-1,1), random.uniform(-1,1), random.uniform(-1,1))
                a.vel += (direcao.norm() * REBELDE_MULT) / MASSA  # massa influencia

            # corte explícito de velocidade (zera de vez se ficar muito pequena)
            if vp.mag(a.vel) < VEL_CUTOFF:
                a.vel = vp.vector(0,0,0)

            a.pos += a.vel * DT
            a.casca.pos = a.pos
            colisao_casca_caixa(a, MASSA)

    # Exportar XYZ final
    XYZ_FILE = f"simulacao_{SYS}_{NUM_ATOMOS}.xyz"
    with open(XYZ_FILE, "w") as f:
        f.write(f"{NUM_ATOMOS}\n")
        f.write(f"MD final positions | RAIO_ATOMO={RAIO_EFETIVO} Å | RAIO_CASCA={RAIO_CASCA} Å | Massa={MASSA}\n")
        for a in atomos:
            f.write(f"{SYS} {a.pos.x:.4f} {a.pos.y:.4f} {a.pos.z:.4f}\n")

    print(f"Arquivo XYZ gerado: {XYZ_FILE}")
    scene.delete()  # fecha a cena antes da próxima simulação
