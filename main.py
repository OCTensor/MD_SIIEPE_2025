import vpython as vp
import random

# ----------------- CONFIGURAÇÕES -----------------
NUM_ATOMOS = 6
MASSA = 112
RAIO_ATOMO = 1.76
SYS = "Cd"


fator = 0.01

RAIO_ATOMO = RAIO_ATOMO - fator

RAIO_CASCA = 1.5
CAIXA_SIZE = 20

DT = 0.06
NUM_PASSOS = 3000

R1 = 8.0
K = 5.0
K_REPULSAO = 3  # força fraca empurrando para fora de R1

P_COLAGEM = 0.2
K_GLOBAL = 0.1

DAMPING_GLOBAL = 0.96

P_REBELDE = 0.02  # 2% de chance por passo de um átomo ganhar impulso extra
REBELDE_MULT = 0.2  # intensidade do impulso extra
XYZ_FILE = f"simulacao_{SYS}_{NUM_ATOMOS}.xyz"
# -------------------------------------------------

# Criar cena e caixa
scene = vp.canvas(title="MD simplificado", width=800, height=600, center=vp.vector(0,0,0))
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
    atomo = vp.sphere(pos=pos, radius=RAIO_ATOMO, color=vp.color.red)
    casca = vp.sphere(pos=pos, radius=RAIO_CASCA, color=vp.color.yellow, opacity=0.3)
    atomo.vel = vel
    atomo.casca = casca
    atomo.cluster = [atomo]
    atomos.append(atomo)

# ----------------- FUNÇÕES -----------------
def colisao_casca_caixa(a):
    for axis in ['x','y','z']:
        limite = CAIXA_SIZE/2 - RAIO_CASCA
        if getattr(a.pos, axis) > limite:
            setattr(a.pos, axis, limite)
            setattr(a.vel, axis, -getattr(a.vel, axis))
        elif getattr(a.pos, axis) < -limite:
            setattr(a.pos, axis, -limite)
            setattr(a.vel, axis, -getattr(a.vel, axis))

def colisao_atomica(a1, a2):
    r_vec = a2.pos - a1.pos
    dist = vp.mag(r_vec)

    # colisão núcleo-núcleo
    if dist < 2*RAIO_ATOMO:
        n = vp.norm(r_vec)
        v_rel = a1.vel - a2.vel
        a1.vel -= (2*MASSA/(MASSA+MASSA)) * vp.dot(v_rel, n) * n
        a2.vel += (2*MASSA/(MASSA+MASSA)) * vp.dot(v_rel, n) * n
        overlap = 2*RAIO_ATOMO - dist
        a1.pos -= n * overlap/2
        a2.pos += n * overlap/2
        # quebrar cluster se houver
        a1.cluster = [a1]
        a2.cluster = [a2]
        return

    # colagem probabilística na casca (faixa núcleo+0.05 até casca)
    min_col = RAIO_ATOMO + 0.05
    max_col = RAIO_CASCA
    if min_col < dist < max_col:
        if not any(x in a2.cluster for x in a1.cluster) and random.random() < P_COLAGEM:
            novo_cluster = list(set(a1.cluster + a2.cluster))
            for x in novo_cluster:
                x.cluster = novo_cluster

# ----------------- LOOP DE SIMULAÇÃO -----------------
energia_total = []

for step in range(NUM_PASSOS):
    vp.rate(100)
    
    # fator de desligamento gradual da força externa
    factor_off = 1.0 if step < 0.9*NUM_PASSOS else max(0.0, (NUM_PASSOS-step)/(0.1*NUM_PASSOS))

    # Colisões e clusters
    for i in range(NUM_ATOMOS):
        for j in range(i+1, NUM_ATOMOS):
            colisao_atomica(atomos[i], atomos[j])

    E_kin = 0.0
    E_pot = 0.0
    
    for a in atomos:
        # força central empurrando para fora de R1
        r_vec = a.pos - esfera_central.pos
        dist_centro = vp.mag(r_vec)
        if dist_centro < R1:
            F_repulsa = K_REPULSAO * r_vec * factor_off
            a.vel += F_repulsa / MASSA * DT
        
        # força global
        F_global = K_GLOBAL * (vp.vector(0,0,0) - a.pos)
        a.vel += F_global * DT

        # perda gradual de velocidade
        a.vel *= DAMPING_GLOBAL

        # impulso “rebelde” ocasional
        if random.random() < P_REBELDE:
            direcao = vp.vector(random.uniform(-1,1), random.uniform(-1,1), random.uniform(-1,1))
            a.vel += direcao.norm() * REBELDE_MULT

        # mover átomo
        a.pos += a.vel * DT
        a.casca.pos = a.pos
        colisao_casca_caixa(a)
        
        # energia
        E_kin += 0.5 * MASSA * vp.mag(a.vel)**2
        if dist_centro < R1:
            E_pot += 0.5 * K * dist_centro**2
        E_pot += 0.5 * K_GLOBAL * vp.mag(a.pos)**2

    energia_total.append((E_kin, E_pot, E_kin + E_pot))

# ----------------- EXPORTAR XYZ -----------------
with open(XYZ_FILE, "w") as f:
    f.write(f"{NUM_ATOMOS}\n")
    f.write(f"MD simulation final positions | RAIO_ATOMO={RAIO_ATOMO} Å | RAIO_CASCA={RAIO_CASCA} Å\n")
    for a in atomos:
        f.write(f"{SYS} {a.pos.x:.4f} {a.pos.y:.4f} {a.pos.z:.4f}\n")

print(f"Arquivo XYZ gerado: {XYZ_FILE}")

# ----------------- PRINT ENERGIA FINAL -----------------
for i, (Ek, Ep, Et) in enumerate(energia_total):
    print(f"Passo {i}: Ecin={Ek:.2f}, Epot={Ep:.2f}, Etotal={Et:.2f}")
