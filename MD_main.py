from vpython import *
import numpy as np

# ================== CONFIGURAÇÕES ==================
CAIXA      = 22
NUM_ATOMOS = 16

RAIO       = 0.5
CASCA      = 1.4

MASSA      = 195

dt = 1

K_atom_atom = -30  # força de atração quando as cascas se sobrepõem

K_atrativo  = 3
K_repulsivo = 4

redutor_vel = 0.95

# ================== FUNÇÕES ==================
def gerar_posicoes(NUM_ATOMOS, RAIO):
    atomos = []
    for i in range(NUM_ATOMOS):
        atomo = vec(*(np.random.rand(3)*CAIXA - CAIXA/2))
        atomos.append(atomo)

    limite = 2 * RAIO
    sobreposto = True
    while sobreposto:
        sobreposto = False
        for i in range(len(atomos)):
            for j in range(i + 1, len(atomos)):
                dist = mag(atomos[i] - atomos[j])
                if dist < limite:
                    desloc = vec(*(np.random.rand(3) - 0.5) * 0.5)
                    atomos[j] += desloc
                    sobreposto = True
    return atomos

def gerar_velocidades(NUM_ATOMOS):
    velocidades = []
    for i in range(NUM_ATOMOS):
        vel = vec(*(np.random.rand(3) - 0.5)*2)
        velocidades.append(vel)
    return velocidades

def salvar_xyz(atomos, filename="configuracao.xyz"):
    with open(filename, "w") as f:
        f.write(f"{len(atomos)}\n")
        f.write("Configuração atual dos átomos\n")
        for atomo in atomos:
            f.write(f"Pt {atomo.pos.x:.6f} {atomo.pos.y:.6f} {atomo.pos.z:.6f}\n")

# ================== INICIA POSIÇÕES E VELOCIDADES ==================
posicoes_iniciais = gerar_posicoes(NUM_ATOMOS, RAIO)
velocidades = gerar_velocidades(NUM_ATOMOS)

# ================== CENA ==================
scene.width  = 800
scene.height = 600
scene.title = "Simulação Molecular com VPython"

# ================== CAIXA VISUAL ==================
box(pos=vec(0,0,0), length=CAIXA, height=CAIXA, width=CAIXA,
    color=color.cyan, opacity=0.2)

# ================== CRIA ÁTOMOS ==================
atomos = []
for i in range(NUM_ATOMOS):
    atomo = sphere(pos=posicoes_iniciais[i], radius=RAIO,
                   color=color.red, make_trail=False)
    atomos.append(atomo)

cascas = []
for i in range(NUM_ATOMOS):
    casca_c = sphere(pos=posicoes_iniciais[i], radius=CASCA,
                   color=color.blue, make_trail=False, opacity=0.3)
    cascas.append(casca_c)
    
# ================== LOOP PRINCIPAL ==================
passo = 0
while True:
    rate(30)
    passo += 1
    
    # Atualiza forças, acelerações, velocidades e posições
    for i, (atomo, casca_c) in enumerate(zip(atomos, cascas)):
        f1 = -K_atrativo * atomo.pos
        r = mag(atomo.pos) + 1e-6
        f2 = (K_repulsivo / r**2) * norm(atomo.pos)
        forca_total = f1 + f2
        a = forca_total / MASSA
        velocidades[i] += a * dt
        velocidades[i] *= redutor_vel
        atomo.pos += velocidades[i] * dt + 0.5 * a * dt**2
        casca_c.pos = atomo.pos
            
    # Colisões com a caixa        
    for i, atomo in enumerate(atomos):
        if atomo.pos.x + RAIO > CAIXA/2:
            atomo.pos.x = CAIXA/2 - RAIO
            velocidades[i].x *= -1
        elif atomo.pos.x - RAIO < -CAIXA/2:
            atomo.pos.x = -CAIXA/2 + RAIO
            velocidades[i].x *= -1
        if atomo.pos.y + RAIO > CAIXA/2:
            atomo.pos.y = CAIXA/2 - RAIO
            velocidades[i].y *= -1
        elif atomo.pos.y - RAIO < -CAIXA/2:
            atomo.pos.y = -CAIXA/2 + RAIO
            velocidades[i].y *= -1
        if atomo.pos.z + RAIO > CAIXA/2:
            atomo.pos.z = CAIXA/2 - RAIO
            velocidades[i].z *= -1
        elif atomo.pos.z - RAIO < -CAIXA/2:
            atomo.pos.z = -CAIXA/2 + RAIO
            velocidades[i].z *= -1

    # Colisões entre átomos
    for m in range(NUM_ATOMOS):
        for n in range(m+1, NUM_ATOMOS):
            rij = atomos[n].pos - atomos[m].pos
            dist = mag(rij)
            if dist < 2*CASCA:  
                rij_hat = norm(rij)
                overlap = 2*CASCA - dist
                f_atr = K_atom_atom * overlap * rij_hat
                a_m = f_atr / MASSA
                a_n = -f_atr / MASSA
                velocidades[m] += a_m * dt
                velocidades[n] += a_n * dt
            rij = atomos[m].pos - atomos[n].pos
            dist = mag(rij)
            if dist < 2*RAIO:
                rij_hat = norm(rij)
                overlap = 2*RAIO - dist
                atomos[m].pos += 0.5 * overlap * rij_hat
                atomos[n].pos -= 0.5 * overlap * rij_hat
                vi = velocidades[m]
                vj = velocidades[n]
                vi_proj = dot(vi, rij_hat)
                vj_proj = dot(vj, rij_hat)
                vi_new = vi + (vj_proj - vi_proj) * rij_hat
                vj_new = vj + (vi_proj - vi_proj) * rij_hat
                velocidades[m] = vi_new
                velocidades[n] = vj_new

    if passo % 100 == 0:
        salvar_xyz(atomos, filename=f"config_{passo}.xyz")
