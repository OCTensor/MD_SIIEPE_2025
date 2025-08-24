import random
import math

# ----------------- CONFIGURAÇÕES -----------------
NUM_ATOMOS = 2
R1 = 8.0
DT = 0.5
NUM_PASSOS = 100

K1 = 5    # força para R1
K2 = 10   # repulsão dentro de R1
K3 = 40   # força global para o centro

MASSA = 195.0
DAMP_COLL = 0.8
DAMP_CASCA2 = 0.9
DAMP_R1 = 0.5  # diminuição de velocidade dentro de R1

RAIO_CASCA1 = 1.15
RAIO_CASCA2 = 1.5

# ----------------- FUNÇÕES -----------------
def mag(vec):
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def norm(vec):
    m = mag(vec)
    if m == 0:
        return [0,0,0]
    return [v/m for v in vec]

def sub(a, b):
    return [a[i]-b[i] for i in range(3)]

def add(a, b):
    return [a[i]+b[i] for i in range(3)]

def scale(vec, s):
    return [v*s for v in vec]

# ----------------- INICIALIZAÇÃO -----------------
atomos = []
for _ in range(NUM_ATOMOS):
    pos = [random.uniform(-R1,R1) for _ in range(3)]
    vel = [random.uniform(-1,1) for _ in range(3)]
    atomos.append({"pos": pos, "vel": vel})

# ----------------- SIMULAÇÃO -----------------
for step in range(NUM_PASSOS):
    # Atualiza forças e velocidades
    for nucleo in atomos:
        r = mag(nucleo["pos"])
        direcao = norm(nucleo["pos"])

        # 3 forças
        F1 = scale(direcao, -K1*(r-R1))
        if r < R1:
            F2 = scale(direcao, K2*(R1-r))
            F2 = scale(F2, DAMP_R1)
        else:
            F2 = [0,0,0]
        F3 = scale(nucleo["pos"], -K3)
        
        # soma forças
        total_F = [F1[i]+F2[i]+F3[i] for i in range(3)]
        nucleo["vel"] = [nucleo["vel"][i] + total_F[i]/MASSA*DT for i in range(3)]

    # Colisões de casca1
    for i in range(NUM_ATOMOS):
        for j in range(i+1, NUM_ATOMOS):
            r_vec = sub(atomos[j]["pos"], atomos[i]["pos"])
            dist = mag(r_vec)
            min_dist = 2*RAIO_CASCA1
            if dist < min_dist:
                n = norm(r_vec)
                overlap = min_dist - dist
                atomos[i]["pos"] = sub(atomos[i]["pos"], scale(n, overlap/2))
                atomos[j]["pos"] = add(atomos[j]["pos"], scale(n, overlap/2))
                v_rel = sub(atomos[i]["vel"], atomos[j]["vel"])
                atomos[i]["vel"] = sub(atomos[i]["vel"], scale(n, sum(v_rel[k]*n[k] for k in range(3))*DAMP_COLL))
                atomos[j]["vel"] = add(atomos[j]["vel"], scale(n, sum(v_rel[k]*n[k] for k in range(3))*DAMP_COLL))

    # Atualiza posições
    for nucleo in atomos:
        nucleo["pos"] = add(nucleo["pos"], scale(nucleo["vel"], DT))

# ----------------- EXPORTAR XYZ -----------------
with open("platina_sim.xyz", "w") as f:
    f.write(f"{NUM_ATOMOS}\n")
    f.write("Platina simulada\n")
    for nucleo in atomos:
        x, y, z = nucleo["pos"]
        f.write(f"Pt {x:.6f} {y:.6f} {z:.6f}\n")
