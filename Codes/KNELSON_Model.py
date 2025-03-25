import matplotlib.pyplot as plt # Gráficas
import numpy as np # Arreglos Vectoriales
from math import * # Símbolos matemáticos escalares
from sympy import symbols, solve, Eq # Resolución de ecuaciones

global g, k_valve, k_universal, k_elbow, k_check, k_out, D, k_te,  k_contraction1,   k_contraction2
g = 9.81 # [m/s2]
     # DATOS
Ddk = 3  * 2.54 / 100 # [m] diametro de salida del KNELSON
Dt = 1.07 #m Diametro del Tanque de agua
Ddt = 2 * 2.54 / 100 # diametro de la salida del Tanque de agua
Dru1 = 1.5 * 2.54 / 100 # diametro de reduccion 2" a 1.5" (acero)
Dru2 = 1 * 2.54 / 100 # diametro de reduccion 1.5" a 1" (pvc)
Deu1 = 1.5  * 2.54 / 100 # diametro de expansion 1" a 1.5" (pvc)
Dru3 = 0.75 * 2.54 / 100 # diametro de reduccion 1 1/2" a 3/4" (acero)


Ao= Dru3**2 * pi /4 # [m2] Diámetro de descarga
roughness_ac = 0.15e-3 # [m] Rugosidad acero
roughness_pvc = 0.0025e-3 # [m] Rugosidad PVC

 # Relaciones de diametros para reducciones y expansiones
R2_12 = Dru1/Ddt # RED 2" a 1.5"
R12_1 = Dru2/Dru1 # RED 1.5" a 1"
EX1_12 = Dru2/Deu1 # EXP 1" a 1.5"
R12_34 = Dru3/Dru1 # RED 1.5" a 3/4"
print(R2_12, R12_1, R12_34, EX1_12)
#longitudes equivalentes

L_Con_1 = 0.8 * 0.3048 # [m] Longitud equivalente de la Reduccion 1
L_Con_2 = 0.6 * 0.3048 # [m] Longitud equivalente de la Reduccion 1
L_Con_3 = 0.65 * 0.3048 # [m] Longitud equivalente de la Reduccion 1
L_Exp_1 = 0.6 * 0.3048 # [m] Longitud equivalente de la Reduccion 1
L_te_D = 20

# factores de friccion

def friction_f(Re, e, D):
    # se resuelve por iteraciones
    f = 0.01 # f inical
    izq = 1 / sqrt(f)
    # reacalcular 1/sqrt(f)
    izq_nuevo = -2 * log10((e/D) / 3.7 + 2.51 / Re * izq )
    while abs(izq - izq_nuevo) > 0.000000000010:
        izq = izq_nuevo
        izq_nuevo = -2 * log10((e / D) / 3.7 + 2.51 / Re * izq)
    return izq ** -2

f_leq_acer_t1 = friction_f(1e20, roughness_ac, Dru1)
f_leq_pvc_t1 = friction_f(1e20, roughness_pvc, Dru2)
f_leq_pvc_t2 = friction_f(1e20, roughness_pvc, Dru3)
f_leq_acer_t2 = friction_f(1e20, roughness_ac, Deu1)

#constantes de los accesorios
k_out = 1 # Salida del Tk de agua
k_elbow_90 = 0.9 # codo de 90
k_elbow_45 = 16 * f_leq_acer_t2 # codo de 45
k_check = 2.5 # valvula check
k_te = 1.8 # Tee
k_universal = 1.1 # Union universal
k_valve = 10 # valvula de bola
k_contraction1 = L_Con_1 / Dru1 * f_leq_acer_t1 # Reduccion 2" a 1.5"
k_contraction2 = L_Con_2 / Dru2 * f_leq_pvc_t1 # Reduccion 1.5" a 1"
k_expansion1 = L_Exp_1 / Deu1 * f_leq_pvc_t2 # Expansion 1" a 1.5"
k_contraction3 = L_Con_3 / Dru3 * f_leq_acer_t2 # Reduccion 1.5" a 3/4"


def reynolds(v,D):

    rho = 1000 # [kg / m3]
    mu = 1e-3 # viscosidad dinámica [Pa * s]
    Re = (rho * v * D) / mu
    return Re

def h_bomba(Q):
    # El Caudal debe estar en L / min
    h_pump = -0.443 * Q + 31.36 # [m]
    return h_pump

def h_pipe(f, D, v, L):
    hfpipe = f * L / D * v ** 2 / (2 * g)
    return hfpipe

def solve_s(h, v_ini):
    global Dru2, Deu1, g, Dru3
    rho = 1000 #kg/m3
    gamma = rho * g
    dPk = 5 * 101325/14.7 #5 psi KNELSON
    dPf = 2 * 101325/14.7 #2 psi Filter
    v = symbols('v')
    Ao = pi / 4 * Dru3 ** 2
    k_all = 2*k_te + 3 * k_valve + k_out + k_check + 2 * k_universal + k_elbow_45 + 3 * k_elbow_90 + k_contraction1 \
            + k_contraction2
    hfacc = (k_all + 1) * (v**2/(2*g))
    L = (4+8+4+6+14+5+4+7)/100 # Longitud de tuberia (uniones, neplos y tramos de tuberia) [m]
    L2 = .1 # Longitud de acoples de acero del filtro [m]
    hfpipe = h_pipe(friction_f(reynolds(v_ini, Dru2), roughness_pvc, Dru2), Dru2, v, L) + \
             h_pipe(friction_f(reynolds(v_ini, Deu1), roughness_ac, Deu1), Deu1, v, L2)
    exp = Eq(+hfpipe + hfacc - h + .27 + v ** 2 / (2 * g) - h_bomba(v * Ao * 1000 * 60) - dPk/gamma + dPf/gamma, 0)
    v = solve(exp)

    if v[0] > 0:
        v = v[0]
    else:
        v = v[1]
    return v

def bernoulli(h):
    global g, k_valve, k_universal, k_elbow, k_check, k_out, k_contraction1, k_contraction2, k_contraction3\
            , k_expansion1

    v_ini =1 #velocidad inicial
    #se recalcula la velocidad
    v = solve_s(h, v_ini)
    while abs(v-v_ini) > 0.00001:
        v_ini = v
        v = solve_s(h, v_ini)

    return float(v) * Ao * 60000 #[L/min]



#f_max = friction_f(reynolds(45 / (1000 * Ao * 60)), 0.0015e-3)
h = [0.33] # altura de llenado del tanque [m]
t = [0] # tiempo inicial de operacion [min]
dt = 0.1
At = 1.07**2 * pi/4
Ao = Dru3**2 * np.pi / 4
print(Ao)
print(dt)
while float(h[-1]) > 0.08:
   if h[-1] > 0:
        h.append(h[-1] - bernoulli(h[-1])/1000 / At * dt)
        t.append(t[-1] + dt)
   else:
       break

# DETERMINACIÓN DE LAS VELOCIDADES
v = np.array([bernoulli(i) for i in h])
# DATOS OBTENIDOS EXPERIMENTALMENTE
t_exp = [0, 1, 2, 3, 4, 5]
h_exp = np.array([33.0, 28.2, 23.5, 18.7, 14.0, 9.2]) * 1e-2
#GRÁFICO DE LOS RESULTADOS DEL MODELO REALIZADO
fig, ax = plt.subplots(1, 2,  layout='tight', figsize=[10, 5]) #Crear tanto la figura como el eje de coordenadas
ax[0].plot(t, h, '--k') #Grafico de las alturas simuladas en función del tiempo
ax[0].plot(t_exp, h_exp, '*k')
ax[0].legend(['Modelo', 'Datos Experimentales'])
# Propiedades de las gráficas
ax[0].set_xlabel('Tiempo [min]')
ax[0].set_ylabel('Altura de agua [m]')
ax[0].set_title('a)', position=[0,1]) # Nivel de agua vs tiempo
ax[1].plot(t, v, color='black') #Gráfico de caudales de descarga en función del tiempo
# Propiedades de las gráficas
ax[1].set_xlabel('Tiempo [min]')
ax[1].set_ylabel('Caudal [L/min]')
ax[1].set_title('b)', position=[0,1])# Caudal vs tiempo
# Guardar la figura
# plt.savefig('/Figuras/modelo.png', 
#             dpi=fig.dpi, 
#             bbox_inches='tight' # ajustar el margen exactamente al texto para evitar tener pérdidas de información
#             )
plt.show()

