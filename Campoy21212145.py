"""
Práctica 2: Caos en Sistema Biologiico

Departamento de Ingeniería Eléctrica y Electrónica, Ingeniería Biomédica
Tecnológico Nacional de México [TecNM - Tijuana]
Blvd. Alberto Limón Padilla s/n, C.P. 22454, Tijuana, B.C., México

Nombre del alumno: Marco Antonio Campoy Alegria
Número de control: 21212145
Correo institucional: L21212145@tectijuana.edu.mx

Asignatura: Modelado de Sistemas Fisiológicos
Docente: Dr. Paul Antonio Valle Trujillo; paul.valle@tectijuana.edu.mx
"""

import control as ctrl
import numpy as np
import matplotlib.pyplot as plt # noqa
from scipy import signal
import pandas as pd

u = np.array(pd.read_excel('signal.xlsx',header=None))
x0,t0,tend,dt,w,h = 0,0,15,1e-3,10,5
"""
Práctica 2: Caos en Sistema Biologiico

Departamento de Ingeniería Eléctrica y Electrónica, Ingeniería Biomédica
Tecnológico Nacional de México [TecNM - Tijuana]
Blvd. Alberto Limón Padilla s/n, C.P. 22454, Tijuana, B.C., México

Nombre del alumno: Marco Antonio Campoy Alegria
Número de control: 21212145
Correo institucional: L21212145@tectijuana.edu.mx

Asignatura: Modelado de Sistemas Fisiológicos
Docente: Dr. Paul Antonio Valle Trujillo; paul.valle@tectijuana.edu.mx
"""
#libreria pip install controll
#libreria pip install slycot
# Librerias para calculo numerico y generacon de graficas
import control as ctrl
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import pandas as pd

u = np.array(pd.read_excel('signal.xlsx',header=None))
x0,t0,tend,dt,w,h = 0,0,15,1e-3,10,5
N = round((tend-t0)/dt) + 1
t = np.linspace(t0,tend,N)
u = np.reshape(signal.resample(u,len(t)),-1)

def cardio(Z,C,R,L):
    num = [L*R,R*Z]
    den = [C*L*R*Z,L*R+L*Z,R*Z]
    sys = ctrl.tf(num,den)
    return sys


#Funcion de transferencia: Normotenso
Z,C,R,L = 0.033,1.5,0.95,0.01
sysnormo = cardio(Z,C,R,L)
print(f'Funcion de transferencia del normotenso (control): {sysnormo}')

#Funcion de transferencia: Hipotenso
Z,C,R,L = 0.02,0.25,0.6,0.005
syshipo = cardio(Z,C,R,L)
print(f'Funcion de transferencia del hipotenso (caso 1): {syshipo}')

#Funcion de transferencia: Hipertenso
Z,C,R,L = 0.05,2.5,1.4,0.02
syshiper = cardio(Z,C,R,L)
print(f'Funcion de transferencia del hipertenso (caso 2): {syshiper}')

#Respuestas en lazo abierto
_,Pp0 = ctrl.forced_response(sysnormo,t,u,x0)
_,Pp1 = ctrl.forced_response(syshipo,t,u,x0)
_,Pp2 = ctrl.forced_response(syshiper,t,u,x0)

fg1 = plt.figure()
plt.plot(t,Pp0,'-',linewidth=1,color=[0.37,0.48,0.22],label='Pp(t): Normotenso')
plt.plot(t,Pp1,'-',linewidth=1,color=[0.42,0.11,0.12],label='Pp(t): Hipotenso')
plt.plot(t,Pp2,'-',linewidth=1,color=[0.61,0.16,0.72],label='Pp(t): Hipertenso')
plt.grid(False)
plt.xlim(0,15); plt.xticks(np.arange(0,16,1))
plt.ylim(-0.6,1.4); plt.yticks(np.arange(-0.6,1.6,0.2))
plt.xlabel('t[s]')
plt.ylabel('Pp(t) [V]')
plt.legend(bbox_to_anchor=(0.5,-0.2),loc='center',ncol=3)
plt.show()
fg1.set_size_inches(w,h)
fg1.tight_layout()
fg1.savefig('Cardiovascular lazo abierto python.pdf')

#Controlador PID
def controlador(kP,kI,kD,sys):
    Cr = 1e-6
    Re = 1/(kI*Cr)
    Rr = kP*Re
    Ce = kD/Rr
    numPID = [Re*Rr*Ce*Cr,(Re*Ce + Rr*Cr),1]
    denPID = [Re*Cr,0]
    PID = ctrl.tf(numPID,denPID)
    X = ctrl.series(PID,sys)
    sysPID = ctrl.feedback(X,1,sign=-1)
    return sysPID


hipoPID = controlador(1.494,352.000,0.000491,syshipo)
print(f'funcion de transferencia del hipotenso en lazo cerrado: {hipoPID}')
hiperPID = controlador(12.715,363.894,0.0343,syshiper)
print(f'funcion de transferencia del hipertenso en lazo cerrado: {hiperPID}')

#Respuestas del sistema del control en lazo cerrado
_,PID1 = ctrl.forced_response(hipoPID,t,Pp0,x0)
_,PID2 = ctrl.forced_response(hiperPID,t,Pp0,x0)





fig, ax = plt.subplots(2, 1, sharex=True)


ax[0].plot(t, Pp0, '-', linewidth=1.5,
color=[0.37,0.48,0.22],
label=r'$P_p(t)$: Normotenso')

ax[0].plot(t, Pp1, '-', linewidth=1.5,
color=[0.42,0.11,0.12],
label=r'$P_p(t)$: Hipotenso')

ax[0].plot(t, PID1, ':', linewidth=2,
color=[0.90, 0.60, 0.70],
label=r'$PID(t)$: Hipotenso')

ax[0].set_title('Normotenso vs Hipotenso')
ax[0].set_ylabel(r'$P_p(t)$ [V]')
ax[0].set_xlim(0, 15)
ax[0].grid(False)

# Leyenda ARRIBA (un poco más abajo)
ax[0].legend(loc='lower center',
bbox_to_anchor=(0.5, 1.12),
ncol=3,
frameon=False,
fontsize=11)

# Etiqueta (a)
ax[0].text(0.01, 1.05, '(a)', transform=ax[0].transAxes)

# Señal completa
ax[0].autoscale(enable=True, axis='y')
ax[0].margins(y=0.1)


ax[1].plot(t, Pp0, '-', linewidth=1.5,
color=[0.37,0.48,0.22],
label=r'$P_p(t)$: Normotenso')

ax[1].plot(t, Pp2, '-', linewidth=1.5,
color=[0.42,0.11,0.12],
label=r'$P_p(t)$: Hipertenso')

ax[1].plot(t, PID2, ':', linewidth=2,
color=[0.95, 0.60, 0.75],
label=r'$PID(t)$: Hipertenso')

ax[1].set_title('Normotenso vs Hipertenso')
ax[1].set_xlabel(r'$t$ [s]')
ax[1].set_ylabel(r'$P_p(t)$ [V]')
ax[1].set_xlim(0, 15)
ax[1].grid(False)

# Leyenda ARRIBA (un poco más abajo)
ax[1].legend(loc='lower center',
bbox_to_anchor=(0.5, 1.12),
ncol=3,
frameon=False,
fontsize=11)

# Etiqueta (b)
ax[1].text(0.01, 1.05, '(b)', transform=ax[1].transAxes)

# Señal completa
ax[1].autoscale(enable=True, axis='y')
ax[1].margins(y=0.1)



fig.set_size_inches(w, 2*h)

# Deja espacio para leyendas superiores
fig.tight_layout(rect=[0, 0, 1, 0.95])

plt.show()

# Guardar figura
fig.savefig('Subplots_Normotenso_Hipotenso_Hipertenso_PID.pdf')
