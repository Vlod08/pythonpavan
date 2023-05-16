import math

G = 6.67430e-11  # constante gravitationnelle en m^3/kg/s^2

class Orbit:
    def __init__(self, a, e, M):
        self.a = a  # demi-grand axe en m
        self.e = e  # excentricité
        self.M = M  # masse du corps central en kg
        self.T = 2*math.pi*math.sqrt(a**3/(G*M))  # période orbitale en s

    def position(self, t):
        n = 2*math.pi/self.T
        E = lambda E: E - self.e*math.sin(E) - n*t
        E0 = t/(self.T/2)  # approximation initiale
        E = self._newton_raphson(E, E0)
        x = self.a*(math.cos(E) - self.e)
        y = self.a*math.sqrt(1 - self.e**2)*math.sin(E)
        return x, y

    def _newton_raphson(self, f, x0):
        eps = 1e-12  # tolérance d'erreur
        dx = eps*10
        while abs(dx) > eps:
            df = (f(x0+dx)-f(x0))/dx
            dx = -f(x0)/df
            x0 += dx
        return x0

    def integrate(self, t0, tf, dt):
        positions = []
        t = t0
        while t <= tf:
            positions.append(self.position(t))
            t += dt
        return positions

# Orbite de la Terre
a_T = 149.60e9  # m
e_T = 0.0167
M_S = 1.989e30  # masse du Soleil en kg
orbit_T = Orbit(a_T, e_T, M_S)

# Orbite de Mars
a_M = 227.92e9  # m
e_M = 0.0934
M_M = 6.39e23  # masse de Mars en kg
orbit_M = Orbit(a_M, e_M, M_M)

# Orbite de transfert
a_transfert = (a_T + a_M)/2
e_transfert = 1 - a_T/a_M
orbit_transfert = Orbit(a_transfert, e_transfert, M_S)

# Fenêtres de tir
T_transfert = orbit_transfert.T
T_T = orbit_T.T
T_M = orbit_M.T
t = 0
windows = []
while t < T_T:
    t_transfert_start = t + T_T/4
    t_transfert_end = t + 3*T_T/4
    if (t_transfert_start % T_transfert) < (t_transfert_end % T_transfert):
        windows.append((t_transfert_start, t_transfert_end))
    t += T_T/10

# Affichage des trajectoires et des fenêtres de tir
import matplotlib.pyplot as plt

# Paramètres des trajectoires
num_trajectories = 5
initial_speeds = [25000, 30000, 35000, 40000, 45000]  # en m/s
initial_angles = [0, math.pi/6, math.pi/4, math.pi/3, 2*math.pi/3]

import matplotlib.pyplot as plt

for i in range(num_trajectories):
    fig, ax = plt.subplots()
    positions_T = orbit_T.integrate(0, orbit_T.T, orbit_T.T/1000)
    positions_M = orbit_M.integrate(0, orbit_M.T, orbit_M.T/1000)
    positions_transfert = orbit_transfert.integrate(0, orbit_transfert.T, orbit_transfert.T/1000)

    ax.plot([x for x,y in positions_T], [y for x,y in positions_T], label='Terre')
    ax.plot([x for x,y in positions_M], [y for x,y in positions_M], label='Mars')
    ax.plot([x for x,y in positions_transfert], [y for x,y in positions_transfert], label='Transfert')
    for t_start, t_end in windows:
        ax.axvspan(t_start, t_end, alpha=0.2, color='gray')

    # Paramètres de la trajectoire spécifique
    v0 = initial_speeds[i]
    theta = initial_angles[i]
    x0 = a_T*(math.cos(theta) - e_T)
    y0 = a_T*math.sqrt(1 - e_T**2)*math.sin(theta)
    ax.scatter(x0, y0, label=f'Trajectoire {i+1}', marker='x')
    
    ax.legend()
    plt.show()
