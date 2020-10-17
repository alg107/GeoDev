import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits import mplot3d

def plot_circle(ax, rad):
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = rad*np.cos(u)*np.sin(v)
    y = rad*np.sin(u)*np.sin(v)
    z = rad*np.cos(v)
    ax.plot_surface(x, y, z,color='blue')

r_s = 1
def odeset(x,tau):
    t, td, r, rd, th, thd, ph, phd = x
    if r < 1.01*r_s:
        return 0*x
    omrsor = 1-r_s/r

    tdd = -omrsor**(-1)*r_s/r**2*rd*td

    rdd = 0.5*omrsor**(-1)*r_s/r**2*rd**2 -0.5*omrsor*r_s/r**2*td**2 + r*omrsor*(thd**2 + np.sin(th)**2*phd**2)

    thdd = -2/r*thd*rd+np.sin(th)*np.cos(th)*phd**2

    phdd = -2/r*phd*rd - 2*(1/np.tan(th))*phd*thd

    dxdt = [
        td,
        tdd,
        rd,
        rdd,
        thd,
        thdd,
        phd,
        phdd
            ]
    return dxdt

tau = np.linspace(0,20, 1000)

y0 = [
        [0,1.5,1.7,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,1.75,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,1.8,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,1.85,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,1.9,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,1.95,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,2,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,2.05,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,2.1,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,2.15,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,2.2,0,np.pi/2,0.5,np.pi,0],
        [0,1.5,2.25,0,np.pi/2,0.5,np.pi,0]
]

ax = plt.axes(projection='3d')


for y in y0:
    sol = odeint(odeset, y, tau)
    t_sol = sol[:,0]
    r_sol = sol[:,2]
    th_sol = sol[:,4]
    ph_sol = sol[:,6]
    x_sol = r_sol*np.sin(th_sol)*np.cos(ph_sol)
    y_sol = r_sol*np.sin(th_sol)*np.sin(ph_sol)
    z_sol = r_sol*np.cos(th_sol)

    ax.plot(x_sol, y_sol, z_sol)

plot_circle(ax,r_s)
#plt.axis('square')
plotrange = 10
plt.xlim(-plotrange/2,plotrange/2)
plt.ylim(-plotrange/2,plotrange/2)
ax.set_zlim(-plotrange/2,plotrange/2)

plt.show()
