from numpy import *
from numpy import linalg as la
import matplotlib.pyplot as plt
import sympy as smp
from scipy.integrate import odeint
class Pendel:
    def __init__(self, fig, ax):
        self.ax = ax
        self.canvas = fig.canvas
        self.px = array([])  # x coordinates for control polygon
        self.py = array([])  # y coordinates
        self.nDrag = -1

        self.canvas.mpl_connect('button_press_event', self.ButtonPress)
        self.canvas.mpl_connect('motion_notify_event', self.Move)
        self.canvas.mpl_connect('button_release_event', self.Release)

        self.xlim = (0, 16)
        self.ylim = (0, 12)
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

        #Pendel Value
        self.points = array([])

        #Sympy Stuff
        self.L1 = smp.symbols('L1')
        self.L2 = smp.symbols('L2')
        self.m1 = smp.symbols('m1')
        self.m2 = smp.symbols('m2')
        self.t = smp.symbols('t')
        self.g = smp.symbols('g')
        self.the1, self.the2 = smp.symbols(r'\theta_1, \theta_2', cls=smp.Function)

        self.the1 = self.the1(self.t)
        self.the2 = self.the2(self.t)

        self.the1_d = smp.diff(self.the1, self.t)
        self.the2_d = smp.diff(self.the2, self.t)
        self.the1_dd = smp.diff(self.the1_d, self.t)
        self.the2_dd = smp.diff(self.the2_d, self.t)

        self.x1 = self.L1 * smp.sin(self.the1)
        self.y1 = -self.L1 * smp.cos(self.the1)
        self.x2 = self.L1 * smp.sin(self.the1) + self.L2 * smp.sin(self.the2)
        self.y2 = -self.L1 * smp.cos(self.the1) - self.L2 * smp.cos(self.the2)

        # Kinetic
        T1 = 1 / 2 * self.m1 * (smp.diff(self.x1, self.t) ** 2 + smp.diff(self.y1, self.t) ** 2)
        T2 = 1 / 2 * self.m2 * (smp.diff(self.x2, self.t) ** 2 + smp.diff(self.y2, self.t) ** 2)
        T = T1 + T2
        # Potential
        V1 = self.m1 * self.g * self.y1
        V2 = self.m2 * self.g * self.y2
        V = V1 + V2
        # Lagrangian
        self.L = T - V

        LE1 = smp.diff(self.L, self.the1) - smp.diff(smp.diff(self.L, self.the1_d), self.t).simplify()
        LE2 = smp.diff(self.L, self.the2) - smp.diff(smp.diff(self.L, self.the2_d), self.t).simplify()

        self.sols = smp.solve([LE1, LE2], (self.the1_dd, self.the2_dd),
                         simplify=False, rational=False)

        self.dz1dt_f = smp.lambdify((self.t, self.g, self.m1, self.m2, self.L1, self.L2, self.the1, self.the2, self.the1_d, self.the2_d), self.sols[self.the1_dd])
        self.dz2dt_f = smp.lambdify((self.t, self.g, self.m1,self.m2, self.L1, self.L2, self.the1, self.the2, self.the1_d, self.the2_d), self.sols[self.the2_dd])
        self.dthe1dt_f = smp.lambdify(self.the1_d, self.the1_d)
        self.dthe2dt_f = smp.lambdify(self.the2_d, self.the2_d)

        t = linspace(0, 40, 1001)
        g = 9.81
        m1 = 2
        m2 = 1
        L1 = 2
        L2 = 1
        #TODO y0 = Startwerte einsetzen
        self.ans = odeint(self.dSdt, y0=[1, -3, -1, 5], t=t, args=(g, m1, m2, L1, L2))
        self.x1, self.y1, self.x2, self.y2 = self.get_x1y1x2y2(t, self.ans.T[0], self.ans.T[2], L1, L2)
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.set_facecolor('k')
        ax.get_xaxis().set_ticks([])  # enable this to hide x axis ticks
        ax.get_yaxis().set_ticks([])  # enable this to hide y axis ticks
        self.ln1, = plt.plot([], [], 'ro--', lw=3, markersize=8)
        ax.set_ylim(-4, 4)
        ax.set_xlim(-4, 4)
        self.ani = animation.FuncAnimation(fig, self.animate, frames=1000, interval=50)
        #ani.save('pen.gif', writer='pillow', fps=25)

        plt.show()  # draw empty subplot with axes

    def get_x1y1x2y2(self, t, the1, the2, L1, L2):
        return (L1 * sin(the1),
                -L1 * cos(the1),
                L1 * sin(the1) + L2 * sin(the2),
                -L1 * cos(the1) - L2 * cos(the2))


    def dSdt(self, S, t, g, m1, m2, L1, L2):
        the1, z1, the2, z2 = S
        return [
            self.dthe1dt_f(z1),
            self.dz1dt_f(t, g, m1, m2, L1, L2, the1, the2, z1, z2),
            self.dthe2dt_f(z2),
            self.dz2dt_f(t, g, m1, m2, L1, L2, the1, the2, z1, z2),

        ]

    def animate(self, i):
        self.ln1.set_data([0, self.x1[i], self.x2[i]], [0, self.y1[i], self.y2[i]])




    def ButtonPress(self, event):  # mouse button pressed event
        # print( event.xdata, event.ydata)
        if event.button == 1:  # add new point
            if self.px.size < 3:
                self.px = r_[self.px, event.xdata]
                self.py = r_[self.py, event.ydata]
                #if self.px.size == 3:
                    #self.calcLength()
        elif event.button == 3:  # move closest point
            if size(self.px) > 0:
                dist2 = (self.px - event.xdata) ** 2 + (self.py - event.ydata) ** 2
                self.nDrag = argmin(dist2)
                self.px[self.nDrag] = event.xdata
                self.py[self.nDrag] = event.ydata
        self.Draw(event)

    def Move(self, event):  # drag point when moving with pressed button
        if self.nDrag != -1:
            self.px[self.nDrag] = event.xdata
            self.py[self.nDrag] = event.ydata
            self.Draw(event)

    def Release(self, event):  # mouse button released
        self.nDrag = -1  # no dragging

    def Draw(self, event):  # re-draw canvas
        # print('draw')
        self.ax.cla()
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

        plt.plot(self.px, self.py, '*-b')

        event.canvas.draw()
        # using plt.show() here will cause stack overflow
    # def calcPoints(self):
    #     p1 = c_[self.px[0], self.py[0]]
    #     p2 = c_[self.px[1], self.py[1]]
    #     p3 = c_[self.px[2], self.py[2]]
    #
    #     self.points = r_[p1,p2,p3]
    #
    #     return p1,p2,p3

    # def calcLength(self):
    #     if self.px.size < 3:
    #         return
    #
    #     self.calcPoints()
    #     self.calcangle()
    #
    #     self.l1 = la.norm(self.points[0]-self.points[1])
    #     self.l2 = la.norm(self.points[1]-self.points[2])
    #
    #     return self.l1, self.l2

   # def calcangle(self):
        #if self.px.size < 3:
        #    return
        #self.theta1 = arctan2(self.py[1] - self.py[0], self.px[1] - self.px[0])
        #self.theta2 = arctan2(self.py[2] - self.py[1], self.px[2] - self.px[1])





if __name__ == '__main__':
    print("LAUNCH")
    fig, ax = plt.subplots()
    ax.set_title('Right button to drag')
    pendel = Pendel(fig, ax)