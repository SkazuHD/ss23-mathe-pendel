from numpy import *
from numpy import linalg as la
import matplotlib.pyplot as plt
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
        self.l1 = 0
        self.l2 = 0
        self.theta1 = 0
        self.theta2 = 0
        self.m1 = 1
        self.m2 = 1
        self.gravity = 9.8

        plt.show()  # draw empty subplot with axes

    def ButtonPress(self, event):  # mouse button pressed event
        # print( event.xdata, event.ydata)
        if event.button == 1:  # add new point
            if self.px.size < 3:
                self.px = r_[self.px, event.xdata]
                self.py = r_[self.py, event.ydata]
                if self.px.size == 3:
                    self.calcLength()
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
    def calcPoints(self):
        p1 = c_[self.px[0], self.py[0]]
        p2 = c_[self.px[1], self.py[1]]
        p3 = c_[self.px[2], self.py[2]]

        self.points = r_[p1,p2,p3]

        return p1,p2,p3

    def calcLength(self):
        if self.px.size < 3:
            return

        self.calcPoints()
        self.calcangle()

        self.l1 = la.norm(self.points[0]-self.points[1])
        self.l2 = la.norm(self.points[1]-self.points[2])

        return self.l1, self.l2

    def calcangle(self):
        if self.px.size < 3:
            return
        self.theta1 = arctan2(self.py[1] - self.py[0], self.px[1] - self.px[0])
        self.theta2 = arctan2(self.py[2] - self.py[1], self.px[2] - self.px[1])



if __name__ == '__main__':
    print("LAUNCH")
    fig, ax = plt.subplots()
    ax.set_title('Right button to drag')
    pendel = Pendel(fig, ax)