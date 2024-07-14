import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import G
from matplotlib.animation import FuncAnimation

def peri2ECI(a, e, i, omega, theta, OMEGA, mu):
    i = np.deg2rad(i)
    OMEGA = np.deg2rad(OMEGA)
    omega = np.deg2rad(omega)
    theta = np.deg2rad(theta)

    h = np.sqrt(a * mu * (1 - e ** 2))

    rp = h ** 2 / mu / (1 + e * np.cos(theta)) * np.array([np.cos(theta), np.sin(theta), 0])
    vp = mu / h * np.array([-np.sin(theta), e + np.cos(theta), 0])

    a1 = np.array([[np.cos(omega), np.sin(omega), 0], [-np.sin(omega), np.cos(omega), 0], [0, 0, 1]])
    a2 = np.array([[1, 0, 0], [0, np.cos(i), np.sin(i)], [0, -np.sin(i), np.cos(i)]])
    a3 = np.array([[np.cos(OMEGA), np.sin(OMEGA), 0], [-np.sin(OMEGA), np.cos(OMEGA), 0], [0, 0, 1]])

    A = (a1 @ a2 @ a3).T

    ri = A @ rp
    vi = A @ vp

    return np.concatenate((ri, vi))


def drawsphere(x, y, z, r, ax):
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    a = r * np.cos(u) * np.sin(v) + x
    b = r * np.sin(u) * np.sin(v) + y
    c = r * np.cos(v) + z
    ax.plot_surface(a, b, c, color='b', alpha=0.3)

