""" quite sensetive to parameters n, k and R
    Easily blows up. Unstable. Should not have
    used Euler method for integration"""

import pygame
from pygame.locals import *
import numpy as np

pygame.init()
sc = pygame.display.set_mode((600, 600))
black = (30, 60, 30)
red = (200, 10, 10)
sc.fill(black)
# initial conditions
n = 100
R = 100
l0 = 2 * np.pi * R / n * 0.95  # under tension
m = 1
k_ball_points = 5000  # spring constant for springs between masses
k_ball_ground = 1000  # spring constant for mass-ground collision
r = 5  # radius of the masses
g = -1  # grav.
c_inter = 0.02  # Energy dissipation factor for spring between masses
c_air = 0.0001  # Energy dissipation factor due to air drag
# modeling the ball
theta = np.arange(n) * 2 * np.pi / n
x = R * np.cos(theta)
y = R * np.sin(theta)
# initial velocities and accelerations
vx = np.ones(n) * 5
vy = np.ones(n) * -25
axn = np.zeros(n)
ayn = np.zeros(n)
dt = 0.01


def draw_circle():
    """drawing the masses that make up the sphere"""
    for i in range(n):
        pygame.draw.circle(sc, red, (int(x[i] + 300), int(y[i] + 300)), r)


def particles_acc():
    """Calculating the acceleration of each particle given their
    positions with respect to each other and the spring constant.

    Returns:
        tuple: tuple consisting of ax and ay each an array of length "n".
    """
    # distance of each mass with its left neighbour
    dxl = np.diff(x, prepend=x[-1])
    dyl = np.diff(y, prepend=y[-1])
    dl = np.linalg.norm((dxl, dyl), axis=0)
    # distance of each mass with its right neighbour
    dxr = -np.diff(x, append=x[0])
    dyr = -np.diff(y, append=y[0])
    dr = np.linalg.norm((dxr, dyr), axis=0)
    # accelerations due to spring connected to the masses
    al = -k_ball_points * (dl - l0) / m
    ar = -k_ball_points * (dr - l0) / m
    # vector decomposition
    alx = al * dxl / dl
    aly = al * dyl / dl
    arx = ar * dxr / dr
    ary = ar * dyr / dr
    # finding net acceleration
    ax = arx + alx
    ay = ary + aly
    return ax, ay


# PV = constant for gas inside the ball gives:
# F_net/A * V = cosntant ==> for sphere: F_net=constant/radius
# hence: cosntant = F_net*radius
# for each particle: constant = constant/n
constant = np.linalg.norm(particles_acc(), axis=0).sum() / n * R

# -----------------
# ---------------
draw_circle()
pygame.display.update()
# clock = pygame.time.Clock()
cont = True
j = 0
while cont:
    # position and velocity of ball's center
    xc, yc, vxc, vyc = x.sum() / n, y.sum() / n, vx.sum() / n, vy.sum() / n
    # obtaining acceleration values
    a = particles_acc()
    # distance from masses to the center of the ball
    d = np.linalg.norm((x - xc, y - yc), axis=0)
    # collision with walls
    conditions_x = [x - r <= -300, (x - r > -300) & (x < 300 - r), x >= 300 - r]
    ax = (
        a[0]
        + constant / d * (x - xc) / d
        - c_inter * (vx - vxc) ** 2 / m * np.copysign(1, vx - vxc)
        - c_air * abs(vx) / m * np.copysign(1, vx)
    )
    choices_x = [
        ax - k_ball_ground * 10 * (x - r + 300) / m,
        ax,
        ax - k_ball_ground * (x - (300 - r)) / m,
    ]

    conditions_y = [y - r <= -300, (y - r > -300) & (y < 300 - r), y >= 300 - r]
    ay = (
        a[1]
        + constant / d * (y - yc) / d
        - g
        - c_inter * (vy - vyc) ** 2 / m * np.copysign(1, vy - vyc)
        - c_air * abs(vy) / m * np.copysign(1, vy)
    )
    choices_y = [
        ay - k_ball_ground * (y - r + 300) / m,
        ay,
        ay - k_ball_ground * (y - (300 - r)) / m,
    ]
    ######################
    # x and y componnents of the acceleration at any given moment
    axn = np.select(conditions_x, choices_x)
    ayn = np.select(conditions_y, choices_y)
    # updating velocities and positions
    vx[:] += axn * dt
    vy[:] += ayn * dt
    x[:] += vx * dt
    y[:] += vy * dt
    j += 1
    if j % 10 == 0:  # updating the screen after each 10 loops
        sc.fill(black)
        draw_circle()
        pygame.display.update()
    # event loop to close the program
    for event in pygame.event.get():
        if event.type == KEYDOWN:
            if event.key == pygame.K_ESCAPE or event.key == K_q:
                cont = False
        if event.type == QUIT:
            cont = False

pygame.quit()
