"""This model is made by considering a set 
of springs connected to the center of the circle
Note: needs documentation and improvement
"""
import pygame
from pygame.locals import *
from math import pi, sin, cos
import numpy as np

pygame.init()
sc = pygame.display.set_mode((600, 600))
red = (255, 0, 0)
blue = (100, 100, 255)
yellow = (255, 255, 0)
white = (255, 255, 255)
black = (30, 60, 30)
sc.fill(black)
n = 100
R = 100
l = 2 * pi * R * 0.3
l0 = l / n
m = 1
k = 100
r = 5

theta = 2 * pi / n
x = np.array([R * sin(i * theta) for i in range(n)])
y = np.array([R * cos(i * theta) for i in range(n)])
vx = np.zeros(n)
vy = np.zeros(n) - 25
g = -1
c = 0.01


def draw_circle():
    for i in range(n):
        pygame.draw.circle(sc, blue, (int(x[i] + 300), int(y[i] + 300)), r)


def erase_circle():
    for i in range(n):
        pygame.draw.circle(sc, black, (int(x[i] + 300), int(y[i] + 300)), r)


# -----------------
# ---------------
draw_circle()
pygame.display.update()
dt = 0.01
cont = True
while cont:
    erase_circle()
    dxl = np.diff(x, prepend=x[-1])
    dxr = -np.diff(x, append=x[0])
    dyl = np.diff(y, prepend=y[-1])
    dyr = -np.diff(y, append=y[0])
    dl = np.linalg.norm((dxl, dyl), axis=0)
    dr = np.linalg.norm((dxr, dyr), axis=0)
    al = -k * (dl - l0) / m
    ar = -k * (dr - l0) / m
    alx = al * dxl / dl
    aly = al * dyl / dl
    arx = ar * dxr / dr
    ary = ar * dyr / dr
    xc = sum(x) / n
    yc = sum(y) / n
    vxc = 0
    vyc = sum(vy) / n
    D = np.linalg.norm((x - xc, y - yc), axis=0)
    Ax = -k * (D - R) * (x - xc) / (m * D) - c * 10 * (vx - vxc) ** 2 / m * np.copysign(
        1, vx - vxc
    )
    Ay = -k * (D - R) * (y - yc) / (m * D) - c * 10 * (vy - vyc) ** 2 / m * np.copysign(
        1, vy - vyc
    )

    ax = arx + alx + Ax - c * abs(vx) / m * np.copysign(1, vx)
    ay = ary + aly - g + Ay - c * abs(vy) / m * np.copysign(1, vy)

    condition_bottom = y + vy * dt >= 300 - r
    ay[condition_bottom] += (
        -k * 10 * (y[condition_bottom] + vy[condition_bottom] * dt - 300 + r) / m
    ) - c * abs(vy[condition_bottom] - vyc) * np.copysign(1, vy[condition_bottom] - vyc)

    vx += ax * dt
    vy += ay * dt
    x += vx * dt
    y += vy * dt

    draw_circle()
    pygame.display.update()

    for event in pygame.event.get():
        if event.type == KEYDOWN:
            if event.key == pygame.K_ESCAPE or event.key == K_q:
                cont = False
        if event.type == QUIT:
            cont = False
pygame.quit()
