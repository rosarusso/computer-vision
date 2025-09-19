
import numpy as np
import matplotlib.pyplot as plt
import cv2  # or from PIL import Image
from scipy.io import loadmat

# Normalize function
def normalize(m):
    return m / m[-1, :]

# Load image
def load_image(filename):
    img = cv2.imread(filename)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)  # Convert BGR to RGB for matplotlib
    return img

# Generate sphere points
def generate_sphere(radius=1.8, resolution=30):
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
    return x, y, z

# Animation 1: Moving point along y-axis
def animation_1():
    img = load_image('Image1.jpg')
    data = loadmat('Calib_direct_1.mat')
    P = data['P']

    Mi = np.array([
        [0, 0, 0, 1],
        [0, 0, 8.3, 1],
        [5.3, 0, 8.3, 1],
        [5.3, 0, 0, 1],
        [0, 14.8, 8.3, 1],
        [5.3, 14.8, 8.3, 1]
    ])

    point = np.array([5.3, 0, 8.3, 1])
    m = P @ point
    m = m / m[2]

    plt.imshow(img)
    # plt.plot(m[0], m[1], 'ro', markersize=5, linewidth=2)
    pl, = plt.plot([m[0]], [m[1]], 'ro', markersize=5, linewidth=2)

    for i in range(2, 52):
        point = np.array([5.3, 14.8 * (i - 1) / 50, 8.3, 1])
        m = P @ point
        m = m / m[2]
        pl.set_data([m[0]], [m[1]])
        plt.pause(0.2)

    plt.show()

# Animation 2: Moving sphere along a trajectory
def animation_2():
    img = load_image('Image2.jpg')
    data = loadmat('Calib_direct_2.mat')
    P = data['P']

    X, Y, Z = generate_sphere()

    def center(t):
        x = -0.47 * t + 2.9
        y = 9.3
        z = 2 * t**2 - 1.64 * t - 1.8
        return x, y, z

    for t in np.arange(0, 3.1, 0.1):
        plt.clf()
        plt.imshow(img)
        x, y, z = center(t)
        Sph = np.array([
            1.8 * X.ravel() + x,
            1.8 * Y.ravel() + y,
            1.8 * Z.ravel() + z,
            np.ones(X.size)
        ])
        sph = normalize(P @ Sph)
        plt.scatter(sph[0, :], sph[1, :], c='red', s=1)
        plt.pause(0.01)

    for t in np.arange(3, -0.1, -0.1):
        plt.clf()
        plt.imshow(img)
        x, y, z = center(t)
        Sph = np.array([
            1.8 * X.ravel() + x,
            1.8 * Y.ravel() + y,
            1.8 * Z.ravel() + z,
            np.ones(X.size)
        ])
        sph = normalize(P @ Sph)
        plt.scatter(sph[0, :], sph[1, :], c='red', s=1)
        plt.pause(0.005)

    plt.show()

# Animation 3: Bouncing sphere
def animation_3():
    img = load_image('Image1.jpg')
    data = loadmat('Calib_direct_1.mat')
    P = data['P']

    box_limits_y = [2.6, 14.8]
    pos = np.array([2.5, 7, 4])
    vel = np.array([0, 0.15, 0])
    raggio = 1.8

    X, Y, Z = generate_sphere(radius=raggio, resolution=30)

    plt.figure()
    for frame in range(300):
        pos += vel

        if pos[1] - raggio < box_limits_y[0]:
            vel[1] = abs(vel[1])
        elif pos[1] + raggio > box_limits_y[1]:
            vel[1] = -abs(vel[1])

        Sph = np.array([
            raggio * X.ravel() + pos[0],
            raggio * Y.ravel() + pos[1],
            raggio * Z.ravel() + pos[2],
            np.ones(X.size)
        ])

        proj = P @ Sph
        proj = proj / proj[2, :]

        plt.clf()
        plt.imshow(img)
        plt.scatter(proj[0, :], proj[1, :], c='red', s=10)
        plt.pause(0.02)

    plt.show()

# Run animations
if __name__ == "__main__":
    animation_1()
    animation_2()
    animation_3()

