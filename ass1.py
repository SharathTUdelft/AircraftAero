import matplotlib.pyplot as plt
from declarations import *

x_vort = (x_c[1:] + x_c[0:-1]) * 0.5
y_vort = (foil[1:] + foil[0:-1]) * 0.5
length = ((foil[1:] - foil[0:-1]) ** 2 + (x_c[1:] - x_c[0:-1]) ** 2) ** 0.5

theta = np.arctan2(foil[1:] - foil[0:-1], x_c[1:] - x_c[0:-1])
cosine = np.cos(theta)
sine = np.sin(theta)

for i in range(panels):
    RHS[i] = np.sin(theta[i] - alpha)


class vlm:
    def __init__(self):
        pass

    def solve(self):
        for i in range(panels):
            for j in range(panels):
                if i == j:
                    cn1[i, j] = -1
                    cn2[i, j] = 1
                    ct1[i, j] = 0.5 * np.pi
                    ct2[i, j] = 0.5 * np.pi
                else:
                    A = - (x_vort[i] - x_c[j]) * cosine[j] - (y_vort[i] - foil[j]) * sine[j]
                    B = (x_vort[i] - x_c[j]) ** 2 + (y_vort[i] - foil[j]) ** 2
                    C = np.sin(theta[i] - theta[j])
                    D = np.cos(theta[i] - theta[j])
                    E = (x_vort[i] - x_c[j]) * sine[j] - (y_vort[i] - foil[j]) * cosine[j] # Check
                    F = np.log(1 + (length[j] ** 2 + 2 * length[j] * A) / B)
                    G = np.arctan2(E * length[j], B + (A * length[j])) # Check
                    P = ((x_vort[i] - x_c[j]) * np.sin(theta[i] - 2 * theta[j])) + \
                        ((y_vort[i] - foil[j]) * np.cos(theta[i] - 2 * theta[j])) # Check
                    Q = ((x_vort[i] - x_c[j]) * np.cos(theta[i] - 2 * theta[j])) - \
                        ((y_vort[i] - foil[j]) * np.sin(theta[i] - 2 * theta[j])) # Check

                    cn2[i, j] = D + (0.5 * Q * F / length[j]) - (A * C + D * E) * G / length[j]
                    cn1[i, j] = 0.5 * D * F + C * G - cn2[i, j]
                    ct2[i, j] = C + 0.5 * P * F / length[j] + (A * D - C * E) * G / length[j]
                    ct1[i, j] = 0.5 * C * F - D * G - ct2[i, j]

                for j in range(panels):
                    AN[j, 0] = cn1[j, 0]
                    AN[j, -1] = cn2[j, -1]

                    AT[j, 0] = ct1[j, 0]
                    AT[j, -1] = ct2[j, -1]

                    for k in range(1, panels):
                        AN[j, k] = cn1[j, k] + cn2[j, k-1]
                        AT[j, k] = ct1[j, k] + ct2[j, k-1]
        AN[-1, 0] = 1
        AN[-1, -1] = 1
        for i in range(1, panels):
            AN[-1, i] = 0
        RHS[-1] = 0
        gamma = np.linalg.solve(AN, RHS)
        for i in range(panels):
            V[i] = np.cos(theta[i] - alpha)
            for j in range(nb_points):
                V[i] = V[i] + AT[i, j] * gamma[j]
                CP[i] = 1 - V[i] ** 2



if __name__ == "__main__":

    solver = vlm()
    solver.solve()
    plt.plot(x_vort, CP, "*--")
    plt.show()



