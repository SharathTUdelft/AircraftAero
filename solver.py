import matplotlib.pyplot as plt
import numpy as np
from airfoil import AirfoilGen

#

# airfoilcons = AirfoilGen(airfoil_name=airfoil_name, chord=1, nb_points=8)
# airfoilcons.run()
# x_c = airfoilcons.x_foil
# foil = airfoilcons.y_foil
#
# foil[0] = 0
# foil[-1] = 0
#
# panels = airfoilcons.panels
# nb_points = panels + 1
# alpha = np.deg2rad(5)
# cn1 = np.zeros((panels, panels))
# cn2 = np.zeros((panels, panels))
# ct1 = np.zeros((panels, panels))
# ct2 = np.zeros((panels, panels))
# AT = np.zeros((panels, nb_points))
# AN = np.zeros((nb_points, nb_points))
# gamma = np.zeros((nb_points, 1))
# RHS = np.zeros((nb_points, 1))
# V = np.zeros(panels)
# CP = np.zeros(panels)
#
# x_c = airfoilcons.x_foil
# foil = airfoilcons.y_foil
#
# foil[0] = 0
# foil[-1] = 0
#
# panels = airfoilcons.panels
# nb_points = panels + 1
# alpha = np.deg2rad(5)
# cn1 = np.zeros((panels,panels))
# cn2 = np.zeros((panels,panels))
# ct1 = np.zeros((panels,panels))
# ct2 = np.zeros((panels,panels))
# AT = np.zeros((panels, nb_points))
# AN = np.zeros((nb_points, nb_points))
# gamma = np.zeros((nb_points, 1))
# RHS = np.zeros((nb_points, 1))
# V = np.zeros(panels)
# CP = np.zeros(panels)
#
# x_vort = (x_c[1:] + x_c[0:-1]) * 0.5
# y_vort = (foil[1:] + foil[0:-1]) * 0.5
# length = ((foil[1:] - foil[0:-1]) ** 2 + (x_c[1:] - x_c[0:-1]) ** 2) ** 0.5
#
# theta = np.arctan2(foil[1:] - foil[0:-1], x_c[1:] - x_c[0:-1])
# cosine = np.cos(theta)
# sine = np.sin(theta)






















class vlm:
    def __init__(self, panels, x_foil, y_foil, x_mid, y_mid, alpha=5):
        self.x_c = x_foil
        self.foil = y_foil
        self.x_vort = x_mid
        self.y_vort = y_mid
        self.panels = panels
        self.nb_points = self.panels + 1
        self.alpha = np.deg2rad(alpha)
        self.theta = np.arctan2(self.foil[1:] - self.foil[0:-1], self.x_c[1:] - self.x_c[0:-1])
        self.cosine = np.cos(self.theta)
        self.sine = np.sin(self.theta)
        self.length = ((self.foil[1:] - self.foil[0:-1]) ** 2 + (self.x_c[1:] - self.x_c[0:-1]) ** 2) ** 0.5

        self.cn1 = np.zeros((self.panels, self.panels))
        self.cn2 = np.zeros((self.panels, self.panels))
        self.ct1 = np.zeros((self.panels, self.panels))
        self.ct2 = np.zeros((self.panels, self.panels))
        self.AT = np.zeros((self.panels, self.nb_points))
        self.AN = np.zeros((self.nb_points, self.nb_points))
        self.gamma = np.zeros((self.nb_points, 1))
        self.RHS = np.zeros((self.nb_points, 1))
        self.V = np.zeros(self.panels)
        self.CP = np.zeros(self.panels)
        self.CL = 0
        self.CM = 0


    def solve(self):
        for i in range(panels):
            self.RHS[i] = np.sin(self.theta[i] - self.alpha)

        for i in range(self.panels):
            for j in range(self.panels):
                if i == j:
                    self.cn1[i, j] = -1
                    self.cn2[i, j] = 1
                    self.ct1[i, j] = 0.5 * np.pi
                    self.ct2[i, j] = 0.5 * np.pi
                else:
                    A = - (self.x_vort[i] - self.x_c[j]) * self.cosine[j] - (self.y_vort[i] - self.foil[j]) * self.sine[j]
                    B = (self.x_vort[i] - self.x_c[j]) ** 2 + (self.y_vort[i] - self.foil[j]) ** 2
                    C = np.sin(self.theta[i] - self.theta[j])
                    D = np.cos(self.theta[i] - self.theta[j])
                    E = (self.x_vort[i] - self.x_c[j]) * self.sine[j] - (self.y_vort[i] - self.foil[j]) * self.cosine[j] # Check
                    F = np.log(1 + (self.length[j] ** 2 + 2 * self.length[j] * A) / B)
                    G = np.arctan2(E * self.length[j], B + (A * self.length[j])) # Check
                    P = ((self.x_vort[i] - self.x_c[j]) * np.sin(self.theta[i] - 2 * self.theta[j])) + \
                        ((self.y_vort[i] - self.foil[j]) * np.cos(self.theta[i] - 2 * self.theta[j])) # Check
                    Q = ((self.x_vort[i] - self.x_c[j]) * np.cos(self.theta[i] - 2 * self.theta[j])) - \
                        ((self.y_vort[i] - self.foil[j]) * np.sin(self.theta[i] - 2 * self.theta[j])) # Check

                    self.cn2[i, j] = D + (0.5 * Q * F / self.length[j]) - (A * C + D * E) * G / self.length[j]
                    self.cn1[i, j] = 0.5 * D * F + C * G - self.cn2[i, j]
                    self.ct2[i, j] = C + 0.5 * P * F / self.length[j] + (A * D - C * E) * G / self.length[j]
                    self.ct1[i, j] = 0.5 * C * F - D * G - self.ct2[i, j]

                for j in range(self.panels):
                    self.AN[j, 0] = self.cn1[j, 0]
                    self.AN[j, -1] = self.cn2[j, -1]

                    self.AT[j, 0] = self.ct1[j, 0]
                    self.AT[j, -1] = self.ct2[j, -1]

                    for k in range(1, self.panels):
                        self.AN[j, k] = self.cn1[j, k] + self.cn2[j, k-1]
                        self.AT[j, k] = self.ct1[j, k] + self.ct2[j, k-1]
        self.AN[-1, 0] = 1
        self.AN[-1, -1] = 1
        for i in range(1, self.panels):
            self.AN[-1, i] = 0
        self.RHS[-1] = 0
        self.gamma = np.linalg.solve(self.AN, self.RHS)
        for i in range(self.panels):
            self.V[i] = np.cos(self.theta[i] - self.alpha)
            for j in range(self.nb_points):
                self.V[i] +=  self.AT[i, j] * self.gamma[j]
                self.CP[i] = 1 - self.V[i] ** 2
            self.CL += 2 * (self.gamma[i] + self.gamma[i + 1]) * self.length[i] * np.pi
            self.CM += self.CP[i] * self.length[i] * (((self.x_vort[i] - 0.25) * np.cos(self.theta[i])) + self.y_vort[i] * np.sin(self.theta[i]))

        return self.CP, self.CL, self.CM, self.x_vort

    def cp_plot(self):
        # plt.plot(self.x_vort, self.CP, "--")
        # plt.show()
        return self.x_vort, self.CP


if __name__ == "__main__":

    airfoilcons = AirfoilGen(airfoil_name="NACA0012", chord=1, nb_points=121)
    x_foil, y_foil, x_mid, y_mid, panels = airfoilcons.run()


    foil_solver = vlm(x_foil=x_foil, y_foil=y_foil,x_mid=x_mid, y_mid=y_mid, panels=panels)
    foil_solver.solve()
    _, cp_foil =foil_solver.cp_plot()
    with open("cp_file_0012", "r+") as f:
        x_cord = []
        y_cord = []
        cp_cord = []
        for line in f:

            x, y, cp = line.strip().split()

            x_cord.append(float(x))
            y_cord.append(float(y))
            cp_cord.append(float(cp))



    ax = plt.gca()
    ax.plot(x_cord, cp_cord, "b", x_mid, cp_foil, "-*r" )

    plt.xlabel(r' Normalized chordwise distribution ($\frac{x}{c}$)  $\longrightarrow$ ', fontsize= 16)
    plt.xticks(np.arange(min(x_cord), max(x_cord) , 0.1))

    plt.ylabel(r' Coefficient of Pressure ($C_p) $ $ \longrightarrow$   ' , fontsize= 16)
    plt.title(r'$C_p$ vs $\frac{x}{c}$ for NACA0012 Airfoil', fontdict = {'fontsize' : 18})
    plt.legend([r'NACA0012 $C_p$ from XFOIL', r'NACA0012 $C_p$ from Code'])
    plt.show()
