import numpy as np
import matplotlib.pyplot as plt


class Airfoil:

    def __init__(self, airfoil_name="NACA0012", airfoil_type="symmetric", resolution=10, alpha = 8):
        self.type = airfoil_type
        self.name = airfoil_name
        self.resolution = resolution
        self.alpha = np.deg2rad(alpha)

        self.construct_airfoil()

        self.cn1 = np.zeros((self.panels, self.panels))
        self.cn2 = np.zeros((self.panels, self.panels))
        self.ct1 = np.zeros((self.panels, self.panels))
        self.ct2 = np.zeros((self.panels, self.panels))
        self.AT = np.zeros((self.panels, self.points))
        self.AN = np.zeros((self.points, self.points))
        self.gamma = np.zeros((self.points, 1))
        self.RHS = np.zeros((self.points, 1))
        self.V = np.zeros(self.panels)
        self.CP = np.zeros(self.panels)
        print("fsd")

    def name_scrape(self):
        return self.name[4], self.name[5], self.name[6], self.name[7]

    def construct_airfoil(self):
        if self.type == "symmetric":
            theta = np.linspace(0, np.pi, self.resolution)
            x_c = 0.5 * (1 - np.cos(theta))
            x_c = np.flip(x_c) # flipped so as per fotran program
            a, b, c, d = self.name_scrape()
            t = int(c+d)/100
            top = 5*t*(0.2969*(x_c**0.5) - 0.1260*x_c - 0.3516*(x_c**2) + 0.2843*(x_c**3) - 0.1015*(x_c**4))
            bot = - top
            top[-1] = 0
            x_c_new = np.append(x_c, np.flip(x_c)[1:])
            self.foil = np.append(top, np.flip(bot)[1:])
            self.x_c = x_c_new
            self.points = len(self.foil)
            self.panels = self.points - 1
            pass

    def vortices(self):
        self.construct_airfoil()
        self.x_vort = (self.x_c[1:] + self.x_c[0:-1]) * 0.5
        self.y_vort = (self.foil[1:] + self.foil[0:-1]) * 0.5
        self.length = (self.foil[1:] - self.foil[0:-1]) ** 2 + (self.x_c[1:] - self.x_c[0:-1]) ** 2

    def angles(self):
        self.theta = np.arctan2(self.foil[1:] - self.foil[0:-1], self.x_c[1:], self.x_c[0:-1])
        self.cosine = np.cos(self.theta)
        self.sine = np.sin(self.theta)

    def solve(self):
        self.construct_airfoil()
        self.vortices()
        self.angles()
        for i in range(self.panels):
            self.RHS[i] = self.theta[i] - self.alpha

        for i in range(self.panels):
            for j in range(self.panels):
                if i == j:
                    self.cn1[i,j] = -1
                    self.cn2[i,j] = -1
                    self.ct1[i,j] = 0.5 * np.pi
                    self.ct2[i,j] = 0.5 * np.pi
                else:
                    A = - (self.x_vort[i] - self.x_c[j]) * self.cosine[j] - (self.y_vort[i] - self.foil[j]) * self.sine[j]
                    B = (self.x_vort[i] - self.x_c[j]) ** 2 + (self.y_vort[i] - self.foil[j]) ** 2
                    C = np.sin(self.theta[i] - self.theta[j])
                    D = np.cos(self.theta[i] - self.theta[j])
                    E = (self.x_vort[i] - self.x_c[j]) * self.sine[j] - (self.y_vort[i] - self.foil[j]) * self.cosine[j]
                    F = np.log(1 + self.length[j]* (self.length[j] + 2*A)/B )
                    G = np.arctan2(E*self.length[j], B+(A*self.length[j]) )
                    P = ((self.x_vort[i] - self.x_c[j]) * np.sin(self.theta[i] - 2*self.theta[j])) + \
                        ((self.y_vort[i] - self.foil[j]) * np.cos(self.theta[i] - 2*self.theta[j]))

                    Q = ((self.x_vort[i] - self.x_c[j]) * np.cos(self.theta[i] - 2 * self.theta[j])) - \
                        ((self.y_vort[i] - self.foil[j]) * np.sin(self.theta[i] - 2 * self.theta[j]))

                    self.cn2[i, j] = D + (0.5 * Q * F/self.length[j]) - (A*C + D*E) * G/self.length[j]
                    self.cn1[i, j] = 0.5 * D * F + C*G - self.cn2[i, j]
                    self.ct2[i, j] = C + 0.5*P*F/self.length[j] + (A*D - C*E) * G/self.length[j]
                    self.ct1[i, j] = 0.5*C*F - D*G - self.ct2[i, j]

        for i in range(self.panels):
            self.AN[i, 0] = self.cn1[i, 0]
            self.AN[i, -1] = self.cn2[i, -1]

            self.AT[i, 0] = self.ct1[i, 0]
            self.AT[i, -1] = self.ct2[i, -1]

            for j in range(1, self.panels):
                self.AN[i, j] = self.cn1[i, j] + self.cn2[i, j-1]
                self.AT[i, j] = self.ct1[i, j] + self.ct2[i, j-1]

        self.AN[-1, 0] = 1
        self.AN[-1, -1] = 1
        for j in range(1, self.panels):
            self.AN[-1, j] = 0

        print(self.AN.shape)
        self.gamma = np.linalg.solve(self.AN, self.RHS)
        print(self.gamma)

        for i in range(self.panels):
            self.V[i] = np.cos(self.theta[i] - self.alpha)
            for j in range(self.points):
                self.V[i] = self.V[i] + self.AT[i, j]*self.gamma[j]
                self.CP[i] = 1 - self.V[i]**2



    def plot(self):
        self.construct_airfoil()
        self.vortices()
        plt.plot(self.x_c, self.foil, "*--", self.x_vort, self.y_vort, "*r")
        plt.show()


if __name__ == "__main__":
    a = Airfoil()

    a.solve()

