import numpy as np
import matplotlib.pyplot as plt


class Airfoil:

    def __init__(self, airfoil_name="NACA0012", airfoil_type="symmetric", resolution=10, alpha = 10):
        self.type = airfoil_type
        self.name = airfoil_name
        self.resolution = resolution
        self.alpha = np.deg2rad(alpha)

    def name_scrape(self):
        return self.name[4], self.name[5], self.name[6], self.name[7]

    def construct_airfoil(self):
        if self.type == "symmetric":
            theta = np.linspace(0, np.pi, 50)
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
            print(self.x_c)
            pass

    def vortices(self):
        self.construct_airfoil()
        self.x_vort = (self.x_c[1:] + self.x_c[0:-1]) * 0.5
        self.y_vort = (self.foil[1:] + self.foil[0:-1]) * 0.5
        self.length = (self.foil[1:] - self.foil[0:-1]) ** 2 + (self.x_c[1:] - self.x_c[0:-1]) ** 2
        pass

    def angles(self):
        self.theta = np.arctan2(self.foil[1:] - self.foil[0:-1], self.x_c[1:], self.x_c[0:-1])
        self.cosine = np.cos(self.theta)
        self.sine = np.sin(self.theta)


    def solve(self):
        self.RHS = np.sin(self.theta - self.alpha)



    def plot(self):
        self.construct_airfoil()
        self.vortices()
        plt.plot(self.x_c, self.foil, "*--", self.x_vort, self.y_vort, "*r")
        plt.show()


if __name__ == "__main__":
    a = Airfoil()
    a.vortices()
    a.plot()

