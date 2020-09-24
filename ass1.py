import numpy as np
import matplotlib.pyplot as plt
class airfoil:

    def __init__(self, airfoil_name="NACA0012", airfoil_type="symmetric", resolution=10):
        self.type = airfoil_type
        self.name = airfoil_name
        self.resolution = resolution

    def name_scrape(self):
        return self.name[4], self.name[5], self.name[6], self.name[7]


    def construct_airfoil(self):
        if self.type == "symmetric":
            theta = np.linspace(0, np.pi, 50)
            x_c = 0.5 * (1 - np.cos(theta))
            a, b, c, d =  self.name_scrape()
            t = int(c+d)/100
            top = 5*t*(0.2969*(x_c**0.5) - 0.1260*x_c - 0.3516*(x_c**2) + 0.2843*(x_c**3) - 0.1015*(x_c**4))
            bot = - top
            top[-1] = 0
            x_c_new = np.append(x_c, np.flip(x_c)[1:])
            foil = np.append(top, np.flip(bot)[1:])
            return foil, x_c_new

    def vortices(self):
        foil, x_c_new  = self.construct_airfoil()
        x_vort = (x_c_new[1:] + x_c_new[0:-1]) * 0.5
        y_vort = (foil[1:] + foil[0:-1]) * 0.5
        return x_vort, y_vort

    def plot(self):
        foil, x_c = self.construct_airfoil()
        x_vort, y_vort = self.vortices()
        plt.plot(x_c, foil, "*--", x_vort, y_vort, "*r")
        plt.show()




if __name__ == "__main__":
    a = airfoil()
    a.vortices()
    a.plot()

