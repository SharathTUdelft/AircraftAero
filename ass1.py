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
            return top, bot, x_c

    def plot(self):
        top, bot, x_c = self.construct_airfoil()
        x_c_new = np.append(x_c, np.flip(x_c))
        foil = np.append(top, np.flip(bot))
        plt.plot(x_c_new, foil, "*--")
        plt.show()




if __name__ == "__main__":
    a = airfoil()

    a.plot()
