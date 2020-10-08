import numpy as np
import matplotlib.pyplot as plt


class AirfoilGen:

    def __init__(self, airfoil_name="NACA2412", nb_points=51, alpha=8, chord=1):
        self.name = airfoil_name
        self.nb_points = nb_points
        self.panels = 2 * (self.nb_points - 1)
        self.alpha = alpha
        self.chord = chord
        self.a0 = 0.2969
        self.a1 = -0.1260
        self.a2 = -0.3516
        self.a3 = 0.2843
        self.a4 = -0.1015  #open trailing edge
        self.a4_t = -0.1036  #closed trailing edge
        self.y_cam = np.zeros(self.nb_points)
        self.y_t = np.zeros(self.nb_points)
        self.dy_dx_cam = np.zeros(self.nb_points)
        self.theta = np.zeros(self.nb_points)
        self.x_u = np.zeros(self.nb_points)
        self.x_l = np.zeros(self.nb_points)
        self.y_u = np.zeros(self.nb_points)
        self.y_l = np.zeros(self.nb_points)
        self.x_cord = np.zeros(self.nb_points)
        self.y_cord = np.zeros(self.nb_points)
        self.run()

    def parsefoil(self):
        camber_max = (float(self.name[4]) / 100) * self.chord
        camber_position = self.chord * (float(self.name[5]) / 10)
        thickness = self.chord * (float(self.name[6] + self.name[7]) / 100)
        return camber_max, camber_position, thickness

    def construct_foil(self):
        cam_max, cam_pos, t = self.parsefoil()
        angle = np.linspace(0, np.pi, self.nb_points, endpoint=True)
        self.x_cord = self.chord * ((np.cos(angle) / 2) + 0.5) #TODO: Check chord = 2

        for i in range(self.nb_points):
            term0 = self.a0 * self.x_cord[i] ** 0.5
            term1 = self.a1 * self.x_cord[i]
            term2 = self.a2 * self.x_cord[i] ** 2
            term3 = self.a3 * self.x_cord[i] ** 3
            term4 = self.a4_t * self.x_cord[i] ** 4

            self.y_t[i] = 5 * t * (term0 + term1 + term2 + term3 + term4)

            if (self.x_cord[i] >= 0 and self.x_cord[i] < cam_pos):

                self.y_cam[i] = (cam_max/cam_pos**2) * ((2*cam_pos * self.x_cord[i]) - self.x_cord[i]**2)
                self.dy_dx_cam[i] = ((2*cam_max)/cam_pos**2)*(cam_pos - self.x_cord[i])

            elif self.x_cord[i] >= cam_pos and self.x_cord[i] <= self.chord:

                self.y_cam[i] = (cam_max / (1 - cam_pos) ** 2) * (1 - 2 * cam_pos + (2 * cam_pos * self.x_cord[i]) - self.x_cord[i] ** 2)
                self.dy_dx_cam[i] = ((2 * cam_max) / (1 - cam_pos) ** 2) * (cam_pos - self.x_cord[i])

            self.theta[i] = np.arctan(self.dy_dx_cam[i])

        # upper surface
            self.x_u[i] = self.x_cord[i] - self.y_t[i] * np.sin(self.theta[i])
            self.y_u[i] = self.y_cam[i] + self.y_t[i] * np.cos(self.theta[i])

        # lower surface
            self.x_l[i] = self.x_cord[i] + self.y_t[i] * np.sin(self.theta[i])
            self.y_l[i] = self.y_cam[i] - self.y_t[i] * np.cos(self.theta[i])

    def plot_foil(self):
        self.construct_foil()

        plt.plot(self.x_cord, self.y_cam, "*", self.x_u, self.y_u, "-*b", self.x_l, self.y_l, "-*r")
        plt.axis('equal')
        plt.show()

    def airfoil(self):
        return self.x_u, self.y_u, self.x_l, self.y_l

    def coordinates(self):
        cord_x = np.flip(self.x_u)
        cord_y = np.flip(self.y_u)
        self.x_foil = np.append(self.x_l, cord_x[1:])
        self.y_foil = np.append(self.y_l, cord_y[1:])
        pass

    def run(self):
        self.construct_foil()
        self.coordinates()


if __name__ == "__main__":
    a = AirfoilGen(chord=1)
    print(a.y_foil)
    # x, y = a.coordinates()
    # plt.plot(x , y)
    # plt.axis("equal")
    # plt.show()

