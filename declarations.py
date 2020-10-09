import numpy as np
from airfoil import AirfoilGen

airfoil_name = input("Enter the airfoil name:")

airfoilcons = AirfoilGen(airfoil_name=airfoil_name, chord=1, nb_points=8)
airfoilcons.run()
x_c = airfoilcons.x_foil
foil = airfoilcons.y_foil

foil[0] = 0
foil[-1] = 0

panels = airfoilcons.panels
nb_points = panels + 1
alpha = np.deg2rad(5)
cn1 = np.zeros((panels,panels))
cn2 = np.zeros((panels,panels))
ct1 = np.zeros((panels,panels))
ct2 = np.zeros((panels,panels))
AT = np.zeros((panels, nb_points))
AN = np.zeros((nb_points, nb_points))
gamma = np.zeros((nb_points, 1))
RHS = np.zeros((nb_points, 1))
V = np.zeros(panels)
CP = np.zeros(panels)

