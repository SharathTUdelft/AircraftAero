from solver import *



if __name__ == "__main__":

    solver = vlm()
    solver.solve()
    plt.plot(x_vort, -CP, "*--")
    plt.xlabel(airfoil_name + " x/c")
    plt.show()
