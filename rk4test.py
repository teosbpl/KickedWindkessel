#!/usr/bin/env python

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

def fLorenz(t, Y):
    x, y, z = Y
    return [10*(y-x), x*(28-z)-y, x*y-(8./3.)*z]

def main():
    # create arrays and set initial condition
    tspan = np.arange(0, 100, 0.01)
    Y = np.empty((3, tspan.size))
    Y[:, 0] = [1.5961, 0.1859, 4.7297] # random (x,y,z)

    # create explicit Runge-Kutta integrator of order (4)5
    r = ode(fLorenz).set_integrator('dopri5')
    r.set_initial_value(Y[:, 0], tspan[0])

    # run the integration
    for i, t in enumerate(tspan):
        if not r.successful():
            break
        if i == 0:
            continue # skip the initial position
        r.integrate(t)
        Y[:, i] = r.y

    # plot the result in x-y plane
    plt.plot(Y[0, :], Y[1, :])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

if __name__ == '__main__':
    main()