# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import logging
import unittest
import numpy as np
import matplotlib.pyplot as plt
from Notifiers import SeriesNotifier, FiringTimesNotifier
from IntegrateAndFire import IntegrateAndFire
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorData




class HeartActionForceTest(unittest.TestCase):
    """
    This class tests all features of a complicated IntegrateAndFire class.
    1. Unperturbed rate.
    2. Perturbed rate which visually responds to stimulations.
    3. Perturbed rate which is visually coupled to a different rhythm.
    4. Perturbed rate which may be plugged into Kicked Windkessel model.
    """

    def test_UnperturbedRate(self):
        data = RungeKutta45IntegratorData(7,0.0)
        iaf = IntegrateAndFire()
        fireNotifier = FiringTimesNotifier()
        iaf.Notify = fireNotifier.Notify
        npoints = 1000
        period = 1.0
        seriesNotifier = SeriesNotifier(7,npoints)
        allTimes = np.linspace(0.0,10.0,npoints)
        for t in allTimes:
            data.t = t
            data.y[3] = 10.0 * np.cos(2 * np.pi * t / period)
            data.y[0] = data.y[1] = data.y[2] = 0.0
            #print(data)
            iaf.ApplyDrive(data) # will open or not.
            seriesNotifier.Notify(data)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(2,1,1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes,seriesNotifier.GetVar(3))
        plt.subplot(2,1,2)
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g")
        plt.plot(fireNotifier.firingTimes,spikes,"ro")
        plt.savefig("test_UnperturbedRate.png")
        #plt.show()

 
def main():

    logging.basicConfig(level=logging.INFO)
    unittest.main()

if __name__ == "__main__":
    main()
