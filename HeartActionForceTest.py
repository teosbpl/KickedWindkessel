# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import sys
import logging
import unittest
import numpy as np
import matplotlib.pyplot as plt
from Notifiers import SeriesNotifier, FiringTimesNotifier
from HeartActionForce import RectangularHeartActionForce, RespiratoryDelayedSmearedHeartActionForce
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorData

class HeartActionForceTest(unittest.TestCase):

    def test_RectangularHeartActionForce(self):
        data = RungeKutta45IntegratorData(5,0.0)
        force = RectangularHeartActionForce()
        fireNotifier = FiringTimesNotifier()
        force.Notify = fireNotifier.Notify
        npoints = 1000
        period = 1.0
        seriesNotifier = SeriesNotifier(5,npoints)
        allTimes = np.linspace(0.0,10.0,npoints)
        for t in allTimes:
            data.t = t
            data.y[3] = 10.0 * np.cos(2 * np.pi * t / period)
            data.y[0] = data.y[1] = data.y[2] = 0.0
            print(data)
            force.ApplyDrive(data) # will open or not.
            seriesNotifier.Notify(data)
        allValues4 = seriesNotifier.GetVar(4)
        allValues3 = seriesNotifier.GetVar(3)
        print("Firing times: %s" % str(fireNotifier.firingTimes))
        #fig = plt.figure()
        plt.subplot(2, 1, 1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes, allValues3)
        plt.subplot(2, 1, 2)
        plt.plot(allTimes, allValues4)
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 
        #plt.show()

    def test_RespiratoryDelayedSmearedHeartActionForce(self):
        data = RungeKutta45IntegratorData(5,0.0)
        period = 0.5
        lastFire = - period #force fire at 0.0
        fireNotifier = FiringTimesNotifier()
        force = RespiratoryDelayedSmearedHeartActionForce()
        force.Drive = 0.0 #Current value of heart drive.
        force.CoordinateNumber = 4 #Coordinate number, to  which the state will be written
        force.KickAmplitude = 1.0 #Kick amplitude
        force.DelayTau = 0.13 #Time from firing order to actual kick
        force.SamplingTime = 0.1 # required to normalize delay time.
        force.DecayTau = 15.0 # Time by which the drive decays
        force.Notify = fireNotifier.Notify

        npoints = 1000
        period = 1.0
        seriesNotifier = SeriesNotifier(5,npoints)
        allTimes = np.linspace(0.0,10.0,npoints)
        for t in allTimes:
            data.t = t
            data.y[3] = 10.0 * np.cos(2 * np.pi * t / period)
            data.y[0] = data.y[1] = data.y[2] = 0.0
            if t >= lastFire + period:
                data.y[0] = 1.0
                force.FireOrderTime = t
                lastFire = t
            print(data)
            force.ApplyDrive(data) # will open or not.
            seriesNotifier.Notify(data)

        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(allTimes,seriesNotifier.GetVar(3))
        plt.subplot(2, 1, 2)
        plt.plot(allTimes,seriesNotifier.GetVar(0),"r")
        plt.plot(allTimes,seriesNotifier.GetVar(4))
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 
def main():

    logging.basicConfig(level=logging.DEBUG)
    unittest.main()

if __name__ == "__main__":
    main()
