# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import logging
import unittest
import numpy as np
import matplotlib.pyplot as plt
from Notifiers import NilNotify, SeriesNotifier
from HeartActionForce import RectangularHeartActionForce
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorData

class FiringTimesNotifier:
    def __init__(self): 
        self.firingTimes = []    
    def Notify(self,data):    
        self.firingTimes.append(data.t)    
        print(self.firingTimes)
            
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
        fig = plt.figure()
        plt.subplot(2,1,1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes,allValues3) 
        plt.subplot(2,1,2)
        plt.plot(allTimes,allValues4) 
        plt.savefig("test_RectangularHeartActionForce.png")
        plt.show()
def main():
        
        logging.basicConfig(level=logging.DEBUG)    
        unittest.main()

if __name__ == "__main__":
    main()       