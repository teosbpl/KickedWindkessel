# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import logging
import warnings
#from enum import Enum
import numpy as np
import unittest
#from scipy.integrate import ode
import matplotlib.pyplot as plt
import pdb
from KickedWindkesselModel import KickedWindkesselModel
#from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorParams, RungeKutta45ConstStepIntegrator, RungeKutta45IntegratorData
from Notifiers import SeriesNotifier
from HeartActionForce import RectangularHeartActionForce

class KickedWindkesselModelTest(unittest.TestCase):
    """
    This class tests all features of a complicated KickedWindkesselModel class.
    1. Unperturbed rate.
    2. Response to high amplitude of intrathoracic pressure.
    3. Response to change of respiratory  period
    4. Response to altered heart flow time.
    """
    def test_UnperturbedRate(self):
        settings = KickedWindkesselModel.KickedWindkesselModelSettings()    
        settings.heartActionForce = RectangularHeartActionForce()
        model = KickedWindkesselModel(settings)
        model.param.Npoints = 1000
        series = SeriesNotifier(model.param.dimension,model.param.Npoints+1)
        model.Notify = series.Notify
        model.IterateToNotifiers()
        fname = "overview.png"
        series.PlotSeries((1,2,3,4),fname)
        print("Plot data in %s" % fname)

    def test_HighAmplitude(self):
        pass
    def test_RespiratoryPeriod(self):
        pass
           
        plt.savefig("test_RespiratoryDelayedSmearedHeartActionForce.png")
def main():

    logging.basicConfig(level=logging.DEBUG)
    unittest.main()

if __name__ == "__main__":
    main()   


if __name__ == "__main__":
    BasicProcessor()

