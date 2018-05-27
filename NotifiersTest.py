# -*- coding: utf-8 -*-
"""
Created on Mon May  7 00:22:50 2018

@author: teos
"""
import logging
import unittest
import numpy as np
import matplotlib.pyplot as plt
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorParams, RungeKutta45ConstStepIntegrator, RungeKutta45IntegratorData
from Notifiers import MinMaxNotifier, SeriesNotifier, NotifierChain, StatsNotifier
from KickedWindkesselModel import KickedWindkesselModel
from HeartActionForce import RectangularHeartActionForce


class NotifiersTest(unittest.TestCase):

    def test_min_max_notifier(self):
        data = RungeKutta45IntegratorData(2, 0.0)
        mmNotifier = MinMaxNotifier(0)#first variable
        npoints = 1000
        seriesNotifier = SeriesNotifier(2, npoints)
        notifier = NotifierChain((mmNotifier, seriesNotifier))
        period = 1.0
        allTimes = np.linspace(0.0, 10.0, npoints)
        for t in allTimes:
            data.t = t
            data.y[0] = 10.0 * np.cos(2 * np.pi * t / period)
            data.y[1] = 10.0 * np.sin(2 * np.pi * t / period)

            notifier.Notify(data)
        maxTime,maxValue,maxIterator = mmNotifier.getMaxima()
        minTime,minValue,minIterator = mmNotifier.getMinima()
        allValues = seriesNotifier.GetVar(0)
        fig = plt.figure()
        plt.ylim(-12.0, 12.0)
        plt.plot(allTimes,allValues)
        plt.plot(minTime, minValue, "bo")
        plt.plot(maxTime, maxValue, "r*")
        plt.show()
    def test_stats_notifier(self):
        settings = KickedWindkesselModel.KickedWindkesselModelSettings()
        settings.heartActionForce = RectangularHeartActionForce()
        model = KickedWindkesselModel(settings)
        model.param.Npoints = 1000
        stats = StatsNotifier(model.param.dimension)
        model.Notify = stats.Notify
        model.IterateToNotifiers()
        AV,SD = stats.GetStats()
        print(AV)
        print(SD)

def main():

        logging.basicConfig(level=logging.DEBUG)
        unittest.main()

if __name__ == "__main__":
    main()