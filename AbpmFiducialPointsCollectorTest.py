# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:24:29 2018

@author: teos
"""
import logging
import unittest
import pdb
import matplotlib.pyplot as plt
import numpy as np
from Notifiers import MinMaxNotifier, NotifierChain, SeriesNotifier
from AbpmFiducialPointsCollector import AbpmFiducialPointsCollector
from KickedWindkesselModel import KickedWindkesselModel
from HeartActionForce import RectangularHeartActionForce, RandomHeartActionForce

class AbpmFiducialPointsCollectorTest(unittest.TestCase):
        def test_fiducial_points_collector_random(self): 
            settings = KickedWindkesselModel.KickedWindkesselModelSettings()    
            settings.heartActionForce = RandomHeartActionForce()
            model = KickedWindkesselModel(settings)
            model.param.Npoints = 10000
            
            collector = AbpmFiducialPointsCollector(0)#ABPM = 0. Other compartments have different numbers.
            settings.heartActionForce.Notify = collector.HeartOpenNotifier
            mmNotifier = MinMaxNotifier(0)
            mmNotifier.MinNotifier = collector.AbpmMinNotifier
            mmNotifier.MaxNotifier = collector.AbpmMaxNotifier
            seriesNotifier = SeriesNotifier(model.param.dimension,model.param.Npoints)
            
            chain = NotifierChain((mmNotifier,collector,seriesNotifier))
            model.Notify = chain.Notify
            
            model.IterateToNotifiers()
            allValues = seriesNotifier.GetVar(0)
            allTimes = np.linspace(0.0,model.param.dT*model.param.Npoints,model.param.Npoints)
            allItems = collector.GetFiducialPointsList()
            fig = plt.figure()
            #plt.ylim(-12.0,12.0)
            plt.plot(allTimes,allValues)        
            #pdb.set_trace()
            plt.plot(allItems[:,0],allItems[:,1],"ro") # systolic
            plt.plot(allItems[:,0],allItems[:,2],"go") # diastolic
            plt.plot(allItems[:,0],allItems[:,3],"bo") # mean
            plt.savefig("test_fiducial_points_collector_random.png")
            logging.warning("Result in test_fiducial_points_collector_random.png")
            #plt.show()
            
        def test_fiducial_points_collector_rect(self): 
            settings = KickedWindkesselModel.KickedWindkesselModelSettings()    
            settings.heartActionForce = RectangularHeartActionForce()
            model = KickedWindkesselModel(settings)
            model.param.Npoints = 10000
            
            collector = AbpmFiducialPointsCollector(0)#ABPM = 0. Other compartments have different numbers.
            settings.heartActionForce.Notify = collector.HeartOpenNotifier
            mmNotifier = MinMaxNotifier(0)
            mmNotifier.MinNotifier = collector.AbpmMinNotifier
            mmNotifier.MaxNotifier = collector.AbpmMaxNotifier
            seriesNotifier = SeriesNotifier(model.param.dimension,model.param.Npoints)
            
            chain = NotifierChain((mmNotifier,collector,seriesNotifier))
            model.Notify = chain.Notify
            
            model.IterateToNotifiers()
            allValues = seriesNotifier.GetVar(0)
            allTimes = np.linspace(0.0,model.param.dT*model.param.Npoints,model.param.Npoints)
            allItems = collector.GetFiducialPointsList()
            fig = plt.figure()
            #plt.ylim(-12.0,12.0)
            plt.plot(allTimes,allValues)        
            #pdb.set_trace()
            plt.plot(allItems[:,0],allItems[:,1],"ro") # systolic
            plt.plot(allItems[:,0],allItems[:,2],"go") # diastolic
            plt.plot(allItems[:,0],allItems[:,3],"bo") # mean
            plt.savefig("test_fiducial_points_collector_rect.png")
            logging.warning("Result in test_fiducial_points_collector_rect.png")
            #plt.show()

#==============================================================================
#         data = RungeKutta45IntegratorData(2,0.0)
#         mmNotifier = MinMaxNotifier(0)#first variable        
#         npoints = 1000
#         seriesNotifier = SeriesNotifier(2,npoints)
#         notifier = NotifierChain((mmNotifier,seriesNotifier))
#         period = 1.0
#         allTimes = np.linspace(0.0,10.0,npoints)
#         for t in allTimes:
#             data.t = t
#             data.y[0] = 10.0 * np.cos(2 * np.pi * t / period)
#             data.y[1] = 10.0 * np.sin(2 * np.pi * t / period)
#             
#             notifier.Notify(data)
#         maxTime,maxValue,maxIterator = mmNotifier.getMaxima()
#         minTime,minValue,minIterator = mmNotifier.getMinima()
#         allValues = seriesNotifier.GetVar(0)
#         fig = plt.figure()
#         plt.ylim(-12.0,12.0)
#         plt.plot(allTimes,allValues)        
#         plt.plot(minTime,minValue,"bo")
#         plt.plot(maxTime,maxValue,"r*")
#         plt.show()
#==============================================================================

def main():
        logging.basicConfig(level=logging.WARNING)    
        unittest.main()

if __name__ == "__main__":
    main()            
