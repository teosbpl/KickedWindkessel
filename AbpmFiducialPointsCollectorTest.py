# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:24:29 2018

@author: teos
"""
import logging
import unittest
from Notifiers import MinMaxNotifier, NotifierChain
from AbpmFiducialPointsCollector import AbpmFiducialPointsCollector
from KickedWindkesselModel import KickedWindkesselModel, HeartActionForce

class AbpmFiducialPointsCollectorTest(unittest.TestCase):
        def test_fiducial_points_collector(self): 
            settings = KickedWindkesselModel.KickedWindkesselModelSettings()    
            settings.heartActionForce = HeartActionForce()
            model = KickedWindkesselModel(settings)
            model.param.Npoints = 1000
            
            collector = AbpmFiducialPointsCollector()
            settings.heartActionForce.Notify = collector.HeartOpenNotifier
            mmNotifier = MinMaxNotifier(0)
            mmNotifier.MinNotifier = collector.AbpmMinNotifier
            mmNotifier.MaxNotifier = collector.AbpmMaxNotifier
            
            chain = NotifierChain((mmNotifier,collector))
            model.Notify = chain.Notify
            
            model.IterateToNotifiers()    
    
            
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
        logging.basicConfig(level=logging.INFO)    
        unittest.main()

if __name__ == "__main__":
    main()            
