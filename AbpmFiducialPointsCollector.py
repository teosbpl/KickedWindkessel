# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:24:29 2018

@author: teos
"""
import logging
import copy
import numpy as np
        
class AbpmFiducialPointsCollector:
    """
    This class collects all events related to blood pressure changes in order to calculate proper 
    """
    def __init__(self,varNumber):
        self.varNumber = varNumber
        self.FiducialPointsList = None
        self.meanAbpmCollector = 0.0
        self.meanAbpmSampleCounter = 0
        self.lastMin = None
        self.lastMax = None
        self.lastMinTime = None
        self.lastMaxTime = None
        self.lastHeartOpenTime = 0.0
        self.last_t = 0.0
        pass
    def AbpmMinNotifier(self,data):
        self.lastMin = data.y[self.varNumber]
        self.lastMinTime = data.t
        logging.debug("In min notifier at %lf with %lf" % (data.t,self.lastMin))        
        
    def AbpmMaxNotifier(self,data):
        self.lastMax = data.y[self.varNumber]
        self.lastMaxTime = data.t
        logging.debug("In max notifier at %lf with %lf" % (data.t,self.lastMax))        

    def Finalize(self,data=None):
        if not self.lastMin:
            self.lastMin = 0.0
        if not self.lastMax:
            self.lastMax = 0.0
        _map = self.meanAbpmCollector / self.meanAbpmSampleCounter # MAP
        if data:
            t = data.t
        else:
            t = self.last_t
        _rr = t - self.lastHeartOpenTime
        _sbp = self.lastMax
        if self.lastMinTime < self.lastMaxTime:
            logging.debug("Registering DBP which was before SBP and RR to next: all in one item.")
        # I assume that HeartOpen is launched before MinNotifier, 
        #if it is launched before it carries the DBP which was before SBP            
        #otherwise MinNotifier carries DBP which is after SBP.
        #The latter behavior can be obtained easily: BP value at 
        #HeartOpen event is the DBP following the registered SBP.
        #MAP value is assigned to previous R, not to the following
        #RR - likewise
        _dbp = self.lastMin
        item = [self.lastHeartOpenTime, _sbp, _dbp, _map, _rr]
        if self.FiducialPointsList is None:
           self.FiducialPointsList = [item]
           logging.debug("Entering item")
        else:    
            self.FiducialPointsList.append(item)
            logging.debug("Appending item %s to %s" % (item,self.FiducialPointsList))
        self.meanAbpmSampleCounter = 0
        self.meanAbpmCollector = 0.0
        self.lastMin = None
        self.lastMax = None
        self.lastHeartOpenTime = t
        
    def HeartOpenNotifier(self,data):
        logging.debug("In heart notifier at %lf" % (data.t))        
        #close everything that was collected so far.
        if self.meanAbpmSampleCounter > 0:
            self.Finalize(data)

        
    
    def Notify(self,data):
        """
        Full model notifier is fireed at each step of simulation.
        Its name is Notifier so that it can be used in chain notifier.
        """        
        self.meanAbpmCollector = self.meanAbpmCollector + data.y[self.varNumber]
        self.meanAbpmSampleCounter = self.meanAbpmSampleCounter + 1
        self.last_t = data.t
        logging.debug("In full model notifier at %lf" % (data.t))
    def GetFiducialPointsList(self):
        if self.meanAbpmSampleCounter > 0:
            self.Finalize()        
        return np.array(self.FiducialPointsList)
        
        

