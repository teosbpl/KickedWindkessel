# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:24:29 2018

@author: teos
"""
import logging

class AbpmFiducialPointsCollector:
    """
    This class collects all events related to blood pressure changes in order to calculate proper 
    """
    def __init__(self):
        pass
    def AbpmMinNotifier(self,data):
        logging.info("In min notifier at %lf" % (data.t))        
        pass
    def AbpmMaxNotifier(self,data):
        logging.info("In max notifier at %lf" % (data.t))        
        pass
    def HeartOpenNotifier(self,data):
        logging.info("In heart notifier at %lf" % (data.t))        
        pass
    
    def Notify(self,data):
        """
        Full model notifier is fireed at each step of simulation.
        Its name is Notifier so that it can be used in chain notifier.
        """
        logging.debug("In full model notifier at %lf" % (data.t))        
        pass
        

