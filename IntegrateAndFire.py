# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import numpy as np
from Notifiers import NilNotify

def phaseEfectivenessCurve(phi):        
            aux1 = np.power(1 - phi, 3.0)
            value = np.power(phi, 1.3) * (phi - 0.45) * (aux1 / ((np.power(1 - 0.8, 3.0)) + aux1));
            return value    
            
class IntegrateAndFire:
    """
    This class represents the heart phase.
    
    
    """
    def __init__(self):
        self.Phase = 0.0 #Current value of heart drive.
        self.r = 0.1 # r in sampling time units - natural phase increment
        self.SamplingTime = 0.1 # required to normalize r.
        self.CoordinateNumberForR = 6        
        self.CoordinateNumberForForceInput = 5
        self.CoordinateNumberForOutput = 4 #Coordinate number, to  which the state will be written
        self.KickAmplitude = 1.0 #Kick amplitude        
        self.Notify = NilNotify
    #pre step length zmiana od 0 do 1.25 - przesuwa w fazie.
    
    def ApplyDrive(self,data):
        currentDrive = 0.0
        effectiveR = self.r +  data[self.CoordinateNumberForForceInput] * phaseEfectivenessCurve(self.Phase)
        self.Phase = self.Phase + effectiveR/self.SamplingTime #at each time step
        isOpen = self.Phase >= 1.0        
        if isOpen: # check if it is time to close                
                currentDrive = self.KickAmplitude                
        data[self.CoordinateNumberForOutput] = currentDrive
        data[self.CoordinateNumberForR] = effectiveR
        if isOpen:
            self.Notify(data)
