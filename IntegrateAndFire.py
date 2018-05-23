# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import numpy as np
from Notifiers import NilNotify

def phaseEfectivenessCurveSH(phi):        
            aux1 = np.power(1 - phi, 3.0)
            value = np.power(phi, 1.3) * (phi - 0.45) * (aux1 / ((np.power(1 - 0.8, 3.0)) + aux1));
            return value    

def phaseEfectivenessCurveNull(phi):        
            return 1.0
            
class IntegrateAndFire:
    """
    This class implements an integrate and fire model.
    This model accepts input from outside, adds it to phase velocity r
    and accumulates changes until phase reaches 1.
    Adding to phase velocity instead of phase itself is a short term: 
    it changes only the phase increment at certain time step.
    """
    def __init__(self):        
        self.Phase = 0.0 # Current value of heart drive.
        self.r = 0.001 # r in sampling time units - natural phase increment
        self.SamplingTime = 0.1 # required to normalize r.
        # the rate may be complicated. It has to be separately monitored.
        self.CoordinateNumberForRate = 7         
        self.CoordinateNumberForPhase = 6
        self.CoordinateNumberForForceInput = 5
        # Coordinate number, to  which the state will be written
        self.CoordinateNumberForOutput = 4 
        self.KickAmplitude = 1.0 #Kick amplitude        
        self.Notify = NilNotify
        self.NotifyFunctions = [] # list of notify functions instead of single Notify
        
        self.CoordinateNumberForOutput = 4
        self.phaseEfectivenessCurve = phaseEfectivenessCurveNull
        self.ResetPhaseAfterFire = True
        
    
    def SetPhaseVelocityFromBPM(self,bpm):
        """
        Full period = phase [0,1]
        Full period is in seconds
        BPM is in beats per minute = nbeats per 60 seconds.
        Full period is T = 60 / bpm [seconds].
        So in order to generate such a period
        Phase velocity must be 1.0 / T
        Single phase increment is per single sampling time step.
        r may be determined from proportion:
        1.0 - T
        phaseIncrement - sampling time increment
        
        phaseIncrement * T = 1.0 * SamplingTime
        phaseIncrement = SamplingTime / T                
        """
        self.r = self.SamplingTime * bpm / 60.0
        return self.r
        
    def SetInitialPhase(self,phase0):
        """
        Sets initial phase             
        """
        self.Phase = phase0
        return self.Phase
        
    def ApplyDrive(self,data):
        """
        The drive value is essentially zero. When the oscillator fires, drive is 1.0.
        Additionally the phase and the local phase slope is also reported.
        """
        currentDrive = 0.0
        if self.CoordinateNumberForForceInput >= 0:
#            effectiveR = self.r / (1.0 +  data[self.CoordinateNumberForForceInput] * self.phaseEfectivenessCurve(self.Phase))
            effectiveR = self.r  +  data[self.CoordinateNumberForForceInput] * self.phaseEfectivenessCurve(self.Phase)
        else:
            effectiveR = self.r
        self.Phase = self.Phase + effectiveR #at each time step        
        #self.Phase = self.Phase + effectiveR/self.SamplingTime #at each time step        
        doFire = self.Phase >= 1.0        
        if doFire:
                currentDrive = self.KickAmplitude
                if self.ResetPhaseAfterFire:
                    self.Phase = 0.0
                else:
                    self.Phase = self.Phase - 1.0
        data[self.CoordinateNumberForOutput] = currentDrive
        if self.CoordinateNumberForPhase != -1:
            data[self.CoordinateNumberForPhase] = self.Phase
        if self.CoordinateNumberForRate != -1:
            data[self.CoordinateNumberForRate] = effectiveR
        if doFire:
            self.Notify(data)
            if self.NotifyFunctions:
                for fn in self.NotifyFunctions:
                    fn(data)

