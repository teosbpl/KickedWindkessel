# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import numpy as np
from Notifiers import NilNotify

class RectangularHeartActionForce:
    """
    This class represents the periodic action of the heart.
    Initially the heart is in phase with respiration,
    as both periods are commensurate (3:1) and initial phase of both is 0.
    """
    def __init__(self):
        self.Drive = 0.0 #Current value of heart drive.
        self.CoordinateNumber = 4 #Coordinate number, to  which the state will be written
        self.KickAmplitude = 1.0 #Kick amplitude
        self.StepWidth = 0.01 #Kick time
        self.StepPeriod = 1.2 #Heart period
        self.LastStepStart = 0.0 #Time of last beat [s]
        self.StepShift = 0.0 # Step shift [s] additionally shifts the heart action.
        self.Notify = NilNotify
    #pre step length zmiana od 0 do 1.25 - przesuwa w fazie.
    def ApplyDrive(self,data):
        isOpen = self.Drive == self.KickAmplitude
        justOpened = False
        timeFromLastStepStart = data.t - self.LastStepStart
        if isOpen: # check if it is time to close
            if timeFromLastStepStart > self.StepWidth:
                self.Drive = 0.0
        else: #check if it is time to open
            if timeFromLastStepStart == 0.0 or timeFromLastStepStart > self.StepPeriod: # initial kick
                justOpened = True
                self.Drive = self.KickAmplitude
                self.LastStepStart = data.t
        data[self.CoordinateNumber] = self.Drive
        if justOpened:
            self.Notify(data)

class RandomHeartActionForce:
    """
    This class represents the periodic action of the heart.
    Initially the heart is in phase with respiration,
    as both periods are commensurate (3:1) and initial phase of both is 0.
    """
    def __init__(self):
        self.Drive = 0.0 #Current value of heart drive.
        self.CoordinateNumber = 4 #Coordinate number, to  which the state will be written
        self.KickAmplitude = 1.0 #Kick amplitude
        self.StepWidth = 0.01 #Kick time
        self.StepPeriod = np.random.poisson(5, 1) #Heart period
        self.LastStepStart = 0.0 #Time of last beat [s]
        self.StepShift = 0.0 # Step shift [s] additionally shifts the heart action.
        self.Notify = NilNotify
    #pre step length zmiana od 0 do 1.25 - przesuwa w fazie.
    def ApplyDrive(self,data):
        isOpen = self.Drive == self.KickAmplitude
        justOpened = False
        timeFromLastStepStart = data.t - self.LastStepStart
        if isOpen: # check if it is time to close
            if timeFromLastStepStart > self.StepWidth:
                self.Drive = 0.0 #Close
                self.StepPeriod = np.random.poisson(5, 1) #Select new.
        else: #check if it is time to open
            if timeFromLastStepStart == 0.0 or timeFromLastStepStart > self.StepPeriod: # initial kick
                justOpened = True
                self.Drive = self.KickAmplitude
                self.LastStepStart = data.t
        data[self.CoordinateNumber] = self.Drive
        if justOpened:
            self.Notify(data)

class RespiratoryDelayedSmearedHeartActionForce:
    """    
    This class is a helper to represent the action of the heart driven by external kick.
    The external kick is applied by setting FireOrderTime to current time.
    Then after delay time a kick is applied to drive value.
    """
    def __init__(self):
        self.Drive = 0.0 #Current value of heart drive.
        self.CoordinateNumber = 4 #Coordinate number, to  which the state will be written
        self.CoordinateNumberForInput = 9 #Coordinate number, from which the state will be read
        self.KickAmplitude = 1.0 #Kick amplitude
        self.DelayTau = 0.1 #Time from firing order to actual kick
        self.SamplingTime = 0.1 # required to normalize delay time.
        self.FireOrderTime = None #Time of last beat [s]
        self.DecayTau = 0.2 # Time by which the drive decays
        self.Notify = NilNotify
    #pre step length zmiana od 0 do 1.25 - przesuwa w fazie.
    def ApplyDrive(self,data):
        if data[self.CoordinateNumberForInput] > 0.0:            
                data[self.CoordinateNumberForInput] = 0.0#reset the kick signal
                self.FireOrderTime = t
        justOpened = False
        if self.FireOrderTime: #is not None
            #waiting for delay time
            timeFromFireOrderTime = data.t - self.FireOrderTime
            if timeFromFireOrderTime >= self.DelayTau:
                self.Drive = self.Drive + self.KickAmplitude
                self.FireOrderTime = None #nullify after use, wait for next kick
                justOpened = True
        data[self.CoordinateNumber] = self.Drive
        #decrement drive
        driveDecrement = self.Drive * (self.SamplingTime / self.DecayTau)
        self.Drive = self.Drive - driveDecrement
        if justOpened:
            self.Notify(data)

class HeartActionForceChain:
    """
    All forces are fired in order of appearance.
    """
    def __init__(self,forces = None):
        self.forces = []
        if forces:
            for force in forces:
                self.RegisterForce(force)
                
    def RegisterForce(self,force):
        self.forces.append(force)
        
    def ApplyDrive(self,data):
        for force in self.forces:
            force.ApplyDrive(data)