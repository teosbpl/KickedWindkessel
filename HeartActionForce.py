# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
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
        self.Amplitude = 1.0 #Kick amplitude
        self.StepWidth = 0.01 #Kick time
        self.StepPeriod = 1.0 #Heart period
        self.LastStepStart = 0.0 #Time of last beat [s]
        self.StepShift = 0.0 # Step shift [s] additionally shifts the heart action.
        self.Notify = NilNotify
    #pre step length zmiana od 0 do 1.25 - przesuwa w fazie.
    def ApplyDrive(self,data):
        isOpen = self.Drive == self.Amplitude
        justOpened = False
        timeFromLastStepStart = data.t - self.LastStepStart
        if isOpen: # check if it is time to close
            if timeFromLastStepStart > self.StepWidth:
                self.Drive = 0.0
        else: #check if it is time to open
            if timeFromLastStepStart == 0.0 or timeFromLastStepStart > self.StepPeriod: # initial kick
                justOpened = True
                self.Drive = self.Amplitude
                self.LastStepStart = data.t       
        data[self.CoordinateNumber] = self.Drive
        if justOpened:
            self.Notify(data)
            
class RespiratoryDelayedSmearedHeartActionForce:
    """
    This class represents the periodic action of the heart.
    Initially the heart is in phase with respiration, 
    as both periods are commensurate (3:1) and initial phase of both is 0.
    """
    def __init__(self):
        self.Drive = 0.0 #Current value of heart drive.
        self.CoordinateNumber = 4 #Coordinate number, to  which the state will be written
        self.Amplitude = 1.0 #Kick amplitude
        self.StepWidth = 0.01 #Kick time
        self.StepPeriod = 1.0 #Heart period
        self.LastStepStart = 0.0 #Time of last beat [s]
        self.StepShift = 0.0 # Step shift [s] additionally shifts the heart action.
        self.Notify = NilNotify
    #pre step length zmiana od 0 do 1.25 - przesuwa w fazie.
    def ApplyDrive(self,data):
        isOpen = self.Drive == self.Amplitude
        justOpened = False
        timeFromLastStepStart = data.t - self.LastStepStart
        if isOpen: # check if it is time to close
            if timeFromLastStepStart > self.StepWidth:
                self.Drive = 0.0
        else: #check if it is time to open
            if timeFromLastStepStart == 0.0 or timeFromLastStepStart > self.StepPeriod: # initial kick
                justOpened = True
                self.Drive = self.Amplitude
                self.LastStepStart = data.t       
        data[self.CoordinateNumber] = self.Drive
        if justOpened:
            self.Notify(data)            