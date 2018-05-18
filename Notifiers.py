# -*- coding: utf-8 -*-
"""
Created on Mon May  7 00:22:50 2018

@author: teos
"""
import logging
import numpy as np
import matplotlib.pyplot as plt
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorData
def NilNotify(data):
    pass
    
class SeriesNotifier:
    def __init__(self,dimension,npoints):        
        self.dimension = dimension
        self.states = np.zeros((npoints,self.dimension+1))# all + time
        self.statesIterator = 0            

    def Notify(self,data):
        self.states[self.statesIterator,0] = data.t
        for ii in range(self.dimension):                    
            self.states[self.statesIterator,ii+1] = data.y[ii]
        self.statesIterator = self.statesIterator+1
    def GetVar(self,varNumber):
        """Var is counted as in data vector"""
        return self.states[:,varNumber+1]
    def PlotSeries(self,rangeOfVars,fileName = None):
        t = self.states[:,0]
        fig = plt.figure()
        for varNumber in rangeOfVars:         
            x = self.states[:,varNumber]        
            plt.plot(t,x)
        if fileName:
            plt.savefig(fileName)
        else:
            plt.show()

class StatsNotifier:
    def __init__(self,dimension):
        self.dimension = dimension
        self.sumValues = np.zeros(self.dimension)
        self.sumSquares = np.zeros(self.dimension)
        self.statesIterator = 0            

    def Notify(self,data):
        for ii in range(self.dimension):                    
            self.sumValues[ii] = self.sumValues[ii] + data.y[ii]
            self.sumSquares[ii] = self.sumSquares[ii] + (data.y[ii] * data.y[ii])
        self.statesIterator = self.statesIterator+1
    def GetStats(self):
        """
        Sensitive to catrastrophic cancellation.
        """
        AV = np.zeros(self.dimension)
        SD = np.zeros(self.dimension)
        for ii in range(self.dimension):                    
            AV[ii] = (1.0 * self.sumValues[ii] / self.statesIterator)
            SD[ii] = np.sqrt((1.0 * self.sumSquares[ii] / self.statesIterator) - AV[ii]*AV[ii])
        return AV,SD

class MinMaxNotifier:
    def __init__(self,varNumber,bufferSizeGuess=100):
        self.varNumber = varNumber
        self.currentValue = None
        self.currentTime = None
        self.trendDirection = None
        self.minTime = np.zeros(bufferSizeGuess)
        self.minValue = np.zeros(bufferSizeGuess)
        self.maxTime = np.zeros(bufferSizeGuess)
        self.maxValue = np.zeros(bufferSizeGuess)
        self.maxIterator = 0
        self.minIterator = 0
        self.MinNotifier = NilNotify
        self.MaxNotifier = NilNotify
        self.notifierData = RungeKutta45IntegratorData(1,0.0)
    def registerExtremum(self,times,values,iterator,notifier):
        self.notifierData.t = self.currentTime
        self.notifierData.y[0] = self.currentValue
        notifier(self.notifierData)
        if(iterator < len(times)):
            times[iterator] = self.currentTime
            values[iterator] = self.currentValue
        else:
            np.append(times,(self.currentTime,))
            np.append(values,(self.currentValue,))
        return iterator + 1
    def registerMinimum(self):
        self.minIterator = self.registerExtremum(self.minTime,self.minValue,self.minIterator,self.MinNotifier)
    def registerMaximum(self):
        self.maxIterator = self.registerExtremum(self.maxTime,self.maxValue,self.maxIterator,self.MaxNotifier)
    def Notify(self,data):
        def getTrendDirection(currentValue,newValue):
            if newValue == currentValue:
                return 0
            elif newValue > currentValue:
                return +1
            else:
                return -1
            
        newValue  = data.y[self.varNumber]
        newTime = data.t
        if self.currentValue is None:
            self.currentValue = newValue
            self.currentTime = newTime
            return
        if self.trendDirection is None:
            self.trendDirection = getTrendDirection(self.currentValue,newValue)
            return
        logging.debug("%lf -> %lf at trend: %lf" % (self.currentValue,newValue,self.trendDirection))
        #now the machine is ready to register.
        #it is assumed that the plot is smooth. If it isn't better method 
        #has to be used.
        #new value is larger at negative trend
        if newValue > self.currentValue and self.trendDirection < 0:
            logging.debug("Is minimum: %lf" % self.currentValue)
            self.registerMinimum()
        #new value is smaller at positive trend
        if newValue < self.currentValue and self.trendDirection > 0:
            logging.debug("Is maximum: %lf" % self.currentValue)
            self.registerMaximum()
        self.trendDirection = getTrendDirection(self.currentValue,newValue)                
        self.currentValue = newValue
        self.currentTime = newTime
    def getMaxima(self):
        if(self.maxIterator < len(self.maxTime)):
            self.maxTime = self.maxTime[0:self.maxIterator]
            self.maxValue = self.maxValue[0:self.maxIterator]
        return self.maxTime,self.maxValue,self.maxIterator
    def getMinima(self):
        if(self.minIterator < len(self.minTime)):
            self.minTime = self.minTime[0:self.minIterator]
            self.minValue = self.minValue[0:self.minIterator]
        return self.minTime,self.minValue,self.minIterator        
        def IterateToNotifiers(self):

            self.integrator.Reset(self.param)
            self.Notify(self.data)
class NotifierChain:
    def __init__(self,notifiers = None):
        self.notifiers = []
        if notifiers:
            for notifier in notifiers:
                self.RegisterNotifier(notifier)        
    def RegisterNotifier(self,notifier):
        self.notifiers.append(notifier)
    def Notify(self,data):
        for notifier in self.notifiers:
            notifier.Notify(data)

class FiringTimesNotifier:
    def __init__(self):
        self.firingTimes = []
    def Notify(self,data):
        self.firingTimes.append(data.t)
        logging.debug(self.firingTimes)