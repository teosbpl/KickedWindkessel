# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import sys
import logging
import warnings
#from enum import Enum
import numpy as np
import unittest
#from scipy.integrate import ode
import matplotlib.pyplot as plt
import pdb
from KickedWindkesselModel import KickedWindkesselModel
from AbpmFiducialPointsCollector import AbpmFiducialPointsCollector
from Notifiers import SeriesNotifier, FiringTimesNotifier, NotifierChain
from IntegrateAndFire import IntegrateAndFire, phaseEfectivenessCurveSH
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorData
from HeartActionForce import RectangularHeartActionForce,RespiratoryDelayedSmearedHeartActionForce

class KickedWindkesselModelTest(unittest.TestCase):
    """
    This class tests all features of a complicated KickedWindkesselModel class.
    1. Unperturbed rate.
    2. Response to high amplitude of intrathoracic pressure.
    3. Response to change of respiratory  period
    4. Response to altered heart flow time.
    """
    def test_UnperturbedRate(self):
        settings = KickedWindkesselModel.KickedWindkesselModelSettings()  
        settings.heartActionForce = RectangularHeartActionForce()
        model = KickedWindkesselModel(settings)
        model.param.Npoints = 1000
        series = SeriesNotifier(model.param.dimension,model.param.Npoints+1)
        model.Notify = series.Notify
        model.IterateToNotifiers()
        fname = sys._getframe().f_code.co_name + ".png"
        series.PlotSeries((1,2,3,4),fname)
        print("Plot data in %s" % fname)

    def test_FullCoupling(self):
        """
        This test handles the majority of physiological control of heart 
        rate. Respiratory IAF periodically forces the heart via the 
        RespiratoryDelayedSmearedHeartActionForce. The force is applied to heart
        drive. Initially the heart is in phase with respiration, as both 
        periods are commensurate (3:1) and initial phase of both is 0.
        """
        dimension = 10
        npoints = 120  
        settings = KickedWindkesselModel.KickedWindkesselModelSettings() 
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 5 #Coordinate number, to  which the state will be written
        force.KickAmplitude = 2.0 #Kick amplitude
        force.DelayTau = 0.1 #Time from firing order to actual kick
        force.SamplingTime = 0.1 # required to normalize delay time.                
        force.DecayTau = 0.9 # Time by which the drive decays
        force.Notify = fireNotifier.Notify        

        iafResp = IntegrateAndFire()
        iafResp.SetPhaseVelocityFromBPM(20)
        iafResp.SamplingTime = 0.1
        iafResp.SetInitialPhase(0.0)
        iafResp.KickAmplitude = 0.1
        iafResp.SamplingTime = 0.1
        iafResp.CoordinateNumberForRate = 7         
        iafResp.CoordinateNumberForPhase = 6
        iafResp.CoordinateNumberForForceInput = -1 #not used
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 4         
        fireNotifierResp = FiringTimesNotifier()
        iafResp.Notify = fireNotifierResp.Notify
        

        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.1        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        iafHeart.CoordinateNumberForRate = 0
        iafHeart.CoordinateNumberForPhase = 1
        iafHeart.CoordinateNumberForForceInput = force.CoordinateNumber        
        # Coordinate number, to  which the state will be written        
        iafHeart.CoordinateNumberForOutput = 2 
        fireNotifierHeart = FiringTimesNotifier()
        iafHeart.Notify = fireNotifierHeart.Notify

        collector = AbpmFiducialPointsCollector(0)#ABPM = 0
        seriesNotifier = SeriesNotifier(dimension,npoints)
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
        
        settings.heartActionForce = HeartActionForceChain([iafResp,force,iafHeart])
        model = KickedWindkesselModel(settings,dimension)
        model.param.Npoints = npoints
        
        chain = NotifierChain((collector,seriesNotifier))
        model.Notify = chain.Notify
            
        model.IterateToNotifiers()
            
        for t in allTimes:
            data.t = t            
            data.y[0] = 0.0 # reset only the variable which gets set manually 
            # to prevent constant firing
            iafResp.ApplyDrive(data)
            if data[iafResp.CoordinateNumberForOutput] > 0.0:            
                data.y[0] = 1.0
                force.FireOrderTime = t
            logging.debug(data)
            force.ApplyDrive(data) # will open or not.
            iafHeart.ApplyDrive(data)
            seriesNotifier.Notify(data)
                    
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(5, 1, 1)
        plt.plot(allTimes,seriesNotifier.GetVar(6),"b")
        plt.plot(fireNotifierResp.firingTimes,fireNotifierResp.firingTimesSpikes(),"bo")        
        plt.ylabel("Resp. phase [1/rad]")
        plt.xlabel("Time [s]")
        plt.ylim(0.0,1.2)
        
        plt.subplot(5, 1, 2)
        plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForOutput),"g",linewidth=2)
        plt.plot(fireNotifier.firingTimes,fireNotifier.firingTimesSpikes(),"go")
        plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"g",linewidth=2)      

        plt.xlabel("Time [s]")
        plt.ylabel("Force")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        plt.subplot(5, 1, 3)
        plt.ylabel("Heart phase velocity")                
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForRate),"v")
        
        plt.subplot(5, 1, 4)
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForPhase),"r")
        #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
        plt.plot(fireNotifierHeart.firingTimes,fireNotifierHeart.firingTimesSpikes(),"ro")
        plt.ylim(0.0,1.2)        
        plt.xlabel("Time [s]")
        plt.ylabel("Heart Phase")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)

        plt.subplot(5, 1, 5)
        plt.ylabel("ISI")                
        plt.plot(fireNotifierHeart.firingTimes[1:],fireNotifierHeart.ISI()[1:],"b+-")

        
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname)         
    def test_HighAmplitude(self):
        pass
    def test_RespiratoryPeriod(self):
        pass
           
        plt.savefig("test_RespiratoryDelayedSmearedHeartActionForce.png")
def main():

    logging.basicConfig(level=logging.DEBUG)
    unittest.main()

if __name__ == "__main__":
    main()   


if __name__ == "__main__":
    BasicProcessor()

