# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import sys
import random
import logging
import unittest
import pdb
import numpy as np
import matplotlib.pyplot as plt
from Notifiers import SeriesNotifier, FiringTimesNotifier
from IntegrateAndFire import IntegrateAndFire, phaseEfectivenessCurveSH
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorData
from HeartActionForce import RespiratoryDelayedSmearedHeartActionForce



class IntegrateAndFireTest(unittest.TestCase):
    """
    This class tests all features of a complicated IntegrateAndFire class.
    1. Unperturbed rate.
    2. Perturbed rate which visually responds to stimulations (frequency modulation).
    3. Perturbed rate which is visually coupled to a different rhythm.
    4. Perturbed rate which may be plugged into Kicked Windkessel model.
    """

    def test_UnperturbedRate(self):
        data = RungeKutta45IntegratorData(8,0.0)
        iaf = IntegrateAndFire()
        fireNotifier = FiringTimesNotifier()
        iaf.Notify = fireNotifier.Notify
        npoints = 1000
        period = 1.0
        seriesNotifier = SeriesNotifier(8,npoints)
        allTimes = np.linspace(0.0,10.0,npoints)
        for t in allTimes:
            data.t = t
            # this is to make the plot visually nicer :-)
            data.y[3] = 10.0 * np.cos(2 * np.pi * t / period)
            #data.y[0] = data.y[1] = data.y[2] = 0.0
            iaf.ApplyDrive(data) # will open or not.
            #print(data)
            seriesNotifier.Notify(data)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(2,1,1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes,seriesNotifier.GetVar(3))
        plt.subplot(2,1,2)
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g")
        plt.plot(fireNotifier.firingTimes,spikes,"ro")
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

        #plt.show()
        
    def test_FrequencyModulation(self):
        data = RungeKutta45IntegratorData(8,0.0)
        iaf = IntegrateAndFire()
        fireNotifier = FiringTimesNotifier()
        iaf.Notify = fireNotifier.Notify
        npoints = 1000
        period = 20.0
        seriesNotifier = SeriesNotifier(8,npoints)
        allTimes = np.linspace(0.0,npoints*iaf.SamplingTime,npoints)
        for t in allTimes:
            data.t = t
            data.y[5] = 2.0 * iaf.r * (1.0 + np.cos(2 * np.pi * t / period))            
            iaf.ApplyDrive(data) # will open or not.
            #print(data)
            seriesNotifier.Notify(data)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(2,1,1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes,seriesNotifier.GetVar(5))#effective rate
        plt.subplot(2,1,2)
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g")
        plt.plot(fireNotifier.firingTimes,spikes,"ro")
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 


    def test_PhaseResponseCurve(self):
        data = RungeKutta45IntegratorData(8,0.0)
        iaf = IntegrateAndFire()
        iaf.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        fireNotifier = FiringTimesNotifier()
        iaf.Notify = fireNotifier.Notify
        npoints = 1000
        #allPhasesCircle = np.linspace(0.0,2.0*np.pi,npoints)
        allPhases = np.linspace(0.0,1.0,npoints)
        allPRC = np.linspace(0.0,1.0,npoints)
        for i,t in enumerate(allPhases):
            allPRC[i] = phaseEfectivenessCurveSH(t)
        fig = plt.figure()
        plt.subplot(2,1,1)
        plt.plot(allPhases,allPRC)
        plt.axhline(y=0.0,color='k',linestyle='-')
        ax = plt.subplot(2,1,2, projection='polar')
        ax.plot(allPhases*2.0*np.pi,0.1*np.ones(npoints),'b')
        ax.plot(allPhases*2.0*np.pi,0.1+allPRC,'r')
        ax.set_rmax(0.3)
        ax.set_rticks([0.1, 0.2, 0.3])  # less radial ticks
        ax.set_rlabel_position(22.5)  # get radial labels away from plotted line
        ax.grid(True)           
#        ax.set_aspect('equal', 'datalim')        
#        plt.plot(0.5*np.cos(2.0*np.pi*allPhases), 
#                 0.5*np.sin(2.0*np.pi*allPhases), 
#                    linewidth=1,color='k')
#        plt.plot((0.5 + phaseEfectivenessCurveSH(allPhases)) * np.cos(2.0*np.pi*allPhases), 
#                 (0.5 + phaseEfectivenessCurveSH(allPhases)) * np.sin(2.0*np.pi*allPhases), 
#                    linewidth=1,color='r')
#        plt.axhline(y=0.0,color='k',linestyle='-')
#        plt.axvline(x=0.0,color='k',linestyle='-')
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

 
 
    def test_FrequencyModulationWithPhaseResponseCurve(self):
        data = RungeKutta45IntegratorData(8,0.0)
        iaf = IntegrateAndFire()
        iaf.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        fireNotifier = FiringTimesNotifier()
        iaf.Notify = fireNotifier.Notify
        npoints = 1000
        period = 20.0
        seriesNotifier = SeriesNotifier(8,npoints)
        allTimes = np.linspace(0.0,npoints*iaf.SamplingTime,npoints)
        for t in allTimes:
            data.t = t
            data.y[5] = 2.0 * iaf.r * (1.0 + np.cos(2 * np.pi * t / period))
            data.y[6] = iaf.Phase+phaseEfectivenessCurveSH(iaf.Phase)
            iaf.ApplyDrive(data) # will open or not.
            #print("%lf\t%lf" % (iaf.Phase,data.y[6]))
            #print(data)
            seriesNotifier.Notify(data)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(2,1,1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes,seriesNotifier.GetVar(5))#effective rate
        plt.subplot(2,1,2)
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g")
        plt.plot(allTimes,seriesNotifier.GetVar(6),"r")
        plt.plot(fireNotifier.firingTimes,spikes,"ro")        
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 


    def test_FrequencyModulationWithPRCandKick(self):
        """
        This test is a comparison of two IAFs: one with PRC and one without.
        """
        npoints = 1000
        period = 20.0
        
        data = RungeKutta45IntegratorData(8,0.0)

        iaf = IntegrateAndFire()
        iaf.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        fireNotifier = FiringTimesNotifier()
        iaf.Notify = fireNotifier.Notify
        seriesNotifier = SeriesNotifier(8,npoints)
        
        data2 = RungeKutta45IntegratorData(8,0.0)
        iaf2 = IntegrateAndFire()        
        fireNotifier2 = FiringTimesNotifier()
        iaf2.Notify = fireNotifier2.Notify
        seriesNotifier2 = SeriesNotifier(8,npoints)

        allTimes = np.linspace(0.0,npoints*iaf.SamplingTime,npoints)
        binaryKick = np.zeros(npoints)
        for i in range(npoints):
            if random.randint(0,10) == 1:
                binaryKick[i] = 0.2

        for i,t in enumerate(allTimes):
            data.t = t
            data[5] = binaryKick[i]
            data2.t = t
            data2[5] = binaryKick[i]

            
            #data.y[5] = 2.0 * iaf.r * (1.0 + np.cos(2 * np.pi * t / period))
            
            #data.y[6] = iaf.Phase+phaseEfectivenessCurveSH(iaf.Phase)
            iaf.ApplyDrive(data) # will open or not.
            iaf2.ApplyDrive(data2)
            #print("%lf\t%lf" % (iaf.Phase,data.y[6]))
            #print(data)
            seriesNotifier.Notify(data)
            seriesNotifier2.Notify(data2)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        spikes2 = np.ones(len(fireNotifier2.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(3,1,1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes,seriesNotifier.GetVar(5))#effective rate
        plt.subplot(3,1,2)
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g")
        plt.plot(allTimes,seriesNotifier.GetVar(6),"r")
        plt.plot(fireNotifier.firingTimes,spikes,"ro")        
        plt.subplot(3,1,3)
        plt.plot(allTimes,seriesNotifier2.GetVar(4),"g")
        plt.plot(allTimes,seriesNotifier2.GetVar(6),"r")
        plt.plot(fireNotifier2.firingTimes,spikes2,"ro")        

        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

    def test_TwoIAFOneWayCoupled(self):
        """
        This test is a comparison of two IAFs: one with PRC and one without.
        """
        npoints = 1000
        
        data = RungeKutta45IntegratorData(10,0.0)

        iafResp = IntegrateAndFire()
        iafResp.r = 0.001
        iafResp.KickAmplitude = 0.01
        iafResp.CoordinateNumberForRate = 7         
        iafResp.CoordinateNumberForPhase = 6
        iafResp.CoordinateNumberForForceInput = 5 #not used
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 4         
        fireNotifier = FiringTimesNotifier()
        iafResp.Notify = fireNotifier.Notify
        seriesNotifier = SeriesNotifier(10,npoints)
        
        iafHeart = IntegrateAndFire()
        iafHeart.r = 0.0021
        iafHeart.CoordinateNumberForRate = 0
        iafHeart.CoordinateNumberForPhase = 1
        iafHeart.CoordinateNumberForForceInput = 4 #essential
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafHeart.CoordinateNumberForOutput = 2 
        fireNotifierHeart = FiringTimesNotifier()
        iafHeart.Notify = fireNotifierHeart.Notify

        allTimes = np.linspace(0.0,npoints*iafHeart.SamplingTime,npoints)
        for t in allTimes:
            data.t = t
 
            #should work out of the box.
            iafResp.ApplyDrive(data) # will open or not.
            iafHeart.ApplyDrive(data)
            #print("%lf\t%lf" % (iaf.Phase,data.y[6]))
            #print(data)
            seriesNotifier.Notify(data)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        spikes2 = np.ones(len(fireNotifierHeart.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(3,1,1)
        #plt.ylim(-12.0,12.0)
        #plt.plot(allTimes,seriesNotifier.GetVar(5))#effective rate
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g",linewidth=2)
        plt.subplot(3,1,2)
        
        plt.plot(allTimes,seriesNotifier.GetVar(6),"b")
        plt.plot(fireNotifier.firingTimes,spikes,"bo")        
        plt.subplot(3,1,3)
        plt.plot(allTimes,seriesNotifier.GetVar(1),"r")
        #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
        plt.plot(fireNotifierHeart.firingTimes,spikes2,"ro")        

        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

    def test_FrequencyModulationWithPhaseResponseCurve(self):
        data = RungeKutta45IntegratorData(8,0.0)
        iaf = IntegrateAndFire()
        iaf.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        fireNotifier = FiringTimesNotifier()
        iaf.Notify = fireNotifier.Notify
        npoints = 1000
        period = 20.0
        seriesNotifier = SeriesNotifier(8,npoints)
        allTimes = np.linspace(0.0,npoints*iaf.SamplingTime,npoints)
        for t in allTimes:
            data.t = t
            data.y[5] = 2.0 * iaf.r * (1.0 + np.cos(2 * np.pi * t / period))
            data.y[6] = iaf.Phase+phaseEfectivenessCurveSH(iaf.Phase)
            iaf.ApplyDrive(data) # will open or not.
            #print("%lf\t%lf" % (iaf.Phase,data.y[6]))
            #print(data)
            seriesNotifier.Notify(data)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(2,1,1)
        #plt.ylim(-12.0,12.0)
        plt.plot(allTimes,seriesNotifier.GetVar(5))#effective rate
        plt.subplot(2,1,2)
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g")
        plt.plot(allTimes,seriesNotifier.GetVar(6),"r")
        plt.plot(fireNotifier.firingTimes,spikes,"ro")        
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

    def test_TwoIAFOneWayCoupledInSeconds(self):
        """
        This test is a comparison of two IAFs: one with PRC and one without.
        """
        npoints = 120
        
        data = RungeKutta45IntegratorData(10,0.0)

        iafResp = IntegrateAndFire()
        iafResp.SetPhaseVelocityFromBPM(20)
        iafResp.KickAmplitude = 0.01
        iafResp.CoordinateNumberForRate = 7         
        iafResp.CoordinateNumberForPhase = 6
        iafResp.CoordinateNumberForForceInput = 5 #not used
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 4         
        fireNotifier = FiringTimesNotifier()
        iafResp.Notify = fireNotifier.Notify
        seriesNotifier = SeriesNotifier(10,npoints)
        
        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.1        
        iafHeart.SetPhaseVelocityFromBPM(66)
        iafHeart.CoordinateNumberForRate = 0
        iafHeart.CoordinateNumberForPhase = 1
        iafHeart.CoordinateNumberForForceInput = 4 #essential
        iafHeart.CoordinateNumberForForceInput = 9 # void  #essential
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafHeart.CoordinateNumberForOutput = 2 
        fireNotifierHeart = FiringTimesNotifier()
        iafHeart.Notify = fireNotifierHeart.Notify

        allTimes = np.linspace(0.0,npoints*iafHeart.SamplingTime,npoints)
        for t in allTimes:
            data.t = t
 
            #should work out of the box.
            iafResp.ApplyDrive(data) # will open or not.
            iafHeart.ApplyDrive(data)
            #print("%lf\t%lf" % (iaf.Phase,data.y[6]))
            #print(data)
            seriesNotifier.Notify(data)
        #Print dots at firing times
        spikes = np.ones(len(fireNotifier.firingTimes))
        spikes2 = np.ones(len(fireNotifierHeart.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(3,1,1)
        #plt.ylim(-12.0,12.0)
        #plt.plot(allTimes,seriesNotifier.GetVar(5))#effective rate
        plt.plot(allTimes,seriesNotifier.GetVar(4),"g",linewidth=2)
        plt.subplot(3,1,2)
        #pdb.set_trace()
        
        plt.plot(allTimes,seriesNotifier.GetVar(6),"b")
        plt.plot(fireNotifier.firingTimes,spikes,"bo")
        plt.xlabel("Time [s]")
        plt.ylim(0.0,1.2)
        plt.subplot(3,1,3)
        plt.plot(allTimes,seriesNotifier.GetVar(1),"r")
        #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
        plt.plot(fireNotifierHeart.firingTimes,spikes2,"ro")
        plt.ylim(0.0,1.2)
        plt.xlabel("Time [s]")
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 


    def test_RespiratoryDelayedSmearedHeartActionForceAndPRC(self):
        data = RungeKutta45IntegratorData(5,0.0)
        period = 0.5
        lastFire = - period #force fire at 0.0
        fireNotifier = FiringTimesNotifier()
        force = RespiratoryDelayedSmearedHeartActionForce()
        force.Drive = 0.0 #Current value of heart drive.
        force.CoordinateNumber = 4 #Coordinate number, to  which the state will be written
        force.KickAmplitude = 1.0 #Kick amplitude
        force.DelayTau = 0.13 #Time from firing order to actual kick
        force.SamplingTime = 0.1 # required to normalize delay time.
        force.StepPeriod = 1.0 #Heart period
        force.FireOrderTime = None #Time of last beat [s]
        force.DecayTau = 15.0 # Time by which the drive decays
        force.Notify = fireNotifier.Notify

        npoints = 1000
        period = 1.0
        seriesNotifier = SeriesNotifier(5,npoints)
        allTimes = np.linspace(0.0,10.0,npoints)
        for t in allTimes:
            data.t = t
            data.y[3] = 10.0 * np.cos(2 * np.pi * t / period)
            data.y[0] = data.y[1] = data.y[2] = 0.0
            if t >= lastFire + period:
                data.y[0] = 1.0
                force.FireOrderTime = t
                lastFire = t
            #print(data)
            force.ApplyDrive(data) # will open or not.
            seriesNotifier.Notify(data)

        print("Firing times: %s" % str(fireNotifier.firingTimes))
        #fig = plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(allTimes,seriesNotifier.GetVar(3))
        plt.subplot(2, 1, 2)
        plt.plot(allTimes,seriesNotifier.GetVar(0),"r")
        plt.plot(allTimes,seriesNotifier.GetVar(4))
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

def main():
 
    logging.basicConfig(level=logging.INFO)
    unittest.main()

if __name__ == "__main__":
    main()
