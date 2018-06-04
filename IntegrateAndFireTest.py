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
            logging.debug(data)
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
            logging.debug(data)
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
            logging.debug(data)
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
            logging.debug(data)
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
        npoints = 3100
        
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
            logging.debug(data)
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
            logging.debug(data)
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
        npoints = 240
        
        data = RungeKutta45IntegratorData(10,0.0)

        iafResp = IntegrateAndFire()
        iafResp.SetPhaseVelocityFromBPM(20)
        iafResp.SetInitialPhase(1.0/7.0)
        iafResp.KickAmplitude = 0.1
        iafResp.SamplingTime = 0.01
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
        iafHeart.SamplingTime = 0.01        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.CoordinateNumberForRate = 0
        iafHeart.CoordinateNumberForPhase = 1
        iafHeart.CoordinateNumberForForceInput = 4 #essential
        #iafHeart.CoordinateNumberForForceInput = 9 # void  #essential
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafHeart.CoordinateNumberForOutput = 2 
        fireNotifierHeart = FiringTimesNotifier()
        iafHeart.Notify = fireNotifierHeart.Notify

        allTimes = np.linspace(0.0,npoints*iafHeart.SamplingTime,npoints)
        for t in allTimes:
            data.t = t
 
            iafResp.ApplyDrive(data) # will open or not.
            iafHeart.ApplyDrive(data)
            #print("%lf\t%lf" % (iaf.Phase,data.y[6]))
            logging.debug(data)
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
        plt.xlabel("Time [s]")
        plt.ylabel("Force")
        plt.ylim(0.0,0.1)
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


    def test_RespiratoryDelayedSmearedHeartActionForceForcedByResp(self):
        """
        This test handles the majority of physiological control of heart 
        rate. Respiratory IAF periodically forces the heart via the 
        RespiratoryDelayedSmearedHeartActionForce.
        The test shows how respiratory kick is converted to delayed smeared kick.
        """
        npoints = 240
        data = RungeKutta45IntegratorData(10,0.0)
        

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
        
        
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()
        force.Drive = 0.0 #Current value of heart drive.
        force.CoordinateNumber = 5 #Coordinate number, to  which the state will be written
        force.KickAmplitude = 1.0 #Kick amplitude
        force.DelayTau = 1.53 #Time from firing order to actual kick
        force.SamplingTime = 0.1 # required to normalize delay time.
        force.StepPeriod = 1.0 #Heart period
        #force.FireOrderTime = None #Time of last beat [s]
        force.DecayTau = 15.0 # Time by which the drive decays
        force.Notify = fireNotifier.Notify


        seriesNotifier = SeriesNotifier(data.dimension,npoints)
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
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
            seriesNotifier.Notify(data)
            
        spikes = np.ones(len(fireNotifier.firingTimes))
        spikesResp = np.ones(len(fireNotifierResp.firingTimes))
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(3, 1, 1)
        plt.plot(allTimes,seriesNotifier.GetVar(6),"b")
        plt.plot(fireNotifierResp.firingTimes,spikesResp,"bo")        
        plt.ylabel("Resp. phase [1/rad]")
        plt.xlabel("Time [s]")
        plt.ylim(0.0,1.2)
        
        plt.subplot(3, 1, 2)
        plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForOutput),"g",linewidth=2)
        plt.xlabel("Time [s]")
        plt.ylabel("Force")
        plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        
        plt.subplot(3, 1, 3)
        plt.plot(fireNotifier.firingTimes,spikes,"bo")
        plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"g",linewidth=2)      
        plt.xlabel("Time [s]")
        plt.ylabel("Delayed Force")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

    
    def test_RspDrivesHeartViaSmearedForce(self):
        """
        This test handles the majority of physiological control of heart 
        rate. Respiratory IAF periodically forces the heart via the 
        RespiratoryDelayedSmearedHeartActionForce. The force is applied to heart
        drive. Initially the heart is in phase with respiration, as both 
        periods are commensurate (3:1) and initial phase of both is 0.
        """
        npoints = 160
        data = RungeKutta45IntegratorData(10,0.0)
        

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
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 5 #Coordinate number, to  which the state will be written
        force.KickAmplitude = 0.1 #Kick amplitude
        force.DelayTau = 0.53 #Time from firing order to actual kick
        force.SamplingTime = 0.1 # required to normalize delay time.                
        force.DecayTau = 5.0 # Time by which the drive decays
        force.Notify = fireNotifier.Notify

        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.1        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.CoordinateNumberForRate = 0
        iafHeart.CoordinateNumberForPhase = 1
        iafHeart.CoordinateNumberForForceInput = force.CoordinateNumber        
        # Coordinate number, to  which the state will be written        
        iafHeart.CoordinateNumberForOutput = 2 
        fireNotifierHeart = FiringTimesNotifier()
        iafHeart.Notify = fireNotifierHeart.Notify


        seriesNotifier = SeriesNotifier(data.dimension,npoints)
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
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
        plt.plot(fireNotifierHeart.firingTimes,fireNotifierHeart.ISI(),"v")

        
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 

    def test_RspDrivesHeartViaSmearedForceWithPRC(self):
        """
        This test handles the majority of physiological control of heart 
        rate. Respiratory IAF periodically forces the heart via the 
        RespiratoryDelayedSmearedHeartActionForce. The force is applied to heart
        drive. Initially the heart is in phase with respiration, as both 
        periods are commensurate (3:1) and initial phase of both is 0.
        """
        npoints = 120
        data = RungeKutta45IntegratorData(10,0.0)
        

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
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 5 #Coordinate number, to  which the state will be written
        force.KickAmplitude = 2.0 #Kick amplitude
        force.DelayTau = 0.1 #Time from firing order to actual kick
        force.SamplingTime = 0.1 # required to normalize delay time.                
        force.DecayTau = 0.9 # Time by which the drive decays
        force.Notify = fireNotifier.Notify

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


        seriesNotifier = SeriesNotifier(data.dimension,npoints)
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
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

    def test_DoubleKickViaSmearedForceWithPRC(self):
        """
        This test handles the majority of physiological control of heart 
        rate. Respiratory IAF periodically forces the heart via the 
        RespiratoryDelayedSmearedHeartActionForce. The force is applied to heart
        drive. Initially the heart is in phase with respiration, as both 
        periods are commensurate (3:1) and initial phase of both is 0.
        """
        npoints = 120
        data = RungeKutta45IntegratorData(10,0.0)
        

        iafResp = IntegrateAndFire()
        iafResp.SetPhaseVelocityFromBPM(20)
        iafResp.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        iafResp.SamplingTime = 0.1
        iafResp.SetInitialPhase(0.0)
        iafResp.KickAmplitude = 0.1
        iafResp.SamplingTime = 0.1
        iafResp.CoordinateNumberForRate = 7         
        iafResp.CoordinateNumberForPhase = 6
        iafResp.CoordinateNumberForForceInput = 8 #coupling in opposite direction.
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 4         
        fireNotifierResp = FiringTimesNotifier()
        iafResp.Notify = fireNotifierResp.Notify
        
        forceCoordinateNumber = 5

        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.1        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        iafHeart.CoordinateNumberForRate = 0
        iafHeart.CoordinateNumberForPhase = 1
        iafHeart.CoordinateNumberForForceInput = forceCoordinateNumber        
        # Coordinate number, to  which the state will be written        
        iafHeart.CoordinateNumberForOutput = 2 
        fireNotifierHeart = FiringTimesNotifier()
        iafHeart.Notify = fireNotifierHeart.Notify

        fireNotifier = FiringTimesNotifier()
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = forceCoordinateNumber #Coordinate number, to  which the state will be written
        force.CoordinateNumberForInput = iafResp.CoordinateNumberForOutput
        force.KickAmplitude = -1.0 #Kick amplitude
        force.DelayTau = 0.1 #Time from firing order to actual kick
        force.SamplingTime = 0.1 # required to normalize delay time.                
        force.DecayTau = 0.9 # Time by which the drive decays
        force.Notify = fireNotifier.Notify

        fireNotifierR = FiringTimesNotifier()
        forceR = RespiratoryDelayedSmearedHeartActionForce()        
        forceR.CoordinateNumber = 8 #Coordinate number, to  which the state will be written
        forceR.KickAmplitude = 1.0 #Kick amplitude
        forceR.DelayTau = 0.1 #Time from firing order to actual kick
        forceR.SamplingTime = 0.1 # required to normalize delay time.                
        forceR.DecayTau = 0.9 # Time by which the drive decays
        forceR.CoordinateNumberForInput = iafHeart.CoordinateNumberForOutput
        forceR.Notify = fireNotifierR.Notify




        seriesNotifier = SeriesNotifier(data.dimension,npoints)
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
        for t in allTimes:
            data.t = t            
            data.y[0] = 0.0 # reset only the variable which gets set manually 
            # to prevent constant firing
            iafResp.ApplyDrive(data)
            #if data[iafResp.CoordinateNumberForOutput] > 0.0:            
            #    data.y[0] = 1.0
            #    force.FireOrderTime = t
            logging.debug(data)
            forceR.ApplyDrive(data)
            force.ApplyDrive(data)
            iafHeart.ApplyDrive(data)
            seriesNotifier.Notify(data)
                    
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        plt.subplot(5, 1, 1)
        plt.plot(allTimes,seriesNotifier.GetVar(6),"b")
        plt.plot(fireNotifierResp.firingTimes,fireNotifierResp.firingTimesSpikes(),"bo")        
        plt.ylabel(r"$\Phi(t)$")
	plt.yticks([])
	plt.xticks([])
        plt.ylim(0.0,1.2)
        
        plt.subplot(5, 1, 2)
        plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForOutput),"g",linewidth=2)
        plt.plot(fireNotifier.firingTimes,fireNotifier.firingTimesSpikes(),"go")
        plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"g",linewidth=2)      
	plt.yticks([])
	plt.xticks([])
        plt.ylabel("Force")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        plt.subplot(5, 1, 3)
        plt.ylabel(r"$d\varphi(t)/dt$")                
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForRate),"v")
	plt.yticks([])    
	plt.xticks([])    
        plt.subplot(5, 1, 4)
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForPhase),"r")
        #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
        plt.plot(fireNotifierHeart.firingTimes,fireNotifierHeart.firingTimesSpikes(),"ro")
        plt.ylim(0.0,1.2)        
	plt.yticks([])
	plt.xticks([])
        plt.ylabel(r"$\varphi(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)

        plt.subplot(5, 1, 5)
        plt.ylabel("RR")      
	plt.yticks([])         
        plt.plot(fireNotifierHeart.firingTimes[1:],fireNotifierHeart.ISI()[1:],"b+-")
        plt.xlabel("Time [s]")
        
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname) 


def main():
 
    logging.basicConfig(level=logging.INFO)
    unittest.main()

if __name__ == "__main__":
    main()
