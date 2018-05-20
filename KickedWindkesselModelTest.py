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
import matplotlib.gridspec as gridspec
#from scipy.integrate import ode
import matplotlib.pyplot as plt
import pdb
from KickedWindkesselModel import KickedWindkesselModel
from AbpmFiducialPointsCollector import AbpmFiducialPointsCollector
from Notifiers import SeriesNotifier, FiringTimesNotifier, NotifierChain, MinMaxNotifier
from IntegrateAndFire import IntegrateAndFire, phaseEfectivenessCurveSH
from HeartActionForce import RectangularHeartActionForce,RespiratoryDelayedSmearedHeartActionForce,HeartActionForceChain
from KickedWindkesselModelVisualization import KickedWindkesselModelVisualization

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
        npoints = 1800
        settings = KickedWindkesselModel.KickedWindkesselModelSettings() 
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 9 #Coordinate number, to  which the state will be written
        force.CoordinateNumberForInput = 6
        force.KickAmplitude = 1.0 #Kick amplitude
        force.DelayTau = 0.1 #Time from firing order to actual kick
        force.SamplingTime = 0.01 # required to normalize delay time.                
        force.DecayTau = 0.9 # Time by which the drive decays
        force.Notify = fireNotifier.Notify        

        iafResp = IntegrateAndFire()
        iafResp.SamplingTime = 0.01
        iafResp.SetPhaseVelocityFromBPM(20)
        iafResp.SetInitialPhase(0.0)
        iafResp.KickAmplitude = 0.1        
        iafResp.CoordinateNumberForRate = -1
        iafResp.CoordinateNumberForPhase = 5
        iafResp.CoordinateNumberForForceInput = -1 #not used
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 6         
        fireNotifierResp = FiringTimesNotifier()
        iafResp.Notify = fireNotifierResp.Notify
        

        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.01        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        iafHeart.CoordinateNumberForRate = 8
        iafHeart.CoordinateNumberForPhase = 7
        iafHeart.CoordinateNumberForForceInput = force.CoordinateNumber        
        # Coordinate number, to  which the state will be written        
        iafHeart.CoordinateNumberForOutput = 4 
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
                                
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        gs = gridspec.GridSpec(6,1,
                               height_ratios = [1,2,2,2,2,4]
                               )
        ax = plt.subplot(gs[0])
        plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForPhase),"b")
        plt.plot(fireNotifierResp.firingTimes,fireNotifierResp.firingTimesSpikes(),"bo")
        plt.yticks([0.0,1.0])
        plt.ylabel(r"$\Phi(t)$")
        #plt.xlabel("Time [s]")
        plt.ylim(0.0,1.2)
        plt.setp(ax.get_xticklabels(), visible = False)
        ax = plt.subplot(gs[1])
        #todo: the guy below is all flat.
        #plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForOutput),"g",linewidth=2)
        plt.plot(fireNotifier.firingTimes,fireNotifier.firingTimesSpikes(),"go")
        plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"g",linewidth=2)      
        plt.yticks([])
        plt.setp(ax.get_xticklabels(), visible = False)        
        #plt.xlabel("Time [s]")
        plt.ylabel(r"$r_{n}(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        ax = plt.subplot(gs[2])
        plt.yticks([])
        plt.ylabel(r"$r(t)$")                
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForRate),"v")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[3])
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForPhase),"r")
        #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
        plt.plot(fireNotifierHeart.firingTimes,fireNotifierHeart.firingTimesSpikes(),"ro")
        plt.ylim(0.0,1.2)
        plt.yticks([0.0,1.0])
        #plt.xlabel("Time [s]")
        plt.ylabel(r"$\varphi(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[4])
        plt.ylabel("ISI")                
        plt.plot(fireNotifierHeart.firingTimes[1:],fireNotifierHeart.ISI()[1:],"b+-")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        plt.subplot(gs[5])        
        for varNumber in (0,1,2,3):         
            plt.plot(allTimes,seriesNotifier.GetVar(varNumber))
        plt.ylim(0.0,300.0)
        plt.xlabel("Time [s]")
        plt.ylabel("BP [mmHg]")
        fname = sys._getframe().f_code.co_name + ".png"
        print("Test result in %s" % fname)
        plt.savefig(fname)   
        
    def test_FullCouplingPaperSettings(self):
        """
        This test handles the majority of physiological control of heart 
        rate. Respiratory IAF periodically forces the heart via the 
        RespiratoryDelayedSmearedHeartActionForce. The force is applied to heart
        drive. Initially the heart is in phase with respiration, as both 
        periods are commensurate (3:1) and initial phase of both is 0.
        """
        dimension = 10
        npoints = 2400
        settings = KickedWindkesselModel.KickedWindkesselModelSettings() 
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 9 #Coordinate number, to  which the state will be written
        force.CoordinateNumberForInput = 6
        force.KickAmplitude = 1.0 #Kick amplitude
        force.DelayTau = 0.3 #Time from firing order to actual kick
        force.SamplingTime = 0.01 # required to normalize delay time.                
        force.DecayTau = 0.3 # Time by which the drive decays
        force.Notify = fireNotifier.Notify        

        iafResp = IntegrateAndFire()
        iafResp.SamplingTime = 0.01
        iafResp.SetPhaseVelocityFromBPM(20)
        iafResp.SetInitialPhase(0.0)
        iafResp.KickAmplitude = 0.1        
        iafResp.CoordinateNumberForRate = -1
        iafResp.CoordinateNumberForPhase = 5
        iafResp.CoordinateNumberForForceInput = -1 #not used
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 6         
        fireNotifierResp = FiringTimesNotifier()
        iafResp.Notify = fireNotifierResp.Notify
        

        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.01        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        iafHeart.CoordinateNumberForRate = 8
        iafHeart.CoordinateNumberForPhase = 7
        iafHeart.CoordinateNumberForForceInput = force.CoordinateNumber        
        # Coordinate number, to  which the state will be written        
        iafHeart.CoordinateNumberForOutput = 4 
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
                                
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        gs = gridspec.GridSpec(6,1,
                               height_ratios = [1,2,2,2,2,4]
                               )
                               
        ax = plt.subplot(gs[0])
        plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForPhase),"b")
        plt.plot(fireNotifierResp.firingTimes,fireNotifierResp.firingTimesSpikes(),"bo")
        plt.yticks([0.0,1.0])
        plt.ylabel(r"$\Phi(t)$")
        #plt.xlabel("Time [s]")
        plt.ylim(0.0,1.2)
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[1])
        #todo: the guy below is all flat.
        #plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForOutput),"g",linewidth=2)
        plt.plot(fireNotifier.firingTimes,fireNotifier.firingTimesSpikes(),"go")
        plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"g",linewidth=2)      
        plt.yticks([])
        plt.setp(ax.get_xticklabels(), visible = False)        
        #plt.xlabel("Time [s]")
        plt.ylabel(r"$r_{n}(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        
        
        ax = plt.subplot(gs[2])
        plt.yticks([])
        plt.ylabel(r"$r(t)$")                
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForRate),"v")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[3])
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForPhase),"r")
        #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
        plt.plot(fireNotifierHeart.firingTimes,fireNotifierHeart.firingTimesSpikes(),"ro")
        plt.ylim(0.0,1.2)
        plt.yticks([0.0,1.0])
        #plt.xlabel("Time [s]")
        plt.ylabel(r"$\varphi(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[4])
        plt.ylabel("ISI")
        plt.plot(fireNotifierHeart.firingTimes[1:],fireNotifierHeart.ISI()[1:],"b+-")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        plt.subplot(gs[5])        
        for varNumber in (0,1,2,3):         
            plt.plot(allTimes,seriesNotifier.GetVar(varNumber))
        plt.ylim(0.0,300.0)
        plt.xlabel("Time [s]")
        plt.ylabel("BP [mmHg]")
        fname = sys._getframe().f_code.co_name + ".png"
        logging.info("Test result in %s" % fname)
        plt.savefig(fname)   

        
    def test_RespiratoryPeriod(self):
        """
        Is the influence of respiratory period visible       
        """
        dimension = 10
        npoints = 2400
        settings = KickedWindkesselModel.KickedWindkesselModelSettings() 
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 9 #Coordinate number, to  which the state will be written
        force.CoordinateNumberForInput = 6
        force.KickAmplitude = 1.0 #Kick amplitude
        force.DelayTau = 0.3 #Time from firing order to actual kick
        force.SamplingTime = 0.01 # required to normalize delay time.                
        force.DecayTau = 0.3 # Time by which the drive decays
        force.Notify = fireNotifier.Notify        

        iafResp = IntegrateAndFire()
        iafResp.SamplingTime = 0.01
        iafResp.SetPhaseVelocityFromBPM(12)
        iafResp.SetInitialPhase(0.0)
        iafResp.KickAmplitude = 0.1        
        iafResp.CoordinateNumberForRate = -1
        iafResp.CoordinateNumberForPhase = 5
        iafResp.CoordinateNumberForForceInput = -1 #not used
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 6         
        fireNotifierResp = FiringTimesNotifier()
        iafResp.Notify = fireNotifierResp.Notify
        

        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.01        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        iafHeart.CoordinateNumberForRate = 8
        iafHeart.CoordinateNumberForPhase = 7
        iafHeart.CoordinateNumberForForceInput = force.CoordinateNumber        
        # Coordinate number, to  which the state will be written        
        iafHeart.CoordinateNumberForOutput = 4 
        fireNotifierHeart = FiringTimesNotifier()
        iafHeart.Notify = fireNotifierHeart.Notify

        collector = AbpmFiducialPointsCollector(0)#ABPM = 0
        seriesNotifier = SeriesNotifier(dimension,npoints)
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
        
        settings.heartActionForce = HeartActionForceChain([iafResp,force,iafHeart])
        settings.CoordinateNumberForRespPhase = iafResp.CoordinateNumberForPhase
        model = KickedWindkesselModel(settings,dimension)
        model.param.Npoints = npoints
        
        chain = NotifierChain((collector,seriesNotifier))
        model.Notify = chain.Notify
            
        model.IterateToNotifiers()
                                
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        gs = gridspec.GridSpec(7,1,
                               height_ratios = [1,2,2,2,2,2,4]
                               )
        ax = plt.subplot(gs[0])
        plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForPhase),"b")
        plt.plot(fireNotifierResp.firingTimes,fireNotifierResp.firingTimesSpikes(),"bo")
        plt.yticks([0.0,1.0])
        plt.ylabel(r"$\Phi(t)$")
        #plt.xlabel("Time [s]")
        plt.ylim(0.0,1.2)
        plt.setp(ax.get_xticklabels(), visible = False)
        

        ax = plt.subplot(gs[1])        
        plt.plot(allTimes,seriesNotifier.GetVar(3))
        plt.ylabel(r"$p_I(t)$")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[2])
        #todo: the guy below is all flat.
        #plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForOutput),"g",linewidth=2)
        plt.plot(fireNotifier.firingTimes,fireNotifier.firingTimesSpikes(),"go")
        plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"g",linewidth=2)      
        plt.yticks([])
        plt.setp(ax.get_xticklabels(), visible = False)        
        #plt.xlabel("Time [s]")
        plt.ylabel(r"$r_{n}(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        ax = plt.subplot(gs[3])
        plt.yticks([])
        plt.ylabel(r"$r(t)$")                
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForRate),"v")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[4])
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForPhase),"r")
        #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
        plt.plot(fireNotifierHeart.firingTimes,fireNotifierHeart.firingTimesSpikes(),"ro")
        plt.ylim(0.0,1.2)
        plt.yticks([0.0,1.0])
        #plt.xlabel("Time [s]")
        plt.ylabel(r"$\varphi(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[5])
        plt.ylabel("ISI")                
        plt.plot(fireNotifierHeart.firingTimes[1:],fireNotifierHeart.ISI()[1:],"b+-")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        plt.subplot(gs[6])        
        for varNumber in (0,1,2):         
            plt.plot(allTimes,seriesNotifier.GetVar(varNumber))
        plt.ylim(0.0,300.0)
        plt.xlabel("Time [s]")
        plt.ylabel("BP [mmHg]")
        fname = sys._getframe().f_code.co_name + ".png"
        logging.info("Test result in %s" % fname)
        plt.savefig(fname)
        
    def test_HighAmplitude(self):
        """
        Is the influence of high intrathoracic pressure amplitude visible
        at disabled heart control
        """
        dimension = 10
        npoints = 2400
        settings = KickedWindkesselModel.KickedWindkesselModelSettings() 
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 9 #Coordinate number, to  which the state will be written
        force.CoordinateNumberForInput = 6
        force.KickAmplitude = 0.0 #Kick amplitude
        force.DelayTau = 0.3 #Time from firing order to actual kick
        force.SamplingTime = 0.01 # required to normalize delay time.                
        force.DecayTau = 0.3 # Time by which the drive decays
        force.Notify = fireNotifier.Notify        

        iafResp = IntegrateAndFire()
        iafResp.SamplingTime = 0.01
        iafResp.SetPhaseVelocityFromBPM(6)
        iafResp.SetInitialPhase(0.0)
        iafResp.KickAmplitude = 0.1        
        iafResp.CoordinateNumberForRate = -1
        iafResp.CoordinateNumberForPhase = 5
        iafResp.CoordinateNumberForForceInput = -1 #not used
        # Coordinate number, to  which the state will be written
        #this variable will be used for coupling.
        iafResp.CoordinateNumberForOutput = 6         
        fireNotifierResp = FiringTimesNotifier()
        iafResp.Notify = fireNotifierResp.Notify
        

        iafHeart = IntegrateAndFire()
        iafHeart.SamplingTime = 0.01        
        iafHeart.SetPhaseVelocityFromBPM(67)
        iafHeart.phaseEfectivenessCurve = phaseEfectivenessCurveSH
        iafHeart.CoordinateNumberForRate = 8
        iafHeart.CoordinateNumberForPhase = 7
        iafHeart.CoordinateNumberForForceInput = force.CoordinateNumber        
        # Coordinate number, to  which the state will be written        
        iafHeart.CoordinateNumberForOutput = 4 
        fireNotifierHeart = FiringTimesNotifier()        
        #iafHeart.Notify = fireNotifierHeart.Notify


        collector = AbpmFiducialPointsCollector(0)#ABPM = 0
        iafHeart.NotifyFunctions = [fireNotifierHeart.Notify,collector.HeartOpenNotifier]
        mmNotifier = MinMaxNotifier(0)
        mmNotifier.MinNotifier = collector.AbpmMinNotifier
        mmNotifier.MaxNotifier = collector.AbpmMaxNotifier



        seriesNotifier = SeriesNotifier(dimension,npoints)
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
        
        settings.heartActionForce = HeartActionForceChain([iafResp,force,iafHeart])
        settings.CoordinateNumberForRespPhase = iafResp.CoordinateNumberForPhase
        settings.p_I1 = 3.0 # default is 0.1
        model = KickedWindkesselModel(settings,dimension)
        model.param.Npoints = npoints
        
        chain = NotifierChain((mmNotifier,collector,seriesNotifier))
        model.Notify = chain.Notify
            
        model.IterateToNotifiers()
        allItems = collector.GetFiducialPointsList()
        series = {}
        series["systolicABP"] = np.array(allItems[:,1])
        series["diastolicABP"] = np.array(allItems[:,2] )
        series["meanABP"] = np.array(allItems[:,3])
        for key in series.keys():
            av = np.average(series[key])
            sd = np.std(series[key])
            logging.info("%s: av: %.2lf sd: %.2lf, npoints: %d"%(key,av,sd,len(series[key])))
        
        print(allItems)
        
        fname = sys._getframe().f_code.co_name + ".png"
        KickedWindkesselModelVisualization(fname,allTimes,allItems,seriesNotifier,
                                       fireNotifierResp,
                                       fireNotifier,
                                       fireNotifierHeart,
                                       iafResp.CoordinateNumberForPhase,
                                       force.CoordinateNumber,
                                       iafHeart.CoordinateNumberForRate,
                                       iafHeart.CoordinateNumberForPhase
                                       )
                   
def main():

    logging.basicConfig(level=logging.INFO)
    unittest.main()

if __name__ == "__main__":
    main()   

