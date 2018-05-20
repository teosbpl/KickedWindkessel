# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import copy
import logging
import pdb
import sys
import numpy as np
import matplotlib.pyplot as plt
from KickedWindkesselModel import KickedWindkesselModel
from HeartActionForce import RectangularHeartActionForce,RespiratoryDelayedSmearedHeartActionForce, HeartActionForceChain
from AbpmFiducialPointsCollector import AbpmFiducialPointsCollector
from Notifiers import SeriesNotifier, MinMaxNotifier, NotifierChain, FiringTimesNotifier
from IntegrateAndFire import IntegrateAndFire, phaseEfectivenessCurveSH
from KickedWindkesselModelVisualization import KickedWindkesselModelVisualization

def PhaseShiftProcessorOld():
    logging.basicConfig(level=logging.ERROR)
    # to jest procesor w funkcji fazy dla heart action force prostokątnego
    #przerobić phase shift processor żeby korzystał z FiducialPoints
    #przerobić też żeby miał nomenklaturę zgodną z publikacją i żeby wprowadzał opóźnienie w kick a nie
    #w oddechu dorobić pętlę w funkcji amplitudy, ale hmm gdzie tu jest amplituda?
    stepShiftLinspace = np.linspace(0,1.25,10)
    stepShiftLinspace = np.linspace(0,2.0*np.pi,100)
    Sav = []
    Ssd = []
    Dav = []
    Dsd = []
    Mav = []
    Msd = []

    for stepShift in stepShiftLinspace:
        settings = KickedWindkesselModel.KickedWindkesselModelSettings()
        settings.breathingPhi0 = stepShift
        settings.heartActionForce = RectangularHeartActionForce()
        settings.heartActionForce.StepShift = stepShift
        model = KickedWindkesselModel(settings)
        model.param.Npoints = 1000
        collector = AbpmFiducialPointsCollector(0)#ABPM = 0. Other compartments have different numbers.
        settings.heartActionForce.Notify = collector.HeartOpenNotifier
        mmNotifier = MinMaxNotifier(0)
        mmNotifier.MinNotifier = collector.AbpmMinNotifier
        mmNotifier.MaxNotifier = collector.AbpmMaxNotifier
        seriesNotifier = SeriesNotifier(model.param.dimension,model.param.Npoints)
        
        chain = NotifierChain((mmNotifier,collector,seriesNotifier))
        model.Notify = chain.Notify
            
        model.IterateToNotifiers()    
        
        allItems = collector.GetFiducialPointsList()
        averages = []
        standard_deviations = []

        for channel in [1,2,3]:# 1-systloic 2-diastolic 3-mean
            averages.append(np.mean(allItems[:,channel]))
            standard_deviations.append(np.std(allItems[:,channel]))

        Sav.append(averages[0])
        Ssd.append(standard_deviations[0])
        Dav.append(averages[1])
        Dsd.append(standard_deviations[1])
        Mav.append(averages[2])
        Msd.append(standard_deviations[2])
        print("%f\t%f\t%f\t%f\t%f\t%f\t%f"%(settings.heartActionForce.StepShift,averages[0],standard_deviations[0],averages[1],standard_deviations[1],averages[2],standard_deviations[2]))

    #normalize values 
    basalValue = Ssd[0]
    for i,v in enumerate(Ssd):
        Ssd[i] = 100.0 * (v / basalValue)
    basalValue = Dsd[0]
    for i,v in enumerate(Dsd):
        Dsd[i] = 100.0 * (v / basalValue)
    basalValue = Msd[0]
    for i,v in enumerate(Msd):
        Msd[i] = 100.0 * (v / basalValue)

    #pdb.set_trace()        
    fig,ax = plt.subplots()
    ax.ticklabel_format(useOffset = False)
    #plt.subplot(3,1,1)    
    ax.plot(stepShiftLinspace,Ssd,"r",linestyle="-",linewidth=1) # SAP thin curve    
    #plt.subplot(3,1,2)
    ax.plot(stepShiftLinspace,Dsd,"g",linestyle="--",linewidth=3) # DAP thick dashed
    #plt.subplot(3,1,3)
    ax.plot(stepShiftLinspace,Msd,"b",linestyle="-",linewidth=3) # MAP thick solid
    for y in (100.0,): # 80.0,100.0,120.0):
        plt.axhline(y=y,color='k',linestyle='-')    
    #ax.fill(30, 30, fill=False, hatch='\\')
    #pdb.set_trace()
    ylim = copy.deepcopy(ax.get_ylim())# keep original limits.\
    xlim = copy.deepcopy(ax.get_xlim())# keep original limits.\    
    
    
#==============================================================================
#     ax.fill([0,0,1,1,0],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
#     offsetOfHalf = len(Msd)/2
#     rightHalfMsd = np.array(Msd[offsetOfHalf:])
#     print(rightHalfMsd)
#     rightLargerThan100 = np.where(rightHalfMsd>100)[0][0]
#     rightHatchedEdge = stepShiftLinspace[offsetOfHalf+rightLargerThan100]
#     print(rightHatchedEdge)
#     ax.fill([rightHatchedEdge,rightHatchedEdge,xlim[1],xlim[1],rightHatchedEdge],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
# 
#==============================================================================
    #ax.bar(range(0,2), np.arange(ylim[0],ylim[1]), color='red', edgecolor='black', hatch="/")
    ax.set_ylim([ylim[0],ylim[1]])    
    ax.set_xlim([xlim[0],xlim[1]])    
    plt.savefig("phaseShiftFunction100.png")
    plt.show()

def PhaseShiftProcessor():
    #stepShiftLinspace = np.linspace(-0.5,1.25,20)
    stepShiftLinspace = np.linspace(0,3.0,50)
    #stepShiftLinspace = np.linspace(0,2.0*np.pi,100)
    Sav = []
    Ssd = []
    Dav = []
    Dsd = []
    Mav = []
    Msd = []

    for stepShift in stepShiftLinspace:    
    
        dimension = 20
        npoints = 3600
        npoints = 5000        
        settings = KickedWindkesselModel.KickedWindkesselModelSettings() 
        
        fireNotifier = FiringTimesNotifier()        
        force = RespiratoryDelayedSmearedHeartActionForce()        
        force.CoordinateNumber = 9 #Coordinate number, to  which the state will be written
        force.CoordinateNumberForInput = 6
        force.KickAmplitude = 1.0 #Kick amplitude
        force.DelayTau = stepShift #Time from firing order to actual kick
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
#        series = {}
#        series["systolicABP"] = np.array(allItems[:,1])
#        series["diastolicABP"] = np.array(allItems[:,2] )
#        series["meanABP"] = np.array(allItems[:,3])
#        for key in series.keys():
#            av = np.average(series[key])
#            sd = np.std(series[key])
#            logging.info("%s: av: %.2lf sd: %.2lf, npoints: %d"%(key,av,sd,len(series[key])))
#        
#        print(allItems)    
        averages = []
        standard_deviations = []
        firstAfterWarmup = 20
        skipThisLine = False
        for channel in [1,2,3]:# 1-systloic 2-diastolic 3-mean
            averages.append(np.mean(allItems[firstAfterWarmup:,channel]))
            standard_deviations.append(np.std(allItems[firstAfterWarmup:,channel]))
            if standard_deviations[-1] > 40.0:
                skipThisLine = True
        if skipThisLine:
            fname = sys._getframe().f_code.co_name + "_delay_%f.png" % (stepShift)
            KickedWindkesselModelVisualization(fname,
                                       allTimes,
                                       allItems,
                                       seriesNotifier,
                                       fireNotifierResp,
                                       fireNotifier,
                                       fireNotifierHeart,
                                       iafResp.CoordinateNumberForPhase,
                                       force.CoordinateNumber,
                                       iafHeart.CoordinateNumberForRate,
                                       iafHeart.CoordinateNumberForPhase)
            print("SKIPPING %f\t%f\t%f\t%f\t%f\t%f\t%f\t%d"%(force.DelayTau,averages[0],standard_deviations[0],averages[1],standard_deviations[1],averages[2],standard_deviations[2],len(allItems)))
            #continue
        Sav.append(averages[0])
        Ssd.append(standard_deviations[0])
        Dav.append(averages[1])
        Dsd.append(standard_deviations[1])
        Mav.append(averages[2])
        Msd.append(standard_deviations[2])
        print("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d"%(force.DelayTau,averages[0],standard_deviations[0],averages[1],standard_deviations[1],averages[2],standard_deviations[2],len(allItems)))

    #normalize values 
    basalValue = Ssd[0]
    for i,v in enumerate(Ssd):
        Ssd[i] = 100.0 * (v / basalValue)
    basalValue = Dsd[0]
    for i,v in enumerate(Dsd):
        Dsd[i] = 100.0 * (v / basalValue)
    basalValue = Msd[0]
    for i,v in enumerate(Msd):
        Msd[i] = 100.0 * (v / basalValue)

    #pdb.set_trace()        
    fig,ax = plt.subplots()
    ax.ticklabel_format(useOffset = False)
    #plt.subplot(3,1,1)    
    ax.plot(stepShiftLinspace,Ssd,"+r",linestyle="-",linewidth=1) # SAP thin curve    
    #plt.subplot(3,1,2)
    ax.plot(stepShiftLinspace,Dsd,"og",linestyle="--",linewidth=3) # DAP thick dashed
    #plt.subplot(3,1,3)
    ax.plot(stepShiftLinspace,Msd,"vb",linestyle="-",linewidth=3) # MAP thick solid
    for y in (100.0,): # 80.0,100.0,120.0):
        plt.axhline(y=y,color='k',linestyle='-')    
    
    #ax.fill(30, 30, fill=False, hatch='\\')
    #pdb.set_trace()
    ylim = copy.deepcopy(ax.get_ylim())# keep original limits.\
    xlim = copy.deepcopy(ax.get_xlim())# keep original limits.\    
    
    
#==============================================================================
#     ax.fill([0,0,1,1,0],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
#     offsetOfHalf = len(Msd)/2
#     rightHalfMsd = np.array(Msd[offsetOfHalf:])
#     print(rightHalfMsd)
#     rightLargerThan100 = np.where(rightHalfMsd>100)[0][0]
#     rightHatchedEdge = stepShiftLinspace[offsetOfHalf+rightLargerThan100]
#     print(rightHatchedEdge)
#     ax.fill([rightHatchedEdge,rightHatchedEdge,xlim[1],xlim[1],rightHatchedEdge],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
# 
#==============================================================================
    #ax.bar(range(0,2), np.arange(ylim[0],ylim[1]), color='red', edgecolor='black', hatch="/")
    ax.set_ylim([ylim[0],ylim[1]])    
    ax.set_xlim([xlim[0],xlim[1]])    
    ax.set_ylim(60,140)
    plt.savefig("phaseShiftProcessor.png")
    plt.show()

if __name__ == "__main__":
    PhaseShiftProcessor()