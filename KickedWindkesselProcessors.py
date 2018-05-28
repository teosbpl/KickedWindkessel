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


def ZCCombinator(npoints=4000,kickAmplitude=-0.05,kickDelay=2.0,p_I1=1.5,respBPM =20.0,heartBPM=66.0):
    #stepShiftLinspace = np.linspace(-0.5,1.25,20)
    #stepShiftLinspace = np.linspace(0,2.0*np.pi,100)
    Sav = []
    Ssd = []
    Dav = []
    Dsd = []
    Mav = []
    Msd = []

    ZavLinspace = np.linspace(0.2,2.2,20)
    ZvcLinspace = np.linspace(0.2,2.2,20)
    ZcaLinspace = np.linspace(100.0,200.0,20)


    for i,Zav in enumerate(ZavLinspace):
        for j,Zvc in enumerate(ZvcLinspace):
            for k,Zca in enumerate(ZcaLinspace):

                dimension = 10
                #npoints = 3600
                #npoints = 4000
                settings,model,force,iafResp,iafHeart,collector,seriesNotifier,fireNotifierResp,fireNotifierHeart,fireNotifier = StandardModelSetup(dimension,npoints)
                iafResp.SetPhaseVelocityFromBPM(respBPM)
                force.DelayTau = kickDelay
                force.KickAmplitude = kickAmplitude
                model.settings.p_I1 = p_I1 # default is 2.0
                model.settings.Zca = Zca #1 / 0.006 # 166.[6]
                #model.settings.Zca = 1 / 0.01 # 166.[6]
                model.settings.Zav = Zav #1 / 0.9 # was 0.9 TB 27.05.2018
                model.settings.Zvc = Zvc # 1 / 0.37 # was 0.005 TB 27.05.2018
                iafHeart.SetPhaseVelocityFromBPM(heartBPM) #45 worked OK
                iafResp.SetPhaseVelocityFromBPM(respBPM)
                model.IterateToNotifiers()
                arterial = seriesNotifier.GetVar(0)
                cardiopulmonary = seriesNotifier.GetVar(1)
                venous = seriesNotifier.GetVar(2)
                perc = (np.percentile(arterial, 10),
                np.percentile(arterial, 90),
                np.percentile(cardiopulmonary, 10),
                np.percentile(cardiopulmonary, 90),
                np.percentile(venous, 10),
                np.percentile(venous, 90))

                logging.warning("%d\tZav\t%f\t%d\tZvc\t%f\t%d\tZca\t%f\tA\t%f\t%f\tC\t%f\t%f\tV\t%f\t%f"
                %(i,Zav,j,Zvc,k,Zca,perc[0],perc[1],perc[2],perc[3],perc[4],perc[4]))


def PhaseShiftProcessorOld():
    logging.basicConfig(level=logging.ERROR)
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


    ax.set_ylim([ylim[0],ylim[1]])
    ax.set_xlim([xlim[0],xlim[1]])
    plt.savefig("phaseShiftFunction100.png")
    plt.show()
    return (Sav, Ssd, Dav, Dsd, Mav, Msd), (settings, model, force, iafResp, iafHeart, collector, seriesNotifier, fireNotifierResp, fireNotifierHeart, fireNotifier)

def StandardModelSetup(dimension,npoints):
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
        iafHeart.SetPhaseVelocityFromBPM(66) #45 worked OK
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


        settings.heartActionForce = HeartActionForceChain([iafResp,force,iafHeart])
        settings.CoordinateNumberForRespPhase = iafResp.CoordinateNumberForPhase

        model = KickedWindkesselModel(settings,dimension)
        model.param.Npoints = npoints

        chain = NotifierChain((mmNotifier,collector,seriesNotifier))
        model.Notify = chain.Notify

        return settings, model, force, iafResp, iafHeart, collector, seriesNotifier, fireNotifierResp, fireNotifierHeart, fireNotifier


def PhaseShiftProcessor(npoints=4000,firstAfterWarmup=25,kickAmplitude=0.04,p_I1=2.0,respBPM =20.0,basalValues = [0.0,0.0,0.0],stepShiftLinspace= np.linspace(0,3.0,20)):
    #stepShiftLinspace = np.linspace(-0.5,1.25,20)
    #stepShiftLinspace = np.linspace(0,2.0*np.pi,100)
    Sav = []
    Ssd = []
    Dav = []
    Dsd = []
    Mav = []
    Msd = []

    for n,stepShift in enumerate(stepShiftLinspace):

        dimension = 10
        #npoints = 3600
        #npoints = 4000
        settings,model,force,iafResp,iafHeart,collector,seriesNotifier,fireNotifierResp,fireNotifierHeart,fireNotifier = StandardModelSetup(dimension,npoints)
        force.DelayTau = stepShift
        force.KickAmplitude = kickAmplitude
        model.settings.p_I1 = p_I1 # default is 2.0

        iafResp.SetPhaseVelocityFromBPM(respBPM)
        model.IterateToNotifiers()
        allItems = collector.GetFiducialPointsList()
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
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
        #firstAfterWarmup = 25
        skipThisLine = False
        for channel in [1,2,3]:# 1-systloic 2-diastolic 3-mean
            averages.append(np.mean(allItems[firstAfterWarmup:,channel]))
            standard_deviations.append(np.std(allItems[firstAfterWarmup:,channel]))
            if standard_deviations[-1] > 30.0:
                logging.info("Standard deviation went crazy for data:")
                logging.info("AV:%lf" % (np.mean(allItems[firstAfterWarmup:,channel])))
                logging.info("SD:%lf" % (np.std(allItems[firstAfterWarmup:,channel])))
                logging.info(allItems[firstAfterWarmup:,channel])
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
        print("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d"%(n,force.DelayTau,averages[0],standard_deviations[0],averages[1],standard_deviations[1],averages[2],standard_deviations[2],len(allItems)))

    #normalize values
    if basalValues[0] == 0.0:
    	basalValues[0] = Ssd[0]
    if basalValues[1] == 0.0:
    	basalValues[1] = Dsd[0]
    if basalValues[2] == 0.0:
    	basalValues[2] = Msd[0]

    for i,v in enumerate(Ssd):
        Ssd[i] = 100.0 * (v / basalValues[0])
    for i,v in enumerate(Dsd):
        Dsd[i] = 100.0 * (v / basalValues[1])
    for i,v in enumerate(Msd):
        Msd[i] = 100.0 * (v / basalValues[2])

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
    #ax.set_ylim(60,140)
    plt.savefig("phaseShiftProcessor_r0-%lf_pI1-%lf_resp%lf.png" % (kickAmplitude,model.settings.p_I1,respBPM))
    plt.close(fig)
    #plt.show()
    return basalValues,(Sav,Ssd,Dav,Dsd,Mav,Msd),(settings, model, force, iafResp, iafHeart, collector, seriesNotifier, fireNotifierResp, fireNotifierHeart, fireNotifier)

def KickAmplitudeProcessor(npoints=4000,firstAfterWarmup=25,forceDelayTau=0.04,forceDecayTau=0.3,p_I1=2.0,respBPM =20.0,basalValues = [0.0,0.0,0.0],KickAmplitudeLinspace= np.linspace(0,0.04,20)):
    #stepShiftLinspace = np.linspace(-0.5,1.25,20)
    #KickAmplitudeLinspace = np.linspace(0,0.1,50)
    #stepShiftLinspace = np.linspace(0,2.0*np.pi,100)
    Sav = []
    Ssd = []
    Dav = []
    Dsd = []
    Mav = []
    Msd = []

    for n,kickAmplitude in enumerate(KickAmplitudeLinspace):

        dimension = 10
        settings,model,force,iafResp,iafHeart,collector,seriesNotifier,fireNotifierResp,fireNotifierHeart,fireNotifier = StandardModelSetup(dimension,npoints)

        force.DelayTau = forceDelayTau
        force.DecayTau = forceDecayTau
        force.KickAmplitude = kickAmplitude
        model.settings.p_I1 = p_I1 # default is 2.0

        iafResp.SetPhaseVelocityFromBPM(respBPM)
        model.IterateToNotifiers()
        allItems = collector.GetFiducialPointsList()
        allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
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

        skipThisLine = False
        for channel in [1,2,3]:# 1-systloic 2-diastolic 3-mean
            averages.append(np.mean(allItems[firstAfterWarmup:,channel]))
            standard_deviations.append(np.std(allItems[firstAfterWarmup:,channel]))
            if standard_deviations[-1] > 30.0:
                logging.info("Standard deviation went crazy for data:")
                logging.info("AV:%lf" % (np.mean(allItems[firstAfterWarmup:,channel])))
                logging.info("SD:%lf" % (np.std(allItems[firstAfterWarmup:,channel])))
                logging.info(allItems[firstAfterWarmup:,channel])
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
        print("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d"%(n,force.DelayTau,averages[0],standard_deviations[0],averages[1],standard_deviations[1],averages[2],standard_deviations[2],len(allItems)))

    #normalize values
    if basalValues[0] == 0.0:
    	basalValues[0] = Ssd[0]
    if basalValues[1] == 0.0:
    	basalValues[1] = Dsd[0]
    if basalValues[2] == 0.0:
    	basalValues[2] = Msd[0]

    for i,v in enumerate(Ssd):
        Ssd[i] = 100.0 * (v / basalValues[0])
    for i,v in enumerate(Dsd):
        Dsd[i] = 100.0 * (v / basalValues[1])
    for i,v in enumerate(Msd):
        Msd[i] = 100.0 * (v / basalValues[2])

    #pdb.set_trace()
    fig,ax = plt.subplots()
    ax.ticklabel_format(useOffset = False)
    #plt.subplot(3,1,1)
    ax.plot(KickAmplitudeLinspace,Ssd,"+r",linestyle="-",linewidth=1) # SAP thin curve
    #plt.subplot(3,1,2)
    ax.plot(KickAmplitudeLinspace,Dsd,"og",linestyle="--",linewidth=3) # DAP thick dashed
    #plt.subplot(3,1,3)
    ax.plot(KickAmplitudeLinspace,Msd,"vb",linestyle="-",linewidth=3) # MAP thick solid
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
    #ax.set_ylim(60,140)
    fname = "kickAmplitudeProcessor_r0-%lf_pI1-%lf_resp%lf_decay%lf-delay%lf.png" % (kickAmplitude,model.settings.p_I1,respBPM,forceDecayTau,forceDelayTau)
    print("Writing to: %s" % (fname))
    plt.savefig(fname)
    #plt.show()
    plt.close(fig)
    #plt.show()
    return basalValues,(Sav,Ssd,Dav,Dsd,Mav,Msd),(settings, model, force, iafResp, iafHeart, collector, seriesNotifier, fireNotifierResp, fireNotifierHeart, fireNotifier)

if __name__ == "__main__":
    logging.basicConfig(level=logging.WARN)
    npoints=4000
    firstAfterWarmup=25
    p_I1 = 3.0
    respBPM =20.0
    basalValues = [0.0,0.0,0.0]
    basalValues0 = [0.0,0.0,0.0]
    doFinalShift = False
    doPhaseShift = False
    doDecay = False
    doKickAmplitude = False
    doCombinator = True
    if doCombinator:
        ZCCombinator()
    elif doFinalShift:
        # normalization
        basalValues0,values = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,0.0,p_I1,respBPM,basalValues,np.linspace(0.0,0.0,1)))

        #kickAmplitudeLinspace = np.linspace(1.0,10.0,20)
        kickAmplitudeLinspace = np.linspace(-0.031,-0.001,31)
        stepShiftLinspace = np.linspace(0,3.0,100)
        for i,kickAmplitude in enumerate(kickAmplitudeLinspace):#[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
        	print("============================= kickAmplitude %lf normed by %s ========================= " %(kickAmplitude,basalValues0))
        	basalValues,values = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues0,stepShiftLinspace))
    elif doPhaseShift:
        # normalization
        basalValues0,values = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,0.0,p_I1,respBPM,basalValues,np.linspace(0.0,0.0,1)))

        #kickAmplitudeLinspace = np.linspace(1.0,10.0,20)
        kickAmplitudeLinspace = np.linspace(-0.05,-0.001,20)
        stepShiftLinspace = np.linspace(0,3.0,50)
        for i,kickAmplitude in enumerate(kickAmplitudeLinspace):#[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
        	print("============================= kickAmplitude %lf normed by %s ========================= " %(kickAmplitude,basalValues0))
        	basalValues,values = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues0,stepShiftLinspace))
    elif doDecay:
        KickAmplitudeLinspace = np.linspace(0.0,0.0,1)
        forceDelayTau=0.5 # dummy as there is no input
        forceDecayTau=0.3 # dummy as there is no input
        basalValues0,values = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues,KickAmplitudeLinspace))

        decayTauLinspace = np.linspace(0.1,3.1,100)
        decayTauLinspace = np.linspace(0.05,0.25,100)
        KickAmplitudeLinspace = np.linspace(0.0,0.1,100)
        for i,forceDecayTau in enumerate(decayTauLinspace):#[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            print("============================= forceDecayTau %lf normed by %s ========================= " %(forceDecayTau,basalValues0))
            basalValues,values = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues0,KickAmplitudeLinspace))
    elif doKickAmplitude:
        KickAmplitudeLinspace = np.linspace(0.0,0.0,1)
        forceDelayTau=0.5 # dummy as there is no input
        forceDecayTau=0.3 # dummy as there is no input
        basalValues0,values = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues,KickAmplitudeLinspace))

        #decayTauLinspace = np.linspace(0.1,3.1,100)
        #decayTauLinspace = np.linspace(0.05,0.25,100)
        stepShiftLinspace = np.linspace(0,3.0,20)
        KickAmplitudeLinspace = np.linspace(0.0,-0.5,50)
        for i,forceDelayTau in enumerate(stepShiftLinspace):#[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            print("============================= forceDecayTau %lf normed by %s ========================= " %(forceDecayTau,basalValues0))
            basalValues,values = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues0,KickAmplitudeLinspace))