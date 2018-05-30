# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import copy
import logging
import pdb
import sys
import os.path
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as plt
from KickedWindkesselModel import KickedWindkesselModel
from HeartActionForce import RectangularHeartActionForce,RespiratoryDelayedSmearedHeartActionForce, HeartActionForceChain
from AbpmFiducialPointsCollector import AbpmFiducialPointsCollector
from Notifiers import SeriesNotifier, MinMaxNotifier, NotifierChain, FiringTimesNotifier
from IntegrateAndFire import IntegrateAndFire, phaseEfectivenessCurveSH
from KickedWindkesselModelVisualization import KickedWindkesselModelVisualization
from KickedWindkesselProcessors import StandardModelSetup,PhaseShiftProcessor,KickAmplitudeProcessor

def findOuterLimits(x,y,limitValue,moreThan=True,divider=0.5):
    """
    The method finds limits of the regions in which the condition is met.
    It is assumed, that the regions extend from the right end to some value
    and from the left end to some value.
    Comparison operator and dividing point may be chosen.
    """
    nCalcPoints = len(y)
    if moreThan:
        ww = np.where(y>limitValue)[0]
    else:
        ww = np.where(y<limitValue)[0]

    leftPart = ww[np.where(ww<nCalcPoints*divider)]
    rightPart = ww[np.where(ww>nCalcPoints*divider)]
    leftAreaRightEdge = x[np.max(leftPart)]
    rightAreaLeftEdge = x[np.min(rightPart)]
    return leftAreaRightEdge,rightAreaLeftEdge

def findLimits(x,y,limitValue,moreThan=True):
    """
    The method finds limits of the regions in which the condition is met.
    It is assumed, that the regions extend from the right end to some value
    and from the left end to some value.
    """
    if moreThan:
        ww = np.where(y>limitValue)[0]
    else:
        ww = np.where(y<limitValue)[0]
    rightEdge = x[np.max(ww)]
    leftEdge = x[np.min(ww)]
    return leftEdge,rightEdge

def reduceTicksFrequency(yaxis,rate,offset):
    every_nth = rate
    for n, label in enumerate(yaxis.get_ticklabels()):
        if n % every_nth == offset:
            label.set_visible(False)


def GenerateFig4(resetDumpFile=False):
    #logging.basicConfig(level=logging.INFO)
    logging.info("In %s"%(sys._getframe().f_code.co_name))
    binDumpFileName = "Fig4CalcDump.bin"
    if resetDumpFile and os.path.isfile(binDumpFileName):
        os.unlink(binDumpFileName)

    nCalcPoints = 200
    stepShiftLinspace = np.linspace(0,3.0,nCalcPoints)
    if os.path.isfile(binDumpFileName):
        values = np.fromfile(binDumpFileName)
        Sav = values[0:nCalcPoints]
        Ssd = values[nCalcPoints:2*nCalcPoints]
        Dav = values[2*nCalcPoints:3*nCalcPoints]
        Dsd = values[3*nCalcPoints:4*nCalcPoints]
        Mav = values[4*nCalcPoints:5*nCalcPoints]
        Msd = values[5*nCalcPoints:6*nCalcPoints]

        #pdb.set_trace()
    else:
        npoints=4000
        firstAfterWarmup=25
        p_I1 = 3.0
        p_I1 = 1.5
        respBPM =20.0
        kickAmplitude = -0.017
        basalValues = [0.0,0.0,0.0]
        basalValues0 = [0.0,0.0,0.0]

        basalValues0,values,_ = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,0.0,p_I1,respBPM,basalValues,np.linspace(0.0,0.0,1)))


        logging.info("============================= kickAmplitude %lf normed by %s ========================= " %(kickAmplitude,basalValues0))
        basalValues,values,setOfModelObjects = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues0,stepShiftLinspace))
        np.array(values).tofile(binDumpFileName)
        (Sav,Ssd,Dav,Dsd,Mav,Msd) = values

    fig,ax = plt.subplots()
    ax.ticklabel_format(useOffset = False)
    ax.plot(stepShiftLinspace,Ssd,"k",linestyle="-",linewidth=1) # SAP thin curve
    ax.plot(stepShiftLinspace,Dsd,"k",linestyle="--",linewidth=3) # DAP thick dashed
    ax.plot(stepShiftLinspace,Msd,"k",linestyle="-",linewidth=3) # MAP thick solid
    for y in (80.0,100.0,120.0):
        plt.axhline(y=y,color='k',linestyle='-')
    #ax.fill(30, 30, fill=False, hatch='\\')
    #pdb.set_trace()
    plt.ylabel("Magnitude of fluctuations [%]")
    plt.xlabel(r"Delay time $\tau$ [s]")
    ylim = copy.deepcopy(ax.get_ylim())# keep original limits.\
    xlim = copy.deepcopy(ax.get_xlim())# keep original limits.\
    leftAreaRightEdge,rightAreaLeftEdge = findOuterLimits(stepShiftLinspace,Msd,100,0.65)
    ax.fill([0,0,leftAreaRightEdge,leftAreaRightEdge,0],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
    ax.fill([rightAreaLeftEdge,rightAreaLeftEdge,xlim[1],xlim[1],rightAreaLeftEdge],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
    logging.info("Edges detected at: %lf and %lf" % (leftAreaRightEdge,rightAreaLeftEdge))
    #ax.bar(range(0,2), np.arange(ylim[0],ylim[1]), color='red', edgecolor='black', hatch="/")
    ax.set_ylim([ylim[0],ylim[1]])
    ax.set_xlim([xlim[0],xlim[1]])
    fname = "Fig4.png"
    plt.savefig(fname)
    logging.info("Writing to: %s" % (fname))
    #plt.show()

def GenerateFig5(resetDumpFile=False):
    logging.info("In %s"%(sys._getframe().f_code.co_name))
    binDumpFileName = "Fig5CalcDump.bin"
    if resetDumpFile and os.path.isfile(binDumpFileName):
        os.unlink(binDumpFileName)
    nCalcPoints = 400
    #KickAmplitudeLinspace = np.linspace(0.0,-0.08,nCalcPoints) to dla opóźnienia 2.0 to pierwsza opcja.
    leftKickLimit = -0.8
    rightKickLimit = 0.8
    KickAmplitudeLinspace = np.linspace(rightKickLimit,leftKickLimit,nCalcPoints)
    if os.path.isfile(binDumpFileName):
        values = np.fromfile(binDumpFileName)
        Sav = values[0:nCalcPoints]
        Ssd = values[nCalcPoints:2*nCalcPoints]
        Dav = values[2*nCalcPoints:3*nCalcPoints]
        Dsd = values[3*nCalcPoints:4*nCalcPoints]
        Mav = values[4*nCalcPoints:5*nCalcPoints]
        Msd = values[5*nCalcPoints:6*nCalcPoints]
        #pdb.set_trace()
    else:
        npoints=4000
        npoints=5000
        firstAfterWarmup=25
        firstAfterWarmup=35
        p_I1 = 3.0
        p_I1 = 1.5
        respBPM =20.0
        forceDelayTau = 1.5 # 2.0
        basalValues = [0.0,0.0,0.0]
        basalValues0 = [0.0,0.0,0.0]

        KickAmplitudeLinspace0 = np.linspace(0.0,0.0,1)
        forceDecayTau=0.3 # dummy as there is no input
        basalValues0,values,_ = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues,KickAmplitudeLinspace0))

        logging.info("============================= forceDecayTau %lf normed by %s ========================= " %(forceDecayTau,basalValues0))
        basalValues,values,setOfModelObjects = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues0,KickAmplitudeLinspace))

        np.array(values).tofile(binDumpFileName)
        (Sav,Ssd,Dav,Dsd,Mav,Msd) = values
    #pdb.set_trace()
    fig,ax = plt.subplots()
    ax.ticklabel_format(useOffset = False)
    ax.plot(KickAmplitudeLinspace,Ssd,"k",linestyle="-",linewidth=1) # SAP thin curve
    ax.plot(KickAmplitudeLinspace,Dsd,"k",linestyle="--",linewidth=3) # DAP thick dashed
    ax.plot(KickAmplitudeLinspace,Msd,"k",linestyle="-",linewidth=3) # MAP thick solid
    for y in (80.0,100.0,120.0):
        plt.axhline(y=y,color='k',linestyle='-')
    plt.axvline(x=0,color='k',linestyle='-')
    ax.set_xlim([leftKickLimit,rightKickLimit])
    #pdb.set_trace()
    plt.ylabel("Magnitude of fluctuations [%]")
    plt.xlabel(r"Kick amplitude $r_0$ [1/s]")

    ylim = copy.deepcopy(ax.get_ylim())# keep original limits.\
    xlim = copy.deepcopy(ax.get_xlim())# keep original limits.\
    leftEdge,rightEdge = findLimits(KickAmplitudeLinspace,Dsd,100,False)
    logging.info("Edges detected at: %lf and %lf" % (leftEdge,rightEdge))
    ax.fill([leftEdge,leftEdge,rightEdge,rightEdge,leftEdge],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
    #ax.bar(range(0,2), np.arange(ylim[0],ylim[1]), color='red', edgecolor='black', hatch="/")
    ax.set_ylim([ylim[0],ylim[1]])
    ax.set_xlim([xlim[0],xlim[1]])
    fname = "Fig5.png"
    plt.savefig(fname)
    logging.info("Writing to: %s" % (fname))
    #plt.show()

def GenerateFig3():
    npoints=2500
    firstAfterWarmup=20
    p_I1 = 3.0
    p_I1 = 1.5
    respBPM =20.0
    kickAmplitude = -0.017
    basalValues = [0.0,0.0,0.0]
    basalValues0 = [0.0,0.0,0.0]
    kickAmplitude = -0.017
    basalValues0,values,setOfModelObjects = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues,np.linspace(0.0,0.0,1)))
    settings, model, force, iafResp, iafHeart, collector, seriesNotifier, fireNotifierResp, fireNotifierHeart, fireNotifier = setOfModelObjects
    fname = "Fig3.png"
    allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
    allItems = collector.GetFiducialPointsList()
    fig = plt.figure()
    gs = gridspec.GridSpec(5,1,
                           height_ratios = [2,2,2,2,2]
                           )

    ax = plt.subplot(gs[0])
    plt.plot(allTimes,seriesNotifier.GetVar(0),'k-')
    plt.plot(allItems[1:,0],allItems[1:,1],"k^") # systolic
    plt.plot(allItems[1:,0],allItems[1:,2],"k^") # diastolic
    plt.plot(allItems[1:,0],allItems[1:,3],"k^") # mean
    plt.ylim(60.0,160.0)
    reduceTicksFrequency(ax.yaxis,2,1)
    plt.ylabel("BP [mmHg]")

    plt.setp(ax.get_xticklabels(), visible = False)

    ax = plt.subplot(gs[1])
    plt.plot(allTimes,seriesNotifier.GetVar(2),'k-',linewidth=2)
    plt.plot(allTimes,seriesNotifier.GetVar(1),'k-')
    #plt.plot(fireNotifier.firingTimes,fireNotifier.firingTimesSpikes(),"go")
    #plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"g",linewidth=2)
    #plt.yticks([])
    reduceTicksFrequency(ax.yaxis,2,1)
    plt.setp(ax.get_xticklabels(), visible = False)
    #plt.xlabel("Time [s]")
    plt.ylabel(r"$p_V,p_C$")
    #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)

    ax = plt.subplot(gs[2])
    plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForPhase),color='lightgray')
    plt.fill_between(allTimes,0,seriesNotifier.GetVar(iafResp.CoordinateNumberForPhase),facecolor='lightgray')
    plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForPhase),"k")
    plt.ylim(0.0,1.1)
    plt.yticks([0.0,1.0])
    #plt.xlabel("Time [s]")
    plt.ylabel(r"$\Phi(t),\varphi(t)$")
    #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
    plt.setp(ax.get_xticklabels(), visible = False)

    ax = plt.subplot(gs[3])
    plt.ylabel("RR [ms]")
    x = fireNotifierHeart.firingTimes[1:]
    y = fireNotifierHeart.ISI()[1:]*1000
    plt.plot(x,y,"k+-")
    reduceTicksFrequency(ax.yaxis,2,1)
    #plt.ylim(0.0,2.0)
    plt.setp(ax.get_xticklabels(), visible = False)

    ax = plt.subplot(gs[4])

    plt.plot(allTimes,seriesNotifier.GetVar(3),"k")
    plt.ylabel(r"$p_I(t)$ [mmHg]")
    plt.ylim(-7.0, -3.0)
    reduceTicksFrequency(ax.yaxis,2,1)
    plt.xlabel("Time [s]")
    logging.info("Test result in %s" % fname)
    plt.savefig(fname)


    # default visualization
#==============================================================================
#     allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
#     allItems = collector.GetFiducialPointsList()
#     KickedWindkesselModelVisualization(fname,
#                                            allTimes,
#                                            allItems,
#                                            seriesNotifier,
#                                            fireNotifierResp,
#                                            fireNotifier,
#                                            fireNotifierHeart,
#                                            iafResp.CoordinateNumberForPhase,
#                                            force.CoordinateNumber,
#                                            iafHeart.CoordinateNumberForRate,
#                                            iafHeart.CoordinateNumberForPhase)
#
#==============================================================================

def GenerateFig3A():
    """
    This is Fig3 for a case of largely reduced Msd, to understand what is going on.
    """
    npoints=4000
    firstAfterWarmup=30
    p_I1 = 3.0
    p_I1 = 1.5
    respBPM =20.0
    kickAmplitude = -0.3
    basalValues = [0.0,0.0,0.0]
    basalValues0 = [0.0,0.0,0.0]
    basalValues0,values,setOfModelObjects = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues,np.linspace(0.0,0.0,1)))
    settings, model, force, iafResp, iafHeart, collector, seriesNotifier, fireNotifierResp, fireNotifierHeart, fireNotifier = setOfModelObjects
    fname = "Fig3A.png"
    allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
    allItems = collector.GetFiducialPointsList()
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


def GenerateFig2():
    """
    Pleural pressure:
    inspration: -9 cm H20 = -6.6 mmHg
    expiration: -5 cm H20 = -3.7 mmHg
    1 cmH20 = 0.735 mmHg
    """
    npoints=1000
    firstAfterWarmup=3
    p_I1 = 1.5
    respBPM =20.0
    kickAmplitude = -0.3
    basalValues = [0.0,0.0,0.0]
    basalValues0 = [0.0,0.0,0.0]
    basalValues0,values,setOfModelObjects = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues,np.linspace(0.3,0.3,1)))
    settings, model, force, iafResp, iafHeart, collector, seriesNotifier, fireNotifierResp, fireNotifierHeart, fireNotifier = setOfModelObjects
    fname = "Fig2.png"
    allTimes = np.linspace(0.0,npoints*force.SamplingTime,npoints)
    allItems = collector.GetFiducialPointsList()

    fig = plt.figure()
    gs = gridspec.GridSpec(5,1,
                           height_ratios = [2,2,2,2,2]
                           )
    ax = plt.subplot(gs[0])
    plt.plot(allTimes,seriesNotifier.GetVar(iafResp.CoordinateNumberForPhase),"k")
    plt.plot(fireNotifierResp.firingTimes,fireNotifierResp.firingTimesSpikes(),"ko")
    plt.yticks([0.0, 1.0])
    plt.ylabel(r"$\Phi(t)$")
    #plt.xlabel("Time [s]")
    plt.ylim(0.0, 1.2)
    plt.setp(ax.get_xticklabels(), visible = False)


    ax = plt.subplot(gs[1])
    plt.plot(allTimes,seriesNotifier.GetVar(3),"k")
    plt.ylabel(r"$p_I(t)$ [mmHg]")
    plt.setp(ax.get_xticklabels(), visible = False)
    reduceTicksFrequency(ax.yaxis,2,1)

    ax = plt.subplot(gs[2])
    #todo: the guy below is all flat.
    #plt.plot(fireNotifier.firingTimes,fireNotifier.firingTimesSpikes(),"go")
    plt.plot(allTimes,seriesNotifier.GetVar(force.CoordinateNumber),"k",linewidth=1)
    #plt.yticks([])
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.yticks([-0.3,-0.2,-0.1,0.0])
#    reduceTicksFrequency(ax.yaxis,2,1)

    #plt.xlabel("Time [s]")
    plt.ylabel(r"$r_{n}(t)$")
    #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
    ax = plt.subplot(gs[3])
    #plt.yticks([])
    plt.ylabel(r"$r(t)$")
    plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForRate),"k",linewidth=1)
    plt.setp(ax.get_xticklabels(), visible = False)
    #reduceTicksFrequency(ax.yaxis,2,1)
    plt.yticks([-0.02,-0.01,0.0,0.01,0.02])

    ax = plt.subplot(gs[4])
    plt.plot(allTimes,seriesNotifier.GetVar(iafHeart.CoordinateNumberForPhase),"k")
    #plt.plot(allTimes,seriesNotifier.GetVar(2),"r")
    plt.plot(fireNotifierHeart.firingTimes,fireNotifierHeart.firingTimesSpikes(),"ko")
    plt.ylim(0.0, 1.2)
    plt.yticks([0.0, 1.0])
    #plt.xlabel("Time [s]")
    plt.ylabel(r"$\varphi(t)$")
    #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
    plt.setp(ax.get_xticklabels(), visible = False)
    #plt.ylim(100.0,300.0)
    plt.xlabel("Time [s]")

    logging.info("Test result in %s" % fname)
    plt.savefig(fname)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    GenerateFig2()
    GenerateFig3()
    GenerateFig3A()
    GenerateFig4()
    GenerateFig5()