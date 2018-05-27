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
import numpy as np
import matplotlib.pyplot as plt
from KickedWindkesselModel import KickedWindkesselModel
from HeartActionForce import RectangularHeartActionForce,RespiratoryDelayedSmearedHeartActionForce, HeartActionForceChain
from AbpmFiducialPointsCollector import AbpmFiducialPointsCollector
from Notifiers import SeriesNotifier, MinMaxNotifier, NotifierChain, FiringTimesNotifier
from IntegrateAndFire import IntegrateAndFire, phaseEfectivenessCurveSH
from KickedWindkesselModelVisualization import KickedWindkesselModelVisualization
from KickedWindkesselProcessors import StandardModelSetup,PhaseShiftProcessor,KickAmplitudeProcessor

def GenerateFig4():
    #logging.basicConfig(level=logging.INFO)
    binDumpFileName = "Fig4CalcDump.bin"
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
        respBPM =20.0
        kickAmplitude = -0.017
        basalValues = [0.0,0.0,0.0]
        basalValues0 = [0.0,0.0,0.0]
   
        basalValues0,values,_ = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,0.0,p_I1,respBPM,basalValues,np.linspace(0.0,0.0,1)))
            
        
        print("============================= kickAmplitude %lf normed by %s ========================= " %(kickAmplitude,basalValues0))
        basalValues,values,setOfModelObjects = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues0,stepShiftLinspace))
        np.array(values).tofile(binDumpFileName)
        (Sav,Ssd,Dav,Dsd,Mav,Msd) = values
    #pdb.set_trace()        
    fig,ax = plt.subplots()
    ax.ticklabel_format(useOffset = False)
    #plt.subplot(3,1,1)    
    ax.plot(stepShiftLinspace,Ssd,"k",linestyle="-",linewidth=1) # SAP thin curve    
    #plt.subplot(3,1,2)
    ax.plot(stepShiftLinspace,Dsd,"k",linestyle="--",linewidth=3) # DAP thick dashed
    #plt.subplot(3,1,3)
    ax.plot(stepShiftLinspace,Msd,"k",linestyle="-",linewidth=3) # MAP thick solid
    for y in (80.0,100.0,120.0):
        plt.axhline(y=y,color='k',linestyle='-')    
    #ax.fill(30, 30, fill=False, hatch='\\')
    #pdb.set_trace()
    plt.ylabel("Magnitude of fluctuations (%)")
    plt.xlabel("Delay time (s)")    
    ylim = copy.deepcopy(ax.get_ylim())# keep original limits.\
    xlim = copy.deepcopy(ax.get_xlim())# keep original limits.\    
    leftAreaRightEdge = 0.9
    rightAreaLeftEdge = 2.23
    
#==============================================================================
    ax.fill([0,0,leftAreaRightEdge,leftAreaRightEdge,0],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
#     offsetOfHalf = len(Msd)/2
#     rightHalfMsd = np.array(Msd[offsetOfHalf:])
#     print(rightHalfMsd)
#     rightLargerThan100 = np.where(rightHalfMsd>100)[0][0]
#     rightHatchedEdge = stepShiftLinspace[offsetOfHalf+rightLargerThan100]
#     print(rightHatchedEdge)
    ax.fill([rightAreaLeftEdge,rightAreaLeftEdge,xlim[1],xlim[1],rightAreaLeftEdge],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
# 
#==============================================================================
    #ax.bar(range(0,2), np.arange(ylim[0],ylim[1]), color='red', edgecolor='black', hatch="/")
    ax.set_ylim([ylim[0],ylim[1]])    
    ax.set_xlim([xlim[0],xlim[1]])
    fname = "Fig4.png"
    plt.savefig(fname)
    print("Writing to: %s" % (fname))
    #plt.show()

def GenerateFig5(resetDumpFile=False):
    #logging.basicConfig(level=logging.INFO)
    binDumpFileName = "Fig5CalcDump.bin"
    if resetDumpFile and os.path.isfile(binDumpFileName):
        os.unlink(binDumpFileName)
    nCalcPoints = 200
    #KickAmplitudeLinspace = np.linspace(0.0,-0.08,nCalcPoints) to dla opóźnienia 2.0 to pierwsza opcja.
    KickAmplitudeLinspace = np.linspace(0.0,-0.8,nCalcPoints)
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
        respBPM =20.0        
        forceDelayTau = 1.5 # 2.0
        basalValues = [0.0,0.0,0.0]
        basalValues0 = [0.0,0.0,0.0]

        KickAmplitudeLinspace0 = np.linspace(0.0,0.0,1)
        forceDecayTau=0.3 # dummy as there is no input
        basalValues0,values,_ = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues,KickAmplitudeLinspace0))
        
        print("============================= forceDecayTau %lf normed by %s ========================= " %(forceDecayTau,basalValues0))        
        basalValues,values,setOfModelObjects = copy.deepcopy(KickAmplitudeProcessor(npoints,firstAfterWarmup,forceDelayTau,forceDecayTau,p_I1,respBPM,basalValues0,KickAmplitudeLinspace))
#################################################################        
        
        np.array(values).tofile(binDumpFileName)
        (Sav,Ssd,Dav,Dsd,Mav,Msd) = values
    #pdb.set_trace()        
    fig,ax = plt.subplots()
    ax.ticklabel_format(useOffset = False)
    #plt.subplot(3,1,1)    
    ax.plot(KickAmplitudeLinspace,Ssd,"k",linestyle="-",linewidth=1) # SAP thin curve    
    #plt.subplot(3,1,2)
    ax.plot(KickAmplitudeLinspace,Dsd,"k",linestyle="--",linewidth=3) # DAP thick dashed
    #plt.subplot(3,1,3)
    ax.plot(KickAmplitudeLinspace,Msd,"k",linestyle="-",linewidth=3) # MAP thick solid
    for y in (80.0,100.0,120.0):
        plt.axhline(y=y,color='k',linestyle='-')    
    #ax.fill(30, 30, fill=False, hatch='\\')
    #pdb.set_trace()
    plt.ylabel("Magnitude of fluctuations (%)")
    plt.xlabel("Kick amplitude (1/s)")

    ylim = copy.deepcopy(ax.get_ylim())# keep original limits.\
    xlim = copy.deepcopy(ax.get_xlim())# keep original limits.\    
    leftAreaRightEdge = 0.9
    rightAreaLeftEdge = 2.23
    
#==============================================================================
    #ax.fill([0,0,leftAreaRightEdge,leftAreaRightEdge,0],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
#     offsetOfHalf = len(Msd)/2
#     rightHalfMsd = np.array(Msd[offsetOfHalf:])
#     print(rightHalfMsd)
#     rightLargerThan100 = np.where(rightHalfMsd>100)[0][0]
#     rightHatchedEdge = stepShiftLinspace[offsetOfHalf+rightLargerThan100]
#     print(rightHatchedEdge)
    #ax.fill([rightAreaLeftEdge,rightAreaLeftEdge,xlim[1],xlim[1],rightAreaLeftEdge],[ylim[0],ylim[1],ylim[1],ylim[0],ylim[0]],color='black',linewidth=1,edgecolor='black',linestyle='solid',hatch='/',fill=False)
# 
#==============================================================================
    #ax.bar(range(0,2), np.arange(ylim[0],ylim[1]), color='red', edgecolor='black', hatch="/")
    ax.set_ylim([ylim[0],ylim[1]])    
    ax.set_xlim([xlim[0],xlim[1]])
    fname = "Fig5.png"
    plt.savefig(fname)
    print("Writing to: %s" % (fname))
    #plt.show()
    
def GenerateFig3():
    npoints=4000
    firstAfterWarmup=25
    p_I1 = 3.0
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
        
def GenerateFig3A():
    npoints=4000
    firstAfterWarmup=30
    p_I1 = 3.0
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

if __name__ == "__main__":
    GenerateFig3()
    GenerateFig3A()
    GenerateFig4()
    GenerateFig5()