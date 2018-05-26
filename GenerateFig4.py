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
from KickedWindkesselProcessors import StandardModelSetup,PhaseShiftProcessor,KickAmplitudeProcessor

def GenerateFig4():
    #logging.basicConfig(level=logging.INFO)
    npoints=4000
    firstAfterWarmup=25
    p_I1 = 3.0
    respBPM =20.0
    kickAmplitude = -0.17
    basalValues = [0.0,0.0,0.0]
    basalValues0 = [0.0,0.0,0.0]
    
#            fname = sys._getframe().f_code.co_name + "_delay_%f.png" % (stepShift)
#            KickedWindkesselModelVisualization(fname,
#                                       allTimes,
#                                       allItems,
#                                       seriesNotifier,
#                                       fireNotifierResp,
#                                       fireNotifier,
#                                       fireNotifierHeart,
#                                       iafResp.CoordinateNumberForPhase,
#                                       force.CoordinateNumber,
#                                       iafHeart.CoordinateNumberForRate,
#                                       iafHeart.CoordinateNumberForPhase)    
    basalValues0,values,_ = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,0.0,p_I1,respBPM,basalValues,np.linspace(0.0,0.0,1)))
        
    stepShiftLinspace = np.linspace(0,3.0,10)    
    print("============================= kickAmplitude %lf normed by %s ========================= " %(kickAmplitude,basalValues0))
    basalValues,values,setOfModelObjects = copy.deepcopy(PhaseShiftProcessor(npoints,firstAfterWarmup,kickAmplitude,p_I1,respBPM,basalValues0,stepShiftLinspace))    
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
    leftAreaRightEdge = 0.88
    rightAreaLeftEdge = 2.88
    
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

def GenerateFig5():
    pass

if __name__ == "__main__":
    GenerateFig4()
    GenerateFig5()