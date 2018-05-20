# -*- coding: utf-8 -*-
"""
Created on Sun May 20 21:30:54 2018

@author: teos
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import logging

def KickedWindkesselModelVisualization(fname,
                                       allTimes,
                                       allItems,
                                       seriesNotifier,
                                       fireNotifierResp,
                                       fireNotifier,
                                       fireNotifierHeart,
                                       iafRespCoordinateNumberForPhase,
                                       forceCoordinateNumber,
                                       iafHeartCoordinateNumberForRate,
                                       iafHeartCoordinateNumberForPhase):
        #print("Firing times: %s" % str(fireNotifier.firingTimes))
        fig = plt.figure()
        gs = gridspec.GridSpec(7,1,
                               height_ratios = [1,2,2,2,2,2,4]
                               )
        ax = plt.subplot(gs[0])
        plt.plot(allTimes,seriesNotifier.GetVar(iafRespCoordinateNumberForPhase),"b")
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
        plt.plot(allTimes,seriesNotifier.GetVar(forceCoordinateNumber),"g",linewidth=2)      
        plt.yticks([])
        plt.setp(ax.get_xticklabels(), visible = False)        
        #plt.xlabel("Time [s]")
        plt.ylabel(r"$r_{n}(t)$")
        #plt.ylim(0.0,5.0 * iafResp.KickAmplitude)
        ax = plt.subplot(gs[3])
        plt.yticks([])
        plt.ylabel(r"$r(t)$")                
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeartCoordinateNumberForRate),"v")
        plt.setp(ax.get_xticklabels(), visible = False)
        
        ax = plt.subplot(gs[4])
        plt.plot(allTimes,seriesNotifier.GetVar(iafHeartCoordinateNumberForPhase),"r")
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
        x = fireNotifierHeart.firingTimes[1:]
        y = fireNotifierHeart.ISI()[1:]                        
        logging.debug(x)
        logging.debug(y)
        plt.plot(x,y,"b+-")
        plt.ylim(0.0,2.0)
        plt.setp(ax.get_xticklabels(), visible = False)
        
        plt.subplot(gs[6])        
        plt.plot(allTimes,seriesNotifier.GetVar(0),'b-')        
        plt.plot(allItems[:,0],allItems[:,1],"ro") # systolic
        plt.plot(allItems[:,0],allItems[:,2],"go") # diastolic
        plt.plot(allItems[:,0],allItems[:,3],"bo") # mean
        
        plt.ylim(100.0,300.0)
        plt.xlabel("Time [s]")
        plt.ylabel("BP [mmHg]")        
        logging.info("Test result in %s" % fname)
        plt.savefig(fname)   