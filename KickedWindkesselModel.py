# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import logging
import warnings
from enum import Enum
import numpy as np
#from scipy.integrate import ode
import matplotlib.pyplot as plt
import pdb
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorParams, RungeKutta45ConstStepIntegrator, RungeKutta45IntegratorData
from Notifiers import SeriesNotifier, StatsNotifier, NotifierChain, NilNotify
from HeartActionForce import RectangularHeartActionForce


""" 
  C_a{dp_a}/{dt}=\sum_{i} \delta(t-t_i)Z_{ca}(p_c-p_I)	- Z_{av}p_{av}
  C_c{dp_c}/{dt}=max(Z_{vc}p_{vc},0) - \sum_{i} \delta(t-t_i)Z_{ca}(p_c-p_I)	
  C_v{dp_v}/{dt}= Z_{av}p_{av} - max(Z_{vc}p_{vc},0)
  p_I(t)=p_{I0}+p_{I1}(1+cos\:2\pi\Phi(t))  
\end{equation}	
\begin{equation}	
	d \phi/dt=1/T_0+\sum_{j} \delta(t-t_j-\tau)F(\tilde{\varphi})/T_{RSA}
\end{equation}	
\begin{equation}	
	d \Phi/dt=1/T_{R0}
\end{equation}
\begin{equation}	
	t_i:\varphi(t_i)=i, \tilde{\varphi} = \varphi\:mod\:1
\end{equation}	
\begin{equation}	
	t_j:\Phi(t_j)=j
\end{equation}
\begin{equation}	
	F(\varphi)=\varphi^{1.3}(\varphi-0.45)\frac{(1-\varphi)^3}{(1-0.8)^3+(1-\varphi)^3}
\end{equation}	
     """                   
    
def kickedWindkesselRHS(t,y,Freturn,settings):
    """
    RHS of kicked Windkessel model.
    It still relies on external parameters and changes their value:
    settings.heartFlowBeginTime.
    Should be a non-autonomous part of model variables
    As both phases should.    
    """

    #p_a, p_c, p_v, p_I, startKick = y
    p_a = y[0]
    p_c = y[1]
    p_v = y[2] 
    p_I = y[3] 
    startKick = y[4]
    # If the Windkessel is kicked start the heart flow
    #the heart flow is active in a rectangle:
    #from heartFlowBeginTime for heartFlowTimespan
    if startKick > 0.0:     
        settings.heartFlowBeginTime = t
    # Positive heartFlowBeginTime = open flow                        
    heartFlowIsOpen = (settings.heartFlowBeginTime >= 0.0)
    # Check closing criteria.
    if heartFlowIsOpen and (t > (settings.heartFlowTimespan + settings.heartFlowBeginTime)):
        # by convention: negative heartFlowBeginTime means, that
        #there is no heart flow.
        settings.heartFlowBeginTime = -1.0
        heartFlowIsOpen = False
    
    #//Fs = 125; time step = 1/125 s
    q_ca = settings.Zca * (p_c - p_I)
    if not heartFlowIsOpen: # all or nothing
        q_ca = 0.0
    q_vc = np.max((settings.Zvc * (p_v - p_c), 0.0))
    q_av = settings.Zav * (p_a - p_v)
    #print("Q:%f\t%f\t%f"%(q_ca,q_vc,q_av))
    Freturn[0] = (q_ca - q_av) / settings.Ca#; //p_a : //C_a{dp_a}/{dt}=\sum_{i} \delta(t-t_i)Z_{ca}(p_c-p_I)	- Z_{av}p_{av}
    Freturn[1] = (q_vc - q_ca) / settings.Cc#;//p_c : //  C_c{dp_c}/{dt}=max(Z_{vc}p_{vc},0) - \sum_{i} \delta(t-t_i)Z_{ca}(p_c-p_I)	
    Freturn[2] = (q_av - q_vc) / settings.Cv#;//p_v // //  C_v{dp_v}/{dt}= Z_{av}p_{av} - max(Z_{vc}p_{vc},0)

    #to nie moze dzialac bo przeciez tu sie nie liczy wartosci zmiennej tylko RHS rownania.            
    #Freturn[3] = settings.p_I0 + settings.p_I1 * (1 + np.cos(settings.breathingPhi0 + 2 * np.pi * t / settings.breathingPeriod)) - p_I#;//p_I : //  p_I(t)=p_{I0}+p_{I1}(1+cos\:2\pi\Phi(t))                          
    #print("t:%f\tF:%f\t%f\t%f\t%f"%(t,Freturn[0],Freturn[1],Freturn[2],Freturn[3]))
    #print("p_a %lf p_c %lf p_v %lf p_I %lf" % (p_a, p_c, p_v, p_I))
    #print("q_ca %lf q_vc %lf q_av %lf" % (q_ca,q_vc,q_av))            
    return Freturn

    
class KickedWindkesselModel:
        class ModelVariable(Enum):
            P_a = 1
            P_c = 2
            P_v = 3
            P_I = 4
            Heart_Drive = 5 #         

        class KickedWindkesselModelSettings:
 
            def __init__(self):            
                #bp
                #Rav=0.9;%th(16);%1
                #Rbv=0.005;%0.0005;%0.01;%0.01;%0.01;%th(20);%0.01
                #Rba=0.006;%0.0029;%0.0029;%0.003;%th(15);%0.006

                #Cb=th(7);%4.3
                #Ca=th(3);%1.6
                #Cv=th(4);%100
                self.Ca = 1.6
                self.Cc = 4.3
                self.Cv = 100.0
                self.Zca = 1 / 0.006 # 166.[6]
                self.Zav = 1 / 0.9
                self.Zvc = 1 / 0.005
                #breathing
                self.p_I0 = -4.0
                self.p_I1 = 2.0
                self.breathingPeriod = 3.0
                #self.breathingPhi0 = 0.0
                self.heartPhase = 0.0
                #self.heartActionForce = NullDrivingForce()
               
                self.heartFlowTimespan = 0.01#;//0.1 sec                
                self.heartFlowBeginTime = -1
                self.CoordinateNumberForRespPhase = -1
        
        @property
        def ModelDimension(self):
                return 5
                
        def __init__(self,settings,dimension=0):                    
            self.settings = settings#;// new KickedWindkesselModelSettings(settings);
            self.param = RungeKutta45IntegratorParams()
            if dimension == 0:
                self.param.dimension = self.ModelDimension
            else:
                self.param.dimension = dimension #may be extended to andle extra data.
            self.param.dT = 0.01#0.01#;//1/100 or 1/125 //ol = 2 * 1E-7;                    //error control tolerance
            self.param.Tmin = 0.0#                         //startpoint
            self.param.Npoints = 10000#
            self.data = RungeKutta45IntegratorData(self.param.dimension,self.param.dT)
            #//initialize
            self.data[0] = 90.0#90.0#;//P(2) = 90;%Pa 90
            self.data[1] = 35.0#;//P(1) = 35;%25;%95;%Pb
            self.data[2] = 45.0#;//P(3) = 45;%35;%70;%Pv % bylo 30 i zle asymtotyczni    
            self.data[3] = settings.p_I0#; 
            self.data[4] = 0.0
            
            self.integrator = RungeKutta45ConstStepIntegrator(self.param,self.data,kickedWindkesselRHS)#,lambda t,y: f(t,y,self.settings),self.data)
                 

        def IterateToNotifiers(self):
            """
            Iterate to notifiers is a main interation loop.
            In this main loop all non-autonomous parameters are set 
            to prevent multiple execution of code in kickedWindkesselRHS
            which is called several times during Runge-Kutta integration.
            The heart action force sets 4-th variable to a nonzero value
            if a kick has to be initiated.
            """
            self.integrator.Reset(self.param)
            self.Notify(self.data)
            logging.info("Iterating from T0=0.0 to Tmax=%f, dT=%f, npoints=%d" % (self.integrator.Tmax(),self.integrator.param.dT,self.integrator.param.Npoints))
            while True:
                self.settings.heartActionForce.ApplyDrive(self.data)
                if self.settings.CoordinateNumberForRespPhase == -1:
                    self.data[3] = self.settings.p_I0 + self.settings.p_I1 * (1 + np.cos(2 * np.pi * self.data.t / self.settings.breathingPeriod))#;//p_I : //  p_I(t)=p_{I0}+p_{I1}(1+cos\:2\pi\Phi(t))
                else:
                    self.data[3] = self.settings.p_I0 + self.settings.p_I1 * (1 + np.cos(2 * np.pi * self.data[self.settings.CoordinateNumberForRespPhase]))#;//p_I : //  p_I(t)=p_{I0}+p_{I1}(1+cos\:2\pi\Phi(t))
                #self.settings.sineOfBreathingPhase = np.sin(2 * np.pi * self.data.t / self.settings.breathingPeriod)
                
                if not self.integrator.Iterate(self.settings):
                    logging.info("Iterations completed.")
                    break
                logging.debug(self.data)
                self.Notify(self.data)
                

def NotifyPlainPrint(data):
    print("NPP:"+str(data))


def DummyPrint(data):
    pass


def BasicProcessor():
    settings = KickedWindkesselModel.KickedWindkesselModelSettings()    
    settings.heartActionForce = RectangularHeartActionForce()
    model = KickedWindkesselModel(settings)
    model.param.Npoints = 1000
    series = SeriesNotifier(model.param.dimension,model.param.Npoints+1)
    model.Notify = series.Notify
    model.IterateToNotifiers()
    fname = "overview.png"
    series.PlotSeries((1,2,3,4),fname)
    print("Plat data in %s" % fname)

if __name__ == "__main__":
    BasicProcessor()

