# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 21:13:26 2017

@author: teos
"""
import numpy as np
from scipy.integrate import ode
import pdb
from RungeKutta45ConstStepIntegrator import RungeKutta45IntegratorParams, RungeKutta45ConstStepIntegrator, RungeKutta45IntegratorData

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
class HeartActionForce:
    def __init__(self):
        self.Drive = 0.0
        self.CoordinateNumber = 4
        self.Amplitude = 1.0
        self.StepWidth = 0.01
        self.StepPeriod = 1.0
        self.LastStepStart = 0.0
   
    def ApplyDrive(self,data):
        isOpen = self.Drive == self.Amplitude
        timeFromLastStepStart = data.t - self.LastStepStart
        if isOpen: # check if it is time to close
            if timeFromLastStepStart > self.StepWidth:
                self.Drive = 0.0
        else: #check if it is time to open
            if timeFromLastStepStart == 0.0 or timeFromLastStepStart > self.StepPeriod: # initial kick
                self.Drive = self.Amplitude
                self.LastStepStart = data.t       
        data[self.CoordinateNumber] = self.Drive
"""        


class RungeKutta45ConstStepIntegrator:

    def __init__(self,param,function,data):
        self.Force = 0.0
        self.CoordinateNumber = 0
        self.param = param
        self.function = function
        self.endTime = self.param.Tmin  + self.param.Tdelta * self.param.Npoints
        self.data = data        
        
    def Reset(self,param):        
        self.param = param
        self.endTime = self.param.Tmin  + self.param.Tdelta * self.param.Npoints        
        # create explicit Runge-Kutta integrator of order (4)5
        self.r = ode(self.function).set_integrator('dopri5')
        #pdb.set_trace()
        self.r.set_initial_value(self.data.y, self.param.Tmin)
        

    def Iterate(self):
        if not self.r.successful():
            return False
        self.data.t = self.data.t + self.param.Tdelta
        print("Przed %s" % self.data.y)
        self.data.y = self.r.integrate(self.data.t)
        print("Po %s" % self.data.y)
        return self.data.t < self.endTime
"""        
               
    
        # Place your FUNCTION HERE
def f(t,y,Freturn,settings):
            """definition of equation-example
            damped harmonic oscillator with harmonic driver
            x" + 100*x + 2*x'= 10*sin 3*t
            we define y[0] = x, y[1] = x'
            You need to transform the second order differential equation
            in a system of two differenti5al equatios of first order:
            f[0] = x'= y[1]
            f[1] = x" = -100*x-2*x' = -100 *y[0] -2*y[1] + 10*sin 3*t
            You may enter your own equation here!"""
            #pdb.set_trace()
            #print(y)
            #Freturn = y
            #print(len(y))
            p_a, p_c, p_v, p_I,dummy = y
            
            if settings.openHeartFlow > 0.0: # //open the heart flow            
                if settings.heartFlowBeginTime == -1:                
                    settings.lastDiastolicBp = p_a #; //register diastolic on opening.                
                settings.heartFlowBeginTime = t
                settings.lastHeartBeatTime = t
                                

            heartFlowIsOpen = (settings.heartFlowBeginTime >= 0.0)#;//ok, it holds the real time: the flow is open.
            #//check if it should not be closed
            openTime = t - settings.heartFlowBeginTime
            if heartFlowIsOpen and (openTime > settings.heartFlowTimespan):
                settings.heartFlowBeginTime = -1.0#; //self parameter is set to current time when flow opens. When the flow is closed, it is set to a negative value.
                #//register systolic on closing.
                #//settings.lastSystolicBp = p_a; //close the heart: register systolic now
                heartFlowIsOpen = False
            
            #//Fs = 125; time step = 1/125 s
            q_ca = settings.Zca * (p_c - p_I)
            if not heartFlowIsOpen:
                q_ca = 0.0
            q_vc = np.max(settings.Zvc * (p_v - p_c), 0.0)
            q_av = settings.Zav * (p_a - p_v)
            #print("Q:%f\t%f\t%f"%(q_ca,q_vc,q_av))
            Freturn[0] = (q_ca - q_av) / settings.Ca#; //p_a : //C_a{dp_a}/{dt}=\sum_{i} \delta(t-t_i)Z_{ca}(p_c-p_I)	- Z_{av}p_{av}
            Freturn[1] = (q_vc - q_ca) / settings.Cc#;//p_c : //  C_c{dp_c}/{dt}=max(Z_{vc}p_{vc},0) - \sum_{i} \delta(t-t_i)Z_{ca}(p_c-p_I)	
            Freturn[2] = (q_av - q_vc) / settings.Cv#;//p_v // //  C_v{dp_v}/{dt}= Z_{av}p_{av} - max(Z_{vc}p_{vc},0)
            #//non-autonomous
            Freturn[3] = settings.p_I0 + settings.p_I1 * (1 + np.cos(2 * np.pi * t / settings.breathingPeriod)) - p_I#;//p_I : //  p_I(t)=p_{I0}+p_{I1}(1+cos\:2\pi\Phi(t))                          
            print("t:%f\tF:%f\t%f\t%f\t%f"%(t,Freturn[0],Freturn[1],Freturn[2],Freturn[3]))
            #y = Freturn
            #print("p_a %lf p_c %lf p_v %lf p_I %lf" % (p_a, p_c, p_v, p_I))
            #print("q_ca %lf q_vc %lf q_av %lf" % (q_ca,q_vc,q_av))
            # pdb.set_trace()
            
            return Freturn
        #//p_a = y[0]
        #//p_c = y[1]
        #//p_v = y[2]
        #//p_I = y[3]
        #//F(\varphi)=\varphi^{1.3}(\varphi-0.45)\frac{(1-\varphi)^3}{(1-0.8)^3+(1-\varphi)^3}
    
class KickedWindkesselModel:
        class KickedWindkesselModelSettings:
 
            def __init__(self):            
                self.lastHeartBeatTime = 0.0;
                self.lastDiastolicBp = 0.0;
                self.lastSystolicBp = 0.0;
            
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
                self.p_I1 = 0.1 #2.0
                self.breathingPeriod = 3.0
                self.heartPhase = 0.0
                #self.heartActionForce = NullDrivingForce()
               
                self.heartFlowTimespan = 0.01#;//0.1 sec
                self.throwAmplitudeDeathException = False
                self.heartFlowBeginTime = -1
                self.openHeartFlow = 1.0
            
            """def __init__(self,settings):   
                self.lastHeartBeatTime = 0.0;
                self.lastDiastolicBp = 0.0;
                self.lastSystolicBp = 0.0;
            
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
                self.Zca = 1 / 0.006
                self.Zav = 1 / 0.9
                self.Zvc = 1 / 0.005
                #breathing
                self.p_I0 = -4.0
                self.p_I1 = 0.1
                self.breathingPeriod = 3.0
                self.heartPhase = 0.0
              
                self.heartFlowTimespan = 0.01#;//0.1 sec
                self.throwAmplitudeDeathException = False
                self.heartFlowBeginTime = -1
                self.openHeartFlow = 0.0
                
                self.breathingPeriod = settings.breathingPeriod
                self.Ca = settings.Ca
                self.Cc = settings.Cc
                self.Cv = settings.Cv
                self.heartPhase = settings.heartPhase
                self.openHeartFlow = settings.openHeartFlow
                self.p_I0 = settings.p_I0
                self.p_I1 = settings.p_I1
                self.Zav = settings.Zav
                self.Zca = settings.Zca
                self.Zvc = settings.Zvc
                self.throwAmplitudeDeathException = settings.throwAmplitudeDeathException
                self.heartFlowBeginTime = settings.heartFlowBeginTime
                self.heartFlowTimespan = settings.heartFlowTimespan
                self.heartActionForce = settings.heartActionForce"""

            

        def phaseEfectivenessCurve(phi):        
            aux1 = np.power(1 - phi, 3.0)
            value = np.power(phi, 1.3) * (phi - 0.45) * (aux1 / ((np.power(1 - 0.8, 3.0)) + aux1));
            return value
        
        @property
        def ModelDimension(self):
                return 5
                
        def __init__(self,settings):                    
            self.settings = settings#;// new KickedWindkesselModelSettings(settings);
            self.param = RungeKutta45IntegratorParams()
            self.param.dimension = self.ModelDimension
            self.param.Tdelta = 0.01#0.01#;//1/100 or 1/125 //ol = 2 * 1E-7;                    //error control tolerance
            self.param.Tmin = 0.0#                         //startpoint
            self.param.Npoints = 10000#
            self.data = RungeKutta45IntegratorData(self.param)
            #//initialize
            self.data[0] = 90.0#90.0#;//P(2) = 90;%Pa 90
            self.data[1] = 35.0#;//P(1) = 35;%25;%95;%Pb
            self.data[2] = 45.0#;//P(3) = 45;%35;%70;%Pv % bylo 30 i zle asymtotyczni    
            self.data[3] = settings.p_I0#;
            self.data[4] = 0#;
            
            self.integrator = RungeKutta45ConstStepIntegrator(self.param,self.data,f)#,lambda t,y: f(t,y,self.settings),self.data)
                 

        def IterateToNotifiers(self):

            self.integrator.Reset(self.param)
            self.Notify(self.data)
            while True:            
                self.settings.heartActionForce.ApplyDrive(self.data)                
                self.settings.openHeartFlow = self.settings.heartActionForce.Drive
                isHeartOpen = (self.settings.heartFlowBeginTime != -1.0)
                self.settings.sineOfBreathingPhase = np.sin(2 * np.pi * self.data.t / self.settings.breathingPeriod)
                if isHeartOpen:                
                    abc = 1
                
                if self.integrator.Iterate(self.settings):
                    print("yawn")
                    break
#                for i in range(len(self.data.y)):
#                    if self.data.y[i] < 0:
#                        self.data.y[i] = 0
#                isHeartClosedAfter = (self.settings.heartFlowBeginTime == -1.0)
#                if  isHeartOpen and isHeartClosedAfter: #//it has just closed.                
#                    self.settings.lastSystolicBp = self.data[0]                
#                if self.settings.throwAmplitudeDeathException:
#                    if self.data[0] < 0.0:
#                        raise Exception("Amplitude death occurred")
                self.Notify(self.data)

def NotifyPlainPrint(data):
    print("NPP:"+str(data))

def DummyPrint(data):
    pass
if __name__ == "__main__":
    settings = KickedWindkesselModel.KickedWindkesselModelSettings()    
    settings.heartActionForce = HeartActionForce()
    model = KickedWindkesselModel(settings)
    model.param.Npoints = 100
    model.Notify = NotifyPlainPrint
    #model.Notify = DummyPrint
    model.IterateToNotifiers()

