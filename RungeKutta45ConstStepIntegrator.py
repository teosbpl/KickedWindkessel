# -*- coding: utf-8 -*-
"""
Created on Thu May  3 19:13:26 2018

@author: teos

Rk45 Ordinary Differential Equations Solver (ODES)
Solve the differential equation  using Runge-Kutta-Fehlberg Method
with variable step size (Adjust "Tol" for error tolerance!)
ref: J.Mathews-"Numerical Methods for Mathematics, Science and Engineering"
Prentice Hall, 1992
copywrite 2002
       Cristian Constantin Bordeianu, Romania and R Landau, Oregon State U
e-mail:cbord@warpnet.ro
 Translated to C# by T.Buchner 2015
e-mail: buchner@if.pw.edu.pl
Partially based on http://csharpcomputing.com/Tutorials/Lesson16.htm Copyright 2001-2006 Aleksey Nudelman
 Translated to Python by T.Buchner
e-mail: buchner@if.pw.edu.pl

"""
import numpy as np
import unittest

class RungeKutta45IntegratorParams:
    def __init__(self, orig=None):
        if orig is None:
            self.non_copy_constructor()
        else:
            self.copy_constructor(orig)
    def non_copy_constructor(self):       
        self.dimension = 0
        self.Tmin = 0.0
        self.Npoints = 0
        self.dT = 0.0

    def copy_constructor(self,src):
        self.dimension = src.dimension
        self.Tmin = src.Tmin
        self.dT = src.dT
        self.Npoints = src.Npoints

class  RungeKutta45IntegratorData:
            def __init__(self,param):                
                self.dimension = param.dimension                
                self.y =  np.zeros(self.dimension)#array([0.0 for x in range(self.dimension)])
                self.t = param.Tmin
            def __str__(self):
                return "%lf,%s"%(self.t,','.join([str(x) for x in self.y]))
            def __getitem__(self, index):
                result = self.y.__getitem__(index)
                return result          
            def __setitem__(self, index,value):
                result = self.y.__setitem__(index,value)
                return result   

class RungeKutta45ConstStepIntegrator:            

        def Dimension(self):
            return self.param.dimension

        def Tmax(self):
            return self.tmax

        def __init__(self,param,data,function):
            self.data = data
            self.Reset(param)
            self.f = function
            self.states = np.zeros((self.param.Npoints,self.param.dimension+1))# all + time
            self.statesIterator = 0
            
        #RungeKutta45IntegratorParams 
        def Reset(self,param):        
            self.param = RungeKutta45IntegratorParams(param)
            self.F = np.zeros(self.param.dimension)
            self.ydumb = np.zeros(self.param.dimension)
            self.k1 = np.zeros(self.param.dimension)
            self.k2 = np.zeros(self.param.dimension)
            self.k3 = np.zeros(self.param.dimension)
            self.k4 = np.zeros(self.param.dimension)
            self.dt = self.param.dT;
            self.tmax = self.param.Tmin + self.param.Npoints * self.dt
            
        def Iterate(self,settings=None):            
            return self.RK4(settings)

        def IterateSeries(self,settings=None):            
            while self.RK4(settings):
                self.states[self.statesIterator,0] = self.data.t
                for ii in range(self.param.dimension):                    
                    self.states[self.statesIterator,ii+1] = self.data.y[ii]
                self.statesIterator = self.statesIterator+1
            return self.states
                                    
        
        def RK4(self,settings):
            t = self.data.t
            maxi = self.param.dimension
            if t > self.tmax:
                return False
            #//printout
            if (t + self.dt) > self.tmax:
                self.dt = self.tmax - t # // the last step
            """/*
             * Based on 
             * http://www.particle.kth.se/~lindsey/JavaCourse/Book/Part1/Physics/Chapter04/RungeKutta4th.html
             * A further refinement of the Runge-Kutta approach uses derivatives 
             * at the beginning of the interval, at the end of the interval and 
             * twice at the midpoint.
             * Using the conventional k# variable names, we obtain the following 
             * increments in the variable xn with :
             * k1=f'(xn,tn)*dt;
             * k2=f'(xn+k1/2,tn+dt/2)*dt;
             * k3=f'(xn+k2/2,tn+dt/2)*dt;
             * k4=f'(xn+k3,tn+dt)*dt;
             * Then the new increment in variable x is calculated as a weighted 
             * average of these estimated increments with k2 and k3 , 
             * the two midpoint values, given double weights.
             * xn+1=nx+(k1+2k2+2k3+k4)/6                          
             */"""
            self.F = self.f(t, self.data.y,self.F,settings) #;  // evaluate both RHS's and return in F
            for i in range(maxi):
                self.k1[i] = self.dt * self.F[i]
                self.ydumb[i] = self.data.y[i] + self.k1[i] / 2.0

            self.F = self.f(t + self.dt / 2, self.ydumb,self.F,settings)
            for i in range(maxi):
                self.k2[i] = self.dt * self.F[i]
                self.ydumb[i] = self.data.y[i] + self.k2[i] / 2.0

            self.F = self.f(t + self.dt / 2, self.ydumb,self.F,settings)
            for i in range(maxi):
                self.k3[i] = self.dt * self.F[i]
                self.ydumb[i] = self.data.y[i] + self.k3[i]

            self.F = self.f(t + self.dt, self.ydumb,self.F,settings)
            for i in range(maxi):
                self.k4[i] = self.dt * self.F[i]

            #if (true)//((err[0] < param.Tol) || (err[1] < param.Tol) || (dt <= 2 * hmin))//accept the approximation
            #{
                #// xn+1=nx+(k1+2k2+2k3+k4)/6                          
            for i in range(maxi):
                self.data.y[i] = self.data.y[i] + (self.k1[i] + self.k4[i]) / 6.0 + (self.k2[i] + self.k3[i]) / 3.0                    
            t = t + self.dt
            self.data.t = t
            #}

            #//if ((err[0] == 0) || (err[1] == 0)) s = 0; //trap division by 0
            #//else s = 0.84 * Math.Pow(param.Tol * dt / err[0], 0.25);//step size scalar

            #//if ((s < 0.75) && (dt > 2 * hmin)) dt /= 2.0;  //reduce step
            #//else if ((s > 1.5) && (2 * dt < hmax)) dt *= 2.0;//increase step
            #for (i = 0; i < maxi; i++)
            #{
            #    if (Double.IsNaN(y[i])) throw new NotFiniteNumberException();
            #}            
            return (t < self.tmax)#;//false if exceeded range

class RungeKutta45IntegratorTest(unittest.TestCase):
    def test_params_copy_constructor(self):        
        params = RungeKutta45IntegratorParams()
        params.dimension = 7
        params.Tmin = 0.01
        params.dT = 0.04
        params.Npoints = 10000
        copy_params = RungeKutta45IntegratorParams(params)
        self.assertEqual(params.dimension, copy_params.dimension)
        self.assertEqual(params.Tmin, copy_params.Tmin)
        self.assertEqual(params.dT, copy_params.dT)
        self.assertEqual(params.Npoints, copy_params.Npoints)

    def testSimpleHarmonicOScillator(self):
        import matplotlib.pyplot as plt
        params = RungeKutta45IntegratorParams()
        params.dimension = 2
        params.Tmin = 0.0
        params.dT = 0.001
        params.Npoints = 100000
        data = RungeKutta45IntegratorData(params)
        data.y[0] = 1.0
        data.y[1] = 2.0
        integrator = RungeKutta45ConstStepIntegrator(params,data,testfunction)
        states = integrator.IterateSeries()
        t = states[:,0]                
        x = states[:,1]        
        y = states[:,2]
        fig = plt.figure()
        plt.plot(x,y)
        plt.show()
        
        print(states)
        

def testfunction(t,y,Freturn,settings):
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
            Freturn[0] = y[1]
            Freturn[1] = -0.1*y[0]
            #print(len(y))
            return Freturn

def main():

    
        unittest.main()

if __name__ == "__main__":
    main()