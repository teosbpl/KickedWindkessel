// Rk45 Ordinary Differential Equations Solver (ODES)
//Solve the differential equation  using Runge-Kutta-Fehlberg Method
//with variable step size (Adjust "Tol" for error tolerance!)
//ref: J.Mathews-"Numerical Methods for Mathematics, Science and Engineering"
//Prentice Hall, 1992
//copywrite 2002
//       Cristian Constantin Bordeianu, Romania and R Landau, Oregon State U
//e-mail:cbord@warpnet.ro
// Translated to C# by T.Buchner
//e-mail: buchner@if.pw.edu.pl
//Partially based on http://csharpcomputing.com/Tutorials/Lesson16.htm Copyright 2001-2006 Aleksey Nudelman
using System;
using UserMath.ModelUtils;
namespace UserMath.Integrators
{
    /// <summary>
    /// Runge - Kutta 45(?) variable step integrator
    /// </summary>
    public class RungeKutta45ConstStepIntegrator
    {
        public delegate void Function(double t, double[] y, ref double[] yReturn); //declare a delegate    
        public class RungeKutta45IntegratorParams
        {
            public int dimension;
            public double Tmin;                         //endpoints
            private double dT;
            public int Npoints;
            public RungeKutta45IntegratorParams() 
            {
                this.dimension = 0;
                this.Tmin = 0.0d;                         //endpoints
                this.Npoints = 0;
                this.dT = 0.0;
            }
            public RungeKutta45IntegratorParams(RungeKutta45IntegratorParams data)
            {
                this.dimension = data.dimension;
                //this.Tol = data.Tol;
                this.Tmin = data.Tmin;
                this.dT = data.dT;
                //this.Tmax = data.Tmax;
                this.Npoints = data.Npoints;
            }
            
            public double Tdelta
            {
                get { return dT; }
                set { dT = value; }
            }
	
        }
        private RungeKutta45IntegratorParams param;

        public class RungeKutta45IntegratorData: IntegratorData
        {
            public RungeKutta45IntegratorData(RungeKutta45IntegratorParams param):
                base(param.dimension,param.Tmin)
            {                                
            }
        }

        public int Dimension
        {
            get { return param.dimension; }
        }
        private double dt;        
        private double[] F;
        private double[] ydumb;
        private double[] k1;
        private double[] k2;
        private double[] k3;
        private double[] k4;
        private double tmax;

        public double Tmax
        {
            get { return tmax; }            
        }
        //private double[] err;

        public RungeKutta45ConstStepIntegrator(RungeKutta45IntegratorParams param)
        {
            Reset(param);
        }
        public void Reset(RungeKutta45IntegratorParams param)
        {
            this.param = new RungeKutta45IntegratorParams(param);
            this.F = new double[this.param.dimension];
            this.ydumb = new double[this.param.dimension];
            this.k1 = new double[this.param.dimension];
            this.k2 = new double[this.param.dimension];
            this.k3 = new double[this.param.dimension];
            this.k4 = new double[this.param.dimension];
            this.dt = this.param.Tdelta;
            this.tmax = this.param.Tmin + this.param.Npoints * this.dt;
        }
        //double
        public bool Iterate(ref double t, ref double[] y, Function f)
        {
            int i;
            int maxi = this.param.dimension;
            if (t > this.tmax) return false;
            //printout
            if ((t + dt) > this.tmax) dt = this.tmax - t; // the last step
            /*
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
             */
            f(t, y, ref F);  // evaluate both RHS's and return in F
            for (i = 0; i < maxi; i++)//was <= i
            {
                k1[i] = dt * F[i];
                ydumb[i] = y[i] + k1[i] / 2.0;
            }

            f(t + dt / 2, ydumb, ref F);
            for (i = 0; i < maxi; i++)
            {
                k2[i] = dt * F[i];
                ydumb[i] = y[i] + k2[i] / 2.0;
            }

            f(t + dt / 2, ydumb, ref F);
            for (i = 0; i < maxi; i++)
            {
                k3[i] = dt * F[i];
                ydumb[i] = y[i] + k3[i];
            }

            f(t + dt, ydumb, ref F);
            for (i = 0; i < maxi; i++)
            {
                k4[i] = dt * F[i];
            }

            if (true)//((err[0] < param.Tol) || (err[1] < param.Tol) || (dt <= 2 * hmin))//accept the approximation
            {
                // xn+1=nx+(k1+2k2+2k3+k4)/6                          
                for (i = 0; i < maxi; i++)
                    y[i] = y[i] + (k1[i] + k4[i]) / 6.0 + (k2[i] + k3[i]) / 3.0;
                /*
                if (y.Length > 2)//space for drive
                {
                    f(t, ydumb, ref F);
                    y[2] = F[2];
                }
                 */ 
                t = t + dt;
            }

            //if ((err[0] == 0) || (err[1] == 0)) s = 0; //trap division by 0
            //else s = 0.84 * Math.Pow(param.Tol * dt / err[0], 0.25);//step size scalar

            //if ((s < 0.75) && (dt > 2 * hmin)) dt /= 2.0;  //reduce step
            //else if ((s > 1.5) && (2 * dt < hmax)) dt *= 2.0;//increase step
            for (i = 0; i < maxi; i++)
            {
                if (Double.IsNaN(y[i])) throw new NotFiniteNumberException();
            }            
            return (t < this.tmax);//false if exceeded range
        }
    }

}
