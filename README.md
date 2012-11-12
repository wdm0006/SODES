SODES
=====

MATLAB Symbolic Differential Equation Solver

SODES is a MATLAB function which takes in systems of 2nd order differential equations in standard MATLAB symbolic form, 
re-writes them into C as MEX function, and then integrates them using any of a number of MATLAB integrators.  SODES will
also plot, print and do other nice little things.  SODES is still in developement, and may include C (MEX) integrators
in the future.

Inputs

  eqs:                {2nd order diffeq, 2nd order diffeq, ....}

  free_var:           'x' or 'dx' or 'velocity', name of free variable.

  t_init:             initial time, in seconds

  t_final:            final time, in seconds

  initial conditions: intial conditions, state vector, in columns

  

  optional arguments:

  'Integrator': options: 'ode45', 'ode23', 'ode15s', 'mexRKF2','mexRKF4'  (default ode45)

  'TimeStep' : options: time step in seconds (default none)

  'Plotting'  : options: ON or OFF (default OFF)

  'Terminal'  : options: ON or OFF (default OFF)

  'StateNames': options: Cell array of state names (default y1,y2,...,yn)

  'StateUnits': options: Cell array of state units (default units)

  'Options'   : send an odeset options object



Current Math Limitations:

Not all math will convert right at this time.  The main source of problems centers around the differences in how C and MATLAB handle exponentials.  C uses ^ as the bitwise XOR operator, MATLAB as pow.  To get around this, there is some code to find any ^'s in the MATLAB code, and back track to find the clause that it is raising to whatever power. Right now the code will not convert correctly if:
Anything is raised to a power that isnt -1<x<10 (one digit only)
Any function (abs(x)^2, inv(x)^2, ect.) with the exception of:

*sin

*cos

*tan

*sec

*csc

*cot

*asin

*acos

*atan

*asec

*acsc

*acot

*sinh

*cosh

*tanh

*sech

*csch

*coth

*asinh

*acosh

*atanh

*asech

*acsch

*acoth

*exp

*log

*sqrt

Work will continue to account for more built in functions to becorrectly converted from MATLAB to C. 

Current Integrator Work:

mexRKF2 is a 2nd order Runge-Kutta implimentation, and seems to be working properly.  mexRKF4 is a 4th order Runge-Kutta implimentation and seems to have some problems, I am unsure of the issue right now.

Both of the mex functions assume a system of 2 2nd order diffeqs right now, and work needs to be done to generalize that so that they will work with any number of 2nd order diffeqs.