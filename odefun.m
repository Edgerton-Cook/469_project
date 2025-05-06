function dXdt = odefun(t,X,params)
%{
This function is to be read in by ODE45 to simulate the system.
- NOTE: May not be usable due to coupeling of the theta's

Inputs:
    t - time series we want to simulate over
    X - vector of the states 
    parms - additional necessary parameteres, including masses, lengths,
    and gravity
Outputs:
    dXdt - X' vector
%}
end