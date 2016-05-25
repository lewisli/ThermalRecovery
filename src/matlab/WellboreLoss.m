% Wellbore Losses
% Authors: Lewis Li, Daniel Brito
%
% Computes the losses incurred within a wellbore during steam injection.
% Uses example from (Prats, M. 1982)
%

clear all;
close all;
clc

%% Input data
% Temperature of steam (F)
Tb = 600;
% Ambient temperature of subsurface
Ta = 100;
% Time (days)
time = 21;
% Time function f(td) for the radiation boundary condition model
% Table 10.1 (Prats, M. 1982)
load('Table_10pt1_Prats.mat');

tableVertical = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0];
tableHorizontal = [100.0, 50.0, 20.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1,...
    0.05, 0.02, 0.01, 0.0];
% Well depth (ft)
depth = 1000;
% Tubing outer radius (ft)
ro = 0.1458;
% Insulation radius (ft)
rIns = 0.2292;
% Casing inner radius(ft)
rci = 0.3556;
% Casing outer diameter (ft)
rco = 0.401;
% Wellbore radius (ft)
rw = 0.5;
% Thermal diffusivity of the earth (sq-ft/D)
alphaE = 0.96;
% Emissivity of insulation surface
epsilonIns = 0.9;
% Emissivity of casing surface
epsilonCi = 0.9;
% Thermal conductivity of earth (Btu/ft-D-F)
lambdaE = 24;
% Thermal conductivity of cement (Btu/ft-D-F)
lambdaCem = 12;
% Thermal conductivity of insulation (Btu/ft-D-F)
lambdaIns = 0.96;
% Air viscosity @ average annulus temperature (Figure B.41)
viscosityAir = 0.023;
% Thermal conductivity of air (Btu/ft-D-F) (Figure B.72)
lambdaAir = 0.45;

%% Step 1: Initial guess for specific thermal resistance (Btu/ft-D-F)^-1.
% Sum of all thermal resistances is twice that due to insulation
RhGuess = (1/pi)*(log(rIns/ro))/lambdaIns;
Rprime = RhGuess;
error = Inf;
ConvergenceCriterion = 0.001;

% Check for convergence
while (error > ConvergenceCriterion)
    % Step 2: Calculate f(td) (Ramey function)
    % Calculate dimensionless time
    tD = alphaE*time/(rw*rw);
    
    if (tD > 100.0)
        % Eqn 10.10 for tD < 100
        func = 0.5*log(tD) + 0.403;
    else
        lookupVal = 2*pi*Rprime*lambdaE;
        % Interpolate table 10.1
        func = interp2(tableHorizontal, tableVertical,...
            Table_10pt1, lookupVal, tD);
    end

    % Step 3: Calculate Tci (temperature @ inner surface of casing)
    % Eqn B.68 (neglecting 1st and 3rd terms in parenthesis)
    Tci = Ta + ( (Tb - Ta)/(2*pi*Rprime) ) * ...
        (log(rw/rco)/lambdaCem + func/lambdaE);
    
    % Step 4: Calculate Tins (temperature @ outer face of insulation)
    % Eqn B.70 (only 5th term in parenthesis)
    Tins = Tb - ((Tb - Ta)/(2*pi*Rprime))*...
        (log(rIns/ro)/lambdaIns);
    
    %Step 5: Calculate hRcAn coefficient of annular heat transfer due
    %          to radiation and convection
    % Average temperature @ annulus
    TavgAnnulus = 0.5*(Tins + Tci);
    densityAir = 0.076*(460 + 60)/(460 + TavgAnnulus);
    
    % Isobaric thermal coefficient of volume 
    % expansion for gas (Assuming that air is 
    % an ideal gas.)
    betaG = 1.0/(460 + TavgAnnulus);
    % Eqn B.66: Grashof number
    Grashof = 7.12e7*((rci - rIns)^3*densityAir^2*betaG*(Tins - ...
        Tci))/(viscosityAir*viscosityAir);
    % Prandtl number (assumption)
    Prandtl = 0.92;
    % Eqn B.65: Apparent thermal conductivity of air in the annulus
    lambdaAAn = 0.049*lambdaAir*(Grashof^0.333)*(Prandtl^0.407);
    % Eqn B.64: Radiation temperature function
    F = ((460 + Tins)^2 + (460 + Tci)^2)*(920 + Tins + Tci);
    % Eqn B.63: Coefficient of annular heat transfer
    hRcAn = (4.11e-8 / (1/epsilonIns + (rIns/rci)*...
        (1/epsilonCi - 1)) ) * F+(1/rIns)*...
        lambdaAAn/log(rci/rIns);
    
    % Step 6: Calculate Rh using Eqn 10.6
    % Heat loss across insulation
    Rh1 = log(rIns/ro) / lambdaIns;
    % Heat loss through radiation and convection across annulus
    Rh2 = 1.0 / (hRcAn * rIns );
    % Heat loss across cement
    Rh3 = log(rw/rco) / lambdaCem;
    % Heat loss related to thermal resistance of earth
    Rh4 = func / lambdaE;
    % Eqn 10.6: Overall coefficient of heat loss
    Rh = (1./(2*pi))* (Rh1 + Rh2 + Rh3 + Rh4);
 
    % Convergence check
    error = abs(Rh - Rprime);
    Rprime = Rh;
end

%% Calculate heat loss
%  Eqn 10.1: Heat loss per unit depth of the well (Btu/ft-D)
Qls = (Tb - Ta) / Rh;
% Heat loss for the given depth of the well and the given time period
Qt = Qls * depth;

display(['Temperature of casing is ' num2str(Tci) ' degrees F']);
display(['Heat loss rate from ' num2str(depth) ...
    ' ft. deep well after ' num2str(time) ' days  is ' ...
    num2str(Qt,'%4.2e') ' BTU/Day']);
