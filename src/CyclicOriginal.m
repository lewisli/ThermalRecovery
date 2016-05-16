

%CyclicSteamInjection.m

clear all
%close all
clc

% Cyclic Steam Injection Example (Note: analytical model and data from SPE13037)
% 
%

%% Input data for field and base cases
kRes = 1.5;                             % Reservoir permeability (D)
kSteam = 2 * kRes;                      % Permeability to steam (D)
phi = 0.32;                             % Porosity 
Swi = 0.25;                             % Initial water saturation
rw = 0.31;                              % Well radius (ft)
SorSteam = 0.05;                        % Residual oil saturation to steam 
SorWater = 0.25;                        % Residual oil saturation to water 
Tr = 110;                               % Initial reservoir temperature (F)
Tstd = 60;                              % Standard temperature (F) 
lambda = 24;                            % Reservoir thermal conductivity (Btu/ft day F)
alpha = 0.48;                           % Reservoir thermal diffusivity (ft^2/day)
rhoOilStd = 61.8;                       % Oil density @ standard condictions (lb/ft^3)
TimeStep = 1.0;                         % Time step size (days)
PatternArea=0.588;                      % (acres) = well spacing of 160ft
PayThickness = 80;                      % Pay thickness (ft)
g = 32.17;                              % Gravitational acceleration (ft/s^2)
%%

%% Operational data
SteamInjectionRate = [647 905 953 972 954 1042 1067];      % Steam injection rate (bbl/day)
InjectionTime = [6 9 8 6 6 7 7];                           % Injection time (days)
SoakTime = [5 2 10 12 9 10 11];                            % Soak time (days)
DownHoleSteamTemp = [360 330 330 300 300 300 300];         % Downhole steam temperature (F)
DownHoleSteamQuality = [0.7 0.6 0.5 0.5 0.5 0.5 0.5];      % Downhole steam quality 
ProductionTime = [55 146 79 98 135 136 148];               % Production time (days)
WIP = 50000;                                               % Amount of mobile water in place at beginning of cycle (Assumption)
qw = 0;                                                    % water production rate (BPD) (Assumption)
CumulativeOil = [];                                        % Cumulative oil production (STB) 
CumulativeTime = [];                                       % Cumulative time (days) 
cycleLength = InjectionTime + SoakTime + ProductionTime;   % Total cycle time (days)
cumulativeCycleLength = cumsum([1;cycleLength']);          % Cumulative cycle time (days) 
%%

%% Setting up plot
figure('Color','w');
set(gcf, 'Position', get(0,'Screensize'));
%%
 
%% Starting cycles
 
for i=1:length(SteamInjectionRate) 
        
    t=0;                                                                   % Time at beginning of each cycle                                          
    hwTr = 68 * ( Tr/100 )^1.24;                                           % Water enthalpy correlation for reservoir temperature (eq 16)
    hwSteam = 68 * ( DownHoleSteamTemp(i) / 100 )^1.24;                    % Water enthalpy correlation for steam temperature (eq 16)
    Cw = ( hwSteam - hwTr ) / ( DownHoleSteamTemp(i) - Tr );               % Specific heat of water of Jones (eq 15)         
    Lvdh = 94 * ( 705 - DownHoleSteamTemp(i) )^0.38;                       % Steam latent heat correlation of Farouq Ali (eq 17)
    Qi = Cw * ( DownHoleSteamTemp(i) - Tr ) +...
        Lvdh * DownHoleSteamQuality(i);                                    % Amount of heat injected per unit mass of steam (eq 14)
    cp = 32.5 + (4.6 * phi^0.32 - 2) * (10 * Swi - 1.5);                   % Bulk volumetric heat capacity of Jones (eq 19)  
    pSteam = ( DownHoleSteamTemp(i) / 115.95 )^4.4543;                     % Steam pressure approximation (eq 7)
    rhoSteam = pSteam^( 0.9588 ) / 363.9;                                  % Steam density (eq 10)
    muSteam = 1e-4 * ( 0.2 * DownHoleSteamTemp(i) + 82 );                  % Steam viscosity (eq 11)
    rhoOil = rhoOilStd - 0.0214 * ( Tr - Tstd );                           % Oil density approximation (eq 32)
    rhoWater = 62.4 - 11 * log( ( 705 - Tstd ) / ( 705 - Tr ) );           % Water density approximation (eq 33)     
    ARD = sqrt( ( 350 * 144 * SteamInjectionRate(i) * muSteam ) / ...      % Dimensionless group for scaling the radial steam zone (eq 9) 
          ( 6.328 * pi * ( rhoOil - rhoSteam ) * PayThickness^2 * ...
          kSteam * rhoSteam ) ); 
    hSt = 0.5 * PayThickness * ARD;                                        % Average steam zone thickness by Van Lookeren (eq 8)  
    Wp = 0.0;                                                              % Cumulative water production at beginning of cycle
    
    if(i~=1)    
        Hlast = VSteam * cp * (Taverage - Tr);                             % Amount of heat in the reservoir for subsequent cycles (eq 20)
    else
        Hlast = 0;                                                         % Amount of heat in the reservoir before first cycle starts 
    end
    
    VSteam = ( SteamInjectionRate(i) * InjectionTime(i) * rhoWater * Qi... % Steam zone volume estimation (eq 13)
        + Hlast ) / ( cp * ( DownHoleSteamTemp(i) - Tr ) );         
    RhSteam = sqrt( VSteam / ( pi * hSt ) );                               % Steam zone radius (eq 12)
    
    qoPreviousCycle = zeros(cycleLength(i),1);                             % Oil rate vector     
    TaveragePreviousCycle = zeros(cycleLength(i),1);                       % Average temperature vector
    
    for t = 1:cycleLength(i)                                               % Cycling through InjectionTime + SoakTime + ProductionTime
                
        if t <= (InjectionTime(i))                                         % Injection Interval 
            
            Taverage = DownHoleSteamTemp(i);                               %Average temperature during injection = downhole steam temperature
            qo = 0;
            
        end
                
        if t > InjectionTime(i) && t <= (InjectionTime(i) + SoakTime(i))   % Soak interval
            
           tDH = alpha * ( t - InjectionTime(i) ) / RhSteam^2;             % Vertical loss (eq 23)
           fHD = 1 / ( 1 + 5 * tDH );                                      % Radial loss (eq 22)                           
           tDV = 4 * alpha * ( t - InjectionTime(i) ) / PayThickness^2;    % (eq 25)       
           fVD = 1 / sqrt( 1 + 5 * tDV );                                  % (eq 24)
           fPD = 0;                                                        % Energy removed with produced fluids during soaking phase
           qo = 0;                                                         % Oil rate during soaking phase
           Taverage = Tr + ( DownHoleSteamTemp(i) ...                      % Average temperature at any time by Boberg and Lantz (eq 21)
                    - Tr ) * ( fHD * fVD * ( 1 - fPD ) - fPD );            
        end
                
        if t > (InjectionTime(i) + SoakTime(i))                            % Production Interval         
            
           tDH  = alpha * ( t - InjectionTime(i) ) / RhSteam^2;            % Vertical loss (eq 23)        
           fHD = 1 / ( 1 + 5 * tDH );                                      % Radial loss (eq 22)  
           tDV = 4 * alpha * ( t - InjectionTime(i) ) / PayThickness^2;    % (eq 25)            
           fVD = 1 / sqrt( 1 + 5 * tDV );                                  % (eq 24)     
           TpreviousTimeStep = Taverage;                                   % Registering average temperature from previous time step
           Taverage = Tr + ( DownHoleSteamTemp(i) ...                      % Average temperature at any time by Boberg and Lantz (eq 21)
                    - Tr ) * ( fHD * fVD * ( 1 - fPD ) - fPD );            
           rhoOil = rhoOilStd - 0.0214 * ( Tr - Tstd );                    % Oil density approximation (eq 32)
           Moil = ( 3.065 + 0.00355 * Taverage ) * sqrt( rhoOil );         % Volumetric heat capacity of oil (eq 30)
           rhoWater = 62.4 - 11 * log( ( 705 - Tstd ) / ( 705 - Tr ) );    % Water density approximation (eq 33)   
           Mwater = Cw * rhoWater;                                         % Volumetric heat capacity of water (eq 31)
           Qp = 5.615 * ( qo * Moil + qw * Mwater )*( Taverage - Tr );     % Rate of heat removal from the reservoir with produced fluids (eq 29)
           HeatInjected = 350 * Qi * SteamInjectionRate(i)...              % Amount of heat injected (eq 28)
               * InjectionTime(i);     
           Qmax = HeatInjected + Hlast - pi * RhSteam^2 * lambda...        % Maximum amount of heat supplied to the reservoir (eq 27) 
               * (DownHoleSteamTemp(i) - Tr) * sqrt(SoakTime(i)/(pi*alpha));      
           deltafPD = 5.615 * ( qo * Moil + qw * Mwater )...               % (eq 35)
                * ( TpreviousTimeStep - Tr ) * TimeStep / ( 2 * Qmax );             
           fPD  = fPD + deltafPD;                                          % Energy removed with produced fluids (eq 34)
           muOil = ( 2.698e-5 ) * exp( ( 1.066e+4 ) /( Taverage + 460 ));  % Oil viscosity (eq 36)               
           Rx = sqrt( RhSteam^2 + PayThickness^2 );                        % Radial distance along the hot oil zone (eq 2)
           sinTheta = PayThickness / Rx;                                   % theta = angle between steam-oil interface and reservoir bed (eq 5)       
           deltaH= PayThickness - hSt;                                     % Difference between height of reservoir and steam zone thickness, (eq 6)
           pwf = 0.6 * pSteam ;                                            % Bottom hole flowing pressure (assumption) 
           deltaPhi = deltaH * g * sinTheta + ( ( ( pSteam - pwf )...      % (eq 4)
               * 6895 )/ ( rhoOil * 16.02 ) ) * 10.76;
           SwBar = 1 - SorWater;                                           % Mobile water around the well (eq 38)
           Sw = SwBar - (SwBar - Swi) * Wp / WIP;                          % Water saturation  (eq 39)
           SwStar = (Sw - Swi)/(1 - Swi - SorWater);                       % Normalized water saturation (eq 40)
           krw = -0.002167 * SwStar + 0.024167 * SwStar^2;                 % Water relative permeability (eq 41)
           
           if SwStar <= 0.2
                kro = 1.0;                                                 % Oil relative permeability (eq 43)
           else
                kro = -0.9416 + 1.0808 / SwStar - 0.13858 / SwStar^2;      % Oil relative permeability (eq 42)
           end
           
           deltaSo = ( 1 - Swi ) - SorSteam;                               % Change in oil saturation (eq 3)
           nuAverage = muOil / rhoOil;                                     % Kinematic viscosity of the oil = oil viscosity/oil density     
           qo = 1.87 * Rx * sqrt ( ( kro * kRes * phi * deltaSo...         % Oil rate (eq 1)
               * alpha * deltaPhi )/( 2.0  *nuAverage...
               * ( log( Rx / rw ) - 0.5 ) ) );          
        end
        
        qoPreviousCycle(t) = qo; 
        TaveragePreviousCycle(t) = Taverage;                           
                                
    end        
        
    CumulativeOil = [ CumulativeOil ; qoPreviousCycle ];
    time_range = linspace( cumulativeCycleLength(i) , cumulativeCycleLength(i+1) , (cumulativeCycleLength(i+1) - cumulativeCycleLength(i)) );        
    CumulativeTime = [ CumulativeTime ; time_range' ];

    %Plotting production
    subplot(2,1,2)
    plot([cumulativeCycleLength(i):1:(cumulativeCycleLength(i+1)')-1],qoPreviousCycle,'-','Color',[0 0.5 0],'LineWidth',2);
    hold on;
    xlabel('Days','FontSize',14,'FontWeight','bold');
    ylabel('q_o (STB/day)','FontSize',14,'FontWeight','bold');
    set(gca,'FontSize',14,'FontWeight','bold');
    axis tight;
    
    % Plotting temperature
    subplot(2,1,1);
    plot([cumulativeCycleLength(i):1:(cumulativeCycleLength(i+1)')-1],TaveragePreviousCycle,'r-','LineWidth',2);
    hold on;
    xlabel('Days','FontSize',14,'FontWeight','bold')
    ylabel('T_{avg}','FontSize',14,'FontWeight','bold')        
    set(gca,'FontSize',14,'FontWeight','bold');
    axis tight;

    clear Tavg_saved qo_saved;    
    
end

% Plotting cumulative oil production
% figure('Color','w');
% set(gcf, 'Position', get(0,'Screensize'));
% plot(CumulativeTime,cumsum(CumulativeOil),'-','Color','k','LineWidth',2);
% xlabel('Days','FontSize',14,'FontWeight','bold','Color','k')
% ylabel('Cumulative oil produced (STB)','FontSize',14,'FontWeight','bold','Color','k')
% set(gca,'FontSize',14,'FontWeight','bold');
% 






%%

