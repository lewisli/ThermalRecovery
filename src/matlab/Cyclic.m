clear all
close all
clc

% Cyclic Steam Injection Example (Note: analytical model and data from SPE13037)
% Input data for field and base cases
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


%% Setting up plot
figure('Color','w');
set(gcf, 'Position', get(0,'Screensize'));

%% Evaluate constants
% Bulk volumetric heat capacity of Jones (eq 19)
cp = 32.5 + (4.6 * phi^0.32 - 2) * (10 * Swi - 1.5);    

% Water density approximation (eq 33)
rhoWater = 62.4 - 11 * log( ( 705 - Tstd ) / ( 705 - Tr ) );    

%% Starting cycles

for CurrentCycleNumber=1:length(SteamInjectionRate)
    
    % Time at beginning of each cycle
    t=0;          
    
    % Figure out how much heat is already in the reservoir
    if(CurrentCycleNumber==1)
        Hlast = 0;        
    else
        % Amount of heat in the reservoir for subsequent cycles (eq 20)
        Hlast = VSteam * cp * (Taverage - Tr);                            
    end
    
    %% Evaluate water and steam thermodynamic properties
    % Water enthalpy correlation for reservoir temperature (eq 16)
    hwTr = 68 * ( Tr/100 )^1.24;                              
    
    % Water enthalpy correlation for steam temperature (eq 16)
    hwSteam = 68 * ( DownHoleSteamTemp(CurrentCycleNumber) / 100 )^1.24;                    
    
    % Specific heat of water of Jones (eq 15)
    Cw = ( hwSteam - hwTr ) / ( DownHoleSteamTemp(CurrentCycleNumber) - Tr );        

    % Steam latent heat correlation of Farouq Ali (eq 17)
    Lvdh = 94 * ( 705 - DownHoleSteamTemp(CurrentCycleNumber) )^0.38;                      
    
    % Amount of heat injected per unit mass of steam (eq 14)
    Qi = Cw * ( DownHoleSteamTemp(CurrentCycleNumber) - Tr ) +...
        Lvdh * DownHoleSteamQuality(CurrentCycleNumber);                      
    
    % Amount of heat injected (eq 28)
    HeatInjected = 350 * Qi * SteamInjectionRate(CurrentCycleNumber)...
        * InjectionTime(CurrentCycleNumber);
        
    % Steam pressure approximation (eq 7)
    pSteam = ( DownHoleSteamTemp(CurrentCycleNumber) / 115.95 )^4.4543;                 
    
    % Steam density (eq 10)
    rhoSteam = pSteam^( 0.9588 ) / 363.9;                                  
    
    % Steam viscosity (eq 11)
    muSteam = 1e-4 * ( 0.2 * DownHoleSteamTemp(CurrentCycleNumber) + 82 );                  
    
    %% Evaluating the steam zone size
    % Oil density approximation (eq 32)
    rhoOil = rhoOilStd - 0.0214 * ( Tr - Tstd );          
  
    % Dimensionless group for scaling the radial steam zone (eq 9)
    ARD = sqrt( ( 350 * 144 * SteamInjectionRate(CurrentCycleNumber) * muSteam ) / ...     
        ( 6.328 * pi * ( rhoOil - rhoSteam ) * PayThickness^2 * ...
        kSteam * rhoSteam ) );
    
    % Average steam zone thickness by Van Lookeren (eq 8)
    hSt = 0.5 * PayThickness * ARD;             
    
    % Steam zone volume estimation (eq 13)
    VSteam = ( SteamInjectionRate(CurrentCycleNumber) * ...
        InjectionTime(CurrentCycleNumber) * rhoWater * Qi...
        + Hlast ) / ( cp * ( DownHoleSteamTemp(CurrentCycleNumber) - Tr ) );
    % Steam zone radius (eq 12)
    RhSteam = sqrt( VSteam / ( pi * hSt ) );       
    
    % Radial distance along the hot oil zone (eq 2)
    Rx = sqrt( RhSteam^2 + PayThickness^2 );       
    
    % theta = angle between steam-oil interface and reservoir bed (eq 5)
    sinTheta = PayThickness / Rx;                                   
    
    % Difference between height of reservoir and steam zone thickness, (eq 6)
    deltaH= PayThickness - hSt;  
    
    % Bottom hole flowing pressure (assumption)
    pwf = 0.6 * pSteam ;                       
    
    % Change in enthalpy (equation 4)
    deltaPhi = deltaH * g * sinTheta + ( ( ( pSteam - pwf )...      
        * 6895 )/ ( rhoOil * 16.02 ) ) * 10.76;
    
    % Change in oil saturation (eq 3)
    deltaSo = ( 1 - Swi ) - SorSteam; 
    
    %% Estimate relative permeability of water and oil
    % Cumulative water production at beginning of cycle
    Wp = 0.0;                                                              
    % Volumetric heat capacity of water (eq 31)
    Mwater = Cw * rhoWater;                                         
    % Mobile water around the well (eq 38)
    SwBar = 1 - SorWater;                                           
    % Water saturation  (eq 39)
    Sw = SwBar - (SwBar - Swi) * Wp / WIP;                          
    % Normalized water saturation (eq 40)
    SwStar = (Sw - Swi)/(1 - Swi - SorWater);                       
    % Water relative permeability (eq 41)
    krw = -0.002167 * SwStar + 0.024167 * SwStar^2;           
    
    % Oil relative permeability (eq 43)
    if SwStar <= 0.2
        kro = 1.0;                                                 
    else
        % Oil relative permeability (eq 42)
        kro = -0.9416 + 1.0808 / SwStar - 0.13858 / SwStar^2;     
    end
    
    % Oil rate vector
    qoPreviousCycle = zeros(cycleLength(CurrentCycleNumber),1);                            
    
    % Average temperature vector
    TaveragePreviousCycle = zeros(cycleLength(CurrentCycleNumber),1);                       
    
    % Cycling through InjectionTime + SoakTime + ProductionTime
    for t = 1:cycleLength(CurrentCycleNumber)                                               
        
        % Injection Interval
        if t <= (InjectionTime(CurrentCycleNumber))                                         
            %Average temperature during injection = downhole steam temperature
            Taverage = DownHoleSteamTemp(CurrentCycleNumber);                               
            qo = 0;
            
        end
        
        % Soaking interval
        if t > InjectionTime(CurrentCycleNumber) && ...
                t <= (InjectionTime(CurrentCycleNumber) + ...
                SoakTime(CurrentCycleNumber))  
            
            % Radial loss (eq 22)
            tDH = alpha * ( t - InjectionTime(CurrentCycleNumber) )...
                / RhSteam^2;             
            fHD = 1 / ( 1 + 5 * tDH );                                      
            
            % Vertical loss (eq 23)
            tDV = 4 * alpha * ( t - InjectionTime(CurrentCycleNumber) ) ...
                / PayThickness^2;    
            fVD = 1 / sqrt( 1 + 5 * tDV );                                  
            
            % Energy removed with produced fluids during soaking phase
            fPD = 0;                       
            
            % Oil rate during soaking phase
            qo = 0;                                                         
            
            % Average temperature at any time by Boberg and Lantz (eq 21)
            TpreviousTimeStep = Taverage;
            Taverage = Tr + ( DownHoleSteamTemp(CurrentCycleNumber) ...                      
                - Tr ) * ( fHD * fVD * ( 1 - fPD ) - fPD );
        end
        
        % Production Interval
        if t > (InjectionTime(CurrentCycleNumber) + SoakTime(CurrentCycleNumber))                            

            %Radial loss (eq 22)
            tDH  = alpha * ( t - InjectionTime(CurrentCycleNumber) )...
                / RhSteam^2;           
            fHD = 1 / ( 1 + 5 * tDH );                                     
            
            % Vertical loss (eq 23)
            tDV = 4 * alpha * ( t - InjectionTime(CurrentCycleNumber) )...
                / PayThickness^2;   
            fVD = 1 / sqrt( 1 + 5 * tDV );                                  
            
            % Registering average temperature from previous time step
            TpreviousTimeStep = Taverage;                                  
            
            % Average temperature at any time by Boberg and Lantz (eq 21)
            Taverage = Tr + ( DownHoleSteamTemp(CurrentCycleNumber) ...                      
                - Tr ) * ( fHD * fVD * ( 1 - fPD ) - fPD );
            
            % Volumetric heat capacity of oil (eq 30)
            Moil = ( 3.065 + 0.00355 * Taverage ) * sqrt( rhoOil );
            
            % Rate of heat removal from the reservoir with produced fluids (eq 29)
            Qp = 5.615 * ( qo * Moil + qw * Mwater )*( Taverage - Tr );
            
            % Maximum amount of heat supplied to the reservoir (eq 27)
            Qmax = HeatInjected + Hlast - pi * RhSteam^2 * lambda...
                * (DownHoleSteamTemp(CurrentCycleNumber) - Tr) * ...
                sqrt(SoakTime(CurrentCycleNumber)/(pi*alpha));
                       
            % Energy removed with produced fluids (eq 34)
            fPD  = fPD + 5.615 * ( qo * Moil + qw * Mwater )...
                * ( TpreviousTimeStep - Tr ) * TimeStep / ( 2 * Qmax );;
            
            % Oil viscosity (eq 36)
            muOil = ( 2.698e-5 ) * exp( ( 1.066e+4 ) /( Taverage + 460 ));
            
            % Kinematic viscosity of the oil = oil viscosity/oil density
            nuAverage = muOil / rhoOil;
            
            % Oil rate (eq 1)
            qo = 1.87 * Rx * sqrt ( ( kro * kRes * phi * deltaSo...
                * alpha * deltaPhi )/( 2.0  *nuAverage...
                * ( log( Rx / rw ) - 0.5 ) ) );
        end
        
        qoPreviousCycle(t) = qo;
        TaveragePreviousCycle(t) = Taverage;
        
    end
    
    CumulativeOil = [ CumulativeOil ; qoPreviousCycle ];
    time_range = linspace( cumulativeCycleLength(CurrentCycleNumber) ,...
        cumulativeCycleLength(CurrentCycleNumber+1) , ...
        (cumulativeCycleLength(CurrentCycleNumber+1) - ...
        cumulativeCycleLength(CurrentCycleNumber)) );
    CumulativeTime = [ CumulativeTime ; time_range' ];
    
    %Plotting production
    subplot(2,1,2)
    plot([cumulativeCycleLength(CurrentCycleNumber):1:(...
        cumulativeCycleLength(CurrentCycleNumber+1)')-1],...
        qoPreviousCycle,'-','Color',[0 0.5 0],'LineWidth',2);
    hold on;
    xlabel('Days','FontSize',14,'FontWeight','bold');
    ylabel('q_o (STB/day)','FontSize',14,'FontWeight','bold');
    set(gca,'FontSize',14,'FontWeight','bold');
    axis tight;
    
    % Plotting temperature
    subplot(2,1,1);
    plot([cumulativeCycleLength(CurrentCycleNumber):1:(...
        cumulativeCycleLength(CurrentCycleNumber+1))-1],...
        TaveragePreviousCycle,'r-','LineWidth',2);
    hold on;
    xlabel('Days','FontSize',14,'FontWeight','bold')
    ylabel('T_{avg}','FontSize',14,'FontWeight','bold')
    set(gca,'FontSize',14,'FontWeight','bold');
    axis tight;
    
    clear Tavg_saved qo_saved;
    
end

% Plotting cumulative oil production
figure('Color','w');
set(gcf, 'Position', get(0,'Screensize'));
plot(CumulativeTime,cumsum(CumulativeOil),'-','Color','k','LineWidth',2);
xlabel('Days','FontSize',14,'FontWeight','bold','Color','k')
ylabel('Cumulative oil produced (STB)','FontSize',14,'FontWeight','bold','Color','k')
set(gca,'FontSize',14,'FontWeight','bold');

