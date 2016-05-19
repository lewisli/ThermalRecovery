%% SAGD.m
% Estimates the production rate and percent recovery of original oil in
% place as a function of time in a SAGD project.
%
%
% Authors: Lewis Li, Daniel Brito
% Date: May 18th 2016

clear all; close all;

% Reservoir temperature (C)
Tr = 15;                    
% Steam temperature (C)
Ts = 188;   
% Kinematic viscosity exponent (m)
m_exp = 3.4;                
% Bitumen kinematic  viscosity @ 188C (cs)
mu_Bitumen188 = 7.8;      
% Reservoir thickness (m)
resThickness = 20;   
% Thermal diffusivity (m^2/D)
alpha = 0.07;    
% Porosity
phi = 0.33;      
% Initial oil saturation
S_o = 0.75;
% Residual oil saturation
S_or = 0.13;              
% Effective permeability for oil flow (Darcy)
K_eff = 0.4;            
% Distance between reservoir's base and producers (m)
baseWellDistance = 2.5;  
% Spacing between wells (m)
w = 75;   
% Evaluation period (years)
TotalTime = 7;     
% Gravity (m/s^2)
g = 9.81;                 

%% Conversion Factors
DarcyToM2 = 9.869233e-13;
DaysToSeconds = 24*60*60;
YearsToDays = 365;
centistokesToM2S= 1e-6;

%% Step 1: Calculate rates for depletion   

n = 100;

%Time discretization (s)
t = linspace(0,TotalTime*YearsToDays*DaysToSeconds,n);    

%Time discretization (years)
tYears = t./(YearsToDays*DaysToSeconds);        

% Height of steam zone
h = resThickness - baseWellDistance;    

% Dimensionless time
tStarConstant = (2.0/w) * sqrt (((K_eff * DarcyToM2)*g*...
    alpha/(DaysToSeconds) ) / ( phi * (S_o - S_or) * ...
    m_exp * mu_Bitumen188*centistokesToM2S * h)); 

tStar = tStarConstant*t;

% Dimensionless flow rate
qStar = sqrt(3/2) - sqrt(2/3)*tStar.^2;
               
FFactor = sqrt( (m_exp*mu_Bitumen188*centistokesToM2S) /...
    ( (K_eff * DarcyToM2) * g *(alpha/(DaysToSeconds)) * ...
    h * phi * (S_o - S_or) ) );                

q = 2*qStar*DaysToSeconds./FFactor;

% Analytical integral of q* polynomial
Recovery = sqrt(3/2)*tStar - tStar.^3*sqrt(2/3)/3;

DepletionResults = [tYears' tStar' qStar' Recovery' q'];
ColumnNames = {'Time, years', 't*', 'q*','Recovery','q,m^3/(m day)'}; 

DepletionResults = [ColumnNames; num2cell(DepletionResults)];
display(DepletionResults)

%% Step 2: Calculate rates for the rising steam chamber
Coef1 = (((K_eff * DarcyToM2)*g*alpha/(DaysToSeconds))/...
    ( m_exp * mu_Bitumen188*centistokesToM2S ) )^(2/3);
Coef2 = (phi * (S_o - S_or) )^(1/3);

qCumRise = 2.25*Coef1*Coef2*t.^(4/3)*DaysToSeconds;
qRise = 3*Coef1*Coef2*t.^(1/3)*DaysToSeconds;
RecoveryRise = qCumRise./(h*phi*(S_o - S_or)*w*DaysToSeconds);

SteamResults = [tYears' qRise' RecoveryRise'];
ColumnNames = {'Time, years', 'q,m^3/(m day)','Recovery',}; 

EndIndex = sum(RecoveryRise<=1)+1;
Results = [ColumnNames; num2cell(SteamResults)];
display(Results(1:EndIndex,:))

%% Step 3: Find the time changeover point
addpath('../matlab');
[RecoveryIntersection,qIntersection] = ...
    intersections(Recovery',q,RecoveryRise',qRise,1); 

% Determine time step where intersect occurs for depletion
RecoveryRow = sum(Recovery<=RecoveryIntersection);

% Determine time step where intersect occurs for depletion
RecoveryRiseRow = sum(RecoveryRise<=RecoveryIntersection);

timeCalcDepletion = YearsToDays * tYears(RecoveryRow-1);
timeCalcSteam = YearsToDays * tYears(RecoveryRiseRow-1);
deltaTime = timeCalcSteam - timeCalcDepletion;

display(['The changeover point occurs at ' num2str(qIntersection) ...
    ' m^3/m of production and recovery of ' ...
    num2str(RecoveryIntersection)]);
display(['Corresponds to ' num2str(timeCalcSteam) ...
    ' days production from steam']);
display(['Corresponds to ' num2str(timeCalcDepletion) ...
    ' days production from depletion']);
display(['Means a difference of ' num2str(deltaTime) ' days']);

%%  Step 4: Calculate the actual rates
% The rates for production are from steam until the changeover point, and
% then the rates are from depeletion.
ChangeOverTime = timeCalcSteam/YearsToDays;
StepSize = 0.25;

% Steam-based production (computed as function of time in seconds)

timeSteam = (0:StepSize:1.5)*ChangeOverTime*YearsToDays*DaysToSeconds;
qCumSteam = (2.25 * Coef1 * Coef2 * (timeSteam).^(4/3)) * DaysToSeconds;
qSteam = (3 * Coef1 * Coef2 * (timeSteam).^(1/3)) * DaysToSeconds;        
RecoverySteam = qCumSteam ./ ( h*phi*(S_o - S_or)*w*DaysToSeconds );

% Set up time for depeletion (including time correction)
StepSize = 0.05;
timeDepletion = (0:(ChangeOverTime*StepSize):(TotalTime+2) - ...
    deltaTime/YearsToDays) * YearsToDays*DaysToSeconds; 
tStarDepletion = tStarConstant * timeDepletion; 
qStarDepletion = sqrt(1.5) - (tStarDepletion.^2)*sqrt(2/3);
q = 2 * qStarDepletion * (DaysToSeconds)/( FFactor ); 
recoveryFactor = sqrt(3/2)*tStarDepletion - ...
    (1/3)*(tStarDepletion.^3)*sqrt(2/3); 

figure;
FontSize = 15;
hold on;
plot(RecoverySteam,qSteam);
plot(recoveryFactor,q);
axis tight;
ylabel('Production rate, m^3/day','Fontsize', FontSize) % label left y-axis
xlabel('Recovery, fraction','Fontsize', FontSize) % label right y-axis
legend('Steam','Depletion');
set(gca,'FontSize',FontSize);

%% Step 5: Tabulate and plot results
% Convert time from seconds to years and remove time correction
timeSteamYears = (1/(YearsToDays*DaysToSeconds)) * timeSteam; 
timeDepletionYears = (1/(YearsToDays*DaysToSeconds)) * ...
    timeDepletion + deltaTime/YearsToDays; 

% Find when steam time ends and depeletion starts
SteamEndTime = sum(timeSteamYears<=ChangeOverTime);
DepletionStartTime = find(timeDepletionYears>=ChangeOverTime,1);

qTotal = [(qSteam(:,1:SteamEndTime))';(q(:,DepletionStartTime:end))'];
timeTotal = [(timeSteamYears(:,1:SteamEndTime))' ; ...
    (timeDepletionYears(:,DepletionStartTime:end))'];

qStarTotal = [ zeros(SteamEndTime,1); ...
    (qStarDepletion(:,DepletionStartTime:end))'];
tStarTotal = [ zeros(SteamEndTime,1); ...
    (tStarDepletion(:,DepletionStartTime:end))'];
RecoveryTotal = [ (RecoverySteam(:,1:SteamEndTime))'; ...
    (recoveryFactor(:,DepletionStartTime:end))'];

figure
[haxes,hline1,hline2] = plotyy(timeTotal,qTotal,timeTotal, ...
    RecoveryTotal,'plot','plot');
ylabel(haxes(1),'Production rate, m^3/day','Fontsize', ...
    FontSize) % label left y-axis
ylabel(haxes(2),'Recovery, fraction','Fontsize', ...
    FontSize) % label right y-axis
xlabel(haxes(2),'Time (years)', 'Fontsize', ...
    FontSize) % label x-axis

Nticks = 5;
FontSize = 12;
set(haxes,'XLim',[0 TotalTime])
set(haxes(1),'YLim',[0 max(qTotal)]);
set(haxes(1),'ytick',linspace(0, max(qTotal), Nticks));
set(haxes(2),'YLim',[0 max(RecoveryTotal)]);
set(haxes(2),'ytick',linspace(0, max(RecoveryTotal), Nticks));
set(gcf,'color','w');
set(haxes(1),'FontSize',FontSize);
set(haxes(2),'FontSize',FontSize);

%% Display results
TotalResults = [timeTotal qTotal tStarTotal qStarTotal RecoveryTotal];
ColumnNames = {'Time, years', 'q,m^3/(m day)','t*','q*','Recovery'}; 
TotalResults = [ColumnNames; num2cell(TotalResults)];

EndIndex = sum(timeTotal<=TotalTime)+1;
display(TotalResults(1:EndIndex,:))


