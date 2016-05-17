clear all
close all
clc

% Typical properties:Example: For Kern River, assume a project life of 20
% years and plot thermal efficiency, area of steam zone, and cumulative oil 
% production (Marx-Langenheim Calculation)

%Typical properties
Tr = 90;                % Reservoir temperature (F)
phi = 0.32;             % Porosity 
h = 55;                 % Net pay thickness(ft)
Sor = 0.15;             % Residual oil saturation           
deltaSo = 0.4;          % Oil saturation variation
Ms = 42;                % BTU/cu.ft-F
Mt = 35;                % BTU/cu.ft-F
Ks = 1.2;               % BTU/cu.ft-F-hr

%Injection conditions
Pbh = 100;              % Bottom hole injection pressure (psig)
rate = 360;             % bbl/day CWE (i.e., the rate on a condensed water basis)
Xbh = 0.5;              % Steam quality
acresPerWell = 2.5;     % Drainage area 

%Steam properties (from steam tables)
Ts = 338;                       % Steam temperature (F)
deltaHv = 880;                  % delta_H for vapour
HwTs = 310;                     % BTU/lb @ steam temp
HwTr = 58;                      % BTU/lb @ res temp 
Cw = (HwTs - HwTr) / (Ts - Tr); % (BTU/lb-F)

%Mass injection rate
mi = 360 * 5.615 * 62.4 * 365;  % Rate of injection (lb/year)

%Thermal diffusivity
alphaSs = 1.2*24/42;                    % thermal diffusivity (ft^2/day)

%Calculating hd



%Time discretization
t = [0.02098296, 0.10491482, 0.20982964, 1.34290969, 2.09829639, 4.19659278,...
    6.29488917, 8.39318557, 10.491482, 20.9829639, 41.9659278, 62.9488917];
t = t * 365;

% Computing critical time (CHECK IF REQUIRED)

%tcd = fsolve(@(tcD) exp(tcD).*erfc(sqrt(tcD)) -1/(1+hd),[1 4], optimset('Display','iter'));  
%tcd = tcd(length(tcd));
%tcrit = tcd * (h^2) * (Mt/Ms)^2 / (4*alphaSs); % Critical time (s)

% Calculating td and Qi
td = 4 * t * ( Ms / Mt )^2 * ( alphaSs / h^2 );
ti = t / 365;
Qi = ( mi * Cw * (Ts - Tr)  + Xbh * mi * deltaHv ) * ti;

% Calculating Ehs (Marx-Langenheim):
N = length(td);
Ehs = zeros (1,N);

for i = 1:N     
    Ehs(i) = (2.0 * sqrt ( td(i) / pi ) - 1 + exp( td(i) ).* erfc( sqrt(td(i)) ))...
        /td(i);   
end

% Computing steam volume, cumulative production and oil/steam ratio:
Vs = Qi.*Ehs / (Mt * (Ts - Tr));
Np = Vs .* phi * deltaSo / 5.615;
OSR= Np./(mi*(t/365)/(62.4*5.615));
Acres = Vs ./ h;

%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot

figure1 = figure('Color','w');
set(gcf, 'Position', get(0,'Screensize'));

%NP
subplot(2,1,2);
plot(t,Vs,'r-','LineWidth',2);
hold on;
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('Vs','FontSize',14,'FontWeight','bold')        
set(gca,'FontSize',14,'FontWeight','bold');

%Production
subplot(2,1,1)
plot(t,OSR,'k-','LineWidth',2);
hold on;
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('OSR','FontSize',14,'FontWeight','bold')        
set(gca,'FontSize',14,'FontWeight','bold');


Results = [td' t' Ehs' Qi' Vs' Acres' Np' OSR'];
format shortEng
format compact

fprintf('\n Results are as follows: \n\n');
disp('           tD              t              Ehs             Qi              Vs              Acres              Np              OSR  ')
Results


