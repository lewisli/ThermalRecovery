%% SteamDrive.m
% Estimates the production rate and percent recovery of original oil in
% place as a function of time in a SAGD project.
%
%
% Authors: Lewis Li, Daniel Brito
% Date: May 18th 2016

% Typical properties:Example: For Kern River, assume a project life of 20
% years and plot thermal efficiency, area of steam zone, and cumulative oil 
% production (Marx-Langenheim Calculation)

clear all;

%Typical properties
Tr = 90;                % Reservoir temperature (F)
phi = 0.32;             % Porosity 
h = 55;                 % Net pay thickness(ft)
Sor = 0.15;             % Residual oil saturation           
deltaSo = 0.4;          % Oil saturation variation
Ms = 42;                % BTU/cu.ft-F
Mt = 35;                % BTU/cu.ft-F
Ks = 1.2;               % BTU/cu.ft-F-hr

% Injection properties
Pbh = 100;                      % Bottom hole injection pressure (psig)
rate = 360;                     % bbl/day CWE (i.e., the rate on a condensed water basis)
Xbh = 0.5;                      % Steam quality
acresPerWell = 2.5;             % Drainage area 
Ts = 338;                       % Steam temperature (F)
deltaHv = 880;                  % delta_H for vapour
HwTs = 310;                     % BTU/lb @ steam temp
HwTr = 58;                      % BTU/lb @ res temp 

%Mass injection rate
DaysInYear = 365;
VolInBarrel = 5.615; %ft3/B
Density = 62.4; %lb/ft^3
mi = rate * VolInBarrel * Density * DaysInYear;  % Rate of injection (lb/year)
Cw = (HwTs - HwTr) / (Ts - Tr); % (BTU/lb-F)

% Thermal diffusivity (ft^2/day)
alphaSs = 1.2*24/42;           


% Descritize time
t = [0.02098296, 0.10491482, 0.20982964, 1.34290969, 2.09829639, 4.19659278,...
    6.29488917, 8.39318557, 10.491482, 20.9829639, 41.9659278, 62.9488917];
t = t * DaysInYear;

% Calculate dimensionless time
td = 4 * t * ( Ms / Mt )^2 * ( alphaSs / h^2 );
ti = t / DaysInYear;

% Calculate injection rate
Qi = ( mi * Cw * (Ts - Tr)  + Xbh * mi * deltaHv ) * ti;

%% Marx-Langenheim
EhMarx = (2*sqrt(td/pi) - 1 + exp(td).* erfc( sqrt(td)))./td;
VsMarx = Qi.*EhMarx / (Mt * (Ts - Tr));
NpMarx = VsMarx .* phi * deltaSo / VolInBarrel;
OSRMarx= NpMarx./(mi*(t/DaysInYear)/(Density*VolInBarrel));
AcresMarx = VsMarx ./ h;

%% Myhill-Stegemeier
hd = Xbh * deltaHv / (Cw *(Ts - Tr) );  % latent heat / sensible heat
tcd = fsolve(@(tcD) exp(tcD).*erfc(sqrt(tcD)) -1/(1+hd),[1], ...
optimset('Algorithm','levenberg-marquardt','Display','iter'));  
tcd = tcd(end);
G = 2 * sqrt ( td / pi ) - 1 + exp( td ).* erfc( sqrt(td));

N = length(td);
EhMyHill = zeros (1,N);

for i = 1:N     
    errorFun = @(u) exp(u).*erfc(sqrt(u))./sqrt(td(i) - u); 
    EhMyHill(i) = 1/td(i) * ( G(i) + ( (td(i) > tcd)/( sqrt(pi)*(1+hd) ) )  ...
        * ( 2*sqrt(td(i)) - ( 2*sqrt(td(i) - tcd) )/( 1 + hd )...
        - integral(errorFun, 0, tcd) - sqrt(pi) * G(i) ) );    
end

% Computing steam volume, cumulative production, oil/steam ratio and Area:
VsMyHill = Qi.*EhMyHill / (Mt * (Ts - Tr));    % Steam volume (ft^3)
NpMyHill = VsMyHill .* phi * deltaSo / VolInBarrel;   % Np (bbl)
OSRMyHill= NpMyHill./(mi*(t/DaysInYear)/(Density*VolInBarrel)); % oil/steam ratio
AreaMyHill = VsMyHill ./ h;                     % Area (acres)

%% Plot results
%plot inline -s 2000,1000
figure(1);
subplot(1,2,1);
plot(t,VsMarx,'r-','LineWidth',2);
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('Vs','FontSize',14,'FontWeight','bold')        
title({'Steam Volume Using Marx-Langenheim';''},'FontSize',14);
set(gca,'FontSize',14,'FontWeight','bold');
axis tight; axis square;

subplot(1,2,2);
plot(t,VsMyHill,'r-','LineWidth',2);
hold on;
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('Vs','FontSize',14,'FontWeight','bold')        
set(gca,'FontSize',14,'FontWeight','bold');
title({'Steam Volume Using Myhill-Stegemeier';''},'FontSize',14);
axis tight; axis square;

figure(2);
%plot inline -s 2000,1000
subplot(1,2,1);
plot(t,OSRMarx,'r-','LineWidth',2);
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('OSR','FontSize',14,'FontWeight','bold')        
title({'OSR Using Marx-Langenheim';''},'FontSize',14);
set(gca,'FontSize',14,'FontWeight','bold');
axis tight; axis square;

subplot(1,2,2);
plot(t,OSRMyHill,'r-','LineWidth',2);
hold on;
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('OSR','FontSize',14,'FontWeight','bold')        
set(gca,'FontSize',14,'FontWeight','bold');
title({'OSR Using Myhill-Stegemeier';''},'FontSize',14);
axis tight; axis square;

figure(3);
%plot inline -s 2000,1000
subplot(1,2,1);
plot(t,AcresMarx,'r-','LineWidth',2);
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('Area (acres)','FontSize',14,'FontWeight','bold')        
title({'Area of Steam Zone Using Marx-Langenheim';''},'FontSize',14);
set(gca,'FontSize',14,'FontWeight','bold');
axis tight; axis square;

subplot(1,2,2);
plot(t,AreaMyHill,'r-','LineWidth',2);
hold on;
xlabel('Days','FontSize',14,'FontWeight','bold')
ylabel('Area (acres)','FontSize',14,'FontWeight','bold')        
set(gca,'FontSize',14,'FontWeight','bold');
title({'Area of Steam Zone Using Myhill-Stegemeier';''},'FontSize',14);
axis tight; axis square;

%% Display results
Results = [td' t' EhMarx' Qi' VsMarx' AcresMarx' NpMarx' OSRMarx'];
format shortEng
format compact

display('Results for Marx-Langenheim are as follows:');
display('              tD               t             Ehs              Qi');
display(Results(:,1:4));
display('              Vs           Acres              Np             OSR');
display(Results(:,5:end));

Results = [td' t' EhMyHill' Qi' VsMyHill' AreaMyHill' NpMyHill' OSRMyHill'];
format shortEng
format compact

display('Results for Myhill-Stegemeier are as follows:');
display('              tD               t             Ehs              Qi');
display(Results(:,1:4));
display('              Vs           Acres              Np             OSR');
display(Results(:,5:end));