clear all; close all;
Tr = 15;                    % Reservoir temperature (C)
Ts = 188;                   % Steam temperature (C)
KmuExp = 3.4;               % Kinematic viscosity exponent (m)
MuOil = 100000;             % Oil viscosity @ Tr (cs)
rhoBitumen = 1.0;           % g/cm^3
KmuBitumen100 = 80;         % Bitumen kinematic viscosity @ 100C (cs)
KmuBitumen188 = 7.8;        % Bitumen kinematic  viscosity @ 188C (cs)
muBitumen100 = 80;          % Bitumen viscosity @ 100C (cs)
resThickness = 20;          % Reservoir thickness (m)
alpha = 0.07;               % Thermal diffusivity (m^2/D)
phi = 0.33;                 % Porosity
So = 0.75;                  % Initial oil saturation
Sor = 0.13;                 % Residual oil saturation
Keff = 0.4;                 % Effective permeability for oil flow (Darcy)
pSteam = 1.2;               % Steam pressure (MPa)
baseWellDistance = 2.5;     % Distance between reservoir's base and producers (m)
w = 75;                     % Spacing between wells (m)
evalPeriod = 7;             % Evaluation period (years)
g = 9.81;                   % Gravity (m/s^2)

% Step 1: Calculate rates for depletion   
h = resThickness - baseWellDistance;    % Height

tStarConstant = (2.0/w) * sqrt ( ( (Keff * 9.869233e-13) * g * alpha/(24*60*60) ) / ...
                   ( phi * (So - Sor) * KmuExp * KmuBitumen188*1e-6 * h)); 
               
FFactor = sqrt( (KmuExp*KmuBitumen188*1e-6) /...
    ( (Keff * 9.869233e-13) * g *(alpha/(24*60*60)) * h * phi * (So - Sor) ) );                

t = linspace(0,7*365*24*60*60,100);     %Time discretization (s)
tYears = t./(365*24*60*60);             %Time discretization (years)
n = length(t);
tStar = zeros(n,1);
qStar = zeros(n,1);
q = zeros(n,1);
Recovery = zeros(n,1);

for i = 1:n       
    tStar(i) = tStarConstant * t(i); 
    qStar(i) = sqrt(1.5) - (tStar(i)^2)*sqrt(2/3);
    q(i) = 2*qStar(i)*60*60*24 / FFactor;    
    Recovery(i) = sqrt(3/2)*tStar(i) - (tStar(i)^3)*sqrt(2/3)/3;               
end

%% Step 2: Calculate rates for the rising steam chamber
qCumRise = zeros(n,1); 
qRise = zeros(n,1); 
RecoveryRise = zeros(n,1);
Coef1 = ( ( (Keff * 9.869233e-13) * g * alpha/(24*60*60) ) /...
    ( KmuExp * KmuBitumen188*1e-6 ) )^(2/3);
Coef2 = (phi * (So - Sor) )^(1/3);

for i = 1:n
    qCumRise(i) = 2.25 * Coef1 * Coef2 * (t(i))^(4/3) * (24*60*60);    
    qRise(i) = 3 * Coef1 * Coef2 * (t(i))^(1/3) * (24*60*60);
    RecoveryRise(i) = qCumRise(i) / ( h * phi * (So - Sor) * w * (24*60*60) );
end
%% Step 3: Find the time changeover point
[RecoveryIntersection,qIntersection] = intersections(Recovery,q,RecoveryRise,qRise,1); % Call to matlab function

RecoveryRow=1;
while Recovery(RecoveryRow,1)<=RecoveryIntersection %determining time for intersection for depletion
    RecoveryRow=RecoveryRow+1;
end

RecoveryRiseRow=1;
while RecoveryRise(RecoveryRiseRow,1)<RecoveryIntersection %determining time for intersection for steam
    RecoveryRiseRow=RecoveryRiseRow+1;
end

timeCalcDepletion = 365 * tYears(RecoveryRow-1);
timeCalcSteam = 365 * tYears(RecoveryRiseRow-1);
deltaTime = timeCalcSteam - timeCalcDepletion;
fprintf('The changeover point occurs at %4.3f m^3/m day, and recovery of %4.2f. \n',qIntersection,RecoveryIntersection)    
fprintf('This point corresponds to %.0f days for production from steam.\n',timeCalcSteam);
fprintf('This point corresponds to %.0f days for production from depletion.\n',timeCalcDepletion);
fprintf('(%.0f days difference)\n',deltaTime)
%%

%%  Step 4: Calculate the actual rates
timeSteam = [0, 0.5, 1, timeCalcSteam/365, 2] * 365*24*60*60; % seconds
timeDepletion = ([0, 2, 3, 4, 5, 6, 7, 8] - deltaTime/365) * 365*24*60*60; % seconds (including time correction)

% depletion-based production
tStarDepletion = tStarConstant * timeDepletion; 
qStarDepletion = sqrt(1.5) - (tStarDepletion.^2)*sqrt(2/3);
q = 2 * qStarDepletion * (24*60*60)/( FFactor ); 
recoveryFactor = sqrt(3/2)*tStarDepletion - (1/3)*(tStarDepletion.^3)*sqrt(2/3); 

% steam-based production:
qCumSteam = (2.25 * Coef1 * Coef2 * (timeSteam).^(4/3)) * (24*60*60) ;
qSteam = (3 * Coef1 * Coef2 * (timeSteam).^(1/3)) * (24*60*60);        
RecoverySteam = qCumSteam ./ ( h*phi*(So - Sor)*w*24*60*60 );


%% Step 5: Tabulate and plot results

timeSteam = (1/(365*24*60*60)) * timeSteam; % years
timeDepletion = (1/(365*24*60*60)) * timeDepletion + deltaTime/365; % (removing time correction)
qTotal = [(qSteam(:,1:4))' ; (q(:,2:7))'];
timeTotal = [(timeSteam(:,1:4))' ; (timeDepletion(:,2:7))'];
qStarTotal = [ zeros((length(timeSteam)-1),1); (qStarDepletion(:,2:7))' ];
tStarTotal = [ zeros((length(timeSteam)-1),1); (tStarDepletion(:,2:7))' ];
RecoveryTotal = [ (RecoverySteam(:,1:4))'; (recoveryFactor(:,2:7))'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = zeros(length(timeTotal),5);
results(:,1) = timeTotal;
results(:,2) = qTotal;
results(:,3) = tStarTotal;
results(:,4) = qStarTotal;
results(:,5) = RecoveryTotal;

fprintf('\n Results are as follows: \n\n');

disp('  time       q (m^3/m day)       t_star       q_star       Recovery  ')
disp(results)

figure
plot(recoveryFactor,q,RecoverySteam,qSteam);
axis([0 1 0 0.145]);
ylabel('Production rate, m^3/day','Fontsize', 16, 'FontWeight', 'bold','Color','b') % label left y-axis
xlabel('Recovery, fraction','Fontsize', 16, 'FontWeight', 'bold','Color',[0 0.5 0]) % label right y-axis

figure
[haxes,hline1,hline2] = plotyy(timeTotal,qTotal,timeTotal, RecoveryTotal,'plot','plot');
ylabel(haxes(1),'Production rate, m^3/day','Fontsize', 16, 'FontWeight', 'bold','Color','b') % label left y-axis
ylabel(haxes(2),'Recovery, fraction','Fontsize', 16, 'FontWeight', 'bold','Color',[0 0.5 0]) % label right y-axis
xlabel(haxes(2),'Time (years)', 'Fontsize', 16, 'FontWeight', 'bold','Color','k') % label x-axis
set(hline1,'LineWidth',4,'Color','b');
set(hline2,'LineWidth',4,'Color',[0 0.5 0]);
set(haxes,'FontSize',16,'FontWeight','bold');

