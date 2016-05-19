%% HeatLosses.m
% Estimates steady state heat losses per year per 100ft when steam is
% injected into a pipe that can be insulated with calcium silicate or not.
%
% Authors: Lewis Li, Daniel Brito
% Date: May 18th 2016

ri = 0.1478;     % ft
ro = 0.1667;     % ft
r_ins = 0.4167;  % ft
lamda_I = 0.96;  % Btu/ft-D-F 
lamda_P = 600;   % BTU/ft-D-degF
hf = 48000;      % Btu/sq ft-D-F
hfc = 154;       % Btu/sq ft-D-F
hpi = inf;       % Btu/sq ft-D-F
hpo = 48000;     % Btu/sq ft-D-F
Tb = 550;        % F
Ta = 60;         % F

%% With insulation
v_w = 20;      
hfc = 18 * (v_w^0.6) * (ro^0.6) / ro;

hrc = 110; % Btu/sq ft-D-F 
h_bare = hfc + hrc;  

R_h = (1 / (hf * ri) + 1 / (hpi * ri) + ...
(1 / lamda_P) * log(ro/ri) + 1 / ( (h_bare) * ro ))/(2*pi);

Qls = (Tb - Ta) / R_h; %Btu/ft-D

L = 100; % ft
t = 365; % days

Ql = Qls * L * t; % BTU

display('For an uninsulated pipe:');
display(['Specific thermal resistance was ' ...
num2str(R_h,3) ' BTU/ft-D'])
display(['Heat loss per unit length was ' ...
num2str(Qls,3) ' BTU/ft-D'])
display(['Heat loss from ' num2str(L)  ...
' of pipe over ' num2str(t) ' days was ' ...
num2str(Ql,3) ' BTU'])

%% Without insulation
R_h = (1 / (hf * ri) + 1 / (hpi * ri) + ...
(1 / lamda_P) * log(ro/ri) + 1 / ( (hpo) * ro ) + ...
(1 / lamda_I)*log(r_ins/ro) + 1 / (hfc * r_ins) )/(2*pi);

Qls = (Tb - Ta) / R_h; %Btu/ft-D
Ql = Qls * L * t; % BTU

display('For an insulated pipe:');
display(['Specific thermal resistance was ' ...
num2str(R_h,3) ' BTU/ft-D'])
display(['Heat loss per unit length was '...
num2str(Qls,3) ' BTU/ft-D'])
display(['Heat loss from ' num2str(L)  ' feet of pipe over ' ...
num2str(t) ' days was ' num2str(Ql,3) ' BTU'])