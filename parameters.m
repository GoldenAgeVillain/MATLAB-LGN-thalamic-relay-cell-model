%
% parameters of the model
%

% duration and time step
DT  = 130000;       % length of simulation [ms]
DTA = 125000;       % length of simulation after discarding the initial 5 sec [ms]
dt  = 0.05;         % time step [ms]
deadTime = 5/(dt*1e-03);  
                    % first 5 discarded seconds [time steps]

% capacitance and reversal potentials
Cm  = 1;            % membrane capacitance [µF/cm2]
EL  = -70;          % reversal potential [mV]                               
ENa = +90;          % sodium reversal potential [mV]                        % OK               
EK  = -105;         % potassium reversal potential [mV]                     % OK
Eh  = -43;          % Ih reversal potential [mV]                            % from Huguenard and McCormick, J Neurophysiol 1992, page 1380
EAMPA = 0;          % reversal potential for AMPARs [mV]

% conductances
gL  = 0.025;        % leak conductance [mS/cm2]                             % OK       
gT  = 1.8;          % low threshold calcium conductance [mS/cm2]
gNa = 4.4;          % sodium conductance [mS/cm2]                  
gK  = 3.3;          % potassium conductance [mS/cm2]              
gKL = 0.025;        % potassium leak conductance [mS/cm2]                   % OK                  
gA  = 3;            % potassium A conductance [mS/cm2]
gmax= 2*0.0127;       % Ih conductance [mS/cm2]                             % OK

% physical constants
R   = 8.31441;      % gas constant [J/(mol*K)]
T   = 309.15;       % temperature [K]
F   = 96489;        % Faraday constant [C/mol]

% calcium dynamics
Ca0 = 2;            % extracellular calcium [mM]
Cai = 2.4e-04;      % baseline intracellular calcium [mM]
Cac = 0.002;        % [mM]
tauCa=5;            % calcium extrusion/buffering time constant [ms]
A   = 5.18e-5;      % calcium entry constant [mM*cm2/(ms*µA)]

% exponents for gating variables
M_iA = 4;           % power for m gating variable in iA current
N_iA = 1;           % power for h gating variable in iA current
M_iT = 2;           % power for m gating variable in iT current
N_iT = 1;           % power for h gating variable in iT current
k    = 2;           % appears in the definition of the ih current

% calcium dependency of ih (removed)
% k2   = 0.0004;      % [ms-1]
% k4   = 0.001;       % [ms-1]
% Pc   = 0.01;

% parameters from Julia's experiments
CS      = 1.52e-4;  % cell surface [cm^2]                                   % OK                                                                                          
Vrest   = -77.4;    % resting membrane potential [mV]                       % OK     
iHold   = -2.05;    % [µA/cm^2]
Kshift  = -62.5;
Nashift = -60.1;