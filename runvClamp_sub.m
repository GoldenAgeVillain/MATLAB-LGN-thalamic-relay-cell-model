%
% run vClamp
%

parameters 

% vClamp to test iNa and iK amplitude
vClamp.Vbaseline= Vrest;
vClamp.Vclamp   = vClamp.Vbaseline+40;
vClamp.tstart   = 500;
vClamp.tend     = vClamp.tstart+200;

[time,s] = BazhenovModel(1,0,1,'vClamp',vClamp);

% regenerate V
[a,b]   = size(s);
V       = zeros(a,1);
for i = 1:1:a
    V(i) = vClamp.Vbaseline;
    if i*dt > vClamp.tstart && i*dt < vClamp.tend
        V(i) = vClamp.Vclamp;
    end
end

% lookup table
m       = s(:,1);
h       = s(:,2);
n       = s(:,3);
m_iA    = s(:,4);
h_iA    = s(:,5);
m_iT    = s(:,6);
h_iT    = s(:,7);
Ca      = s(:,8);
O       = s(:,9);

% reversals
ET      = 1e03*(R*T/(2*F)).*log(Ca0./Ca);

% currents
iL      = gL.*(V-EL);                           % iL current
iKL     = gKL.*(V-EK);                          % iKL current
iA      = gA*m_iA.^M_iA.*h_iA.^N_iA.*(V-EK);    % iA current
iT      = gT*m_iT.^M_iT.*h_iT.^N_iT.*(V-ET);    % iT current
ih      = gmax*O.*(V-Eh);                       % ih current
iNa     = gNa*m.^3.*h.*(V-ENa);                 % iNa current
iK      = gK*n.^4.*(V-EK);                      % iK current

iTot_1  = 1e06*(iL+iKL+iA+iT+ih+iNa+iK)*CS;

% vClamp to test iNa and iK amplitude
vClamp.Vbaseline= Vrest;
vClamp.Vclamp   = vClamp.Vbaseline-40;
vClamp.tstart   = 500;
vClamp.tend     = vClamp.tstart+200;

[time,s] = BazhenovModel(1,0,1,'vClamp',vClamp);

% regenerate V
[a,b]   = size(s);
V       = zeros(a,1);
for i = 1:1:a
    V(i) = vClamp.Vbaseline;
    if i*dt > vClamp.tstart && i*dt < vClamp.tend
        V(i) = vClamp.Vclamp;
    end
end

% lookup table
m       = s(:,1);
h       = s(:,2);
n       = s(:,3);
m_iA    = s(:,4);
h_iA    = s(:,5);
m_iT    = s(:,6);
h_iT    = s(:,7);
Ca      = s(:,8);
O       = s(:,9);

% reversals
ET      = 1e03*(R*T/(2*F)).*log(Ca0./Ca);

% currents
iL      = gL.*(V-EL);                           % iL current
iKL     = gKL.*(V-EK);                          % iKL current
iA      = gA*m_iA.^M_iA.*h_iA.^N_iA.*(V-EK);    % iA current
iT      = gT*m_iT.^M_iT.*h_iT.^N_iT.*(V-ET);    % iT current
ih      = gmax*O.*(V-Eh);                       % ih current
iNa     = gNa*m.^3.*h.*(V-ENa);                 % iNa current
iK      = gK*n.^4.*(V-EK);                      % iK current

iTot_2  = 1e06*(iL+iKL+iA+iT+ih+iNa+iK)*CS;

% figure
figure
subplot(211);
plot(time, V, '-k'); 
ylabel('V [mV]');
subplot(212);
plot(time, iTot_1-iTot_2, '-k');
xlabel('time [ms]');
ylabel('[pA]');