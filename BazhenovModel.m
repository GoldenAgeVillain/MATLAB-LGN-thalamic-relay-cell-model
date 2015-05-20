function [time,S] = BazhenovModel(gain,switchPrint,switchPlot,condition,vClamp)

%% parameters

global tstart DT DTA dt Kshift Nashift iHold;
DT = []; DTA = []; dt = []; Cm = []; EL = []; ENa = []; EK = []; Eh = [];...
    EAMPA = []; gL = []; gT = []; gNa = []; gK = []; gKL = []; gA  = [];...
    gmax = []; R = []; T = []; F = []; Ca0 = []; Cai = []; Cac = [];...
    tauCa=[]; A = []; M_iA = []; N_iA = []; M_iT = []; N_iT = [];...
    k = []; k2 = []; k4 = []; Pc = []; Vrest = []; CS = []; Kshift = [];...
    Nashift = []; deadTime = []; ampa = []; nmda = []; iHold = [];
parameters

%% generate a random conductance

% seed    = 62;
% rng(seed);    % seed the RNG
% SynCond = 0.001+0.12*randn(DT/dt+1,1);

%% load the experimentally recorded conductances

load('ampa conductance');               % nS
load('nmda conductance');               % nS
ampa = [0;ampa];
nmda = [0;nmda];
gAMPA = gain*1e-06/CS*ampa;             % mS/cm2
gNMDA = gain*1e-06/CS*nmda;             % mS/cm2
gTot  = gAMPA + gNMDA;                  % mS/cm2
clear ampa nmda;

%% initialize the variables and run

if strcmp(condition,'rest')
    DT = 1000; dt = 0.5; DTA = 0;
elseif strcmp(condition,'vClamp')
    DT = 1000; dt = 0.05; DTA = 0;
end
s0      = [0.0398 0.9849 0.0833 0.6342 0.0223 0.6441 0.0010 0.0003 0.0301 -55.3223];
%s0      = [0.0006 0.7565 0.0029 0.1147 0.4739 0.0491 0.1967 0.0003 0.6062 Vrest];
options = odeset('MaxStep',dt,'RelTol',1e-03,'OutputFcn',@myfun);
tspan   = 0:dt:DT;
tstart  = clock;
[time,S]= ode15s(@fxn,tspan,s0,options);

%% print data

if switchPrint == 1
    filename = ['output/results_' condition '_' num2str(gain,'%05.2f') '.mat'];
    timeClipped     = time(deadTime:end)-time(deadTime);
    SClipped        = S(deadTime:end,:);
    save(filename,'timeClipped','SClipped');
end

%% figure

if switchPlot == 1
    figure;
    lbl     = {'m','h','n','m_{iA}','h_{iA}','m_{iT}','h_{iT}','Ca^{2+} [mM]','O','V [mV]'};
    for i = 1:1:10
        subplot(2,5,i);
        p = plot(time,S(:,i),'-k');
        set(p,'LineWidth',2);
        set(p,'Color',[1 0 0]);
        if i == 8
            set(p,'Color','k');
        elseif i == 10
            set(p,'Color','k');
        end
        ylabel(lbl(i));
    end
end

%% calculate and display frequency

if strcmp(condition,'CTRL') || strcmp(condition,'NMDA')
    [pks,locs]  = findpeaks(S(deadTime:end,end),'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',3/dt);
    if exist('pks')
        frequency   = 1e03*length(pks)/DTA;
    else
        frequency = 0;
    end
    fprintf(['\n']);
    fprintf(['output frequency = ' num2str(frequency,'%3.2f') ' (Hz)\n\n']);
    finish      = (clock-tstart)*[0 0 24*60*60 60*60 60 1]';
    fprintf(['finished after = ' num2str(finish/60,'%3.2f') ' min\n\n']);
end

%% main nested subfunction

    function ds = fxn(t,s)
        
        % boundaries
        s(1:9)  = max(s(1:9),0);
        
        % lookup table
        m       = s(1);
        h       = s(2);
        n       = s(3);
        m_iA    = s(4);
        h_iA    = s(5);
        m_iT    = s(6);
        h_iT    = s(7);
        Ca      = s(8);
        O       = s(9);
        V       = s(10);
        if strcmp(condition,'vClamp')
            if t > vClamp.tstart && t < vClamp.tend
                V = vClamp.Vclamp;
            else
                V = vClamp.Vbaseline;
            end
        end
        
        % calcium reversal potential
        ET      = 1e03*(R*T/(2*F))*log(Ca0/Ca);
        
        % injected current
        if strcmp(condition,'vClamp')
            iHold   = 0;
            iInj    = 0;
        elseif strcmp(condition,'rest')
            iInj    = 0;                            % µA/cm2
        elseif strcmp(condition,'CTRL')
            iInj    = ...
                gTot(floor(t/dt)+1)*(V-EAMPA);
        elseif strcmp(condition,'NMDA')
            iInj    =...
                gAMPA(floor(t/dt)+1)*(V-EAMPA)+...
                gNMDA(floor(t/dt)+1)*(V-EAMPA)*9.69/(1+0.1688*exp(-0.0717*V));
        end
        %         iInj    = SynCond(floor(t/dt)+1)*(V-EAMPA);
        
        % currents
        iL      = gL*(V-EL);                        % iL current
        iKL     = gKL*(V-EK);                       % iKL current
        iA      = gA*m_iA^M_iA*h_iA^N_iA*(V-EK);    % iA current
        iT      = gT*m_iT^M_iT*h_iT^N_iT*(V-ET);    % iT current
        ih      = gmax*O*(V-Eh);                    % ih current
        iNa     = gNa*m^3*h*(V-ENa);                % iNa current
        iK      = gK*n^4*(V-EK);                    % iK current
        
        % total intrinsic currents
        iInt    = iL+iKL+iA+iT+ih+iNa+iK;
        
        % main equation
        ds(1)   = mGate(m,V);                       % dm/dt
        ds(2)   = hGate(h,V);                       % dh/dt
        ds(3)   = nGate(n,V);                       % dn/dt
        ds(4)   = mGateiA(m_iA,V);                  % dm_iA/dt
        ds(5)   = hGateiA(h_iA,V);                  % dh_iA/dt
        ds(6)   = mGateiT(m_iT,V);                  % dm_iT/dt
        ds(7)   = hGateiT(h_iT,V);                  % dh_iT/dt
        ds(8)   = -1/tauCa*(Ca-Cai)-A*iT;           % dCa2+/dt
        ds(9)   = openState(O,V);                   % dO/dt
        ds(10)  = -1/Cm*(iInt+iInj+iHold);          % dV/dt
        
        ds      = ds';                              % transpose the vector of derivatives
        ds(isnan(ds)) = 0;                          % avoids NaN in the vector of derivatives
        ds(isinf(ds)) = 0;                          % avoids Inf in the vector of derivatives
        if strcmp(condition,'vClamp')
            ds(10) = 0;
        end
        
    end

end

%% subfunctions

% iNa current
function dmdt = mGate(m,V)
global Nashift;
V       = V-Nashift;        % convention in the Traub model
alpham  = 0.32*(13.1-V)/(exp((13.1-V)/4)-1);
betam   = 0.28*(V-40.1)/(exp((V-40.1)/5)-1);
dmdt    = alpham*(1-m)-betam*m;
end

function dhdt = hGate(h,V)
global Nashift;
V       = V-Nashift;        % convention in the Traub model
alphah  = 0.128*exp((17-V)/18);
betah   = 4/(1+exp((40-V)/5));
dhdt    = alphah*(1-h)-betah*h;
end

% iK current
function dndt = nGate(n,V)
global Kshift;
V       = V-Kshift;         % convention in the Traub model
alphan  = 0.032*(15-V)/(exp((15-V)/5)-1);
betan   = 0.5*exp((10-V)/40);
dndt    = alphan*(1-n)-betan*n;
end

% iA current
function dmdt = mGateiA(m_iA,V)
minf = 1/(1+exp(-(V+60)/8.5));
taum = 0.27/(exp((V+35.8)/19.7)+exp(-(V+79.7)/12.7))+0.1;
dmdt = 1/taum*(minf-m_iA);
end

function dhdt = hGateiA(h_iA,V)
hinf = 1.0/(1+exp((V+78)/6));
if V < -63
    tauh = 0.27/(exp((V+46)/5)+exp(-(V+238)/37.5));
else
    tauh = 5.1;
end
dhdt = 1/tauh*(hinf-h_iA);
end

% iT current
function dmdt = mGateiT(m_iT,V)
minf = 1/(1+exp(-(V+57)/6.2));
taum = (0.22/(exp(-(V+132)/16.7)+exp((V+16.8)/18.2))+0.13);
dmdt = 1/taum*(minf-m_iT);
end

function dhdt = hGateiT(h_iT,V)
hinf = 1/(1+exp((V+83)/4));
tauh = (8.2+(56.6+0.27*exp((V+115.2)/5))/(1+exp((V+86)/3.2)));
dhdt = 1/tauh*(hinf-h_iT);
end

% ih current voltage dependency (from Huguenard and McCormick, J Neurophysiol 1992)
function dOdt = openState(O,V)
minf = 1/(1+exp((V+75)/5.5));                       % their Equ. 5 and parameters in page 1380
taum = 1/(exp(-14.59-0.086*V)+exp(-1.87+0.0701*V)); % their Equ. 13
dOdt = 1/taum*(minf-O);
end

% subfunction for output
function status = myfun(t,s,flag)
global tstart DT;
eta = (clock-tstart)*[0 0 24*60*60 60*60 60 1]';
fprintf([...
    't = ' num2str(t,'%0.2f') ' ms || ' num2str(100*t/DT,'%0.2f')...
    '%% completed || ETA = ' num2str(eta*(DT-t)./(t*60),'%03.2f') ' min\n']);
status = 0;
end
