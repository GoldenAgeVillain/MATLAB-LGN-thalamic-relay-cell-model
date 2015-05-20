%
% plot steady-state IVs
%

parameters
V       = -120:1:100;
VsK     = V-Kshift;
VsNa    = V-Nashift;
figure;

% iNa
alpham  = 0.32*(13.1-VsNa)./(exp((13.1-VsNa)./4)-1);
betam   = 0.28*(VsNa-40.1)./(exp((VsNa-40.1)./5)-1);
minf    = alpham./(alpham+betam);
alphah  = 0.128.*exp((17-VsNa)./18);
betah   = 4./(1+exp((40-VsNa)./5));
hinf    = alphah./(alphah+betah);
s(1)    = subplot(511);
iNagate = minf.^3.*hinf;
iNa     = gNa.*iNagate.*(V-ENa); 
p1      = plot(V,iNagate,'-r'); hold on;
p1p     = plot(V,iNa,'-r');
m(1)    = max(minf.^3.*hinf);
ylabel('iNa');

% iK
alphan  = 0.032.*(15-VsK)./(exp((15-VsK)./5)-1); 
betan   = 0.5.*exp((10-VsK)./40);
ninf    = alphan./(alphan+betan);
s(2)    = subplot(512);
iKgate  = ninf.^4;
iK      = gK.*iKgate.*(V-EK);
p2      = plot(V,iKgate,'-k'); hold on;
p2p     = plot(V,iK,'-k');
m(2)    = max(ninf.^4);
ylabel('iK');

% iA
minf    = 1./(1+exp(-(V+60)./8.5));
hinf    = 1.0./(1+exp((V+78)./6));
s(3)    = subplot(513);
iAgate  = minf.^M_iA.*hinf.^N_iA;
iA      = gA.*iAgate.*(V-EK);
p3      = plot(V,iAgate,'-g'); hold on;
p3p     = plot(V,iA,'-g');
m(3)    = max(minf.^M_iA.*hinf.^N_iA);
ylabel('iA');

% iT
minf    = 1./(1+exp(-(V+57)./6.2));
hinf    = 1./(1+exp((V+83)./4));
s(4)    = subplot(514);
iTgate  = minf.^M_iT.*hinf.^N_iT;
iT      = gT.*iTgate.*(V-120);
p4      = plot(V,iTgate,'-g'); hold on;
p4p     = plot(V,iT,'-g');
m(4)    = max(minf.^M_iT.*hinf.^N_iT);
ylabel('iT');

% ih
minf    = 1./(1+exp((V+75)./5.5));                     
s(5)    = subplot(515);
ihgate  = minf;
ih      = gmax.*ihgate.*(V-Eh);
p5      = plot(V,ihgate,'-c'); hold on;
p5p     = plot(V,ih,'-c');
m(5)    = 1;
xlabel('voltage [mV]');
ylabel('ih');

% leak currents;
iL      = gL*(V-EL);
iKL     = gKL*(V-EK);

% others
set([p1 p2 p3 p4 p5],'LineWidth',3);
for i = 1:1:5
   axes(s(i));
   l = line([Vrest Vrest],[0 1.3]);
   set(l,'LineStyle','-','Color','k');
   l = line([-55 -55],[0 1.3]);
   set(l,'LineStyle','-','Color','r');
   axis([V(1) V(end) 0 1.2*m(i)]);
end

fig = figure;
axe = axes();
set([p1p p2p p3p p4p p5p],'Parent',axe);
xlabel('voltage [mV]');
ylabel('current [mS/cm^2]');

figure;
load('averageIVs');
avgiv   = plot(output(:,1),output(:,2),'ok'); hold on;
p       = plot(V,iNa+iK+iT+iA+ih+iL+iKL,'-r');
plot([V(1) V(end)],[0 0],'-k');
set(p,'LineWidth',3);
set(avgiv,'MarkerFaceColor','k');
plot([Vrest Vrest],[-1000 1000],'--k');
plot([-55 -55],[-1000 1000],'-g');
axis([V(1) V(end) -10 10]);
xlabel('voltage [mV]');
ylabel('current [µA/cm^2]');
axis([-120 0 -10 10]);