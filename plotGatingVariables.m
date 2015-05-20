%
% plot gating variables
%

V       = -120:1:100;
parameters              
VsK     = V-Kshift;
VsNa    = V-Nashift;

alpham  = 0.32*(13.1-VsNa)./(exp((13.1-VsNa)./4)-1);
betam   = 0.28*(VsNa-40.1)./(exp((VsNa-40.1)./5)-1);
taum    = 1./(alpham+betam);
minf    = alpham./(alpham+betam);

alphah  = 0.128*exp((17-VsNa)./18);
betah   = 4./(1+exp((40-VsNa)./5));
tauh    = 1./(alphah+betah);
hinf    = alphah./(alphah+betah);

alphan  = 0.032*(15-VsK)./(exp((15-VsK)./5)-1);
betan   = 0.5*exp((10-VsK)./40);
taun    = 1./(alphan+betan);
ninf    = alphan./(alphan+betan);

subplot(211);
p.minf = plot(V,minf,'-r'); hold all;
p.hinf = plot(V,hinf); set(p.hinf,'Color',[1 .5 .5])
p.ninf = plot(V,ninf,'-k');
set([p.minf p.hinf p.ninf],'LineWidth',2);
ylabel('x_{\rm inf}');
    
subplot(212);
p.minf = plot(V,taum,'-r'); hold all;
p.hinf = plot(V,tauh); set(p.hinf,'Color',[1 .5 .5])
p.ninf = plot(V,taun,'-k');
set([p.minf p.hinf p.ninf],'LineWidth',2);
xlabel('Voltage (mV)');
ylabel('\tau_{\rm inf} (ms)');

legend({'m','h','n'});