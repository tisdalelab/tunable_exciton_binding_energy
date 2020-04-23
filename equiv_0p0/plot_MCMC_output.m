function [] = plot_MCMC_output(xdata,ydata,pfull,logP,call_model,func_struct,burnin,sample,production_run)

%% Extract parameters
% Time progression
final_time = size(pfull,3);
time = 1:1:final_time;
Nwalkers = size(pfull,2);

%% Burn-in plot

figure();
for i = 1:size(logP,2)
    logPwalker = squeeze(logP(2,i,:));
    plot(time,logPwalker)
    hold on
end
h = gca;
if production_run
    plot(round(burnin*final_time)*[1 1],h.YLim,'--k')
    titlespec = 'Production run';
else
    titlespec = 'Preliminary run';
end   
xlabel('t (every 10 steps per walkers)')
ylabel('Log-likelihood')
title(titlespec)

%% Fit to data

[~,MLEind] = max(logP(2,:,end));

% Plot data
figure(); plot(xdata,ydata,'or');
hold on

% Plot last position of all walkers
for i = 1:Nwalkers
    p = pfull(:,i,end);
    yfit = call_model(p,xdata);
    plot(xdata,yfit,'-','Color',[0 0 1 0.15],'LineWidth',2)
end

% Plot best performer
pMLE = pfull(:,MLEind,end);
yMLE = call_model(pMLE,xdata);
plot(xdata,yMLE,'-c','LineWidth',2)
xlabel('Energy (eV)')
ylabel('Absorbance (a.u.)')
title(titlespec)

% Plot comprising functions for the MLE
if sample == 3  
    
    fitL = func_struct.first_ex_L(pMLE,xdata);
    fitGL = func_struct.first_ex_GL(pMLE,xdata);
    fitHx2 = func_struct.Hx2(pMLE,xdata);
    fitHx3 = func_struct.Hx3(pMLE,xdata);
    fitHx4 = func_struct.Hx4(pMLE,xdata);
    fitS = func_struct.BAND_step(pMLE,xdata);
    fitfine1 = func_struct.BAND_fine1(pMLE,xdata);
    fitfine2 = func_struct.BAND_fine2(pMLE,xdata);
    fitfine3 = func_struct.BAND_fine3(pMLE,xdata);
    
    plot(xdata,fitL,'--k')
    plot(xdata,fitGL,':k')
    plot(xdata,fitHx2,'--','Color',[1 0.5 0])
    plot(xdata,fitHx3,'-.','Color',[1 0.5 0])
    plot(xdata,fitHx4,':','Color',[1 0.5 0])
    plot(xdata,fitS,'--b')
    plot(xdata,fitfine1,'--m')
    plot(xdata,fitfine2,'-.m')
    plot(xdata,fitfine3,':m')
    
elseif sample == 4
    fitGL = func_struct.doped_1s(pMLE,xdata);
    fitHx2 = func_struct.Hx2(pMLE,xdata);
    fitHx3 = func_struct.Hx3(pMLE,xdata);
    fitHx4 = func_struct.Hx4(pMLE,xdata);
    fitS = func_struct.BAND_step(pMLE,xdata);
    fitfine1 = func_struct.BAND_fine1(pMLE,xdata);
    fitfine2 = func_struct.BAND_fine2(pMLE,xdata);
    fitfine3 = func_struct.BAND_fine3(pMLE,xdata);
    
    plot(xdata,fitGL,'--k')
    plot(xdata,fitHx2,'--','Color',[1 0.5 0])
    plot(xdata,fitHx3,'-.','Color',[1 0.5 0])
    plot(xdata,fitHx4,':','Color',[1 0.5 0])
    plot(xdata,fitS,'--b')
    plot(xdata,fitfine1,'--m')
    plot(xdata,fitfine2,'-.m')
    plot(xdata,fitfine3,':m')
    
elseif sample == 5
    
    fitGL = func_struct.doped_1s(pMLE,xdata);
    fitHx2 = func_struct.Hx2(pMLE,xdata);
    fitHx3 = func_struct.Hx3(pMLE,xdata);
    fitHx4 = func_struct.Hx4(pMLE,xdata);
    fitS = func_struct.BAND_step(pMLE,xdata);
    fitfine1 = func_struct.BAND_fine1(pMLE,xdata);
    fitfine2 = func_struct.BAND_fine2(pMLE,xdata);
    fitfine3 = func_struct.BAND_fine3(pMLE,xdata);
    
    plot(xdata,fitGL,'--k')
    plot(xdata,fitHx2,'--','Color',[1 0.5 0])
    plot(xdata,fitHx3,'-.','Color',[1 0.5 0])
    plot(xdata,fitHx4,':','Color',[1 0.5 0])
    plot(xdata,fitS,'--b')
    plot(xdata,fitfine1,'--m')
    plot(xdata,fitfine2,'-.m')
    plot(xdata,fitfine3,':m')
    
end

end

