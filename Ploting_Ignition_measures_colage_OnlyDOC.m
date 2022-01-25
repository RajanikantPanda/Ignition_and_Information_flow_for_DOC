clc
clear
% load the data
CNT = load('E:\CSG\fMRI\FrontoParietalDOC\GT_DFC\IntensicIgnition\4TR_SS_Ignition_CASE1_Phases_WB1.mat');
MCS = load('E:\CSG\fMRI\FrontoParietalDOC\GT_DFC\IntensicIgnition\4TR_SS_Ignition_CASE2_Phases_WB1.mat');
UWS = load('E:\CSG\fMRI\FrontoParietalDOC\GT_DFC\IntensicIgnition\4TR_SS_Ignition_CASE3_Phases_WB1.mat');
load ('E:\CSG\fMRI\FrontoParietalDOC\GT_DFC\Metastability\DFCresult\DFC_DOC_plt');
load('E:\MyDell_Laptop\Desktop\DynEC\DynamicCommunicability_MG2\Communic_res\tau_results\BOLD_tau_DOC_0.01_0.15hz.mat');
BOLD_tau_c = BOLD_tau(1:33,:,:);  
BOLD_tau_m = BOLD_tau(34:59,:,:); 
BOLD_tau_u = BOLD_tau(60:73,:,:);
%%
% DOC mean DFC
figure(); 
% subplot(2,5,1);
% notBoxPlot(Grp_dFC_DOC,Grp_DOC,0.5,'patch',ones(length(Grp_dFC_DOC),1));
% %title('Whole Brain DFC', 'FontSize', 24);
% ylabel('Mean DFC', 'FontSize', 24);
% xticklabels({'Control','MCS','UWS'})
% set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
% Ignition of DOC
Grp_con = ones(length(CNT.mignitionAN),1); Grp_mcs = 2*ones(length(MCS.mignitionAN),1); Grp_uws = 3*ones(length(UWS.mignitionAN),1);
Grp = [Grp_con; Grp_mcs; Grp_uws]; %Grp = Gzarp';
clear Grp_con Grp_mcs Grp_uws
Grp_mignitionAN = [CNT.mignitionAN'; MCS.mignitionAN'; UWS.mignitionAN'];
% Bar plot of standard daviation of evoked ignition
b1= squeeze(mean(CNT.stdevokedintegSS,2));
b2= squeeze(mean(MCS.stdevokedintegSS,2));
b3= squeeze(mean(UWS.stdevokedintegSS,2));
Grp_mevokedinteg = [b1; b2; b3];
figure()%subplot(2,3,1);
notBoxPlot(Grp_mignitionAN,Grp,0.5,'patch',ones(length(Grp_mignitionAN),1));
%title('Ignition-Driven Mean Intigration (IDMI)', 'FontSize', 24);
ylabel('Mean IDMN', 'FontSize', 24);
xticklabels({'Control','MCS','UWS'})
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
figure()%subplot(2,3,2);
notBoxPlot(Grp_mevokedinteg,Grp,0.5,'patch',ones(length(Grp_mevokedinteg),1));
%title('Event Integration Variability', 'FontSize', 24);
xticklabels({'Control','MCS','UWS'})
ylabel('BOLD Event Variability', 'FontSize', 24);
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%%
figure()%subplot(2,3,4);
CI = 1.9; % 95% = 1.96 / 99% = 2.576 confidental Interval
a1= squeeze(mean(CNT.mevokedintegSS,3));
a1s= CI*squeeze(std(CNT.mevokedintegSS,[],3))/sqrt (size(CNT.mevokedintegSS,3));
a2= squeeze(mean(MCS.mevokedintegSS,3));
a2s= CI*squeeze(std(MCS.mevokedintegSS,[],3))/sqrt (size(MCS.mevokedintegSS,3));
a3= squeeze(mean(UWS.mevokedintegSS,3));
a3s= CI*squeeze(std(UWS.mevokedintegSS,[],3))/sqrt (size(UWS.mevokedintegSS,3));
errorbar (a1,a1s, 'k', 'LineWidth',2); grid on;
hold on
errorbar (a2,a2s, 'b', 'LineWidth',2); grid on;
errorbar (a3,a3s, 'r', 'LineWidth',2); grid on;
%title("Brain Regional Ignition-Driven integration", 'FontSize', 18);
legend('Control', 'MCS', 'UWS', 'FontWeight', 'bold')
xlabel('No of ROIs', 'FontSize', 24);
ylabel('IDMN', 'FontSize', 24);
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')

%% % Taw
figure()%subplot(2,3,3);
BOLD_tau_l1= [squeeze(BOLD_tau_c(:,:,2));squeeze(BOLD_tau_m(:,:,2));squeeze(BOLD_tau_u(:,:,2))];
BOLD_tau_l1_M = mean(BOLD_tau_l1,2);
Grp=[ones(size(BOLD_tau_c,1),1);2*ones(size(BOLD_tau_m,1),1);3*ones(size(BOLD_tau_u,1),1)];
notBoxPlot(BOLD_tau_l1_M,Grp,0.5,'patch',ones(length(BOLD_tau_l1_M),1));
xticklabels({'Control','MCS','UWS'})
ylabel('Tau', 'FontSize', 24);
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')




%% %%% Secontime rearenged colage 31-05-2021
%%% Mean intrinsic ignition
figure(); subplot(1,2,1);
Grp = [ones(length(CNT.mignitionAN),1); 3*ones(length(MCS.mignitionAN),1); 2*ones(length(UWS.mignitionAN),1)];
Grp_mignitionAN = [CNT.mignitionAN'; MCS.mignitionAN'; UWS.mignitionAN'];
notBoxPlot(Grp_mignitionAN,Grp,0.5,'patch',ones(length(Grp_mignitionAN),1));
%title('Ignition-Driven Mean Intigration (IDMI)', 'FontSize', 24);
ylabel('Mean intrinsic ignition', 'FontSize', 24);
xticklabels({'Control','UWS','MCS'})
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')

%%% Temporal Bold responce (Taw)
%figure(); 
subplot(1,2,2);
BOLD_tau_l1= [squeeze(BOLD_tau_c(:,:,2));squeeze(BOLD_tau_m(:,:,2));squeeze(BOLD_tau_u(:,:,2))];
BOLD_tau_l1_M = mean(BOLD_tau_l1,2);
Grp=[ones(size(BOLD_tau_c,1),1);3*ones(size(BOLD_tau_m,1),1);2*ones(size(BOLD_tau_u,1),1)];
notBoxPlot(BOLD_tau_l1_M,Grp,0.5,'patch',ones(length(BOLD_tau_l1_M),1));
xticklabels({'Control','UWS','MCS'})
ylabel('Mean memory depth', 'FontSize', 24);
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')

%% stats of figure 1- mean IDMI
clc
c=1:33; m=34:59; u=60:73;
[h_ig_cu,p_ig_cu,~,stats_ig_cu] =ttest2(Grp_mignitionAN(c),Grp_mignitionAN(u))
[h_ig_cm,p_ig_cm,~,stats_ig_cm] =ttest2(Grp_mignitionAN(c),Grp_mignitionAN(m))
[h_ig_mu,p_ig_mu,~,stats_ig_mu] =ttest2(Grp_mignitionAN(m),Grp_mignitionAN(u),0.016,'right')

%% stats of figure 1- mean Tau
clc
[h_tau_cu,p_tau_cu,~,stats_tau_cu] =ttest2(BOLD_tau_l1_M(c),BOLD_tau_l1_M(u))
[h_tau_cm,p_tau_cm,~,stats_tau_cm] =ttest2(BOLD_tau_l1_M(c),BOLD_tau_l1_M(m))
[h_tau_mu,p_tau_mu,~,stats_tau_mu] =ttest2(BOLD_tau_l1_M(m),BOLD_tau_l1_M(u),0.016)

%% Even MEan and verience
clc
CNT_Event = CNT.SubjEvents_std;
MCS_Event = MCS.SubjEvents_std;
UWS_Event = UWS.SubjEvents_std;

CNT_Event_mean = mean(mean(CNT_Event,2))
CNT_Event_std = std(mean(CNT_Event,2))

MCS_Event_mean = mean(mean(MCS_Event,2))
MCS_Event_std = std(mean(MCS_Event,2))

UWS_Event_mean = mean(mean(UWS_Event,2))
UWS_Event_std = std(mean(UWS_Event,2))

[h_event_cu,p_event_cu,~,stats_event_cu] =ttest2(mean(CNT_Event,2),mean(UWS_Event,2))
[h_event_cm,p_event_cm,~,stats_event_cm] =ttest2(mean(CNT_Event,2),mean(MCS_Event,2))
[h_event_mu,p_event_mu,~,stats_event_mu] =ttest2(mean(MCS_Event,2),mean(UWS_Event,2))

figure(); 
Grp = [ones(length(mean(CNT_Event,2)),1); 3*ones(length(mean(MCS_Event,2)),1); 2*ones(length(mean(UWS_Event,2)),1)];
Grp_mean_evnt= [mean(CNT_Event,2); mean(MCS_Event,2); mean(UWS_Event,2)];
notBoxPlot(Grp_mean_evnt,Grp,0.5,'patch',ones(length(Grp_mean_evnt),1));
ylabel('Mean BOLD Events', 'FontSize', 24);
xticklabels({'Control','UWS','MCS'})