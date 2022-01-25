%%%%%%%%%%%%%% fMRI_Intrensic_Ignition_Measure
%%%%%%%%%%%% %%%%%%% Developed By: Gustavo Deco
%%%%%%% Modified and Adopted By: Rajanikant Panda 
%%%%%%% Date of Modification: 1st September 2019 
%%%%%%% Supervised: Gustavo Deco (UPF), Steven Laureys(Uliege), Jitka
%%%%%%% Annen(Uliege)  and Ane A López-González (UPF)

%Input Bold Time Series: 3D mat file in the format of No_os_Subj x ROIs x TimePoi
%output
    %mean(m) and variability(std) of ignition across the events of all subjects (GL = Group level)
       %'mevGL', 'stdevGL'
    %mean and variability of ignition across the events for each single subject (SS = single subject)
       %'mevokedintegSS', 'stdevokedintegSS',  
    % Spontanoues events for each subject
        %'SubjEvents'
    %mean and variability of ignition across nodes for each single subject (AN = across nodes)
        %'mignitionAN', 'stdignitionAN'
    % view this imagesc(phasematrix), imagesc(events)
clear
clc

CNT = load('F:\CSG\fMRI\FrontoParietalDOC\GT_DFC\Data_NxTxB\timeseries_all_CNT_04092019.mat');
CNT=CNT.ts_all; CNT=[CNT(1:22,:,:);CNT(24:32,:,:);CNT(34:35,:,:)]; CNT = num2cell(CNT,[2,3])';
for i=1:size(CNT,2)
    CNT{i} = squeeze(CNT{1,i});
end
%MCS = load('F:\CSG\fMRI\FrontoParietalDOC\GT_DFC\Data_NxTxB\New_150120202\data_MCS_all_2.mat')
MCS =load('F:\CSG\fMRI\FrontoParietalDOC\GT_DFC\IntensicIgnition\MCS.mat')
MCS=MCS.ts_all;  
MCS = num2cell(MCS,[2,3])';
for i=1:size(MCS,2)
    MCS{i} = squeeze(MCS{1,i});
end
%UWS = load('F:\CSG\fMRI\FrontoParietalDOC\GT_DFC\Data_NxTxB\New_150120202\data_UWS_all_2_fMRI.mat')
UWS=load('F:\CSG\fMRI\FrontoParietalDOC\GT_DFC\IntensicIgnition\UWS.mat')
UWS=UWS.ts_all; UWS = num2cell(UWS,[2,3])';
for i=1:size(UWS,2)
    UWS{i} = squeeze(UWS{1,i});
end

timeseries{1,1}=CNT;
timeseries{1,2}=MCS;
timeseries{1,3}=UWS;
%% %% Ignition Measures
for cond=1:3   %to run in Intrensic ignition measures for each group
    % count the number of subjects in each group/condition  
    numb_subj(1)=size(timeseries{1,1},2); %numb_subj(1)=size(timeseries{1,1},1);
    numb_subj(2)=size(timeseries{1,2},2); %numb_subj(2)=size(timeseries{1,2},1);
    numb_subj(3)=size(timeseries{1,3},2); %numb_subj(3)=size(timeseries{1,3},1);
    % To take data of one group/condition
     if cond==1
        TC=timeseries{1};
    elseif cond==2
        TC=timeseries{2};
    elseif cond==3
        TC=timeseries{3};
    end
  %% adapt parameters to YOUR data
  
  TR= 2; %sampling interval(TR)  
  %Tmax=217; %timepoints (computed later by size)
  N = size(TC{1},1); %246; %214;%758;%90;  regions
  NSUB = numb_subj(cond); %subjects
  Isubdiag = find(tril(ones(N),-1));
  nTRs = 5; % nTRs-1 == TRs to compute ignition after spontanous events (Ideal 5)
  Tmax=297; % 295;  % No of dynamics/time point
  CASE = cond;

  
%% %%%%%%%%%%% To filter the data
% Basic filtering parameters
  nevents = zeros(1,N);

  FC = zeros(NSUB,N,N);

  flp = 0.03;   %0.04     % lowpass frequency of filter
  fhi = 0.08;   %0.07     % highpass
  delt = TR;              % sampling interval
  k = 2;                  % 2nd order butterworth filter
  fnq = 1/(2*delt);       % Nyquist frequency
  Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
  [bfilt2,afilt2] = butter(k,Wn);   % construct the filter
   
  ttotal = 1;

  %% compute ignition for each subject
  for nsub = 1:NSUB 
      %xs=UWS{1}
    xs=TC{nsub};  %xs=squeeze(TC(nsub,:,5:end));     
    Tmax = size(xs,2); %time points
    T = 1:Tmax;  % T = 10:Tmax-10; 
    clear  x timeseriedata events Phases ev1 ev2
    
    % obtain events for each seed for individual subject (look - imagesc(events))
    for seed = 1:N
      x = demean(detrend(xs(seed,:)));
      timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
      Xanalytic = hilbert(demean(timeseriedata(seed,:)));
      tise = detrend(demean(timeseriedata(seed,T)));
      ev1 = tise>std(tise)+mean(tise);
      ev2 = [0 ev1(1:end-1)];
      events(seed,:) = (ev1-ev2)>0; 
      Phases(seed,:) = angle(Xanalytic);
      events_std(seed) = std(tise);
      events_mean_std(seed) = std(tise)+mean(tise);
    end
   SubjEvents{nsub} = events;
   SubjEvents_std(nsub,:) = events_std;
   SubjEvents_mean_std(nsub,:) = events_mean_std;
   
%% % Computing Integration for each subject (look- imagesc(phasematrix))
    
    %  Time to Phase transfermation using Hilbert Transfer
    %  T=1:1:size(Phases,2);

    %obtain 'events connectivity matrix' and integration value (integ)
    %for each time point
    
    
    %%%% to change DMN intigration for betwen network of FP (140-175)
    
    
    for t = T
      for i = 1:N %nodes from the DMN
        for j = 1:N %nodes dmn
          phasematrix(i,j)=exp(-3*adif(Phases(i,t),Phases(j,t)));
          %phasematrix(i,j) = events(i,t)*events(j,t); % phasematrix(i,j) = events(i,t-9)*events(j,t-9); 
        end
      end 
      cc = phasematrix; %*Cbin;
      cc = cc-eye(N);   %'events connectivity matrix'
%%% computing intigration values for whole time serise  
%       [comps csize] = get_components(cc);
%       integ(t) = max(csize)/N;  %  integ(t-9) = max(csize)/N;
  
%  % Sparsity based thresholding and event obtain for each regions 



  pp=1;
  PR=0:0.01:0.99;
  for p=PR
   A=abs(cc)>p;
   [comps csize]=get_components(A);
   cs(pp)=max(csize);
   %[M modularity]=community_louvain(A);    %%%% chaned for new intigration
   %cs(pp)=max(M);
   pp=pp+1;
  end
  integ(t)=sum(cs)*0.01/N;  %  integ(t-9)=sum(cs)*0.01/N;

  
  integt(ttotal)=sum(cs)*0.01/N;
  %[M modularity(t)]=community_louvain(phasematrix);  %  [M modularity(t-9)]=community_louvain(phasematrix);

 ttotal=ttotal+1;

end %end obtain integ
%%%%% 
 %% %
 %%%% event trigger
     %%%% to change DMN intigration for betwen network for FP 
     
    nevents2 = zeros(1,N);
    % save events and integration values for nTRs after the event
    for seed = 1:N %nodes from the FP
      flag = 0;
      for t = T
        %detect first event (nevents = matrix with 1xnode and number of events in each cell)
        if events(seed,t) == 1 && flag == 0  % if events(seed,t-9) == 1 && flag == 0
          flag = 1;
          %real events accumulated over subjects
          nevents(seed) = nevents(seed)+1;
          %real events for each subject
          nevents2(seed) = nevents2(seed)+1;
        end
        %save integration value for nTRs after the first event (nodesx(nTR-1)xevents)
        if flag > 0
          %integration accumulated over subjects
          IntegStim(seed,flag,nevents(seed)) = integ(t); %  IntegStim(seed,flag,nevents(seed)) = integ(t-9);
          %integration for each subject
          IntegStim2(seed,flag,nevents2(seed)) = integ(t); %  IntegStim2(seed,flag,nevents2(seed)) = integ(t-9);
          %IntegStimQ(seed,flag,nevents(seed))=modularity(t-9);
          flag = flag+1;
        end
        %after nTRs, set flag to 0 and wait for the next event (then, integ saved for (nTRs-1) events)
        if flag == nTRs
          flag = 0;
        end
       %Plot integration and events
       
       
      end
      index_events{seed}=find(events(seed,:)==1);
    end
%     plot(integt,'Linewidth',2)
% hold on
%     for seed=1:N
%  
% plot(index_events{seed},integt(index_events{seed}),'o')
%     end
  %mean and std of the max ignition in the nTRs for each subject and for each node  
    for seed = 1:N
      %mevokedinteg2OLD(seed)=max(mean(squeeze(IntegStim2(seed,:,1:nevents2(seed))),2));
      mevokedinteg2(seed) = mean(max(squeeze(IntegStim2(seed,:,1:nevents2(seed)))));
      stdevokedinteg2(seed) = std(max(squeeze(IntegStim2(seed,:,1:nevents2(seed)))));
    end
    
    %mean and std ignition across events for each subject in each node(Single Subject, SS)
    mevokedintegSS(:,:,nsub) = mevokedinteg2;
    stdevokedintegSS(:,:,nsub) = stdevokedinteg2;
    
    %mean and std ignition for each subject across nodes(AN)
    mignitionAN(nsub) = mean(mevokedinteg2);
    stdignitionAN(nsub) = std(mevokedinteg2);

       clearvars  mevGL stdevGL mevokedinteg % SubjEvents 
    
  end %end loop compute ignition for each subject
  
%   mignitionon=mean(mignitionon2);
%   stdignitionon=mean(stdignitionon2);

  %minteg=mean(integt);
  
  %mean and std ignition computed for all events in all subjects(Group level)
  for seed = 1:N
     mevokedinteg(seed,:) = mean(squeeze(IntegStim(seed,:,1:nevents(seed))),2);
     mevGL(seed) = mean(squeeze(max(squeeze(IntegStim(seed,:,1:nevents(seed))),[],1)));
     stdevGL(seed) = std(squeeze(max(squeeze(IntegStim(seed,:,1:nevents(seed))),[],1)));
     %mevokedQ(seed,:)=mean(squeeze(IntegStimQ(seed,:,1:nevents(seed))),2);
  end
  
  cd ('F:\CSG\fMRI\FrontoParietalDOC\GT_DFC\IntensicIgnition\ReCal_forEvent')

  save (sprintf('%dTR_SS_Ignition_CASE%d_Phases', nTRs-1, CASE), 'mevGL', 'stdevGL', 'mevokedinteg',...
      'mevokedintegSS', 'stdevokedintegSS', 'mignitionAN', 'stdignitionAN', 'SubjEvents', 'SubjEvents_std', 'SubjEvents_mean_std');
      % figure(cond)
      % mevokedinteg2OLD(seed)=max(mean(squeeze(IntegStim2(seed,:,1:nevents2(seed))),2));
      mevokedinteg2(seed) = mean(max(squeeze(IntegStim2(seed,:,1:nevents2(seed)))));
      stdevokedinteg2(seed) = std(max(squeeze(IntegStim2(seed,:,1:nevents2(seed)))));
 clearvars 'mevGL' 'stdevGL' 'mevokedinteg' 'mevokedintegSS' 'stdevokedintegSS' 'mignitionAN' 'stdignitionAN' 'SubjEvents' 'SubjEvents_std'
cd ..
end
   
