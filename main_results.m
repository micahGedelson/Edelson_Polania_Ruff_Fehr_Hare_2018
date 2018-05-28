% MAIN CODE

clear; close all
%% general 
ori_replic=3;% 1=original group; 2=repli group; 3=all subs 
dependent_variable=1; %1=leadership; 2=RA 
use_FM=0;%use field measure (only with all subs)
test_control_taking=1;%Examine control taking in each condition separately 1=self 2= group test 0=dont use
normalize=1;
plot_w_CI=1;

%% load variables and subs
all_dependent_variables=get_L_variables;  all_dependent_variables2=all_dependent_variables;
econ_params=get_econ_params; 
normalized_baseline_preferences=get_baseline_preferences;
econ_params2=group_preferences(normalized_baseline_preferences);
[RA_Regg,RA_Regg2,L_regg, subs,AVG_RA,SD_RA,p_sign_RA,stats_sign_RA,perc_control_Self,perc_control_Group,AVG_L_Group,AVG_L_Self,stats_L_Group,stats_L_Self]=get_reggs(all_dependent_variables,ori_replic,dependent_variable);
  
%% create addi reggs
econ_params=econ_params(:,subs);
all_dependent_variables=all_dependent_variables(:,subs);
normalized_baseline_preferences=normalized_baseline_preferences(:,subs);
econ_params2=econ_params2(:,subs);
if ori_replic==3
    L_regg=[zscore(all_dependent_variables(1,subs<39),[],2) zscore(all_dependent_variables(2,subs>38),[],2)];
    RA_Regg(:,1)=[zscore(RA_Regg(subs<39)); RA_Regg(subs>38)];
else
    L_regg=L_regg(subs);
end

%% normalize prams? 
index_temp=find(subs==55);index_temp2=1:length(all_dependent_variables);index2=index_temp2;index2(index_temp)=[];% sub 55 in didnt fill in in-group affilation
if normalize==1
    econ_params([1:2 4:end],:)=zscore(econ_params([1:2 4:end],:),[],2);
    econ_params(3,index2)=zscore(econ_params(3,index2),[],2); %for in in-group affilation normalize without NaN
    if ori_replic ~=3;L_regg=zscore(L_regg);RA_Regg=zscore(RA_Regg);end 
end

%% reggresion leadership as dependent
if dependent_variable==1
    %parametric
    [b_regg_L,dev,stats_regg_L] =glmfit([econ_params([1:4 6:7],(index2))' econ_params(5,(index2))' normalized_baseline_preferences(:,index2)' econ_params2(:,index2)' RA_Regg(index2)],L_regg(index2));
    %non parmaetric
    [beta_regg_L,tstats_regg_L,F_stat_regg_L]=nonparamreg(L_regg(index2)',[econ_params([1:4 6:7],(index2))' econ_params(5,(index2))' normalized_baseline_preferences(:,index2)' econ_params2(:,index2)' RA_Regg(index2)]);

% corrolation with L Q
[r_L,p_L]=corr(RA_Regg,L_regg','type', 'spearman');
end


%%  reggresion RA as dependent
if dependent_variable==2
    %parametric 
    if (ori_replic==1 || ori_replic==3)
    index3=index2;index3(index3==28)=[]; %one extreme sub (very high RA==>3.6) has to be excluded since this is a parametric regg
    [b_regg_RA,dev,stats_regg_RA] =glmfit([econ_params([1:4 6:7],(index3))' econ_params(5,(index3))' normalized_baseline_preferences(:,index3)' econ_params2(:,index3)' ],RA_Regg(index3)');   
    [beta_regg_RA,tstats_regg_RA,F_stat_regg_RA]= nonparamreg(RA_Regg(index3),[econ_params([1:4 6:7],(index3))' econ_params(5,(index3))' normalized_baseline_preferences(:,index3)' econ_params2(:,index3)']);
    else    
    [b_regg_RA,dev,stats_regg_RA] = glmfit([econ_params([1:4 6:7],(index2))' econ_params(5,(index2))' normalized_baseline_preferences(:,index2)' econ_params2(:,index2)' ],RA_Regg(index2)');
    [beta_regg_RA,tstats_regg_RA,F_stat_regg_RA]= nonparamreg(RA_Regg(index2),[econ_params([1:4 6:7],(index2))' econ_params(5,(index2))' normalized_baseline_preferences(:,index2)' econ_params2(:,index2)']);
    end

%% corrolation with RA
[r_RA,p_RA]=corr([econ_params([1:4 6:7],index2)' econ_params(5,index2)' normalized_baseline_preferences(:,index2)' econ_params2(:,index2)'],RA_Regg(index2),'type', 'spearman'); 
[r_RA(3),p_RA(3)]=corr(econ_params(3,(index2))',RA_Regg(index2),'type', 'spearman'); %fill in the value for affi In group
end

%% RA with field measure(Not enough subs to combine measures- use raw RA)
if (ori_replic==3 && use_FM==1)
subs2=find(~isnan(all_dependent_variables2(3,:)));
L_regg2=all_dependent_variables2(3,subs2); 
[r_FM,p_FM]=corr(L_regg2',RA_Regg2(subs2),'type', 'spearman','rows','complete');
end

 
%% plot corrolation with L
if plot_w_CI==1
if dependent_variable==1
    figure;
    width=15;height=15;
    set(gcf,'units','centimeters','position',[0,0,width,height])
    
    
    %Use appropriate values for requested test 
    
    if use_FM==1;RA_Regg=zscore(RA_Regg2(subs2));L_regg=L_regg2;
    elseif test_control_taking==1;RA_Regg=zscore(perc_control_Self(subs))';
    elseif test_control_taking==2;RA_Regg=zscore(perc_control_Group(subs))';
    end
    

    fitresult=fit(RA_Regg,L_regg','poly1');
    [ci,y] = predint(fitresult,RA_Regg,0.95,'functional');
    temp=[RA_Regg y ci];temp=sortrows(temp,1);
    [l,p] = boundedline_git(temp(:,1)', temp(:,2)', temp(:,3:4));
    set(l,'LineWidth',2);set(p,'FaceColor'); Marker_Color='b';
    if test_control_taking~=0; set(p,'FaceColor',[0.75,0.75,0.75]);  set(l,'color',[0.4,0.4,0.4]);Marker_Color=[0.5,0.5,0.5];xlim([-3.1 2]);end 

    if (ori_replic==3 && use_FM~=1)
    scatter(RA_Regg(1:18),L_regg(1:18),72,'filled', '^','MarkerFaceColor',Marker_Color,'LineWidth',2);
    scatter(RA_Regg(19:end),L_regg(19:end),72,'s','filled','MarkerFaceColor',Marker_Color,'LineWidth',2);
    else
    scatter(RA_Regg,L_regg,72,'filled','MarkerFaceColor','b','LineWidth',2);
    end
    

    [r,p]=corr(RA_Regg,L_regg','type', 'spearman','rows','complete');
    title(['rho= ' num2str(r), '  p= '  num2str(p)]);
end    
end

% FUNCTIONS CODE

%% extract data
function [B_mes,L_variables] = get_L_variables()
B_mes=[NaN(1,38),25,28,41,35,38,33,42,24,34,31,22,34,48,31,23,38,41,34,29,41,43,18,23,31,20,34,35,30,39,40,43,31,36,49,32,39,41,35,NaN,21,39,37,NaN; 
NaN(1,38),10.3,8.3,31.5,30.3,28,10.5,37.3,13.3,21.5,8.8,2.5,24,26.3,16.5,6.5,23.5,33.3,30,13,36,37,1.5,12.5,10.5,4.75,16,24.3,15.3,32.8,36.3,35.3,13,24.3,41.50,20.8,31,36,17.8,NaN,3.5,29.3,22.8,NaN;
NaN(1,81)];
B_mes(1,[2,4,6,8,11,13,17,19,21:25,27,32,34,36,38])=[27,37,29,34,30,21,45,28,29,23,38,38,37,41,43,35,38,35]; %half ori group was assigned to relationship oriented leadership Q - other half filled in LDBQ see SOM
B_mes(3,[3,9,23,24,29,34,38,42,43,44,46,48,51,53,54,59,60,66,76,77,81])=[-0.78,-0.67,0.22,-0.012,1.27,1.27,2.24,1.6,0.34,-0.67,-0.67,-0.33,-0.65,-0.16,-0.67,-0.33,-0.67,-0.67,0.30,-0.95,0.30]; %normalized field measure available only for S with military rank or scout experience
L_variables={'LDBQ','composite','field measure normalized'};
end

function [econ_params,econ_names] = get_econ_params()
econ_params=[3.5,2,0,3,7,5,7,0,0,2,2,3,3.5,4,0,3.5,0,5,3,3.5,3,3.5,3,4,3.5,2,4,4.5,3,3.5,4,4,6,0,7,0,3.5,0,3.5,0,4,3.5,1,3,2,2,3.5,3.5,3.5,0,3,3,4,3,3,7,4,0,3.5,4,4.5,3.5,2,1,3,3.5,3,3.5,2,0,3,0,3,3.5,4,3.5,3,3,3,2,3.5
0,1,0,2,0,0,7,0,0,2,0,3,3.5,4,0,2,0,2,3,0,1,3.5,0,3,3.5,0,3,4.5,0,0,3,0,2,0,7,0,3.5,0,2,0,4,3.5,1,0,2,0,3.5,3.5,0,0,0,3,3,0,2,7,4,0,3.5,3,4.5,0,0,0,2,2,2,0,0,0,3,0,0,0,0,3.5,3,3,0,0,3.5
7,6,6,9,8,7,7,9,8,5,6,7.5,8,8,6,6,6,7,8,7,4,6,7,8,5,8,9,7,7,8,6,9,4,8,7,9,4,8,6,3,8,8,2,3,3.5,7,4,5,4,1,5,3,2,9,NaN,8,1,3,3,3,5,3,2,0,6,3,1,3,1,2,4,0,8,4,4,2,0.5,3,8,0,7
0,0,3,3,0,3,1,3,5,7,4,8.5,1,6,1,3,2,7,7,3,4,1,3,8,0,8,3,4,0,0,5,2,4,1,1,9,5,7,0,0,0,1,1,1,0,1,0,1,1,0,0,3,1,2,0,6,0,0,0,0,2,0,0,0,2,0,0,0,0,1,0,0,0,1,1,1,0.5,2,3,0,3
23,13,22,22,21,24,22,21,21,22,21,22,22,12,21,22,21,20,25,22,22,22,25,11,22,21,21,21,21,20,26,22,27,32,19,22,22,21,21,21,22,21,21,16,14,21,21,16,17,26,22,21,20,15,21,21,21,21,26,13,22,22,21,13,22,17,16,11,18,22,22,21,19,21,20,22,21,21,22,11,18
10,9.4,9.7,13.1,12.9,11.2,9.4,13.7,14,8.6,13.6,10.9,11.3,9.5,11.7,12,12.1,14.4,12.8,14.8,16,9.3,14.9,15.4,13,14.7,13.8,11.9,10.9,15.9,12.2,15.5,13.3,17.8,16.3,19.7,14.4,19.7,14.5,14.1,15.4,13.1,17.8,17.2,10.4,13.6,12.2,12.2,14.5,13.8,7.8,13.2,13.5,8.5,14.4,12,18.9,12.7,10,15.4,6.7,6.6,13.3,10.1,9.5,12.4,9.4,13,14.6,8.5,15.3,16.7,14.1,15.5,12.7,10,12.6,10.9,12.6,8.7,10.4;
15,10.5,12.1,13.7,8.4,7.9,6.5,8.9,10.7,11.1,11,12.4,12,11,12.3,15.9,15.9,9.4,13.6,14.6,11.5,10.9,14.1,13.7,15.7,16.2,18.7,18.5,9.6,11.5,10.8,18.6,11.4,18.4,19.6,20.9,13.9,10.2,11.2,10.8,17.4,15.7,11.3,15.3,14.9,12.4,12.3,9.5,13.5,10.5,14.7,14.5,11.7,12.8,8.3,12.3,11.6,17.3,11.1,12,8,10.6,8.1,6.5,7.8,12.9,11.3,11.8,11.4,8.7,11.5,12.9,16.8,13.9,14.7,15,13.8,7.9,11.5,13.8,15.6];
econ_names={'dictator ingroup','dictator outgroup','affiliation ingroup','affiliation outgroup','ambiguity indifference','outcome baseline','outcome delegation'};
end
function [normalized_baseline_preferences,addi_T1_params_names] = get_baseline_preferences()
normalized_baseline_preferences=[-0.92,-0.76,0.25,-0.48,-0.83,0.34,-0.34,-0.51,-0.81,-0.37,1.47,-0.48,-0.68,-1.09,-0.46,2.44,-1.05,-0.5,-0.38,3.18,-0.13,-0.56,-0.71,-0.42,-0.04,1.77,0.33,1.31,-0.85,0.12,-0.98,-0.23,0.58,0.24,-0.74,0,1.54,0.74,1.05,-1.25,-1.23,-0.29,-0.47,0.05,1.81,2.02,-0.38,-0.31,-0.47,-0.94,0.83,-1.29,0.25,-0.47,-1.04,-1.3,-0.94,2.12,0.21,1.13,-0.77,1.17,1.08,0.17,0.42,-0.21,-1.02,-1.24,-0.04,0.44,0.3,-0.7,1.23,0.44,-0.29,0.17,0.04,1.71,-1.52,-1.41,0.97;
-0.3,0.87,0.21,-0.39,0.77,0.33,-0.2,-0.47,0.08,2.05,-0.3,0.18,-0.28,-0.11,-0.69,0,1.63,-0.95,-1.46,2.12,-0.5,2.41,0.13,0.06,-0.09,-0.07,-0.49,-2.67,1.32,0.17,0.09,-0.02,0.74,-0.48,-0.41,-1.12,-1.31,-0.87,-0.49,-0.16,0.22,-0.05,-0.67,-1.3,-1.25,-0.98,0.02,0.4,-0.22,-0.02,-0.17,-0.2,-0.83,0.66,0.07,1.5,-0.21,-1.26,1.32,-0.6,2.44,0.83,0.31,-0.27,2.73,-0.51,-1.2,-0.51,0.04,2.76,-0.75,-0.31,-0.19,-1.11,-0.67,1.06,0.25,-0.74,-0.35,-0.66,1.06];
addi_T1_params_names={'risk','loss'};
end

% avg other group members risk preference
function [econ_params2] = group_preferences(normalized_baseline_preferences)
G{1}=(1:4); G{2}=(5:8); G{3}=(9:11); G{4}=(12:15); G{5}=(16:19); G{6}=(20:23); G{7}=(24:27); G{8}=(28:30); G{9}=(31:34); G{10}=(35:38); G{11}=(39:42); G{12}=(43:46); G{13}=(47:50); G{14}=(51:54); G{15}=(55:57); G{16}=(58:61); G{17}=(62:65); G{18}=(66:69); G{19}=(70:73);G{20}=(74:77); G{21}=(78:81); %subs in each group
tt=1;
for g=1:length(G)
    for ind=1:length(G{g})
        temp_G=G{g};temp_G(ind)=[];
        risk_others(tt)=mean(normalized_baseline_preferences(1,temp_G));
        loss_others(tt)=mean(normalized_baseline_preferences(2,temp_G));
        tt=tt+1;
    end
end
econ_params2=[risk_others;loss_others];
end

%% get appropriate RA, L, and subs
function [RA_Regg,RA_Regg2,L_regg, subs,AVG_RA,SD_RA,p_sign_RA,stats_sign_RA,perc_control_Self,perc_control_Group,AVG_L_Group,AVG_L_Self,stats_L_Group,stats_L_Self]=get_reggs(all_dependent_variables,ori_replic,dependent_variable)
if  ori_replic==1
    L_regg=all_dependent_variables(1,:);
elseif ori_replic==2
    L_regg=all_dependent_variables(2,:);
elseif ori_replic==3
    L_regg=[];
end

% define appropriate subs 
if dependent_variable==1 %leadership
    if ori_replic==1;subs=find(~isnan(all_dependent_variables(1,1:38)));
    elseif ori_replic==2;subs=find(~isnan(all_dependent_variables(2,39:end)));subs=subs+38;
    else
        subs=find(~isnan(all_dependent_variables(1,:)));
    end
else %RA
    if ori_replic==1;subs=1:38;
    elseif ori_replic==2;subs=39:81;
    else
        subs=1:81;
    end
end

bound_shift_normalized=[0.887,1.182,0.607,-0.445,-0.338,-0.266,-0.934,-0.098,-0.304,1.216,-0.864,0.324,1.778,-2.057,-1.42,-0.097,-0.771,0.582,0.788,0.443,-0.558,-0.023,-0.799,-1.279,0.591,0.393,-1.139,2.431,-0.136,1.887,-0.329,-0.976,0.264,0.438,1.448,-1.042,-0.328,-1.054,1.225,-1.228,0.502,-1.212,-0.789,1.988,-0.728,-0.098,-0.488,1.382,-0.188,0.285,0.02,1.315,0.256,-0.632,1.327,-0.152,1.123,-2.036,0.055,-0.315,1.331,-0.503,-0.498,0.36,0.882,0.004,-2.493,-1.231,-1.005,0.238,0.646,-0.527,0.084,-1.074,1.007,0.502,1.903,0.503,-0.126,-0.606,-0.948];
for s=1:81
        file_name = fullfile('beh_data',sprintf('Result_matrix_sub_%d.mat',s));
        load(file_name); %col 13 is condition (1= group 2=self) col 1 is answer (1=risky 2= safe 3= defer)
        perc_D_Group=sum(Result_matrix(Result_matrix(:,13)==1,1)==3)/length(Result_matrix(Result_matrix(:,13)==1));
        perc_control_Group(s)=1-perc_D_Group;
        perc_D_Self=sum(Result_matrix(Result_matrix(:,13)==2,1)==3)/length(Result_matrix(Result_matrix(:,13)==2));
        perc_control_Self(s)=1-perc_D_Self;
        RA_Regg(s,1)=(perc_D_Group-perc_D_Self)/perc_D_Self;
        clear Result_matrix
end
% Average values
AVG_RA=[];SD_RA=[];p_sign_RA=[];stats_sign_RA=[];AVG_L_Group=[];AVG_L_Self=[];stats_L_Group=[];stats_L_Self=[];
if (dependent_variable==2 && ori_replic==3)
    AVG_RA=mean(RA_Regg(subs));
    SD_RA=std(RA_Regg(subs));
    [p_sign_RA,h,stats_sign_RA] = signrank(RA_Regg(subs));
    
    AVG_L_Group=mean(perc_control_Group(subs));
    AVG_L_Self=mean(perc_control_Self(subs));
    [p_sign_L_Group,h,stats_L_Group] = signrank(perc_control_Group(subs),0.5);
    [p_sign_L_Self,h,stats_L_Self] = signrank(perc_control_Self(subs),0.5);

end   

RA_Regg2=RA_Regg;
%For replication group take bound size directly as RA
    RA_Regg(39:end,1)=bound_shift_normalized(39:end);

    RA_Regg=RA_Regg(subs);

    
end



