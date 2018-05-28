% MAIN CODE
clear;close all

folder_name= fullfile('mri_data','DCM'); 

%set up
dependent_variable=1; %1=leadership; 2=RA
normalize=1;
calc_exp_var=1; 
LOO=1; %leave one out
stepwise=1;
plot_figs=1;



%% load variables and subs
dependent_variables = get_D_variables;
econ_params=get_econ_params;
normalized_baseline_preferences=get_baseline_preferences;
[econ_params2,normalized_baseline_preferences]=group_preferences(normalized_baseline_preferences);
reggression_params=[econ_params; normalized_baseline_preferences; econ_params2];
DV_mat=dependent_variables(dependent_variable,:);
subs=find(~isnan(DV_mat));
if dependent_variable==1; DV_mat=DV_mat(subs); end

%% normalize prams?
index_temp=find(subs==17);index2=subs;index2(index_temp)=[];% sub 17 in didnt fill in in-group affilation

if normalize==1
    econ_params([1:2 4:end],:)=zscore(econ_params([1:2 4:end],:),[],2);
    econ_params(3,index2)=zscore(econ_params(3,index2),[],2); %for in in-group affilation normalize without NaN
    reggression_params=[econ_params; normalized_baseline_preferences; econ_params2];
end

%% subs loop
for s=1:length(subs)
    runs=(1:4);if subs(s)==5; runs=(2:4);elseif subs(s)==29; runs=(1:3); elseif subs(s)==30; runs=[1:2 4];end %Runs with repeated movement over 3 mm
   
    
    %% variance is explained by the model 
    if calc_exp_var==1
        for r=runs
           
            load(fullfile(folder_name,sprintf('DCM_file_sub_%02d_run%d',subs(s),r)));
            for roi=1:4
                PSS   = sum(sum(DCM.y(:,roi).^2));
                RSS   = sum(sum(DCM.R(:,roi).^2));
                var_exp(r,roi)  = 100*PSS/(PSS + RSS);
            end
        end
        AVG_ex_var_sub(s,1:4)=mean(var_exp);AVG_var_exp_sub_all_roi(s)=mean(AVG_ex_var_sub(s,:));
        var_exp=[];
    end
    
    %% load DCM prams
    A_mat=zeros(4,4);B_mat=zeros(4,4,10);C_mat=zeros(4,10);
    load(fullfile(folder_name,sprintf('DCM_avg_sub_%02d',subs(s))));
    %% All possible connections
    Connections(s,:)=reshape(DCM.Ep.A(eye(length(DCM.Ep.A))'~=1),1,12).'; %First 3 are connections from ROI1 (to ROI 2 3 ect), 3:6 is ROI2 etc
    
    %% inputs
    SEV(s)=DCM.Ep.C(1,10);Gr(s)=DCM.Ep.C(2,1);self(s)=DCM.Ep.C(2,4);
    PD_G(s)=DCM.Ep.C(3,3);PD_S(s)=DCM.Ep.C(3,6);
    lead(s)=DCM.Ep.C(4,7);defer(s)=DCM.Ep.C(4,8);
    inputs(s,:)=[SEV(s); Gr(s); self(s); PD_G(s); PD_S(s); lead(s); defer(s)];
    
    %% Modulations of experimental factors on connections
    %group
    Gr_on_conn_temp=DCM.Ep.B(:,:,1);
    Gr_on_conn(s,:)=reshape(Gr_on_conn_temp(eye(length(Gr_on_conn_temp))'~=1),1,12).';
    %self
    S_on_conn_temp=DCM.Ep.B(:,:,4);
    S_on_conn(s,:)=reshape(S_on_conn_temp(eye(length(S_on_conn_temp))'~=1),1,12).';
    
    %Differences in modulations
    Gr_vs_S_on_conn(s,:)=Gr_on_conn(s,:)- S_on_conn(s,:);
    
    %% Possible nonlinear Connections (Multiplied by relevant main effect)
    nonlinear_temp_for_regg(s,:)=[DCM.Ep.D(1,3,2); DCM.Ep.D(3,1,2);DCM.Ep.D(1,3,4); DCM.Ep.D(3,1,4)];
    ba21_on_ins_to_mpfc(s)=DCM.Ep.D(1,3,2)*(Gr(s)-self(s)); ba21_on_mpfc_to_ins(s)=DCM.Ep.D(3,1,2)*(Gr(s)-self(s));
    TPJ_on_ins_to_mpfc(s)=DCM.Ep.D(1,3,4)*(defer(s)-lead(s)); TPJ_on_mpfc_to_ins(s)=DCM.Ep.D(3,1,4)*(defer(s)-lead(s));
    nonlinear(s,:)=[ba21_on_ins_to_mpfc(s); ba21_on_mpfc_to_ins(s); TPJ_on_ins_to_mpfc(s); TPJ_on_mpfc_to_ins(s)];
    DCM=[];
end

%% var exp
if calc_exp_var==1;mean_var_exp_percent=mean(AVG_var_exp_sub_all_roi);end


%% Regression
Total=[ones(size(Connections,1),1) Connections Gr_vs_S_on_conn inputs nonlinear nonlinear_temp_for_regg];
if normalize==1
    Total(:,2:end)=zscore(Total(:,2:end));
end
% Regression with RA
[b2,dev2,stats2] = glmfit(Total(:,2:end),DV_mat');

% Holms MC correction
real_p_value3=stats2.p;  t_values=stats2.t;
for ii=1:length(stats2.p)
    holms_p_value3(ii)= 0.01/(length(stats2.p)-ii+1);
end
real_p_value4=[real_p_value3'; 1:length(real_p_value3);t_values';stats2.beta']'; %col 1 real p; 2 pram num;3 tvalue; 4 beta;5 holms  p value -if col 1 is lower then 5, this pram is sig
real_p_value4=sortrows(real_p_value4,1);real_p_value4(:,5)=holms_p_value3';
sig_holms2=find (real_p_value4(:,1)<holms_p_value3');

%% leave-one-out for RA using DCM prams only
if (LOO==1 && dependent_variable==2)
    for s2=1:length(subs)
        Total_temp=Total;DV_mat_temp=DV_mat;
        Total_temp(s2,:)=[]; DV_mat_temp(s2)=[]; %remove target subject 
        Total_temp(:,1)=[]; %remove constant
        [betas_LOO, finfo2] = lasso(Total_temp, DV_mat_temp, 'alpha', 0.3); %Lasso reggresion for all the group except the current sub
        [~, Bindex]=min(finfo2.MSE);
        betas_LOO=[finfo2.Intercept(Bindex); betas_LOO(:,Bindex)]; %get betas 
        f= [1 Total(s2,2:end)]; %real values for sub 
        SP.predk(s2)=f*betas_LOO; % multiply these values by the weights 
        SP.predbin(s2)=SP.predk(s2) > median(DV_mat); %  above or below the median?
        SP.truep(s2)=DV_mat(s2) > median(DV_mat); %is the real value above the median
        SP.crt(s2)=SP.predbin(s2)==SP.truep(s2);%was this prediction correct
    end
    
    [r_lasso,p_lasso]=corr(SP.predk',DV_mat','type','Spearman');
    
    mean_bin_pred=mean(SP.crt);[h,p_ttest_bin_pred]=ttest(double(SP.crt),0.5);SE_bin_pred=std(SP.crt)/sqrt(length(SP.crt));
end

%% leave-one-out for leadership using imag+beh prams
if (LOO==1 && dependent_variable==1)
    for s2=1:length(subs)
        
        if s2==17; reggression_params_temp=reggression_params(:,[1:2 4:end]);else;reggression_params_temp=reggression_params;end
        
        Total_temp=[Total reggression_params_temp(:,subs)' dependent_variables(2,subs)'];
        DV_mat_temp=DV_mat;
        Total_temp(s2,:)=[]; DV_mat_temp(s2)=[]; 
        Total_temp(:,1)=[]; %remove constant
        [betas_LOO, finfo2] = lasso(Total_temp, DV_mat_temp, 'alpha', 0.3); %Lasso reggresion for all the group except the current sub
        [~, Bindex]=min(finfo2.MSE); 
        betas_LOO=[finfo2.Intercept(Bindex); betas_LOO(:,Bindex)]; %get betas 
        f= [1 Total(s2,2:end) reggression_params_temp(:,s2)' dependent_variables(2,s2)']; %values of the new fmri subject 
        SP.predk(s2)=f*betas_LOO; % multiply by the weights
        SP.predbin(s2)=SP.predk(s2) > median(DV_mat); % predicting above or below the median
        SP.truep(s2)=DV_mat(s2) > median(DV_mat); %is the real value above the median
        SP.crt(s2)=SP.predbin(s2)==SP.truep(s2);%was prediction correct?
    end
    
    [r_lasso,p_lasso]=corr(SP.predk',DV_mat','type','Spearman');
    
    mean_bin_pred=mean(SP.crt);[h,p_ttest_bin_pred]=ttest(double(SP.crt),0.5);SE_bin_pred=std(SP.crt)/sqrt(length(SP.crt));
    diff_from_prdic=SP.predk'-DV_mat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%  
    if plot_figs==1
        figure;
        width=15;height=15;
        set(gcf,'units','centimeters','position',[0,0,width,height])
        
        fitresult=fit(zscore(SP.predk)',zscore(DV_mat)','poly1');
        [ci,y] = predint(fitresult,zscore(SP.predk),0.95,'functional');
        temp=[zscore(SP.predk)' y ci];temp=sortrows(temp,1);
        [l,p] = boundedline_git(temp(:,1)', temp(:,2)', temp(:,3:4));
        set(l,'LineWidth',2);set(p,'FaceColor',[0,1,0]);  set(l,'color',[0,0.7,0]);Marker_Color=[0,0.8,0];alpha(p,0.2); %xlim([-2.2 2.5]);
        
        scatter(zscore(SP.predk),zscore(DV_mat),72,'filled','MarkerFaceColor',Marker_Color,'LineWidth',2);

        [r,p]=corr(zscore(SP.predk)',zscore(DV_mat)','type', 'spearman','rows','complete');
        title(['rho= ' num2str(r), '  p= '  num2str(p)]);
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%
%% stepwise reggression predicting the Q data from the CV / CV+brain
if (stepwise==1 && dependent_variable==1)

    %% lasso
    %with brain
    Total_temp= Total; Total_temp(:,1)=[];%remove constant
    reggression_params=reggression_params(:,subs)'; RA_regg=dependent_variables(2,subs)';
    [betas_LOO1, finfo1] = lasso([Total_temp reggression_params RA_regg], DV_mat, 'alpha', 0.3); %Lasso reggresion
    [~, Bindex2]=min(finfo1.MSE);
    betas_LOO1=[finfo1.Intercept(Bindex2); betas_LOO1(:,Bindex2)];
    
    %calc BIC/AIC
    BIC_brain_and_beh_lasso= (length(subs)*log(min(finfo1.MSE))) + (length(betas_LOO1) *log(length(subs))); 
    AIC_brain_and_beh_lasso=  (length(subs)*log(min(finfo1.MSE))) + 2*(length(betas_LOO1)); 
    
    %without brain
    Total_temp=[RA_regg reggression_params];
    [betas_LOO2, finfo2] = lasso(Total_temp, DV_mat, 'alpha', 0.3); %Lasso reggresion
    [~, Bindex2]=min(finfo2.MSE);
    betas_LOO2=[finfo2.Intercept(Bindex2); betas_LOO2(:,Bindex2)]; %add intersept so we can predict using a reggrasion that has the constant as well
    
    BIC_only_beh_lasso= (length(subs)*log(min(finfo2.MSE))) + (length(betas_LOO2) *log(length(subs))); %penalize for all
    AIC_only_beh_lasso=  (length(subs)*log(min(finfo2.MSE))) + 2*(length(betas_LOO2)); %penalize for all
    
    %diff
    BIC_diff=BIC_brain_and_beh_lasso-BIC_only_beh_lasso;
    AIC_diff=AIC_brain_and_beh_lasso-AIC_only_beh_lasso;
end

% FUNCTIONS CODE 

%% extract data
function [dependent_variables,D_variables] = get_D_variables()
dependent_variables=[10.3,8.3,31.5,30.3,28,10.5,37.3,13.3,21.5,8.8,2.5,24,26.3,16.5,6.5,23.5,33.3,30,13,36,37,1.5,12.5,10.5,4.75,16,24.3,15.3,32.8,36.3,35.3,13,24.3,41.50,20.8,31,36,17.8,NaN,3.5,29.3,22.8,NaN
    1.225,-1.228,0.502,-1.212,-0.789,1.988,-0.728,-0.0980,-0.488,1.382,-0.188,0.285,0.0200,1.315,0.256,-0.632,1.327,-0.152,1.123,-2.036,0.0550,-0.315,1.331,-0.503,-0.498,0.360,0.882,0.00400,-2.493,-1.231,-1.005,0.238,0.646,-0.527,0.0840,-1.074,1.007,0.502,1.903,0.503,-0.126,-0.606,-0.948];
D_variables={'composite','bound_shift_normalized'};
end


function [econ_params,econ_names] = get_econ_params()
econ_params=[3.5,2,0,3,7,5,7,0,0,2,2,3,3.5,4,0,3.5,0,5,3,3.5,3,3.5,3,4,3.5,2,4,4.5,3,3.5,4,4,6,0,7,0,3.5,0,3.5,0,4,3.5,1,3,2,2,3.5,3.5,3.5,0,3,3,4,3,3,7,4,0,3.5,4,4.5,3.5,2,1,3,3.5,3,3.5,2,0,3,0,3,3.5,4,3.5,3,3,3,2,3.5
    0,1,0,2,0,0,7,0,0,2,0,3,3.5,4,0,2,0,2,3,0,1,3.5,0,3,3.5,0,3,4.5,0,0,3,0,2,0,7,0,3.5,0,2,0,4,3.5,1,0,2,0,3.5,3.5,0,0,0,3,3,0,2,7,4,0,3.5,3,4.5,0,0,0,2,2,2,0,0,0,3,0,0,0,0,3.5,3,3,0,0,3.5
    7,6,6,9,8,7,7,9,8,5,6,7.5,8,8,6,6,6,7,8,7,4,6,7,8,5,8,9,7,7,8,6,9,4,8,7,9,4,8,6,3,8,8,2,3,3.5,7,4,5,4,1,5,3,2,9,NaN,8,1,3,3,3,5,3,2,0,6,3,1,3,1,2,4,0,8,4,4,2,0.5,3,8,0,7
    0,0,3,3,0,3,1,3,5,7,4,8.5,1,6,1,3,2,7,7,3,4,1,3,8,0,8,3,4,0,0,5,2,4,1,1,9,5,7,0,0,0,1,1,1,0,1,0,1,1,0,0,3,1,2,0,6,0,0,0,0,2,0,0,0,2,0,0,0,0,1,0,0,0,1,1,1,0.5,2,3,0,3
    23,13,22,22,21,24,22,21,21,22,21,22,22,12,21,22,21,20,25,22,22,22,25,11,22,21,21,21,21,20,26,22,27,32,19,22,22,21,21,21,22,21,21,16,14,21,21,16,17,26,22,21,20,15,21,21,21,21,26,13,22,22,21,13,22,17,16,11,18,22,22,21,19,21,20,22,21,21,22,11,18
    10,9.4,9.7,13.1,12.9,11.2,9.4,13.7,14,8.6,13.6,10.9,11.3,9.5,11.7,12,12.1,14.4,12.8,14.8,16,9.3,14.9,15.4,13,14.7,13.8,11.9,10.9,15.9,12.2,15.5,13.3,17.8,16.3,19.7,14.4,19.7,14.5,14.1,15.4,13.1,17.8,17.2,10.4,13.6,12.2,12.2,14.5,13.8,7.8,13.2,13.5,8.5,14.4,12,18.9,12.7,10,15.4,6.7,6.6,13.3,10.1,9.5,12.4,9.4,13,14.6,8.5,15.3,16.7,14.1,15.5,12.7,10,12.6,10.9,12.6,8.7,10.4
    15,10.5,12.1,13.7,8.4,7.9,6.5,8.9,10.7,11.1,11,12.4,12,11,12.3,15.9,15.9,9.4,13.6,14.6,11.5,10.9,14.1,13.7,15.7,16.2,18.7,18.5,9.6,11.5,10.8,18.6,11.4,18.4,19.6,20.9,13.9,10.2,11.2,10.8,17.4,15.7,11.3,15.3,14.9,12.4,12.3,9.5,13.5,10.5,14.7,14.5,11.7,12.8,8.3,12.3,11.6,17.3,11.1,12,8,10.6,8.1,6.5,7.8,12.9,11.3,11.8,11.4,8.7,11.5,12.9,16.8,13.9,14.7,15,13.8,7.9,11.5,13.8,15.6];
econ_params=econ_params(:,39:end);
econ_names={'dictator ingroup','dictator outgroup','affiliation ingroup','affiliation outgroup','ambiguity indifference','outcome baseline','outcome delegation'};
end

function [normalized_baseline_preferences,addi_T1_params_names] = get_baseline_preferences()
normalized_baseline_preferences=[-0.92,-0.76,0.25,-0.48,-0.83,0.34,-0.34,-0.51,-0.81,-0.37,1.47,-0.48,-0.68,-1.09,-0.46,2.44,-1.05,-0.5,-0.38,3.18,-0.13,-0.56,-0.71,-0.42,-0.04,1.77,0.33,1.31,-0.85,0.12,-0.98,-0.23,0.58,0.24,-0.74,0,1.54,0.74,1.05,-1.25,-1.23,-0.29,-0.47,0.05,1.81,2.02,-0.38,-0.31,-0.47,-0.94,0.83,-1.29,0.25,-0.47,-1.04,-1.3,-0.94,2.12,0.21,1.13,-0.77,1.17,1.08,0.17,0.42,-0.21,-1.02,-1.24,-0.04,0.44,0.3,-0.7,1.23,0.44,-0.29,0.17,0.04,1.71,-1.52,-1.41,0.97
    -0.3,0.87,0.21,-0.39,0.77,0.33,-0.2,-0.47,0.08,2.05,-0.3,0.18,-0.28,-0.11,-0.69,0,1.63,-0.95,-1.46,2.12,-0.5,2.41,0.13,0.06,-0.09,-0.07,-0.49,-2.67,1.32,0.17,0.09,-0.02,0.74,-0.48,-0.41,-1.12,-1.31,-0.87,-0.49,-0.16,0.22,-0.05,-0.67,-1.3,-1.25,-0.98,0.02,0.4,-0.22,-0.02,-0.17,-0.2,-0.83,0.66,0.07,1.5,-0.21,-1.26,1.32,-0.6,2.44,0.83,0.31,-0.27,2.73,-0.51,-1.2,-0.51,0.04,2.76,-0.75,-0.31,-0.19,-1.11,-0.67,1.06,0.25,-0.74,-0.35,-0.66,1.06];
addi_T1_params_names={'risk','loss'};
end

% avg other group members risk prefrence
function [econ_params2,normalized_baseline_preferences] = group_preferences(normalized_baseline_preferences)
G{1}=(1:4); G{2}=(5:8); G{3}=(9:11); G{4}=(12:15); G{5}=(16:19); G{6}=(20:23); G{7}=(24:27); G{8}=(28:30); G{9}=(31:34); G{10}=(35:38); G{11}=(39:42); G{12}=(43:46); G{13}=(47:50); G{14}=(51:54); G{15}=(55:57); G{16}=(58:61); G{17}=(62:65); G{18}=(66:69); G{19}=(70:73);G{20}=(74:77); G{21}=(78:81); %subs in each group
tt=1;
for g=1:length(G)
    for ind=1:length(G{g})
        temp_G=G{g};temp_G(ind)=[];
        risk_othres(tt)=mean(normalized_baseline_preferences(1,temp_G));
        loss_othres(tt)=mean(normalized_baseline_preferences(2,temp_G));
        tt=tt+1;
    end
end
econ_params2=[risk_othres;loss_othres];
econ_params2=econ_params2(:,39:end);
normalized_baseline_preferences=normalized_baseline_preferences(:,39:end);
end
