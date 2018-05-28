% MAIN CODE
clear;close all
plot=1;

%% load variables and subs
L_regg=get_L_variables; L_regg_ori=L_regg(1:38);L_regg_replic=L_regg(39:end);
econ_params=get_econ_params; 
normalized_baseline_preferences=get_baseline_preferences;
econ_params2=group_preferences(normalized_baseline_preferences);
RA_Regg=get_reggs;

% divide to appropriate subs
index_ori=find(~isnan(L_regg_ori));
index_replic=find(~isnan(L_regg_replic));
L_regg_ori=L_regg_ori(index_ori);L_regg_replic=L_regg_replic(index_replic);
index1=1:38;index2=39:81;

%% Lasso reggresion Original group 
%using all prams
[betas_LOO, finfo] = lasso([econ_params(:,(index_ori))' normalized_baseline_preferences(:,index_ori)' econ_params2(:,index_ori)' RA_Regg(index_ori)], L_regg_ori', 'alpha', 0.3); %Lasso reggresion RA
[~, Bindex]=min(finfo.MSE); betas_LOO=[finfo.Intercept(Bindex); betas_LOO(:,Bindex)]; %best fitting lambda
[betas_LOO_no_affil,finfo_no_affil] = lasso([econ_params([1:2 4:end],(index_ori))' normalized_baseline_preferences([1:2 4:end],index_ori)' econ_params2([1:2 4:end],index_ori)' RA_Regg(index_ori)], L_regg_ori', 'alpha', 0.3); %Lasso reggresion RA
[~, Bindex_no_affil]=min(finfo_no_affil.MSE); betas_LOO_no_affil=[finfo_no_affil.Intercept(Bindex_no_affil); betas_LOO_no_affil(:,Bindex_no_affil)]; %best fitting lambda


%using all except for RA
[betas_LOO2, finfo2] = lasso([econ_params(:,(index_ori))' normalized_baseline_preferences(:,index_ori)' econ_params2(:,index_ori)'], L_regg_ori', 'alpha', 0.3); %Lasso reggresion RA
[~, Bindex]=min(finfo2.MSE); betas_LOO2=[finfo2.Intercept(Bindex); betas_LOO2(:,Bindex)]; %best fitting lambda
[betas_LOO2_no_affil,finfo2_no_affil] = lasso([econ_params([1:2 4:end],(index_ori))' normalized_baseline_preferences([1:2 4:end],index_ori)' econ_params2([1:2 4:end],index_ori)'], L_regg_ori', 'alpha', 0.3); %Lasso reggresion RA
[~, Bindex_no_affil2]=min(finfo2_no_affil.MSE); betas_LOO2_no_affil=[finfo2_no_affil.Intercept(Bindex_no_affil2); betas_LOO2_no_affil(:,Bindex_no_affil2)]; %best fitting lambda

%% Lasso reggresion Replication group 
%using all prams
prediction_with_RA=betas_LOO(1)+([econ_params(:,(index_replic+38))' normalized_baseline_preferences(:,index_replic+38)' econ_params2(:,index_replic+38)' RA_Regg(index_replic+38)]*betas_LOO(2:end)); 
%reggression without ingroup afil for sub with missing data
prediction_with_RA(17)=betas_LOO_no_affil(1)+([econ_params([1:2 4:end],55)' normalized_baseline_preferences(:,55)' econ_params2(:,55)' RA_Regg(55)]*betas_LOO_no_affil((2:end))); 

%all except RA
prediction_w_no_RA=betas_LOO2(1)+([econ_params(:,(index_replic+38))' normalized_baseline_preferences(:,index_replic+38)' econ_params2(:,index_replic+38)' ]*betas_LOO2(2:end)); 
%reggression without ingroup afil for sub with missing data
prediction_w_no_RA(17)=betas_LOO2_no_affil(1)+([econ_params([1:2 4:end],55)' normalized_baseline_preferences(:,55)' econ_params2(:,55)' ]*betas_LOO2_no_affil((2:end))); 

%% Correlation predicted and observed 
%all prams
[r_using_RA,p_using_RA]=corr(prediction_with_RA,L_regg_replic','type', 'spearman', 'rows','complete'); 
%all except RA
[r_w_no_RA,p_w_no_RA]=corr(prediction_w_no_RA,L_regg_replic','type', 'spearman', 'rows','complete'); 

%directly compare models
%f-test of nested model
SSE_full_model=sum((prediction_with_RA-L_regg_replic').^2);
SSE_reduced_model=sum((prediction_w_no_RA-L_regg_replic').^2);
df1=1;df2=length(prediction_w_no_RA)-(length(betas_LOO2)+2);
f_value= ((SSE_reduced_model-SSE_full_model)/1)/(SSE_full_model/df2); %www.public.iastate.edu/~alicia/stat328/Multiple%20regression%20-%20nested%20models.pdf
p=1-fcdf(f_value,df1,df2);

%% plot
if plot==1
    figure;
    width=15;height=15;
    set(gcf,'units','centimeters','position',[0,0,width,height])

    fitresult=fit(prediction_with_RA,L_regg_replic','poly1');
    [ci,y] = predint(fitresult,prediction_with_RA,0.95,'functional');
    temp=[prediction_with_RA y ci];temp=sortrows(temp,1);
    [l,p] = boundedline_git(temp(:,1)', temp(:,2)', temp(:,3:4));
    set(l,'LineWidth',2);set(p,'FaceColor'); Marker_Color='b';
    %set(p,'FaceColor',[0.75,0.75,0.75]);  set(l,'color',[0.4,0.4,0.4]);Marker_Color=[0.5,0.5,0.5];xlim([-3.1 2]);

    scatter(prediction_with_RA,L_regg_replic,72,'filled','MarkerFaceColor','b','LineWidth',2);
    [r,p]=corr(prediction_with_RA,L_regg_replic','type', 'spearman','rows','complete');
    title(['rho= ' num2str(r), '  p= '  num2str(p)]);
end

% FUNCTIONS CODE

%% supp
%% extract data
function [B_mes,L_variables] = get_L_variables()
B_mes=[NaN(1,38),25,28,41,35,38,33,42,24,34,31,22,34,48,31,23,38,41,34,29,41,43,18,23,31,20,34,35,30,39,40,43,31,36,49,32,39,41,35,NaN,21,39,37,NaN]; 
B_mes(1,[2,4,6,8,11,13,17,19,21:25,27,32,34,36,38])=[27,37,29,34,30,21,45,28,29,23,38,38,37,41,43,35,38,35]; %half ori group was assigned to relationship oriented leadership Q - other half filled in LDBQ see SOM
L_variables={'LDBQ'};
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

% avg other group members risk prefrence
function [econ_params2] = group_preferences(normalized_baseline_preferences)
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
end

%% get appropriate RA, L, and subs
function [RA_Regg]=get_reggs()
for s=1:81
        file_name = fullfile('beh_data',sprintf('Result_matrix_sub_%d.mat',s));
        load(file_name); %col 13 is condition (1= group 2=self) col 1 is answer (1=risky 2= safe 3= defer)
        perc_D_Group=sum(Result_matrix(Result_matrix(:,13)==1,1)==3)/length(Result_matrix(Result_matrix(:,13)==1));
        perc_D_Self=sum(Result_matrix(Result_matrix(:,13)==2,1)==3)/length(Result_matrix(Result_matrix(:,13)==2));
        RA_Regg(s,1)=(perc_D_Group-perc_D_Self)/perc_D_Self;
        clear Result_matrix
end
end