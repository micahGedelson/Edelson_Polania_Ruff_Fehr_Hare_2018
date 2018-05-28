clear;%close all
NOB=10; %number of bins

% sub loop
for s=(1:81)
    [choice_matrix]=get_reggs(s);    
    choice_matrix=sortrows(choice_matrix,3);%sort by SVd
        n_per_bin=length(choice_matrix)/NOB;
        RT_mat=reshape(choice_matrix(:,2),n_per_bin,[]);
        C_mat=reshape(choice_matrix(:,1),n_per_bin,[]);
        
 mean_RT_mat_per_sub(s,1:NOB)=mean(RT_mat);
 mean_C_mat_per_sub(s,1:NOB)=mean(C_mat);

end 
mean_RT=mean(mean_RT_mat_per_sub); SE_RT=std(mean_RT_mat_per_sub)/sqrt(length(mean_RT_mat_per_sub));
mean_C=mean(mean_C_mat_per_sub);   SE_C=std(mean_C_mat_per_sub)/sqrt(length(mean_C_mat_per_sub));

figure;errorbar(mean_C,SE_C,'color','r')
title('choice')
figure;errorbar(mean_RT,SE_RT,'color','r')
title('RT')

%% get appropriate params
function [choice_matrix]=get_reggs(s)
%get RT and choice
clear choice RT Svd_sub choice_matrix
file_name = fullfile('beh_data',sprintf('Result_matrix_sub_%d.mat',s));
load(file_name); %col 13 is condition (1= group 2=self) col 1 is answer (1=risky 2= safe 3= defer)

choice=Result_matrix(:,1);
choice(choice==1 | choice==2)=0; choice(choice==3)=1;% convert to 0/1;

RT=Result_matrix(:,4);

self_condition_index=find(Result_matrix(:,13)==2);

%get Svd
Svd = csvread('Svd_reduced.csv');
start_pos_sub=min(find(Svd(:,1)==s));end_pos_sub=max(find(Svd(:,1)==s));Svd_sub=Svd(start_pos_sub:end_pos_sub,2);%get SVd

% use self trials
choice=choice(self_condition_index);
RT=RT(self_condition_index);
Svd_sub=Svd_sub(self_condition_index);

choice_matrix=[choice RT Svd_sub];

if s==78; choice_matrix(61:64,:)=[];end %for sub with 124 trials remove 4 center trials to have 12 in each bin; sub was not used for article plot

clear Result_matrix
end