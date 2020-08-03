% This analysis code was written by Annika Boldt, 2015-2019
% Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
% modulates exploration and exploitation in value-based learning. Neuroscience
% of Consciousness, 2019(1), 1?12. https://doi.org/10.1093/nc/niz004

% This script fits three reinforcement learning (RL) models based on a
% Random Search Algorithm (RSA)

clear all
close all
clc

load('') % put your data path and the '.mat' file here (you can use the output from the SMC here)

allblocks = [9, 19, 7, 14]; % blocks we wish to fit the data to
allsubs = 1:21;   % participant list
npp = length(allsubs); % number of participants
nblocks = length(allblocks); % number of blocks

% lower and upper bounds for the parameters
lrate_list = [0 1];
confmu_list = [0 1.2];
confsig_list = [0 0.4];
confm_list = [0 0.04];
confn_list = [0 0.7];
confq_list = [0 0.002];

npars = 6;
num_psets = 5000; % number of parameter sets

param_set_list=[];

costfun = inline('(1./sqrt(2*pi)).*exp(1).^(-1./2*(z-zhat).^2)');

% RANDOM NUMBER GENERATOR
resetrn = sum(100*clock);
rand('state',resetrn);

% create the list of parameter sets
for k=1:num_psets
    if length(lrate_list)>1
        lrate = (lrate_list(2)-lrate_list(1)).*rand + lrate_list(1);
    else
        lrate = lrate_list;
    end
    if length(confmu_list)>1
        confmu = (confmu_list(2)-confmu_list(1)).*rand + confmu_list(1);
    else
        confmu = confmu_list;
    end
    if length(confsig_list)>1
        confsig = (confsig_list(2)-confsig_list(1)).*rand + confsig_list(1);
    else
        confsig = confmu_list;
    end    
    if length(confm_list)>1
        confm = (confm_list(2)-confm_list(1)).*rand + confm_list(1);
    else
        confm = confm_list;
    end    
    if length(confn_list)>1
        confn = (confn_list(2)-confn_list(1)).*rand + confn_list(1);
    else
        confn = confn_list;
    end  
    if length(confq_list)>1
        confq = (confq_list(2)-confq_list(1)).*rand + confq_list(1);
    else
        confq = confq_list;
    end       
    c_psetting = [lrate confmu confsig confm confn confq];
    param_set_list = [param_set_list; c_psetting];
end



cost_blockarmtrial = NaN(3,num_psets,npp);
for ipar=1:num_psets
    lrate = param_set_list(ipar,1);
    confmu = param_set_list(ipar,2);
    confsig = param_set_list(ipar,3);
    confm = param_set_list(ipar,4);
    confn = param_set_list(ipar,5);
    confq = param_set_list(ipar,6);
    for isub=1:npp
        cost_armtrial = NaN(3,nblocks);
        for iblock=1:nblocks
            newvalestimates = NaN(2,size(themodeldata(isub).block(iblock).rating,1)+1);
            newvalestimates(:,1) = 50;
            newconf = NaN(3,2,size(themodeldata(isub).block(iblock).rating,1));
            trating=NaN(2,length(newvalestimates)-1);
            trating(1,themodeldata(isub).block(iblock).whichobs==1) = themodeldata(isub).block(iblock).rating(themodeldata(isub).block(iblock).whichobs==1);
            trating(2,themodeldata(isub).block(iblock).whichobs==2) = themodeldata(isub).block(iblock).rating(themodeldata(isub).block(iblock).whichobs==2);
            tvalconf=NaN(2,size(newconf,3));
            tvalconf(1,themodeldata(isub).block(iblock).whichobs==1) = themodeldata(isub).block(iblock).valcj(themodeldata(isub).block(iblock).whichobs==1);
            tvalconf(2,themodeldata(isub).block(iblock).whichobs==2) = themodeldata(isub).block(iblock).valcj(themodeldata(isub).block(iblock).whichobs==2);
            ntrials = size(themodeldata(isub).block(iblock).rating,1);
            for itrial=1:ntrials
                if themodeldata(isub).block(iblock).whichobs(itrial)==1 || themodeldata(isub).block(iblock).resp(itrial)==1
                    newvalestimates(1,itrial+1) = newvalestimates(1,itrial) + lrate * (themodeldata(isub).block(iblock).outcome(itrial) - newvalestimates(1,itrial));
                    newvalestimates(2,itrial+1) = newvalestimates(2,itrial);
                elseif themodeldata(isub).block(iblock).whichobs(itrial)==2 || themodeldata(isub).block(iblock).resp(itrial)==2
                    newvalestimates(1,itrial+1) = newvalestimates(1,itrial);
                    newvalestimates(2,itrial+1) = newvalestimates(2,itrial) + lrate * (themodeldata(isub).block(iblock).outcome(itrial) - newvalestimates(2,itrial));
                end
                newconf(1,:,itrial) = normrnd(confmu,confsig,2,1);
                newconf(2,:,itrial) = confm * itrial + confn;
                newconf(3,:,itrial) = confq * trating(:,itrial).^2 - confq * 100 * trating(:,itrial) + confq * 10000 * 0.25;
            end
            cost_trial = NaN(3,2);
            for imod=1:3
                for iarm=1:2
                    sim_data = [newvalestimates(iarm,2:end) (squeeze(newconf(imod,iarm,:))*100)'];
                    emp_data = [trating(iarm,:) tvalconf(iarm,:)*100];
                    sim_data = sim_data(~isnan(emp_data));
                    emp_data = emp_data(~isnan(emp_data));
                    tcost = costfun(emp_data,sim_data);
                    tcost(tcost==0) = realmin;
                    cost_trial(imod,iarm) = 1/length(emp_data) * sum(-2 * log(tcost));
                end
                cost_armtrial(imod,iblock) = 1/2 * sum(cost_trial(imod,:));
            end
        end
        for imod=1:3
            cost_blockarmtrial(imod,ipar,isub) = 1/nblocks * sum(cost_armtrial(imod,:));
        end
    end
end

disp('best_group_parameters:minimum cost function (value and confidence 1)')
[i1a,i2a]=min(mean(squeeze(cost_blockarmtrial(1,:,:)),2))
disp('groupfit:the best parameters based on cost function was (value and confidence 1)')
param_set_list(i2a,[1 2:3])

disp('best_group_parameters:minimum cost function (value and confidence 2)')
[i1b,i2b]=min(mean(squeeze(cost_blockarmtrial(2,:,:)),2))
disp('groupfit:the best parameters based on cost function was (value and confidence 2)')
param_set_list(i2b,[1 4:5])

disp('best_group_parameters:minimum cost function (value and confidence 3)')
[i1c,i2c]=min(mean(squeeze(cost_blockarmtrial(3,:,:)),2))
disp('groupfit:the best parameters based on cost function was (value and confidence 3)')
param_set_list(i2c,[1 6])

the_filename=['RLmodels_' datestr(now) '_sim.mat'];
save(the_filename,'cost_blockarmtrial','param_set_list');