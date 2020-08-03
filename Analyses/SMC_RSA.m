% This analysis code was written by Annika Boldt and Charles Blundell, 2015-2019
% With thanks to Dharshan Kumaran for an earlier version of the SMC script.
% Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
% modulates exploration and exploitation in value-based learning. Neuroscience
% of Consciousness, 2019(1), 1?12. https://doi.org/10.1093/nc/niz004

function SMC_RSA(SMCopt)

% The function performs Sequential Monte Carlo (SMC) fitting based on a
% Random Search Algorithm (RSA) and should be called with a struct as
% argument, which has the with the following fields:

% SMCopt.dpath: specifies the '.mat' file containing the empirical data
% SMCopt.curr_exp: whether we currently fit Exp1 or 2
% SMCopt.costf: which cost function to use (1 or 2)
% SMCopt.noise: which noise to use (1 = gaussian, 2 = cauchy)
% SMCopt.curr_subs: a vector containing all participant numbers
% SMCopt.curr_blocks: a vector containing the blocks we wish to fit the data to
% SMCopt.num_runs: the number of runs (here: 1)
% SMCopt.n_part: number of particles (e.g. 10000)
% SMCopt.var_ratio_list: lower and upper bound for the variance ratio (i.e. init_var/innov_var)
% SMCopt.likelihoodnoise_list: lower and upper bound for the likelihoodnoise
% SMCopt.var_part_init_list: lower and upper bound for the particle variance
% SMCopt.alpha_list: lower and upper bound for alpha
% SMCopt.save_var: 1 = save var to .mat file
% SMCopt.N: number of parameter combinations to try

load(SMCopt.dpath)

fn = fieldnames(datanow);
for i=1:numel(fn)
    tdatanow.(fn{i}) = datanow.(fn{i})(ismember(datanow.block,SMCopt.curr_blocks));
end
datanow=tdatanow;
clear tdatanow

fn = fieldnames(datanow);
for i=1:numel(fn)
    tdatanow.(fn{i}) = datanow.(fn{i})(ismember(datanow.sub,SMCopt.curr_subs));
end
datanow=tdatanow;
clear tdatanow


% RANDOM NUMBER GENERATOR
resetrn = sum(100*clock);
rand('state',resetrn);

try
    gcp();
catch
    warning('Couldn''t create parallel pools.');
end

num_runs=SMCopt.num_runs;
param.num_runs=num_runs;

sub_ind=unique(datanow.sub)';
block_ind=unique(datanow.block)';
n_part=SMCopt.n_part; %number of particles

% Generate parameter space
var_ratio_list = SMCopt.var_ratio_list; %i.e. init_var/innov_var
likelihoodnoise_list = SMCopt.likelihoodnoise_list;
var_part_init_list = SMCopt.var_part_init_list;
alpha_list = SMCopt.alpha_list;

%save?
save_var = SMCopt.save_var; %1= save to .mat file

num_subs=length(sub_ind);
num_block=length(block_ind);
n_states=2; %number of items/ hidden states

likelihoodfun = inline('(1./sqrt(2*pi*R)).*exp(1).^(-(z-zhat).^2./(2*R))');

the_filename=[datestr(now) '_sim.mat'];


%__________________________________________________
param=([]);
param.n_part=n_part;
param.n_states=n_states;
param.likelihoodfun=likelihoodfun;

%---------------------------
param.var_ratio_list=var_ratio_list;
param.likelihoodnoise_list=likelihoodnoise_list;
param.var_part_init_list=var_part_init_list;
param.alpha_list=alpha_list;
param.N = SMCopt.N;

param_set_list=[];

for k=1:param.N
    if length(var_ratio_list)>1
        c_var_ratio = (var_ratio_list(2)-var_ratio_list(1)).*rand + var_ratio_list(1);
    else
        c_var_ratio = var_ratio_list;
    end
    if length(likelihoodnoise_list)>1
        c_likelihoodnoise = (likelihoodnoise_list(2)-likelihoodnoise_list(1)).*rand + likelihoodnoise_list(1);
    else
        c_likelihoodnoise = likelihoodnoise_list;
    end
    if length(var_part_init_list)>1
        c_var_part_init = (var_part_init_list(2)-var_part_init_list(1)).*rand + var_part_init_list(1);
    else
        c_var_part_init = var_part_init_list;
    end
    if length(alpha_list)>1
        c_alpha = (alpha_list(2)-alpha_list(1)).*rand + alpha_list(1);
    else
        c_alpha = alpha_list;
    end
    c_psetting = [c_var_ratio c_likelihoodnoise c_var_part_init c_alpha];
    param_set_list = [param_set_list; c_psetting];
end

param.param_set_list=param_set_list; % list of the parameter combinations
num_psets=size(param_set_list,1); % num parameter combinations to be tested

num_modelruns=num_subs*num_block*num_psets*num_runs;

themodeldata=([]);
the_new_log_L=NaN(num_psets,num_subs);
% load trialorder and data
for curr_pset_num=1:num_psets % PARAMETER SET
    
    fprintf('parameter set number: %s \n',num2str(curr_pset_num));
    fprintf('num parameter sets is: %s \n',num2str(num_psets));
        
    var_ratio=param_set_list(curr_pset_num,1);
    likelihoodnoise=param_set_list(curr_pset_num,2);
    var_part_init=param_set_list(curr_pset_num,3);
    alpha=param_set_list(curr_pset_num,4);
    
    
    var_innov=var_part_init*var_ratio; % define innovation variance as a ratio of initial variance
    
    tmp_n_eff = NaN(num_subs,1);
    
    
    
    respartrunblocksub = NaN(num_subs,1);
    for c_sub=1:num_subs    
        curr_sub=sub_ind(c_sub);
      
        resam = 0;
        respartrunblock = NaN(num_block,1);
        for c_block=1:num_block
            curr_block=block_ind(c_block);
            
            allrun_ppinf_logL=zeros(1,num_runs);
            respartrun = NaN(num_runs,1);
            for rr=1:num_runs
                curr_run=rr;

                
                outcome=datanow.outcome(datanow.sub==curr_sub & datanow.block==curr_block);
                rating=datanow.rating(datanow.sub==curr_sub & datanow.block==curr_block);
                valcj=datanow.valcj(datanow.sub==curr_sub & datanow.block==curr_block);
                whichobs=datanow.whichobs(datanow.sub==curr_sub & datanow.block==curr_block);
                value1=datanow.value1(datanow.sub==curr_sub & datanow.block==curr_block);
                value2=datanow.value2(datanow.sub==curr_sub & datanow.block==curr_block);
                resp=datanow.resp(datanow.sub==curr_sub & datanow.block==curr_block);

                n_iter=length(datanow.rating(datanow.sub==curr_sub & datanow.block==curr_block));
                numtrials=n_iter;
                
                % for the observed data
                themodeldata(c_sub).block(c_block).outcome=outcome;
                themodeldata(c_sub).block(c_block).resp=resp;
                themodeldata(c_sub).block(c_block).rating=rating;
                themodeldata(c_sub).block(c_block).whichobs=whichobs;
                themodeldata(c_sub).block(c_block).value1=value1;
                themodeldata(c_sub).block(c_block).value2=value2;
                themodeldata(c_sub).block(c_block).valcj=valcj;
                %%%%%%%%%
                
                
                tr_count={};
                tr_count{1}=0; % rating trials
                tr_count{2}=0; % choice trials
                
                
                %initialize particles - means, weights, likelihood
                %ESS, resampling MATRICES
                
                part=zeros(n_iter+1,n_part,n_states); % n_iter+1 because the first timestep = initial state
                part_postdiff=zeros(n_iter+1,n_part,n_states); % n_iter+1 because the first timestep = initial state
                norm_weight_log=log(ones(n_iter+1, n_part)/n_part); % normalized
                norm_weight_beforeresampling=log(ones(n_iter+1, n_part)/n_part);
                weight=(ones(n_iter+1, n_part)/n_part);
                weight_log=(ones(n_iter+1, n_part)/n_part); % unnormalized
                norm_weight=(ones(n_iter+1, n_part)/n_part); % normalized
                
                mu_part_init=50; % initial mean of particles
                
                if SMCopt.noise==1
                    part(1,:,:) = mu_part_init + randn(1,n_part, n_states)*sqrt(var_part_init);
                else 
                    part(1,:,:) = mu_part_init + cauchyrnd(0,sqrt(var_part_init),1,n_part,n_states);
                end
                
                n_eff=zeros(1,n_iter+1);
                
                meanlikelihood = zeros(1,n_iter);
                
                
                themodeldata(c_sub).block(c_block).meanpart(1,1) = mean(part(1,:,1),2);
                themodeldata(c_sub).block(c_block).meanpart(2,1) = mean(part(1,:,2),2);
                themodeldata(c_sub).block(c_block).medianpart(1,1) = median(part(1,:,1),2);
                themodeldata(c_sub).block(c_block).medianpart(2,1) = median(part(1,:,2),2);
                themodeldata(c_sub).block(c_block).upperpart(1,1) = quantile(part(1,:,1),0.75);
                themodeldata(c_sub).block(c_block).upperpart(2,1) = quantile(part(1,:,2),0.75);
                themodeldata(c_sub).block(c_block).lowerpart(1,1) = quantile(part(1,:,1),0.25);
                themodeldata(c_sub).block(c_block).lowerpart(2,1) = quantile(part(1,:,2),0.25);   
                themodeldata(c_sub).block(c_block).meanlikelihood(1) = mean(norm_weight_beforeresampling(1,:));
                themodeldata(c_sub).block(c_block).n_eff(1)=tmp_n_eff(c_sub); 
                themodeldata(c_sub).block(c_block).varpart(1,1) = var(part(1,:,1));
                themodeldata(c_sub).block(c_block).varpart(2,1) = var(part(1,:,2));
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                
                for itrial=1:n_iter
                    
                    if whichobs(itrial)>0
                        curr_ppinftype=1; % i.e. this is a rating trial
                        tr_count{1}=tr_count{1}+1;
                    else
                        curr_ppinftype=2;% i.e. this is a choice trial
                        tr_count{2}=tr_count{2}+1;
                    end
                    
                    
                    % SMC: propagate particles forward using state equation
                    % only for RATING TRIALS (SMC frozen in choice trial)     
                    if SMCopt.noise==1
                        part(itrial+1,:,:) = part(itrial,:,:) + randn(1,n_part, n_states)*sqrt(var_innov);  %+ innov noise
                    else
                        part(itrial+1,:,:) = part(itrial,:,:) + cauchyrnd(0,sqrt(var_innov),1,n_part,n_states);
                    end
                    part_postdiff(itrial,:,:)=part(itrial+1,:,:);  
                    
                          

                    
                    % Choice models (particles have to be saved before any
                    % transformations have been made)
                    
                    switch curr_ppinftype % This is where the particle update happens.
                        case 1 %rating
                            
                            likelihood = likelihoodfun(likelihoodnoise,outcome(itrial),part(itrial,:,whichobs(itrial)));
                            likelihood(part(itrial,:,whichobs(itrial))<0 | part(itrial,:,whichobs(itrial))>100) = 0; % bounded likelihood
                            meanlikelihood(itrial) = mean(likelihood);
                            
                            weight_log(itrial+1,:)=(norm_weight_log(itrial,:))+log(likelihood); % use log normalized wts for update
                            weight(itrial+1,:)=exp(weight_log(itrial+1,:)); % unlogged unnormalized wts - use this for normalization and ESS calculation below
                            % normalize
                            norm_weight(itrial+1,:)=(weight(itrial+1,:))/((sum(weight(itrial+1,:)))); % unlogged
                            norm_weight_beforeresampling(itrial+1,:)=(weight(itrial+1,:))/((sum(weight(itrial+1,:)))); % unlogged
                            norm_weight_log(itrial+1,:)=log(norm_weight(itrial+1,:));%
                            % effective sample size
                            tmp_n_eff(c_sub)=1/sum(norm_weight(itrial+1,:).^2);
                            n_eff(itrial+1)=tmp_n_eff(c_sub);

                            if sum(weight(itrial+1,:))
                                aa=datasample(1:n_part, n_part, 'Weights', (weight(itrial+1,:))); % sampling with replacement given unnormalized weights equation 9.
                                part(itrial+1, :, :) = part(itrial+1, aa, :); % now replace current particles with the sample particles themselves.
                                themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).resampleworked(itrial) = 1;
                            else
                                themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).resampleworked(itrial) = 0;
                                resam = resam + 1;
                            end
                            norm_weight(itrial+1,:)=(1/n_part);
                            norm_weight_log(itrial+1,:)=log(norm_weight(itrial+1,:)); % logged
                            
                        case 2 % choice

                            likelihood = likelihoodfun(likelihoodnoise,outcome(itrial),part(itrial,:,resp(itrial)));
                            likelihood(part(itrial,:,resp(itrial))<0 | part(itrial,:,resp(itrial))>100) = 0; % bounded likelihood
                            
                            weight_log(itrial+1,:)=(norm_weight_log(itrial,:))+log(likelihood); % use log normalized wts for update
                            weight(itrial+1,:)=exp(weight_log(itrial+1,:)); % unlogged unnormalized wts - use this for normalization and ESS calculation below
                            % normalize
                            norm_weight(itrial+1,:)=(weight(itrial+1,:))/((sum(weight(itrial+1,:)))); % unlogged
                            norm_weight_beforeresampling(itrial+1,:)=(weight(itrial+1,:))/((sum(weight(itrial+1,:)))); % unlogged
                            norm_weight_log(itrial+1,:)=log(norm_weight(itrial+1,:));
                            % effective sample size
                            tmp_n_eff(c_sub)=1/sum(norm_weight(itrial+1,:).^2);
                            n_eff(itrial+1)=tmp_n_eff(c_sub);

                            if sum(weight(itrial+1,:))
                                aa=datasample(1:n_part, n_part, 'Weights', (weight(itrial+1,:))); % sampling with replacement given unnormalized weights
                                part(itrial+1, :, :) = part(itrial+1, aa, :); % now replace current particles with the sample particles themselves.
                                themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).resampleworked(itrial) = 1;
                            else
                                themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).resampleworked(itrial) = 0;
                                resam = resam + 1;
                            end
                            norm_weight(itrial+1,:)=(1/n_part);
                            norm_weight_log(itrial+1,:)=log(norm_weight(itrial+1,:)); % logged
 
                    end                 
                    
                    themodeldata(c_sub).block(c_block).meanpart(1,itrial+1) = mean(part(itrial+1,:,1),2);
                    themodeldata(c_sub).block(c_block).meanpart(2,itrial+1) = mean(part(itrial+1,:,2),2);
                    themodeldata(c_sub).block(c_block).medianpart(1,itrial+1) = median(part(itrial+1,:,1),2);
                    themodeldata(c_sub).block(c_block).medianpart(2,itrial+1) = median(part(itrial+1,:,2),2);
                    themodeldata(c_sub).block(c_block).upperpart(1,itrial+1) = quantile(part(itrial+1,:,1),0.75);
                    themodeldata(c_sub).block(c_block).upperpart(2,itrial+1) = quantile(part(itrial+1,:,2),0.75);
                    themodeldata(c_sub).block(c_block).lowerpart(1,itrial+1) = quantile(part(itrial+1,:,1),0.25);
                    themodeldata(c_sub).block(c_block).lowerpart(2,itrial+1) = quantile(part(itrial+1,:,2),0.25);
                    themodeldata(c_sub).block(c_block).meanlikelihood(itrial+1) = mean(norm_weight_beforeresampling(itrial+1,:));
                    themodeldata(c_sub).block(c_block).n_eff(itrial+1)=tmp_n_eff(c_sub); 
                    themodeldata(c_sub).block(c_block).varpart(1,itrial+1) = var(part(itrial+1,:,1));
                    themodeldata(c_sub).block(c_block).varpart(2,itrial+1) = var(part(itrial+1,:,2));
                    
                end
                
                
                if SMCopt.costf==1
                    
                    idx = sub2ind(size(part), sort(repmat(find(whichobs>0)',1,n_part)), ...
                                              repmat(1:n_part,1,length(whichobs(whichobs>0))), ...
                                              sort(repmat(whichobs(whichobs>0)',1,n_part)));
                    tpart = reshape(part(idx),n_part,length(whichobs(whichobs>0)))';
                    res = NaN(n_part,1);                    
                    for ipart=1:n_part
                        res(ipart) = norm_weight_beforeresampling(itrial+1,ipart) * prod((1./sqrt(2*pi*valcj(whichobs>0)*alpha)).*exp(1).^(-(rating(whichobs>0)-tpart(ipart)).^2./(2*valcj(whichobs>0)*alpha)));
                    end
                    respartrun(rr) = 1/num_runs * sum(res);      
                    
                elseif SMCopt.costf==2
                    
                    res = NaN(length(whichobs(whichobs>0)),1);
                    ttrials = 1:n_iter;
                    ttrials = ttrials(whichobs>0);
                    idx = sub2ind(size(part), sort(repmat(find(whichobs>0)',1,n_part)), ...
                                              repmat(1:n_part,1,length(whichobs(whichobs>0))), ...
                                              sort(repmat(whichobs(whichobs>0)',1,n_part)));
                    tpart = reshape(part(idx),n_part,length(whichobs(whichobs>0)))';
                    for it=1:sum(whichobs>0)
                        tweights = norm_weight_beforeresampling(ttrials(it)+1,:);
                        tvalcj = valcj(ttrials(it));
                        trating = rating(ttrials(it));
                        tp = tpart(it,:);
                        res(it) = sum((tweights * (1./sqrt(2*pi*tvalcj*alpha)).*exp(1).^(-(trating-tp).^2./(2*tvalcj*alpha))));
                    end
                    themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedres = zeros(1,length(res));
                    themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedres(res==0) = 1;
                    res(res==0) = realmin; % padding
                    respartrun(rr) = 1/num_runs * prod(res);
                    if respartrun(rr)==0
                        respartrun(rr) = realmin; % padding
                        themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedrespartrun = 1;
                    else
                        themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedrespartrun = 0;
                    end
                    
                else
                    
                    % use resampled particles instead
                    res = NaN(length(whichobs(whichobs>0)),1);
                    ttrials = 1:n_iter;
                    ttrials = ttrials(whichobs>0);
                    idx = sub2ind(size(part), sort(repmat((find(whichobs>0)+1)',1,n_part)), ... % increased find to refer to itrial+1 particles
                                              repmat(1:n_part,1,length(whichobs(whichobs>0))), ...
                                              sort(repmat(whichobs(whichobs>0)',1,n_part)));
                    tpart = reshape(part(idx),n_part,length(whichobs(whichobs>0)))';
                    for it=1:sum(whichobs>0)
                        tweights = norm_weight_beforeresampling(ttrials(it)+1,:);
                        tvalcj = valcj(ttrials(it));
                        trating = rating(ttrials(it));
                        tp = tpart(it,:);
                        res(it) = sum((tweights * (1./sqrt(2*pi*tvalcj*alpha)).*exp(1).^(-(trating-tp).^2./(2*tvalcj*alpha))));
                    end
                    themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedres = zeros(1,length(res));
                    themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedres(res==0) = 1;
                    res(res==0) = realmin; % padding
                    respartrun(rr) = 1/num_runs * prod(res);
                    if respartrun(rr)==0
                        respartrun(rr) = realmin; % padding
                        themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedrespartrun = 1;
                    else
                        themodeldata(c_sub).block(c_block).run(rr).pset(curr_pset_num).correctedrespartrun = 0;
                    end

                    

                end
                
            end
            respartrunblock(c_block) = 1/num_block * sum(respartrun);    
        
        end
        the_new_log_L(curr_pset_num,c_sub) = -log(sum(respartrunblock));
        if resam
            the_new_log_L(curr_pset_num,c_sub) = NaN;
        end      
        
    end

if save_var==1
    save(the_filename,'themodeldata','param','the_new_log_L','SMCopt');
end

end


for c_sub=1:num_subs
    curr_sub=sub_ind(c_sub);
    [i1,i2]=min(the_new_log_L(:,c_sub));
    themodeldata(c_sub).minCost=i1;
    themodeldata(c_sub).best_param_set=param_set_list(i2,:); % i.e. actual best parameters
end
        
disp('best_group_parameters:minimum cost function')
[i1,i2]=min(mean(the_new_log_L,2))

disp('groupfit:the best parameters based on cost function was')
param.param_set_list(i2,:)


% SAVE
if save_var==1
    save(the_filename,'themodeldata','param','the_new_log_L','SMCopt');
end

try
    delete(gcp)
catch
    warning('Couldn''t close parallel pools.');
end