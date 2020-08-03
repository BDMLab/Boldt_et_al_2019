% written by Annika Boldt, 2015
% Exp2 from Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
% modulates exploration and exploitation in value-based learning. Neuroscience
% of Consciousness, 2019(1), 1–12. https://doi.org/10.1093/nc/niz004

clear all
clc

resetrn = sum(100*clock);
rand('state',resetrn);

setting.probs = [1 2 3 3 3;
                 3 3 3 2 1]';
setting.order = nchoosek(1:size(setting.probs,1),2);

setting.datapath = fullfile(pwd,'data');

setting.ntrialsprac1 = 5; % choice + CJ only
setting.ntrialsprac2 = 12; % mixed practice
setting.nblocksprac1 = 1;
setting.nblocksprac2 = 1;
setting.ntrialsexp = 600;
setting.mintrials = 20;
setting.maxtrials = 60;
setting.minblocks = 10;
setting.maxblocks = 20;
setting.ntrialsall = [repmat(setting.ntrialsprac1,1,setting.nblocksprac1)...
                      repmat(setting.ntrialsprac2,1,setting.nblocksprac2)];

setting.ntotalblocks = 0;

    while setting.ntotalblocks<setting.minblocks || setting.ntotalblocks>setting.maxblocks
        setting.totalN = [];
        while sum(setting.totalN) < (setting.ntrialsexp - setting.mintrials)
            x = round(betarnd(3,4) * 100);
            while x < setting.mintrials || x > setting.maxtrials
                x = round(betarnd(3,4) * 100);
            end
            if (sum(setting.totalN) + x) > setting.ntrialsexp
                x = setting.ntrialsexp - sum(setting.totalN);
            end
            setting.totalN = [setting.totalN x];
        end
        if sum(setting.totalN) < setting.ntrialsexp
            setting.totalN(end) = setting.ntrialsexp - sum(setting.totalN(1:(end-1)));
        end
        if (setting.totalN(end)) > setting.maxtrials
            setting.totalN = [setting.totalN setting.mintrials];
            setting.totalN(end-1) = setting.ntrialsexp - sum(setting.totalN) + ...
                                    setting.totalN(end-1);
        end
        setting.ntotalblocks = length(setting.totalN);
    end

setting.ntrialsall = [setting.ntrialsall setting.totalN];
setting.ntotalblocksexp = length(setting.ntrialsall);
setting.maxtrials = max(setting.ntrialsall);



% We have to find the order of the probabilities for the entire experiment:
setting.banditsprob = repmat(setting.order,fix(setting.ntotalblocks/size(setting.order,1)),1);
ttemp = randperm(size(setting.order,1));
ttemp = ttemp(1:mod(setting.ntotalblocks,size(setting.order,1)));
rest = setting.order;
setting.banditsprob = [setting.banditsprob; rest(ttemp,:)];
setting.ntotalblocks = size(setting.banditsprob,1);
setting.banditsprob = setting.banditsprob(shuffle(1:setting.ntotalblocks),:);
setting.banditsprob = [repmat(setting.banditsprob(:,1),1,2) repmat(setting.banditsprob(:,2),1,2)];
for i=1:setting.ntotalblocks
    setting.banditsprob(i,1:2) = setting.probs(setting.banditsprob(i,1),:);
    setting.banditsprob(i,3:4) = setting.probs(setting.banditsprob(i,3),:);
end

% shuffle sides on which bandits are shown (not balanced)
for i=1:setting.ntotalblocks
    if rand(1)>0.5
        setting.banditsprob(i,1:4) = setting.banditsprob(i,[3,4,1,2]);
    end
end
% Now sort everything, here is the resulting array:
% probability A, probability B, N trials shown, block, blocktype, order within block

% Now add practice trials to this
setting.banditsprob = [rand(2,4)*2+1; setting.banditsprob];



fordata.block = [];
fordata.part = [];
fordata.trial = 1:sum(setting.ntrialsall);
fordata.withinblocktrial = [];
fordata.maxblocktrial = [];
fordata.ntotalblocks = [];
fordata.maxtrials = [];
fordata.oo = [];
fordata.outcomedec1 = [];
fordata.outcomedec2 = [];
fordata.bb_scoredec1 = [];
fordata.bb_scoredec2 = [];
fordata.metatype = [];
fordata.banditsprobs1a = [];
fordata.banditsprobs1b = [];
fordata.banditsprobs2a = [];
fordata.banditsprobs2b = [];
fordata.banditsmean1 = [];
fordata.banditsmean2 = [];
fordata.banditsvar1 = [];
fordata.banditsvar2 = [];
fordata.banditsdiff = [];
fordata.Nobsblock = [];

for block=1:setting.ntotalblocksexp

    % calculate variable for experiment part: 1=adjustment block
    % (always at end), 2=confidence practice, 3=confidence practice
    % and cues, 4=experiment
    if block <= setting.nblocksprac1
        part = 1;
    elseif block == setting.nblocksprac1 + setting.nblocksprac2
        part = 2;
    else part = 3;
    end

    banditsnow = setting.banditsprob(block,:);
    [t1, t2] = betastat(banditsnow(1),banditsnow(2));
    [t3, t4] = betastat(banditsnow(3),banditsnow(4));
    banditsmean = [t1 t3];
    banditsvar = [t2 t4];
    diffs = banditsmean(1)-banditsmean(2);
    banditsnowlong = repmat(banditsnow,setting.ntrialsall(block),1); % banditsp
    diffslong = repmat(diffs,setting.ntrialsall(block),1); % banditspdiff
    banditsmeanlong = repmat(banditsmean,setting.ntrialsall(block),1);
    banditsvarlong = repmat(banditsvar,setting.ntrialsall(block),1);

    % metatype=0 no confidence judgement, metatype=1 first-order
    % confidence, metatype=2 second-order confidence
    switch part
        case 1
            Nobsblock = 0;
            metatype = repmat(2,1,setting.ntrialsall(block));
        case 2
            Nobsblock = 3;
            metatype = [2 2 2 2 1 2 2 2 1 2 2 1];
        case 3
            shuff = randperm(setting.ntrialsall(block));
            Nobsblock = ceil(setting.ntrialsall(block)*0.25);
            ttemp = [ones(1,Nobsblock) repmat(2,1,setting.ntrialsall(block)-Nobsblock)];
            metatype = ttemp(shuff);
            ttemp = metatype;
            ttemp(ttemp==2) = NaN;
            while metatype(1)==1 || length(find(ttemp-[NaN ttemp(1:(end-1))]==0))
                shuff = randperm(setting.ntrialsall(block));
                metatype = metatype(shuff);
                ttemp = metatype;
                ttemp(ttemp==2) = NaN;
            end
    end

    cumdec = NaN(setting.ntrialsall(block),2);
    cumdec(1,:) = [0 0];

    oo = NaN(1,setting.ntrialsall(block));
    for i=1:setting.ntrialsall(block)
        [~,oo(i)] = max(banditsmean);
    end



    fordata.block = [fordata.block repmat(block,1,setting.ntrialsall(block))];
    fordata.part = [fordata.part repmat(part,1,setting.ntrialsall(block))];
    fordata.withinblocktrial = [fordata.withinblocktrial 1:setting.ntrialsall(block)];
    fordata.maxblocktrial = [fordata.maxblocktrial repmat(setting.ntrialsall(block),1,setting.ntrialsall(block))];
    fordata.ntotalblocks = [fordata.ntotalblocks repmat(setting.ntotalblocksexp,1,setting.ntrialsall(block))];
    fordata.maxtrials = [fordata.maxtrials repmat(max(setting.ntrialsall),1,setting.ntrialsall(block))];
    fordata.oo = [fordata.oo oo];

    for itrial=1:setting.ntrialsall(block)
        fordata.outcomedec1 = [fordata.outcomedec1 round(betarnd(banditsnow(1),banditsnow(2))*100)];
        fordata.outcomedec2 = [fordata.outcomedec2 round(betarnd(banditsnow(3),banditsnow(4))*100)];
        fordata.bb_scoredec1 = [fordata.bb_scoredec1 binornd(1,fordata.outcomedec1(end)/100)];
        fordata.bb_scoredec2 = [fordata.bb_scoredec2 binornd(1,fordata.outcomedec2(end)/100)];
    end

    fordata.metatype = [fordata.metatype metatype];
    fordata.banditsprobs1a = [fordata.banditsprobs1a repmat(setting.banditsprob(block,1),1,setting.ntrialsall(block))];
    fordata.banditsprobs1b = [fordata.banditsprobs1b repmat(setting.banditsprob(block,2),1,setting.ntrialsall(block))];
    fordata.banditsprobs2a = [fordata.banditsprobs2a repmat(setting.banditsprob(block,3),1,setting.ntrialsall(block))];
    fordata.banditsprobs2b = [fordata.banditsprobs2b repmat(setting.banditsprob(block,4),1,setting.ntrialsall(block))];
    fordata.banditsmean1 = [fordata.banditsmean1 banditsmeanlong(:,1)'];
    fordata.banditsmean2 = [fordata.banditsmean2 banditsmeanlong(:,2)'];
    fordata.banditsvar1 = [fordata.banditsvar1 banditsvarlong(:,1)'];
    fordata.banditsvar2 = [fordata.banditsvar2 banditsvarlong(:,2)'];
    fordata.banditsdiff = [fordata.banditsdiff diffslong'];
    fordata.Nobsblock = [fordata.Nobsblock repmat(Nobsblock,1,setting.ntrialsall(block))];

end

save(fullfile(setting.datapath,['Exp2_fordata_',datestr(now,'yyyymmdd_HHMMSS'),'.mat']),'fordata','setting');
