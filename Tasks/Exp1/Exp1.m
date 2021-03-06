function Exp1(subNo,restart)
% this experiment was written by Annika Boldt, 2015
% Exp1 from Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
% modulates exploration and exploitation in value-based learning. Neuroscience
% of Consciousness, 2019(1), 1–12. https://doi.org/10.1093/nc/niz004

% To run the experiment, please call function with 2 parameters from the
% Matlab Command Window:
% 1. subNo: participant number
% 2. restart: whether or not you want to restart the experiment after a
% crash (0: run from beginning, 1: run from where it crashed)
% For example, Exp1(99,0) runs the experiment for participant 99 from the
% start.
% To terminate the experiment, please hold down the 'e', 'x', 'i', and 't'
% keys.
% Please note that if you start the experiment from the beginning with
% participant numbers other than '99', the experiment will prevent overwriting.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% any preliminary stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% stop keyboard output in MATLAB window
ListenChar(2);
HideCursor;

% clear Matlab window
clc;

% check for Opengl compatibility, abort otherwise
AssertOpenGL;

% Check if all needed parameters given
if nargin < 2
    error('Must provide required input parameter "subNo" and "restart"!');
end

allsettings;

% initialise variables
data = [];
testlist = [];

% in case experiment was restarted after crash
if (restart)
    load(['data/Exp1_',int2str(subNo),'.mat']);
    restart = 1;
    startblock = block;
    if (trial == setting.ntrialsall(blockorder(block)))
        trial = 1;
        block = data.block(length(data.block))+1;
    end
else startblock = 1;
end

try

% Get the maximum screen number i.e. get an external screen if avaliable
Screen('Preference','SkipSyncTests', 1);
whichScreen = max(Screen('Screens'));
whichScreen = 0;
[w,rect] = Screen('OpenWindow',whichScreen,0,[],32,2);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;

Screen(w,'FillRect',0);
Screen('TextFont',w,'Helvetica');


% define rectangles
bandits = NaN(2,2,2);
bandits(:,:,1) = [rect(3)/2-250, rect(3)/2-50; rect(4)/2-150, rect(4)/2-40];
bandits(:,:,2) = [rect(3)/2+50, rect(3)/2+250; rect(4)/2-150, rect(4)/2-40];
banditscentre = NaN(2,2);
banditscentre(:,1) = [bandits(1,1,1)+(bandits(1,2,1)-bandits(1,1,1))/2; ...
                      bandits(2,1,1)+(bandits(2,2,1)-bandits(2,1,1))/2];
banditscentre(:,2) = [bandits(1,1,2)+(bandits(1,2,2)-bandits(1,1,2))/2; ...
                      bandits(2,1,2)+(bandits(2,2,2)-bandits(2,1,2))/2];
banditsborder = NaN(2,8,2);
banditsborder(:,:,1) = [bandits(1,1,1) bandits(1,2,1) bandits(1,1,1) bandits(1,2,1) ...
                        bandits(1,1,1) bandits(1,1,1) bandits(1,2,1) bandits(1,2,1); ...
                        bandits(2,1,1) bandits(2,1,1) bandits(2,2,1) bandits(2,2,1) ...
                        bandits(2,1,1) bandits(2,2,1) bandits(2,1,1) bandits(2,2,1)];
banditsborder(:,:,2) = [bandits(1,1,2) bandits(1,2,2) bandits(1,1,2) bandits(1,2,2) ...
                        bandits(1,1,2) bandits(1,1,2) bandits(1,2,2) bandits(1,2,2); ...
                        bandits(2,1,2) bandits(2,1,2) bandits(2,2,2) bandits(2,2,2) ...
                        bandits(2,1,2) bandits(2,2,2) bandits(2,1,2) bandits(2,2,2)];
probmeta = NaN(2,2,2);
probmeta(:,:,1) = [rect(3)/2-250, rect(3)/2-50; rect(4)/2+30, rect(4)/2+230];
probmeta(:,:,2) = [rect(3)/2+50, rect(3)/2+250; rect(4)/2+30, rect(4)/2+230];
probmetacentre = NaN(2,2);
probmetacentre(:,1) = [probmeta(1,1,1)+(probmeta(1,2,1)-probmeta(1,1,1))/2; ...
                       probmeta(2,1,1)+(probmeta(2,2,1)-probmeta(2,1,1))/2];
probmetacentre(:,2) = [probmeta(1,1,2)+(probmeta(1,2,2)-probmeta(1,1,2))/2; ...
                       probmeta(2,1,2)+(probmeta(2,2,2)-probmeta(2,1,2))/2];
probmetaborder = NaN(2,12,2);
probmetaborder(:,:,1) = [probmeta(1,1,1) probmeta(1,2,1) probmeta(1,1,1) probmeta(1,2,1) ...
                         probmeta(1,1,1) probmeta(1,1,1) probmeta(1,2,1) probmeta(1,2,1) ...
                         probmetacentre(1) probmetacentre(1) ...
                         probmeta(1,1,1) probmeta(1,2,1); ...
                         probmeta(2,1,1) probmeta(2,1,1) probmeta(2,2,1) probmeta(2,2,1) ...
                         probmeta(2,1,1) probmeta(2,2,1) probmeta(2,1,1) probmeta(2,2,1) ...
                         probmeta(2,1,1) probmeta(2,2,1) ...
                         probmeta(2,1,1)+(probmeta(2,2,1)-probmeta(2,1,1))/2 probmeta(2,1,1)+(probmeta(2,2,1)-probmeta(2,1,1))/2];
probmetaborder(:,:,2) = [probmeta(1,1,2) probmeta(1,2,2) probmeta(1,1,2) probmeta(1,2,2) ...
                         probmeta(1,1,2) probmeta(1,1,2) probmeta(1,2,2) probmeta(1,2,2) ...
                         probmeta(1,1,2)+(probmeta(1,2,2)-probmeta(1,1,2))/2 probmeta(1,1,2)+(probmeta(1,2,2)-probmeta(1,1,2))/2 ...
                         probmeta(1,1,2) probmeta(1,2,2); ...
                         probmeta(2,1,2) probmeta(2,1,2) probmeta(2,2,2) probmeta(2,2,2) ...
                         probmeta(2,1,2) probmeta(2,2,2) probmeta(2,1,2) probmeta(2,2,2) ...
                         probmeta(2,1,2) probmeta(2,2,2) ...
                         probmetacentre(2) probmetacentre(2)];
cjrect = [rect(3)/2-100,rect(3)/2+100;rect(4)/2+100,rect(4)/2+150];
cjrectborder = [cjrect(1,1) cjrect(1,2) cjrect(1,1) cjrect(1,2) ...
                cjrect(1,1) cjrect(1,1) cjrect(1,2) cjrect(1,2); ...
                cjrect(2,1) cjrect(2,1) cjrect(2,2) cjrect(2,2) ...
                cjrect(2,1) cjrect(2,2) cjrect(2,1) cjrect(2,2)];
frame = [rect(3)/2-400, rect(3)/2+400; rect(4)/2-300, rect(4)/2+300];
frameborder = [frame(1,1) frame(1,2) frame(1,1) frame(1,2) ...
               frame(1,1) frame(1,1) frame(1,2) frame(1,2); ...
               frame(2,1) frame(2,1) frame(2,2) frame(2,2) ...
               frame(2,1) frame(2,2) frame(2,1) frame(2,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set priority for script execution to realtime priority:
priorityLevel=MaxPriority(w);
Priority(priorityLevel);

startexp = GetSecs;

block=startblock;
endblock=setting.ntotalblocksexp;


while block<=endblock

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pseudo-randomising trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~restart

    cumdec = NaN(setting.ntrialsall(blockorder(block)),2);
    cumdec(1,:) = [0 0];

    %save randomized trials
    testlist.subNo(block) = subNo;
    testlist.restart(block) = setting.restarted;
    testlist.blockstart(block,:) = datestr(clock);
    testlist.ntrialsprac1(block) = setting.ntrialsprac1;
    testlist.ntrialsprac2(block) = setting.ntrialsprac2;
    testlist.ntrialsprac3(block) = setting.ntrialsprac3;
    testlist.ntrialsexp(block) = setting.ntrialsexp;
    testlist.ntrialsall(block) = setting.ntrialsall(blockorder(block));
    testlist.nblocksprac1(block) = setting.nblocksprac1;
    testlist.nblocksprac2(block) = setting.nblocksprac2;
    testlist.nblocksprac3(block) = setting.nblocksprac3;
    testlist.nblocks(block) = setting.ntotalblocks;
    testlist.ntotalblocks(block) = setting.ntotalblocksexp;
    testlist.mintrials(block) = setting.mintrials;
    testlist.maxtrials(block) = setting.maxtrials;
    testlist.part(block) = unique(fordata.part(fordata.block==blockorder(block)));
    testlist.cjmap(block) = setting.mapcj;

    testlist.Nobsblock(block) = unique(fordata.Nobsblock(fordata.block==blockorder(block)));
    testlist.Nobs(block,:) = [unique(fordata.Nobs1(fordata.block==blockorder(block))) unique(fordata.Nobs2(fordata.block==blockorder(block)))];
    testlist.banditsprobs(block,:) = [unique(fordata.banditsprobs1a(fordata.block==blockorder(block))) unique(fordata.banditsprobs1b(fordata.block==blockorder(block)))...
                                      unique(fordata.banditsprobs2a(fordata.block==blockorder(block))) unique(fordata.banditsprobs2b(fordata.block==blockorder(block)))];
    testlist.banditsmean(block,:) = [unique(fordata.banditsmean1(fordata.block==blockorder(block))) unique(fordata.banditsmean2(fordata.block==blockorder(block)))];
    testlist.banditsdiff(block,:) = unique(fordata.banditsdiff(fordata.block==blockorder(block)));

    metatype = fordata.metatype(fordata.block==blockorder(block));
    whichobs = fordata.whichobs(fordata.block==blockorder(block));
    oo = fordata.oo(fordata.block==blockorder(block));
    maxtrials = unique(fordata.maxtrials);

    testlist.Nobslong(block,:,:) = [fordata.Nobs1(fordata.block==blockorder(block))' fordata.Nobs2(fordata.block==blockorder(block))';
                                    NaN(maxtrials-size(fordata.Nobs2(fordata.block==blockorder(block)),2),2)];
    testlist.banditsproblong(block,:,:) = [fordata.banditsprobs1a(fordata.block==blockorder(block))' fordata.banditsprobs1b(fordata.block==blockorder(block))' ...
                                           fordata.banditsprobs2a(fordata.block==blockorder(block))' fordata.banditsprobs2b(fordata.block==blockorder(block))'; ...
                                           NaN(maxtrials-size(fordata.banditsprobs1a(fordata.block==blockorder(block)),2),4)];
    testlist.banditsmeanlong(block,:,:) = [fordata.banditsmean1(fordata.block==blockorder(block))' fordata.banditsmean2(fordata.block==blockorder(block))';
                                           NaN(maxtrials-size(fordata.banditsmean1(fordata.block==blockorder(block)),2),2)];
    testlist.banditsdifflong(block,:) = [fordata.banditsdiff(fordata.block==blockorder(block)) NaN(1,maxtrials-length(fordata.banditsdiff(fordata.block==blockorder(block))))];

    testlist.oo(block,:) = [oo NaN(1,maxtrials-size(oo,2))];
    testlist.metatype(block,:) = [metatype NaN(1,maxtrials-size(metatype,2))];
    testlist.whichobs(block,:) = [whichobs NaN(1,maxtrials-size(whichobs,2))];
    testlist.cumobs(block,:,:) = [fordata.cumobs1(fordata.block==blockorder(block))' fordata.cumobs2(fordata.block==blockorder(block))';
                                  NaN(maxtrials-size(fordata.cumobs1(fordata.block==blockorder(block)),2),2)];
    if unique(fordata.part(fordata.block==blockorder(block)))==4
        testlist.casinonum(block) = setting.casinos(block-3);
    else testlist.casinonum(block) = NaN;
    end
    mintrials = tabulate(fordata.withinblocktrial(fordata.part==4));
    ttemp = find(mintrials(:,3)==max(mintrials(:,3)));
    mintrials = mintrials(ttemp(end),1);


    starttrial = 1;
else
        starttrial = trial;
        restart = 0;
        setting.restarted = 1;
end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % instruction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % wait a bit between blocks
    WaitSecs(0.600);
    instruction;

    % Clear screen to background color
    Screen('Flip',w);
    WaitSecs(0.600);



    for trial = starttrial:setting.ntrialsall(blockorder(block))

        timestamp_blankstart = NaN;
        timestamp_blankend = NaN;
        timestamp_startobs = NaN;
        timestamp_highlightobs = NaN;
        timestamp_spinobs = NaN;
        timestamp_feedbackobs = NaN;
        timestamp_cj1 = NaN;
        timestamp_cj1fbstart = NaN;
        timestamp_cj1fbend = NaN;
        timestamp_blankobs = NaN;
        timestamp_endobs = NaN;
        timestamp_startdec = NaN;
        timestamp_resp = NaN;
        timestamp_startcj2 = NaN;
        timestamp_cj2 = NaN;
        timestamp_cj2fbstart = NaN;
        timestamp_cj2fbend = NaN;
        timestamp_spindec = NaN;
        timestamp_feedbackdec = NaN;
        timestamp_enddec = NaN;

        outcome = NaN;
        cj1 = NaN;
        cj1rt = NaN;
        prob = NaN;
        resp = NaN;
        rt = NaN;
        cj2 = NaN;
        cj2rt = NaN;

        buttons = NaN;
        xcoord = NaN;
        ycoord = NaN;
        buttonscj = NaN;
        xcoordcj = NaN;
        ycoordcj = NaN;
        buttonsobs = NaN;
        xcoordobs = NaN;
        ycoordobs = NaN;
        jitter = [];
        theX = [];
        theY = [];


        Screen('DrawLines',w,banditsborder(:,:,1),5,255);
        Screen('DrawLines',w,banditsborder(:,:,2),5,255);
        if block>3
            Screen('TextSize',w, 28);
            DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
        end
        [~, timestamp_blankstart] = Screen('Flip', w);
        timestamp_blankend = WaitSecs('UntilTime', timestamp_blankstart+setting.iti.blank);

        if metatype(trial) == 1

            Screen('DrawLines',w,frameborder,7,[255 0 0]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end

            outcome = fordata.outcome(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);

            Screen('DrawLines',w,banditsborder(:,:,1),5,255);
            Screen('DrawLines',w,banditsborder(:,:,2),5,255);
            [~, timestamp_startobs] = Screen('Flip', w);
            WaitSecs('UntilTime', timestamp_startobs+setting.iti.OBS.start_highlight);

            % highlight bandit
            Screen('DrawLines',w,frameborder,7,[255 0 0]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            Screen('DrawLines',w,banditsborder(:,:,abs(3-whichobs(trial))),5,255);
            Screen('DrawLines',w,banditsborder(:,:,whichobs(trial)),5,[255 255 0]);
            [~, timestamp_highlightobs] = Screen('Flip', w);
            WaitSecs('UntilTime', timestamp_highlightobs+setting.iti.OBS.highlight_spin);

            % spin bandit
            timestamp_spinobs = GetSecs;
            jitter = unifrnd(-0.01,0.01);
            for i=[0:10:120 110:-10:0]
                Screen('DrawLines',w,frameborder,7,[255 0 0]);
                if block>3
                    Screen('TextSize',w, 28);
                    DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
                end
                Screen(w,'FillRect',i,bandits(:,:,whichobs(trial)));
                Screen('DrawLines',w,banditsborder(:,:,abs(3-whichobs(trial))),5,255);
                Screen('DrawLines',w,banditsborder(:,:,whichobs(trial)),5,[255 255 0]);
                Screen('Flip', w);
                WaitSecs(setting.iti.spin+jitter);
            end

            Screen('DrawLines',w,frameborder,7,[255 0 0]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            Screen('DrawLines',w,banditsborder(:,:,abs(3-whichobs(trial))),5,255);
            Screen('DrawLines',w,banditsborder(:,:,whichobs(trial)),5,[255 255 0]);
            Screen('DrawLines',w,probmetaborder(:,:,whichobs(trial)),2,255);
            Screen('TextSize',w,14);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{1}),probmeta(1,1,whichobs(trial))-70,probmetacentre(2,whichobs(trial))-70);
            Screen('DrawText',w,setting.confpoles{1},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{2}),probmeta(1,2,whichobs(trial))+70,probmetacentre(2,whichobs(trial))-70);
            Screen('DrawText',w,setting.confpoles{2},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.probpoles{2}),probmetacentre(1,whichobs(trial)),probmeta(2,1,whichobs(trial))-90);
            Screen('DrawText',w,setting.probpoles{2},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.probpoles{1}),probmetacentre(1,whichobs(trial)),probmeta(2,1,whichobs(trial))+150);
            Screen('DrawText',w,setting.probpoles{1},rct(1),rct(2)+70,255);
            Screen('TextSize',w,36);
            rct = CenterRectOnPoint(Screen('TextBounds',w,int2str(outcome)),banditscentre(1,whichobs(trial)),banditscentre(2,whichobs(trial)));
            Screen('DrawText',w,int2str(outcome),rct(1),rct(2),255);
            [~, timestamp_feedbackobs] = Screen('Flip',w);

            ShowCursor('Arrow');
            % reset the cursor to somewhere within the scale
            theX = round(unifrnd(probmeta(1,1,whichobs(trial)),probmeta(1,2,whichobs(trial))));
            theY = round(unifrnd(probmeta(2,1,whichobs(trial)),probmeta(2,2,whichobs(trial))));
            while 1
                SetMouse(theX,theY);
                [checkX,checkY] = GetMouse;
                if (checkX==theX) && (checkY==theY)
                    break;
                end
            end

            press=0;
            while press==0
                [xcoordobs, ycoordobs, buttonsobs] = GetMouse;

                % check escape key
                [~,~,KeyCode] = KbCheck;
                if (KeyCode(KbName('e')) && KeyCode(KbName('x')) && KeyCode(KbName('i')) && KeyCode(KbName('t')))
                    save(fullfile(setting.datapath,['Exp1_' num2str(subNo) '.mat']),'-v6');
                    Screen('CloseAll');
                    ListenChar(0);
                    Priority(0);
                    break
                end

                if buttonsobs(1)==1
                timestamp_cj1 = GetSecs;
                    if (xcoordobs>probmeta(1,1,whichobs(trial))) && (xcoordobs<probmeta(1,2,whichobs(trial))) && (ycoordobs>probmeta(2,1,whichobs(trial))) && (ycoordobs<probmeta(2,2,whichobs(trial)))
                        press = 1;
                        cj1 = (xcoordobs-probmeta(1,1,whichobs(trial)))/(probmeta(1,2,whichobs(trial))-probmeta(1,1,whichobs(trial)));
                        prob = 100 - round((ycoordobs-probmeta(2,1,whichobs(trial)))/(probmeta(2,2,whichobs(trial))-probmeta(2,1,whichobs(trial))) * 100);
                    end
                end

                WaitSecs(0.001);
            end
            HideCursor;

            Screen('DrawLines',w,frameborder,7,[255 0 0]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            Screen('DrawLines',w,banditsborder(:,:,abs(3-whichobs(trial))),5,255);
            Screen('DrawLines',w,banditsborder(:,:,whichobs(trial)),5,[255 255 0]);
            Screen('DrawLines',w,probmetaborder(:,:,whichobs(trial)),2,255);
            Screen('TextSize',w,14);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{1}),probmeta(1,1,whichobs(trial))-70,probmetacentre(2,whichobs(trial))-70);
            Screen('DrawText',w,setting.confpoles{1},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{2}),probmeta(1,2,whichobs(trial))+70,probmetacentre(2,whichobs(trial))-70);
            Screen('DrawText',w,setting.confpoles{2},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.probpoles{2}),probmetacentre(1,whichobs(trial)),probmeta(2,1,whichobs(trial))-90);
            Screen('DrawText',w,setting.probpoles{2},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.probpoles{1}),probmetacentre(1,whichobs(trial)),probmeta(2,1,whichobs(trial))+150);
            Screen('DrawText',w,setting.probpoles{1},rct(1),rct(2)+70,255);
            Screen('TextSize',w,36);
            rct = CenterRectOnPoint(Screen('TextBounds',w,int2str(outcome)),banditscentre(1,whichobs(trial)),banditscentre(2,whichobs(trial)));
            Screen('DrawText',w,int2str(outcome),rct(1),rct(2),255);
            thedot = CenterRectOnPoint([0 0 15 15],xcoordobs,ycoordobs);
            Screen('FillOval',w,[255 255 0],thedot);
            [~, timestamp_cj1fbstart] = Screen('Flip',w);
            timestamp_cj1fbend = WaitSecs('UntilTime', timestamp_cj1fbstart+setting.iti.OBS.cj1fbstart_cj1fbend);


            if setting.mapcj==21
                cj1 = 1-cj1;
            end
            cj1rt = round(1000*(timestamp_cj1-timestamp_feedbackobs));

            Screen('DrawLines',w,frameborder,7,[255 0 0]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            Screen('DrawLines',w,banditsborder(:,:,abs(3-whichobs(trial))),5,255);
            Screen('DrawLines',w,banditsborder(:,:,whichobs(trial)),5,[255 255 0]);
            Screen('TextSize',w,36);
            rct = CenterRectOnPoint(Screen('TextBounds',w,int2str(outcome)),banditscentre(1,whichobs(trial)),banditscentre(2,whichobs(trial)));
            Screen('DrawText',w,int2str(outcome),rct(1),rct(2),255);

            Screen('DrawLines',w,probmetaborder(:,:,whichobs(trial)),2,255);
            Screen('TextSize',w,14);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{1}),probmeta(1,1,whichobs(trial))-70,probmetacentre(2,whichobs(trial))-70);
            Screen('DrawText',w,setting.confpoles{1},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{2}),probmeta(1,2,whichobs(trial))+70,probmetacentre(2,whichobs(trial))-70);
            Screen('DrawText',w,setting.confpoles{2},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.probpoles{2}),probmetacentre(1,whichobs(trial)),probmeta(2,1,whichobs(trial))-90);
            Screen('DrawText',w,setting.probpoles{2},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.probpoles{1}),probmetacentre(1,whichobs(trial)),probmeta(2,1,whichobs(trial))+150);
            Screen('DrawText',w,setting.probpoles{1},rct(1),rct(2)+70,255);
            Screen('TextSize',w,36);
            rct = CenterRectOnPoint(Screen('TextBounds',w,int2str(outcome)),banditscentre(1,whichobs(trial)),banditscentre(2,whichobs(trial)));
            Screen('DrawText',w,int2str(outcome),rct(1),rct(2),255);
            thedot = CenterRectOnPoint([0 0 15 15],xcoordobs,ycoordobs);
            Screen('FillOval',w,[255 255 0],thedot);

            [~, timestamp_blankobs] = Screen('Flip',w);
            timestamp_endobs = WaitSecs('UntilTime', timestamp_blankobs+setting.iti.OBS.blank_end);

            if trial~=1
                cumdec(trial,:) = cumdec(trial-1,:);
            end



        elseif metatype(trial) == 2

            Screen('DrawLines',w,frameborder,7,[0 0 255]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            ShowCursor('Arrow');
            % reset the cursor of the mouse to inbetween the two bandits
            theX = rect(3)/2;
            theY = banditscentre(2,1);
            while 1
                SetMouse(theX,theY);
                [checkX,checkY] = GetMouse;
                if (checkX==theX) && (checkY==theY)
                    break;
                end
            end

            Screen('DrawLines',w,banditsborder(:,:,1),5,255);
            Screen('DrawLines',w,banditsborder(:,:,2),5,255);
            [~, timestamp_startdec] = Screen('Flip', w);

            press=0;
            while press==0
                [xcoord, ycoord, buttons] = GetMouse;

                % check escape key
                [~,~,KeyCode] = KbCheck;
                if (KeyCode(KbName('e')) && KeyCode(KbName('x')) && KeyCode(KbName('i')) && KeyCode(KbName('t')))
                    save(fullfile(setting.datapath,['Exp1_' num2str(subNo) '.mat']),'-v6');
                    Screen('CloseAll');
                    ListenChar(0);
                    Priority(0);
                    break
                end

                if buttons(1)==1
                timestamp_resp = GetSecs;
                    if     (xcoord>bandits(1,1,1)) && (xcoord<bandits(1,2,1)) && (ycoord>bandits(2,1,1)) && (ycoord<bandits(2,2,1)) || ...
                           (xcoord>bandits(1,1,2)) && (xcoord<bandits(1,2,2)) && (ycoord>bandits(2,1,2)) && (ycoord<bandits(2,2,2))
                        press = 1;
                    end
                end

                WaitSecs(0.001);
            end

            if (xcoord>bandits(1,1,1)) && (xcoord<bandits(1,2,1)) && (ycoord>bandits(2,1,1)) && (ycoord<bandits(2,2,1))
                resp = 1;
                outcome = fordata.outcomedec1(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);
            end
            if (xcoord>bandits(1,1,2)) && (xcoord<bandits(1,2,2)) && (ycoord>bandits(2,1,2)) && (ycoord<bandits(2,2,2))
                resp = 2;
                outcome = fordata.outcomedec2(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);
            end
            rt = round(1000*(timestamp_resp-timestamp_startdec));

            Screen('DrawLines',w,frameborder,7,[0 0 255]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            Screen('DrawLines',w,banditsborder(:,:,abs(3-resp)),5,255);
            Screen('DrawLines',w,banditsborder(:,:,resp),5,[255 255 0]);
            Screen('DrawLines',w,cjrectborder,5,255);
            Screen('TextSize',w,14);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.cjq),rect(3)/2,rect(4)/2);
            Screen('DrawText',w,setting.cjq,rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{1}),cjrect(1,1)-70,cjrect(2,1)-50);
            Screen('DrawText',w,setting.confpoles{1},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{2}),cjrect(1,2)+70,cjrect(2,1)-50);
            Screen('DrawText',w,setting.confpoles{2},rct(1),rct(2)+70,255);
            [~, timestamp_startcj2] = Screen('Flip', w);

            press=0;
            while press==0
                [xcoordcj, ycoordcj, buttonscj] = GetMouse;

                % check escape key
                [~,~,KeyCode] = KbCheck;
                if (KeyCode(KbName('e')) && KeyCode(KbName('x')) && KeyCode(KbName('i')) && KeyCode(KbName('t')))
                    save(fullfile(setting.datapath,['Exp1_' num2str(subNo) '.mat']),'-v6');
                    Screen('CloseAll');
                    ListenChar(0);
                    Priority(0);
                    break
                end

                if (buttonscj(1)==1)
                    timestamp_cj2 = GetSecs;
                    if     (xcoordcj>cjrect(1,1)) && (xcoordcj<cjrect(1,2)) && (ycoordcj>cjrect(2,1)) && (ycoordcj<cjrect(2,2))
                        press = 1;
                        cj2 = (xcoordcj-cjrect(1,1))/(cjrect(1,2)-cjrect(1,1));
                    end
                end

                WaitSecs(0.001);
            end
            HideCursor;


            Screen('DrawLines',w,frameborder,7,[0 0 255]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            Screen('DrawLines',w,banditsborder(:,:,abs(3-resp)),5,255);
            Screen('DrawLines',w,banditsborder(:,:,resp),5,[255 255 0]);
            Screen('DrawLines',w,cjrectborder,5,255);
            Screen('TextSize',w,14);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.cjq),rect(3)/2,rect(4)/2);
            Screen('DrawText',w,setting.cjq,rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{1}),cjrect(1,1)-70,cjrect(2,1)-50);
            Screen('DrawText',w,setting.confpoles{1},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{2}),cjrect(1,2)+70,cjrect(2,1)-50);
            Screen('DrawText',w,setting.confpoles{2},rct(1),rct(2)+70,255);
            Screen('DrawLines',w,[xcoordcj, xcoordcj; cjrect(2,:)],5,[255 255 0]);
            [~, timestamp_cj2fbstart] = Screen('Flip',w);
            timestamp_cj2fbend = WaitSecs('UntilTime', timestamp_cj2fbstart+setting.iti.DEC.cj2fbstart_cj2fbend);

            if setting.mapcj==21
                cj2 = 1-cj2;
            end
            cj2rt = round(1000*(timestamp_cj2-timestamp_startcj2));

            % spin bandit
            timestamp_spindec = GetSecs;
            jitter = unifrnd(-0.01,0.01);
            for i=[0:10:120 110:-10:0]
                Screen('DrawLines',w,frameborder,7,[0 0 255]);
                if block>3
                    Screen('TextSize',w, 28);
                    DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
                end
                Screen(w,'FillRect',i,bandits(:,:,resp));
                Screen('DrawLines',w,banditsborder(:,:,abs(3-resp)),5,255);
                Screen('DrawLines',w,banditsborder(:,:,resp),5,[255 255 0]);
                Screen('DrawLines',w,cjrectborder,5,255);
                Screen('TextSize',w,14);
                rct = CenterRectOnPoint(Screen('TextBounds',w,setting.cjq),rect(3)/2,rect(4)/2);
                Screen('DrawText',w,setting.cjq,rct(1),rct(2)+70,255);
                rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{1}),cjrect(1,1)-70,cjrect(2,1)-50);
                Screen('DrawText',w,setting.confpoles{1},rct(1),rct(2)+70,255);
                rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{2}),cjrect(1,2)+70,cjrect(2,1)-50);
                Screen('DrawText',w,setting.confpoles{2},rct(1),rct(2)+70,255);
                Screen('DrawLines',w,[xcoordcj, xcoordcj; cjrect(2,:)],5,[255 255 0]);
                Screen('Flip', w);
                WaitSecs(setting.iti.spin+jitter);
            end

            Screen('DrawLines',w,frameborder,7,[0 0 255]);
            if block>3
                Screen('TextSize',w, 28);
                DrawFormattedText(w, setting.casinonames{setting.casinos(block-3)}, 'center', 50, 255, [], [], [], 1.5);
            end
            Screen('DrawLines',w,banditsborder(:,:,abs(3-resp)),5,255);
            Screen('DrawLines',w,banditsborder(:,:,resp),5,[255 255 0]);
            Screen('DrawLines',w,cjrectborder,5,255);
            Screen('TextSize',w,14);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.cjq),rect(3)/2,rect(4)/2);
            Screen('DrawText',w,setting.cjq,rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{1}),cjrect(1,1)-70,cjrect(2,1)-50);
            Screen('DrawText',w,setting.confpoles{1},rct(1),rct(2)+70,255);
            rct = CenterRectOnPoint(Screen('TextBounds',w,setting.confpoles{2}),cjrect(1,2)+70,cjrect(2,1)-50);
            Screen('DrawText',w,setting.confpoles{2},rct(1),rct(2)+70,255);
            Screen('DrawLines',w,[xcoordcj, xcoordcj; cjrect(2,:)],5,[255 255 0]);
            Screen('TextSize',w,36);
            rct = CenterRectOnPoint(Screen('TextBounds',w,int2str(outcome)),banditscentre(1,resp),banditscentre(2,resp));
            Screen('DrawText',w,int2str(outcome),rct(1),rct(2),255);
            [~, timestamp_feedbackdec] = Screen('Flip',w);
            timestamp_enddec = WaitSecs('UntilTime', timestamp_feedbackdec+setting.iti.DEC.feedback_end);

            if trial==1
                cumdec(trial,resp) = 1;
            else
                cumdec(trial,resp) = cumdec(trial-1,resp) + 1;
                cumdec(trial,3-resp) = cumdec(trial-1,3-resp);
            end

        end

        % log data
        t = sum(setting.ntrialsall(blockorder(1:(block-1))))+trial;
        data.sub(t) = subNo;
        data.restart(t) = setting.restarted;
        data.mapcj(t) = setting.mapcj;
        data.block(t) = block;
        data.blocktype(t) = blockorder(block);
        data.part(t) = unique(fordata.part(fordata.block==blockorder(block)));
        data.trial(t) = t;
        data.withinblocktrial(t) = trial;
        data.maxblocktrial(t) = setting.ntrialsall(blockorder(block));
        data.ntotalblocks(t) = length(setting.ntrialsall);
        data.mintrials(t) = mintrials;
        data.maxtrials(t) = maxtrials;
        data.resp(t) = resp;
        data.oo(t) = oo(trial);
        if metatype(trial)==2
            data.cor(t) = oo(trial)==resp;
            data.err(t) = 1-data.cor(t);
        else
            data.cor(t) = NaN;
            data.err(t) = NaN;
        end
        data.rt(t) = rt;
        data.outcome(t) = outcome;
        data.cumoutcome(t) = sum(data.outcome);
        data.cumoutcomeblock(t) = sum(data.outcome(data.block==block));
        if metatype(trial)==2
            if resp==1
                data.bb_score(t) = fordata.bb_scoredec1(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);
            elseif resp==2
                data.bb_score(t) = fordata.bb_scoredec2(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);
            end
        else data.bb_score(t) = fordata.bb_score(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);
        end
        if metatype(trial)==2
            data.win(t) = outcome;
        else
            data.win(t) = NaN;
        end
        data.cumwin(t) = sum(data.win(~isnan(data.win) & data.block>3));
        data.cumwinblock(t) = sum(data.win(~isnan(data.win) & data.block==block));
        if metatype(trial)==1
            data.obswin(t) = outcome;
            data.theX(t) = theX;
            data.theY(t) = theY;
        else
            data.obswin(t) = NaN;
            data.theX(t) = NaN;
            data.theY(t) = NaN;
        end
        data.cumobswin(t) = sum(data.obswin(~isnan(data.obswin) & data.block>3));
        data.cumobswinblock(t) = sum(data.obswin(~isnan(data.obswin) & data.block==block));

        data.jitter(t) = jitter;
        data.xcoord(t) = xcoord;
        data.ycoord(t) = ycoord;
        data.cj1(t) = cj1;
        if whichobs(trial)==1
            data.cj11(t) = cj1;
            data.cj12(t) = NaN;
        elseif whichobs(trial)==2
            data.cj11(t) = NaN;
            data.cj12(t) = cj1;
        else
            data.cj11(t) = NaN;
            data.cj12(t) = NaN;
        end
        data.cj1rt(t) = cj1rt;
        if whichobs(trial)==1
            data.cj1rt1(t) = cj1rt;
            data.cj1rt2(t) = NaN;
            data.prevcj11(t) = cj1;
            if trial==1
                data.prevcj12(t) = NaN;
            else
                data.prevcj12(t) = data.prevcj12(t-1);
            end
            data.prevcj1_diff(t) = data.prevcj11(t) - data.prevcj12(t);
            data.prevcj1_mean(t) = mean([data.prevcj11(t) data.prevcj12(t)]);
        elseif whichobs(trial)==2
            data.cj1rt1(t) = NaN;
            data.cj1rt2(t) = cj1rt;
            if trial==1
                data.prevcj11(t) = NaN;
            else
                data.prevcj11(t) = data.prevcj11(t-1);
            end
            data.prevcj12(t) = cj1;
            data.prevcj1_diff(t) = data.prevcj11(t) - data.prevcj12(t);
            data.prevcj1_mean(t) = mean([data.prevcj11(t) data.prevcj12(t)]);
        else
            data.cj1rt1(t) = NaN;
            data.cj1rt2(t) = NaN;
            if trial==1
                data.prevcj11(t) = NaN;
                data.prevcj12(t) = NaN;
            else
                data.prevcj11(t) = data.prevcj11(t-1);
                data.prevcj12(t) = data.prevcj12(t-1);
            end
            data.prevcj1_diff(t) = data.prevcj11(t) - data.prevcj12(t);
            data.prevcj1_mean(t) = mean([data.prevcj11(t) data.prevcj12(t)]);
        end
        data.prob(t) = prob;
        if whichobs(trial)==1
            data.prob1(t) = prob;
            data.prob2(t) = NaN;
            data.prevprob1(t) = prob;
            if trial==1
                data.prevprob2(t) = NaN;
            else
                data.prevprob2(t) = data.prevprob2(t-1);
            end
            data.prevprob_diff(t) = data.prevprob1(t) - data.prevprob2(t);
            data.prevprob_mean(t) = mean([data.prevprob1(t) data.prevprob2(t)]);
        elseif whichobs(trial)==2
            data.prob1(t) = NaN;
            data.prob2(t) = prob;
            if trial==1
                data.prevprob1(t) = NaN;
            else
                data.prevprob1(t) = data.prevprob1(t-1);
            end
            data.prevprob2(t) = prob;
            data.prevprob_diff(t) = data.prevprob1(t) - data.prevprob2(t);
            data.prevprob_mean(t) = mean([data.prevprob1(t) data.prevprob2(t)]);
        else
            data.prob1(t) = NaN;
            data.prob2(t) = NaN;
            if trial==1
                data.prevprob1(t) = NaN;
                data.prevprob2(t) = NaN;
            else
                data.prevprob1(t) = data.prevprob1(t-1);
                data.prevprob2(t) = data.prevprob2(t-1);
            end
            data.prevprob_diff(t) = data.prevprob1(t) - data.prevprob2(t);
            data.prevprob_mean(t) = mean([data.prevprob1(t) data.prevprob2(t)]);
        end

        data.xcoordobs(t) = xcoordobs;
        data.ycoordobs(t) = ycoordobs;
        data.cj2(t) = cj2;
        data.cj2rt(t) = cj2rt;
        data.xcoordcj(t) = xcoordcj;
        data.ycoordcj(t) = ycoordcj;
        data.timestamp_blankstart(t) = timestamp_blankstart;
        data.timestamp_blankend(t) = timestamp_blankend;
        data.timestamp_startobs(t) = timestamp_startobs;
        data.timestamp_highlightobs(t) = timestamp_highlightobs;
        data.timestamp_spinobs(t) = timestamp_spinobs;
        data.timestamp_feedbackobs(t) = timestamp_feedbackobs;
        data.timestamp_cj1(t) = timestamp_cj1;
        data.timestamp_cj1fbstart(t) = timestamp_cj1fbstart;
        data.timestamp_cj1fbend(t) = timestamp_cj1fbend;
        data.timestamp_blankobs(t) = timestamp_blankobs;
        data.timestamp_endobs(t) = timestamp_endobs;
        data.timestamp_startdec(t) = timestamp_startdec;
        data.timestamp_resp(t) = timestamp_resp;
        data.timestamp_startcj2(t) = timestamp_startcj2;
        data.timestamp_cj2(t) = timestamp_cj2;
        data.timestamp_cj2fbstart(t) = timestamp_cj2fbstart;
        data.timestamp_cj2fbend(t) = timestamp_cj2fbend;
        data.timestamp_spindec(t) = timestamp_spindec;
        data.timestamp_feedbackdec(t) = timestamp_feedbackdec;
        data.timestamp_enddec(t) = timestamp_enddec;
        data.metatype(t) = metatype(trial);
        data.whichobs(t) = whichobs(trial);
        data.banditsprobs1a(t) = unique(fordata.banditsprobs1a(fordata.block==blockorder(block)));
        data.banditsprobs1b(t) = unique(fordata.banditsprobs1b(fordata.block==blockorder(block)));
        data.banditsprobs2a(t) = unique(fordata.banditsprobs2a(fordata.block==blockorder(block)));
        data.banditsprobs2b(t) = unique(fordata.banditsprobs2b(fordata.block==blockorder(block)));
        data.banditsmean1(t) = unique(fordata.banditsmean1(fordata.block==blockorder(block)));
        data.banditsmean2(t) = unique(fordata.banditsmean2(fordata.block==blockorder(block)));
        data.banditsvar1(t) = unique(fordata.banditsvar1(fordata.block==blockorder(block)));
        data.banditsvar2(t) = unique(fordata.banditsvar2(fordata.block==blockorder(block)));
        data.banditsdiff(t) = unique(fordata.banditsdiff(fordata.block==blockorder(block)));
        data.banditsoutcomeprobs1a(t) = 1+sum(data.bb_score(data.block==block & data.whichobs==1)==1);
        data.banditsoutcomeprobs1b(t) = 1+length(data.bb_score(data.block==block & data.whichobs==1))-sum(data.bb_score(data.block==block & data.whichobs==1)==1);
        [tm,tv] = betastat(data.banditsoutcomeprobs1a(t),data.banditsoutcomeprobs1b(t));
        data.banditsoutcomemean1(t) = tm;
        data.banditsoutcomevar1(t) = tv;
        data.banditsoutcomeprobs2a(t) = 1+sum(data.bb_score(data.block==block & data.whichobs==2)==1);
        data.banditsoutcomeprobs2b(t) = 1+length(data.bb_score(data.block==block & data.whichobs==2))-sum(data.bb_score(data.block==block & data.whichobs==2)==1);
        [tm,tv] = betastat(data.banditsoutcomeprobs2a(t),data.banditsoutcomeprobs2b(t));
        data.banditsoutcomemean2(t) = tm;
        data.banditsoutcomevar2(t) = tv;
        if whichobs(trial)==1
            data.deltap1(t) = data.banditsoutcomemean1(t) - data.prob1(t);
            data.deltap2(t) = NaN;
        elseif whichobs(trial)==2
            data.deltap1(t) = NaN;
            data.deltap2(t) = data.banditsoutcomemean2(t) - data.prob2(t);
        else
            data.deltap1(t) = NaN;
            data.deltap2(t) = NaN;
        end
        data.banditsoutcomediff(t) = data.banditsoutcomemean1(t)-data.banditsoutcomemean2(t);
        data.Nobsblock(t) = unique(fordata.Nobsblock(fordata.block==blockorder(block)));
        data.Nobs1(t) = unique(fordata.Nobs1(fordata.block==blockorder(block)));
        data.Nobs2(t) = unique(fordata.Nobs2(fordata.block==blockorder(block)));
        data.cumobs1(t) = fordata.cumobs1(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);
        data.cumobs2(t) = fordata.cumobs2(fordata.block==blockorder(block) & fordata.withinblocktrial==trial);
        data.cumdec1(t) = cumdec(trial,1);
        data.cumdec2(t) = cumdec(trial,2);

        data.cumobswinblock1(t) = sum(data.obswin(~isnan(data.obswin) & data.block==block & data.whichobs==1));
        data.cumobswinblock2(t) = sum(data.obswin(~isnan(data.obswin) & data.block==block & data.whichobs==2));
        data.cumdecwinblock1(t) = sum(data.win(~isnan(data.win) & data.block==block & data.resp==1));
        data.cumdecwinblock2(t) = sum(data.win(~isnan(data.win) & data.block==block & data.resp==2));

        if data.resp(t)==1
            data.prob_c(t) = data.prevprob1(t);
            data.prob_uc(t) = data.prevprob2(t);
            data.cj1_c(t) = data.prevcj11(t);
            data.cj1_uc(t) = data.prevcj12(t);
            data.cumdec_c(t) = data.cumdec1(t);
            data.cumdec_uc(t) = data.cumdec2(t);
            data.cumdecwinblock_c(t) = data.cumdecwinblock1(t);
            data.cumdecwinblock_uc(t) = data.cumdecwinblock2(t);
            data.cumobs_c(t) = data.cumobs1(t);
            data.cumobs_uc(t) = data.cumobs2(t);
            data.cumobswinblock_c(t) = data.cumobswinblock1(t);
            data.cumobswinblock_uc(t) = data.cumobswinblock2(t);
        elseif data.resp(t)==2
            data.prob_c(t) = data.prevprob2(t);
            data.prob_uc(t) = data.prevprob1(t);
            data.cj1_c(t) = data.prevcj12(t);
            data.cj1_uc(t) = data.prevcj11(t);
            data.cumdec_c(t) = data.cumdec2(t);
            data.cumdec_uc(t) = data.cumdec1(t);
            data.cumdecwinblock_c(t) = data.cumdecwinblock2(t);
            data.cumdecwinblock_uc(t) = data.cumdecwinblock1(t);
            data.cumobs_c(t) = data.cumobs2(t);
            data.cumobs_uc(t) = data.cumobs1(t);
            data.cumobswinblock_c(t) = data.cumobswinblock2(t);
            data.cumobswinblock_uc(t) = data.cumobswinblock1(t);
        else
            data.prob_c(t) = NaN;
            data.prob_uc(t) = NaN;
            data.cj1_c(t) = NaN;
            data.cj1_uc(t) = NaN;
            data.cumdec_c(t) = NaN;
            data.cumdec_uc(t) = NaN;
            data.cumdecwinblock_c(t) = NaN;
            data.cumdecwinblock_uc(t) = NaN;
            data.cumobs_c(t) = NaN;
            data.cumobs_uc(t) = NaN;
            data.cumobswinblock_c(t) = NaN;
            data.cumobswinblock_uc(t) = NaN;
        end

        if unique(fordata.part(fordata.block==blockorder(block)))==4
            data.casinonum(t) = setting.casinos(block-3);
        else data.casinonum(t) = NaN;
        end

        data.itiought_blankstart_blankend(t) = NaN;
        data.itiought_OBS_start_highlight(t) = NaN;
        data.itiought_OBS_highlight_spin(t) = NaN;
        data.itiought_OBS_spin_feedback(t) = NaN;
        data.itiought_OBS_cj1fbstart_cj1fbend(t) = NaN;
        data.itiought_OBS_blank_end(t) = NaN;
        data.itiought_DEC_cj2fbstart_cj2fbend(t) = NaN;
        data.itiought_DEC_spin_feedback(t) = NaN;
        data.itiought_DEC_feedback_end(t) = NaN;

        data.acttime_blankstart_blankend(t) = NaN;
        data.acttime_OBS_blankend_start(t) = NaN;
        data.acttime_OBS_start_highlight(t) = NaN;
        data.acttime_OBS_highlight_spin(t) = NaN;
        data.acttime_OBS_spin_feedback(t) = NaN;
        data.acttime_OBS_feedback_cj1(t) = NaN;
        data.acttime_OBS_cj1_cj1fbstart(t) = NaN;
        data.acttime_OBS_cj1fbstart_cj1fbend(t) = NaN;
        data.acttime_OBS_cj1fbend_blank(t) = NaN;
        data.acttime_OBS_blank_end(t) = NaN;
        data.acttime_DEC_blankend_start(t) = NaN;
        data.acttime_DEC_start_resp(t) = NaN;
        data.acttime_DEC_resp_startcj2(t) = NaN;
        data.acttime_DEC_startcj2_cj2(t) = NaN;
        data.acttime_DEC_cj2_cj2fbstart(t) = NaN;
        data.acttime_DEC_cj2fbstart_cj2fbend(t) = NaN;
        data.acttime_DEC_cj2fbend_spin(t) = NaN;
        data.acttime_DEC_spin_feedback(t) = NaN;
        data.acttime_DEC_feedback_end(t) = NaN;

        data.itiought_blankstart_blankend(t) = setting.iti.blank;
        data.acttime_blankstart_blankend(t) = round(1000*(timestamp_blankend-timestamp_blankstart));
        if data.metatype(t)==1
            data.itiought_OBS_start_highlight(t) = setting.iti.OBS.start_highlight;
            data.itiought_OBS_highlight_spin(t) = setting.iti.OBS.highlight_spin;
            data.itiought_OBS_spin_feedback(t) = 25 * (setting.iti.spin+jitter);
            data.itiought_OBS_cj1fbstart_cj1fbend(t) = setting.iti.OBS.cj1fbstart_cj1fbend;
            data.itiought_OBS_blank_end(t) = setting.iti.OBS.blank_end;
            data.acttime_OBS_blankend_start(t) = round(1000*(timestamp_startobs-timestamp_blankend));
            data.acttime_OBS_start_highlight(t) = round(1000*(timestamp_highlightobs-timestamp_startobs));
            data.acttime_OBS_highlight_spin(t) = round(1000*(timestamp_spinobs-timestamp_highlightobs));
            data.acttime_OBS_spin_feedback(t) = round(1000*(timestamp_feedbackobs-timestamp_spinobs));
            data.acttime_OBS_feedback_cj1(t) = round(1000*(timestamp_cj1-timestamp_feedbackobs));
            data.acttime_OBS_cj1_cj1fbstart(t) = round(1000*(timestamp_cj1fbstart-timestamp_cj1));
            data.acttime_OBS_cj1fbstart_cj1fbend(t) = round(1000*(timestamp_cj1fbend-timestamp_cj1fbstart));
            data.acttime_OBS_cj1fbend_blank(t) = round(1000*(timestamp_blankobs-timestamp_cj1fbend));
            data.acttime_OBS_blank_end(t) = round(1000*(timestamp_endobs-timestamp_blankobs));
        else
            data.itiought_DEC_cj2fbstart_cj2fbend(t) = setting.iti.DEC.cj2fbstart_cj2fbend;
            data.itiought_DEC_spin_feedback(t) = 25 * (setting.iti.spin+jitter);
            data.itiought_DEC_feedback_end(t) = setting.iti.DEC.feedback_end;
            data.acttime_DEC_blankend_start(t) = round(1000*(timestamp_startdec-timestamp_blankend));
            data.acttime_DEC_start_resp(t) = round(1000*(timestamp_resp-timestamp_startdec));
            data.acttime_DEC_resp_startcj2(t) = round(1000*(timestamp_startcj2-timestamp_resp));
            data.acttime_DEC_startcj2_cj2(t) = round(1000*(timestamp_cj2-timestamp_startcj2));
            data.acttime_DEC_cj2_cj2fbstart(t) = round(1000*(timestamp_cj2fbstart-timestamp_cj2));
            data.acttime_DEC_cj2fbstart_cj2fbend(t) = round(1000*(timestamp_cj2fbend-timestamp_cj2fbstart));
            data.acttime_DEC_cj2fbend_spin(t) = round(1000*(timestamp_spindec-timestamp_cj2fbend));
            data.acttime_DEC_spin_feedback(t) = round(1000*(timestamp_feedbackdec-timestamp_spindec));
            data.acttime_DEC_feedback_end(t) = round(1000*(timestamp_enddec-timestamp_feedbackdec));
        end

        data.timing(t,:) = datestr(clock);
        save(fullfile(setting.datapath,['Exp1_' num2str(subNo) '.mat']),'-v6');

    end



    % calculate stats
    testlist.cumdec(block,:,:) = [cumdec; NaN(maxtrials-size(cumdec,1),2)];
    where = find(sum(setting.ntrialsall(blockorder(1:(block-1))))+1:sum(setting.ntrialsall(blockorder(1:block)))) ...
            + sum(setting.ntrialsall(blockorder(1:(block-1))));
    rt = data.rt(where);
    testlist.rt(block) = round(mean(rt(rt > 0 & ~isnan(rt))));
    resp = data.resp(where);
    testlist.resp(block,:) = [length(find(resp==1)) length(find(resp==2))];
    whichobs = data.whichobs(where);
    testlist.cumwhichobs(block,:) = [length(find(whichobs==1)) length(find(whichobs==2))];
    err = data.err(where);
    testlist.err(block) = mean(err(~isnan(err)));
    cj1 = data.cj1(where);
    testlist.cj1(block,:) = [mean(cj1(~isnan(cj1) & whichobs==1)) ...
                             mean(cj1(~isnan(cj1) & whichobs==2))];
    cj1rt = data.cj1rt(where);
    testlist.cj1rt(block,:) = [round(mean(cj1rt(cj1rt > 0 & ~isnan(cj1rt) & whichobs==1))) ...
                               round(mean(cj1rt(cj1rt > 0 & ~isnan(cj1rt) & whichobs==2)))];
    prob = data.prob(where);
    testlist.prob(block,:) = [mean(prob(~isnan(prob) & whichobs==1)) ...
                              mean(prob(~isnan(prob) & whichobs==2))];
    outcome = data.outcome(where);
    testlist.outcomeprob(block,:) = [sum(outcome(whichobs==1))/length(outcome(whichobs==1)) ...
                                     sum(outcome(whichobs==2))/length(outcome(whichobs==2))];
    cj2 = data.cj2(where);
    testlist.cj2(block) = mean(cj2(~isnan(cj2)));
    cj2rt = data.cj2rt(where);
    testlist.cj2rt(block) = round(mean(cj2rt(cj2rt > 0 & ~isnan(cj2rt))));
    testlist.outcome(block) = data.cumoutcomeblock(end);
    testlist.win(block) = data.cumwinblock(end);
    testlist.obswin(block) = data.cumobswinblock(end);

    testlist.blockend(block,:) = datestr(clock);

    save(fullfile(setting.datapath,['Exp1_' num2str(subNo) '.mat']),'-v6');


    blockstart = datevec(testlist.blockstart(block,:),'dd-mmm-yyyy HH:MM:SS');
    blockend = datevec(testlist.blockend(block,:),'dd-mmm-yyyy HH:MM:SS');
    mins = fix(etime(blockend,blockstart)/60);
    secs = mod(etime(blockend,blockstart),60);

    % Show feedback for the experimenter on screen:
    fprintf('THE PARTICIPANT JUST COMPLETED BLOCK %2.0f (%d min %d sec):\n',block,mins,secs);
    fprintf('LEFT BANDIT:\n');
    fprintf('Real value was %1.2f, observed value %1.2f, and rated value %1.2f.\n',testlist.banditsmean(block,1)*100,testlist.outcomeprob(block,1),testlist.prob(block,1));
    fprintf('This value was judged with an average confidence of %1.2f after %4.0f ms.\n',testlist.cj1(block,1),testlist.cj1rt(block,1));
    fprintf('This bandit was observed %d times.\n',testlist.cumwhichobs(block,1));
    fprintf('RIGHT BANDIT:\n');
    fprintf('Real value was %1.2f, observed value %1.2f, and rated value %1.2f.\n',testlist.banditsmean(block,2)*100,testlist.outcomeprob(block,2),testlist.prob(block,2));
    fprintf('This value was judged with an average confidence of %1.2f after %4.0f ms.\n',testlist.cj1(block,2),testlist.cj1rt(block,2));
    fprintf('This bandit was observed %d times.\n',testlist.cumwhichobs(block,2));
    fprintf('The left bandit was chosen %d time(s), and the right was chosen %d time(s).\n',testlist.resp(block,1),testlist.resp(block,2));
    fprintf('The Participant chose the bandit with the higher value on %1.2f%% of the trials.\n',1-testlist.err(block));
    fprintf('The Participant won %d points and observed %d points won.\n',testlist.win(block),testlist.obswin(block));
    fprintf('Choices were made with an average RT of %4.0f ms.\n',testlist.rt(block));
    fprintf('Choices were made with at average confidence of %1.2f after %4.0f ms.\n\n',testlist.cj2(block),testlist.cj2rt(block));

    block=block+1;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Thanks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Beeper(261.63,.4,1);
Screen('TextSize',w,20);
endscr = 'Thank you for your participation!';
rct = CenterRectOnPoint(Screen('TextBounds',w,endscr),rect(3)/2,rect(4)/2);
Screen('DrawText',w,endscr,rct(1),rct(2),255);
txtline = ['You have won ' int2str(data.cumwin(end)) ' points in this experiment (£' num2str(data.cumwin(end)*0.001,'%1.2f') ').'];
DrawFormattedText(w, txtline, 'center', rect(4)/2-50, 255, [], [], [], 1.5);

% Show stimulus on screen
Screen('Flip',w);
WaitSecs(2.000);
KbWait;

Screen('Flip',w);





%MATLAB can now receive input from the keyboard again
ListenChar(0);
KbWait;
    % Cleanup at end of experiment - Close window, show mouse cursor, close
    % result file, switch Matlab/Octave back to priority 0 -- normal
    % priority:
    Screen('CloseAll');
    Priority(0);
    toc
    % End of experiment:
    return;
catch
    % catch error: This is executed in case something goes wrong in the
    % 'try' part due to programming error etc.:

    % Do same cleanup as at the end of a regular session...
    Screen('CloseAll');
    Priority(0);
    toc
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
end % try ... catch %
