% written by Annika Boldt, 2015
% Exp2 from Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
% modulates exploration and exploitation in value-based learning. Neuroscience
% of Consciousness, 2019(1), 1–12. https://doi.org/10.1093/nc/niz004

fordatafile = 'Exp2_fordata_20150818_120757.mat'; % try creating your own randomisation and replace this file name
setting.datapath=[cd '/data'];
load(['data/',fordatafile]);

% Make sure keyboard mapping is the same on all supported operating systems
% Apple MacOS/X, MS-Windows and GNU/Linux:
KbName('UnifyKeyNames');

resetrn = sum(100*clock);
rand('state',resetrn);

setting.datapath = fullfile(pwd,'data');

blockorder = unique(fordata.block);
blockorder = blockorder(blockorder>2);
blockorder = [1:2 shuffle(blockorder)];
for i=1:length(blockorder)
    if sum(setting.ntrialsall(blockorder(1:i)))>300
        break
    end
end
setting.whenbreak = i+1; % coffee break before this block

% key assigment according to id number (both for colour and confidence)
if (mod(subNo,2))
    setting.confpoles = {'guessing','very confident'}; % very confident that correct is right
    setting.confpoles2 = {'very unconfident','guessing','very confident'}; % very confident that correct is right
    setting.mapcj = 12;
else
    setting.confpoles = {'very confident','guessing'}; % very confident that correct is left
    setting.confpoles2 = {'very confident','guessing','very unconfident'}; % very confident that correct is left
    setting.mapcj = 21;
end
setting.probpoles = {'0','100'};
setting.colours = [255 55];

% check for existing result file to prevent accidentally overwriting
% files from a previous subject/session (except for subject numbers ~= 99):
if subNo ~= 99 && fopen(fullfile(setting.datapath,['Exp2_' num2str(subNo) '.mat'])) ~= -1 && restart ~= 1
    error('Result data file already exists! Choose a different subject number.');
end

% Reseed the random-number generator for each expt.
resetrn = sum(100*clock);
rand('state',resetrn);

setting.restarted = 0;
instr.wait = 'Please wait for the experimenter.';
setting.cjq='How confident are you that you picked the better slot machine?';

setting.iti.spin = 0.02;
setting.iti.blank = 0.8;
setting.iti.OBS.start_grids = 0.5;
setting.iti.OBS.cj1fbstart_cj1fbend = 0.5;
setting.iti.OBS.blank_end = 0.5;
setting.iti.DEC.cj2fbstart_cj2fbend = 0.5;
setting.iti.DEC.feedback_end = 1.0;

setting.casinos = 1:30;
shuff = randperm(30);
setting.casinos = setting.casinos(shuff);

setting.casinonames = {'Aria (Las Vegas, USA)',...
'Atlantis Resort (Bahamas)',...
'Bellagio (Las Vegas, USA)',...
'Borgata (Atlantic City, USA)',...
'Caesars Palace (Las Vegas, USA)',...
'Casino Baden-Baden (Baden-Baden, Germany)',...
'Casino de Monte Carlo (Monaco)',...
'Casino Lisboa (Lisbon, Portugal)',...
'Casino Metropol (Moscow, Russia)',...
'Sofitel Macau at Ponte 16 (Macau, China)',...
'Casino-de-Charlevoix (La Malbaie, Canada)',...
'Cosmopolitan (Las Vegas, USA)',...
'Crown Casino (Melbourne, Australia)',...
'Foxwoods Resort Casino (Ledyard, USA)',...
'Gold Coast Hotel and Casino (Las Vegas, USA)',...
'Grand Lisboa (Macau, China)',...
'Hippodrome Casino (London, UK)',...
'Holland Casino (Amsterdam, The Netherlands)',...
'Mandalay Bay (Las Vegas, USA)',...
'Mandarin Oriental (Las Vegas, USA)',...
'MGM Grand Casino (Las Vegas, USA)',...
'Mirage (Las Vegas, USA)',...
'Red Rock Resort (Las Vegas, USA)',...
'Rio (Las Vegas, USA)',...
'Sands (Macau, China)',...
'The Palace of the Lost City (Sun City Resort, South Africa)',...
'The Casino at the Empire (London, UK)',...
'The Venetian (Macau, China)',...
'Trump Taj Mahal (Atlantic City, USA)',...
'Wynn (Las Vegas, USA)'};
