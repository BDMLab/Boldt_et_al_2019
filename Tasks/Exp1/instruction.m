% written by Annika Boldt, 2015
% Exp1 from Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
% modulates exploration and exploitation in value-based learning. Neuroscience
% of Consciousness, 2019(1), 1–12. https://doi.org/10.1093/nc/niz004

Screen('TextSize', w, 18);

switch block
    case 1
        if ~restart
            insimdata = imread(char('instruction/Exp1_instructions1a.jpg'));
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;

            insimdata = imread(char('instruction/Exp1_instructions1b.jpg'));
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;

            insimdata = imread(char('instruction/Exp1_instructions1c.jpg'));
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;

            insimdata = imread(char('instruction/Exp1_instructions1d.jpg'));
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;

            if setting.mapcj==12
                insimdata = imread(char('instruction/Exp1_instructions1e_12.jpg'));
            else insimdata = imread(char('instruction/Exp1_instructions1e_21.jpg'));
            end
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;

            if setting.mapcj==12
                insimdata = imread(char('instruction/Exp1_instructions1f_12.jpg'));
            else insimdata = imread(char('instruction/Exp1_instructions1f_21.jpg'));
            end
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;

            if setting.mapcj==12
                insimdata = imread(char('instruction/Exp1_instructions1g_12.jpg'));
            else insimdata = imread(char('instruction/Exp1_instructions1g_21.jpg'));
            end
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;

            insimdata = imread(char('instruction/Exp1_instructions1h.jpg'));
            texins = Screen('MakeTexture', w, insimdata);
            Screen('DrawTexture', w, texins);
            txtline = ['Experimenter, please press any mouse button to continue with block ' int2str(block) ' out of ' int2str(setting.ntotalblocksexp) '.'];
            DrawFormattedText(w, txtline, 'center', rect(4)-50, 255, [], [], [], 1.5);
            Screen('Flip',w);
            WaitSecs(2.000);
            GetClicks;
        end

        clear insimdata

    case setting.nblocksprac1+1

        Beeper(261.63,.4,.3);
        insimdata = imread(char('instruction/Exp1_instructions2a.jpg'));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins);
        DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        if setting.mapcj==12
            insimdata = imread(char('instruction/Exp1_instructions2b_12.jpg'));
        else insimdata = imread(char('instruction/Exp1_instructions2b_21.jpg'));
        end
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins);
        DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        insimdata = imread(char('instruction/Exp1_instructions2c.jpg'));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins);
        txtline = ['Experimenter, please press any mouse button to continue with block ' int2str(block) ' out of ' int2str(setting.ntotalblocksexp) '.'];
        DrawFormattedText(w, txtline, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        clear insimdata

    case setting.nblocksprac1+setting.nblocksprac2+1

        Beeper(261.63,.4,.3);
        insimdata = imread(char('instruction/Exp1_instructions3a.jpg'));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins);
        DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        insimdata = imread(char('instruction/Exp1_instructions3b.jpg'));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins);
        txtline = ['Experimenter, please press any mouse button to continue with block ' int2str(block) ' out of ' int2str(setting.ntotalblocksexp) '.'];
        DrawFormattedText(w, txtline, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        clear insimdata

    case setting.nblocksprac1+setting.nblocksprac2+setting.nblocksprac3+1

        Beeper(261.63,.4,.3);
        insimdata = imread(char('instruction/Exp1_instructions4a.jpg'));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins);
        DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        insimdata = imread(char(['stimuli/casino' int2str(setting.casinos(block-3)) '.jpg']));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins, [], CenterRectOnPoint([0 0 720 540],rect(3)/2,rect(4)/2));
        txtline = ['Welcome to the next casino: ' setting.casinonames{setting.casinos(block-3)} '!'];
        DrawFormattedText(w, txtline, 'center', 50, 255, [], [], [], 1.5);
        DrawFormattedText(w, instr.wait, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        insimdata = imread(char('instruction/Exp1_instructions4b.jpg'));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins);
        txtline = ['Experimenter, please press any mouse button to continue with block ' int2str(block) ' out of ' int2str(setting.ntotalblocksexp) '.'];
        DrawFormattedText(w, txtline, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        clear insimdata

    case setting.whenbreak

        Beeper(261.63,.4,.3);
        txtline = ['It''s time to leave ' setting.casinonames{setting.casinos(block-4)} '.'];
        rct = CenterRectOnPoint(Screen('TextBounds',w,txtline),rect(3)/2,rect(4)/2);
        Screen('DrawText',w,txtline,rct(1),rct(2)-50,255);
        txtline = ['In the ' int2str(setting.ntrialsall(blockorder(block-1))-testlist.Nobsblock(block-1)) ' choice trials, you won ' int2str(testlist.win(block-1)) ' points (£' num2str(testlist.win(block-1)*0.001,'%1.2f') ').'];
        rct = CenterRectOnPoint(Screen('TextBounds',w,txtline),rect(3)/2,rect(4)/2);
        Screen('DrawText',w,txtline,rct(1),rct(2)-25,255);
        txtline = ['So far, you have won ' int2str(data.cumwin(end)) ' points in this experiment (£' num2str(data.cumwin(end)*0.001,'%1.2f') ').'];
        rct = CenterRectOnPoint(Screen('TextBounds',w,txtline),rect(3)/2,rect(4)/2);
        Screen('DrawText',w,txtline,rct(1),rct(2),255);
        DrawFormattedText(w, 'Please press any mouse button to continue.', 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        Screen('TextSize', w, 24);
        txtline = 'It''s time for a short break now!';
        DrawFormattedText(w, txtline, 'center', rect(4)/2-50, 255, [], [], [], 1.5);
        txtline = 'Please press the button on the table and wait for the experimenter.';
        DrawFormattedText(w, txtline, 'center', rect(4)/2+50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(30.000);
        KbWait;
        Screen('TextSize', w, 18);

        txtline = 'Are you ready to continue with the remaining blocks?';
        rct = CenterRectOnPoint(Screen('TextBounds',w,txtline),rect(3)/2,rect(4)/2);
        Screen('DrawText',w,txtline,rct(1),rct(2),255);
        DrawFormattedText(w, 'Please press any mouse button to continue to the next casino.', 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        insimdata = imread(char(['stimuli/casino' int2str(setting.casinos(block-3)) '.jpg']));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins, [], CenterRectOnPoint([0 0 720 540],rect(3)/2,rect(4)/2));
        txtline = ['Welcome to the next casino: ' setting.casinonames{setting.casinos(block-3)} '!'];
        DrawFormattedText(w, txtline, 'center', 50, 255, [], [], [], 1.5);
        txtline = ['Please press any mouse button to continue with block ' int2str(block) ' out of ' int2str(setting.ntotalblocksexp) '.'];
        DrawFormattedText(w, txtline, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        clear insimdata

    otherwise

        Beeper(261.63,.4,.3);
        txtline = ['It''s time to leave ' setting.casinonames{setting.casinos(block-4)} '.'];
        rct = CenterRectOnPoint(Screen('TextBounds',w,txtline),rect(3)/2,rect(4)/2);
        Screen('DrawText',w,txtline,rct(1),rct(2)-50,255);
        txtline = ['In the ' int2str(setting.ntrialsall(blockorder(block-1))-testlist.Nobsblock(block-1)) ' choice trials, you won ' int2str(testlist.win(block-1)) ' points (£' num2str(testlist.win(block-1)*0.001,'%1.2f') ').'];
        rct = CenterRectOnPoint(Screen('TextBounds',w,txtline),rect(3)/2,rect(4)/2);
        Screen('DrawText',w,txtline,rct(1),rct(2)-25,255);
        txtline = ['So far, you have won ' int2str(data.cumwin(end)) ' points in this experiment (£' num2str(data.cumwin(end)*0.001,'%1.2f') ').'];
        rct = CenterRectOnPoint(Screen('TextBounds',w,txtline),rect(3)/2,rect(4)/2);
        Screen('DrawText',w,txtline,rct(1),rct(2),255);
        DrawFormattedText(w, 'Please press any mouse button to continue to the next casino.', 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        insimdata = imread(char(['stimuli/casino' int2str(setting.casinos(block-3)) '.jpg']));
        texins = Screen('MakeTexture', w, insimdata);
        Screen('DrawTexture', w, texins, [], CenterRectOnPoint([0 0 720 540],rect(3)/2,rect(4)/2));
        txtline = ['Welcome to the next casino: ' setting.casinonames{setting.casinos(block-3)} '!'];
        DrawFormattedText(w, txtline, 'center', 50, 255, [], [], [], 1.5);
        txtline = ['Please press any mouse button to continue with block ' int2str(block) ' out of ' int2str(setting.ntotalblocksexp) '.'];
        DrawFormattedText(w, txtline, 'center', rect(4)-50, 255, [], [], [], 1.5);
        Screen('Flip',w);
        WaitSecs(2.000);
        GetClicks;

        clear insimdata

end
