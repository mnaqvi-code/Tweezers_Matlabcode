val = get(hfilter,'Value');

figure(1);

if val == 7
    ht  = uicontrol('Style','text','Position',[left3-edge bot9 wid1+wid2+3*edge height],'String',...
        ' ','BackgroundColor',color);
    jh = findjobj(ht);
    jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
    ht  = uicontrol('Style','text','Position',[left3 bot9 (wid1+wid2)/6 height],'String',...
        ' ','BackgroundColor',color);
    jh = findjobj(ht);
    jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
elseif val <= 6 && val > 1
    Nfilter  =   val - 1;
    
    ht  = uicontrol('Style','text','Position',[left3-edge bot9 wid1+wid2+3*edge height],'String',...
        ' ','BackgroundColor',color);
    jh = findjobj(ht);
    jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
    ht  = uicontrol('Style','text','Position',[left3 bot9 (wid1+wid2)/6 height],'String',...
        'f3dB (kHz):','BackgroundColor',color1);
    jh = findjobj(ht);
    jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    
% Default suggested values, in kHz     
    f3dB1 = 50;
    f3dB2 = 80;
    f3dB3 = 80;
    f3dB4 = 80;
    f3dB5 = 80;
   
    for i = 1 : Nfilter
        eval(['hf3dB' num2str(i) '= uicontrol(''Style'',''edit'',' ...
                '''Position'',[left3+(wid1+wid2+edge)/6*i bot9 (wid1+wid2)/6 height],' ...
                '''String'',num2str(f3dB' num2str(i) '),'...
                '''BackgroundColor'',color2,'...
                '''Callback'',''f3dB' num2str(i) ' = str2num(get(hf3dB' num2str(i) ',''''String'''')); '');']);
    end;
    
elseif val == 1
    Nfilter  =   0;
    if prgVr.setupStt == 0
        elec_filters    =   '(1 ./ (1 + (f/22000).^16))';
    elseif prgVr.setupStt == 1
        elec_filters    =   '(1 ./ (1 + (f/22000).^16))';
    elseif prgVr.setupStt == 2
        elec_filters    =   '(1 ./ (1 + (f/22789).^2))';
    end
    ht  = uicontrol('Style','text','Position',[left3-edge bot9 wid1+wid2+3*edge height],'String',' ','BackgroundColor',color);
    jh = findjobj(ht);
    jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
    hfilter_formula  = uicontrol('Style','edit', ...
        'Position',[left3 bot9 wid1+wid2+edge height], ...
        'String',elec_filters, ...
        'Callback','elec_filters = get(hfilter_formula,''String'');', ...
        'BackgroundColor',color2);
    
end;

