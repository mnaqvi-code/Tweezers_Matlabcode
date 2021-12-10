function check_hydro(number,r1,tt,color1)

if (number <= r1)
    bot         =   0.485;
    wid         =  0.5; 
    left1       =   .25;
    height      =   0.045;
    str = ['Value for ' tt ' out of range. Must be larger than ' num2str(r1)];
    OK = uicontrol('Style','text', ...
    'Position',[left1 bot 1.1*wid 3*height], ...
    'String',str,...
     'BackgroundColor',[1 .4 .4],'Fontsize',.3);

    pause(4)
    set(OK,'visible','off')
    error('ErrorTests:convertTest', ...
    ['Error message from program weezercalib: \n Value for ',tt,' out of range. Must be larger than ',num2str(r1)])
end;
