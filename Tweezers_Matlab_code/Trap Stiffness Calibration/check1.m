function flag = check1(number,color1)
%checks if number is positive or not
if (number <= 0)
    bot         =   0.485;
    wid         =  0.6; 
    left1       =   .25;
    height      =   0.045;
    str = {'Number of data points in a block must be larger than zero'};

    OK = uicontrol('Style','text', ...
        'Position',[left1 bot 1.1*wid height], ...
        'String',str,...
        'BackgroundColor',[1 .4 .4],'Fontsize',.3);
    pause(4)
    set(OK,'visible','off')
    error('ErrorTests:convertTest', ...
        ['Error message from program weezercalib: \n Value must be larger than zero '])
end

