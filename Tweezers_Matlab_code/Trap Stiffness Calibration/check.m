function check(number,r1,r2,tt,color1)
%checks if number is between r1 and r2
if (number <= r1  | number > r2 )
    bot         =   0.485;
    wid         =  0.6; 
    left1       =   .25;
    height      =   0.045;
    str = ['Value for ' tt ' out of range. Range is   ' num2str(r1) ' -  ' num2str(r2)];
   OK = uicontrol('Style','text', ...
    'Position',[left1 bot 1.1*wid height], ...
    'String',str,...
     'BackgroundColor',[1 .4 .4],'Fontsize',.3);

    pause(4)
    set(OK,'visible','off')
    
    error('ErrorTests:convertTest', ...
        ['Error message from program tweezercalib: \n Value for ',tt,' out of range. Range is > ',num2str(r1),' to <= ',num2str(r2)])
end;