% SPECTRUM_HISTOGRAM calculates a histogram of the spectrum divided by its
% fit. Should show an exponential distribution.

check1(nblock,color1);  %Checks if nblock is positive nonzero

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%
if ~exist('chi2X') %fit_powerspectrum; end;
htfit  = uicontrol('Style','text', ...
    'Position',[left3+(wid1+edge) bot12 wid2 height], ...
    'String','Activate ''Fit power spectrum'' before clicking this button',...
    'BackgroundColor',color1);
else

  for ixy = IXY % IXY = [1 2] if there exist 2 channels; IXY = 1 if there is only 1 channel
    
    if ixy == 1
        P = Px{jj}; Title = 'X';
    else
        P = Py{jj}; Title = 'Y';
    end;

    if ixy == 1 
        color_now = color_x;
        figure(16); clf; set(gcf,'Numbertitle','off','Name','Phist, X'); hold on; h = title('Final Fit, X');plot_Phist(f,P,scal_fitX,parameters0X,color_now);
        set(h,'FontWeight','Bold')
    else 
        color_now = color_y;
        figure(17); clf; set(gcf,'Numbertitle','off','Name','Phist, Y'); hold on; h = title('Final Fit, Y');plot_Phist(f,P,scal_fitY,parameters0Y,color_now);
        set(h,'FontWeight','Bold')
    end;
end; %for ixy = IXY

end; %if ~exist('parametersX')
