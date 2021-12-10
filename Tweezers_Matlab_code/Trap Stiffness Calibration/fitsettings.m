disp(' ');
addpath('D:\Dropbox\AMOLF PhD\New Tweezers Setup\Data Analysis\Literature\Fitting script\');

cd('D:\Dropbox\AMOLF PhD\New Tweezers Setup\Data Analysis\Literature\Fitting script\');
disp(['pwd = ' pwd]);
disp(' ');


set(0,         'defaultTextColor', [0 0 0],...
              'defaultAxesXColor', [0 0 0],...
              'defaultAxesYColor', [0 0 0],...
              'defaultAxesZColor', [0 0 0],...
          'defaultPatchEdgeColor', [0 0 0],...
    'defaultFigureInvertHardcopy', 'on',...
             'defaultFigureColor', [0.8 0.7 1],...
               'defaultAxesColor', [1 1 1],...
          'defaultFigureColormap', jet(64),...
        'defaultSurfaceEdgeColor', [0 0 0],...
            'defaultAxesPosition', [0.12  0.12  0.8  0.8 ],...
           'defaultAxesXTickMode', 'auto',...
           'defaultAxesYTickMode', 'auto',...
           'defaultAxesZTickMode', 'auto',...
          'defaultFigurePosition', [234 102 841 793],...
         'defaultFigurePaperType', 'A4',...
        'defaultFigurePaperUnits', 'centimeters',...
           'defaultAxesFontUnits', 'normalized',...
            'defaultAxesFontSize', 0.04,...
          'defaultAxesFontWeight', 'bold',...
           'defaultAxesLineWidth', 0.8,...
      'defaultAxesLineStyleOrder', {'-o',    '-s',   '-d',    '--o',    '--s',    '--d',    ':o',    ':s',    ':d'},...
          'defaultAxesColorOrder', [    0.4     0.8     0.1 ;
                                        0.5     0.2     0.8 ;
                                        1.0     0.5     0.3 ;
                                        0.5     0.6     1.0;
                                        0.5     0.3     0.4 ;
                                        1.0     0.0     0.0],...
                 'defaultAxesBox', 'on',...
          'defaultAxesTickLength', [0.0180 0.0250],...
           'defaultTextFontUnits', 'normalized',...
            'defaultTextFontSize', 0.04,...
          'defaultTextFontWeight', 'bold',...
           'defaultLineLineWidth', 1.5,...
               'defaultLineColor', [0.5 0.2 0.8],...
          'defaultLineMarkerSize', 4,...
     'defaultLineMarkerFaceColor', [0.4     0.8     0.1],...
     'defaultLineMarkerEdgeColor', [0 0 0],...
          'defaultPatchFaceColor', [0.5  0.6  1]);

disp(' Ready to fit!'); disp(' ');
