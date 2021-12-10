function [rulerVal] = createRuler( name, kB, T, bp, bpL, pDNA, KDNA, rd, rdL, pProt,  KPROT  )
%UNTITLED2 Summary of this function goes here
%   createRuler('test', datVr.dualCst.T, datVr.dualCst.kB, datVr.basepairs, datVr.baseL, datVr.pDNA, datVr.residues, datVr.resL, 1.65E-9)

x = 0.1:0.001:1.2*((bp*bpL+max(rd)*rdL)*1E-3);
x = x.*1E-6;
%KDNA = 1.2E-9;
LDNA = bp*bpL*1E-9;
pDNA = pDNA*1E-9;
[ force ] = calcEWLC( x, kB, T, pDNA, LDNA, KDNA  );
Lp = max(rd)*rdL*1E-9;
LpCalc = (Lp-(rd.*rdL.*1E-9));
xprotGuess = [0.01:(1.5*Lp*1E6-0.01)/(length(force)-1):1.5*Lp*1E6].*(1E-6);
fMinOptions = optimset('MaxFunEvals', 2000, 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.

h = waitbar(0,'Initializing waitbar...');
set(h,'Name','Creating Ruler');
tTot = 0;
tCount = 0;
%KPROT = 2.0E-9;
for jj = 1:length(rd)
    for ii = 1:length(force)
        tic
        xprot{jj}(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, pProt, Lp-LpCalc(jj), KPROT ), xprotGuess(ii), fMinOptions);
        t = toc;
        tCount = tCount + 1;
        tTot = t + tTot;
        hours = floor((tTot/tCount)*((length(force)*length(rd))-tCount)/3600);
        mins = floor(((tTot/tCount)*((length(force)*length(rd))-tCount)-hours*3600)/60);
        secs = (tTot/tCount)*((length(force)*length(rd))-tCount)-hours*3600-mins*60;
        waitbar((ii+(jj-1)*length(force))/((length(force)*length(rd))),h,[num2str(round(hours)) ':' num2str(round(mins)) ':' num2str(round(secs)) ' remaining.'])
    end
end

constants = zeros(length(force),1);
constants(1) = kB;
constants(2) = T;
constants(3) = bp;
constants(4) = bpL;
constants(5) = pDNA;
constants(6) = KDNA;
constants(7) = rdL;
constants(8) = pProt;
constants(9) = KPROT;
for jj = 1:length(rd)
    constants(9+jj) = rd(jj);
end

force = force';

extension = [x'];
for jj = 1:length(rd)
    extension = [extension x'+xprot{jj}'];
end

rulerVal = [constants force extension];
close(h);
[filename, pathname, filterindex] = uiputfile('*.txt','Save Ruler',[pwd '\CustomRulers\']);
if isequal(filename,0) || isequal(pathname,0)
    disp('User selected Cancel')
else
    
    data.pathname = fullfile(pathname,filename);
    save([data.pathname(1:end-4) '.txt'],'rulerVal','-ascii','-double','-tabs');
    
    disp(['User saved ',fullfile(pathname,filename)])
end

end

