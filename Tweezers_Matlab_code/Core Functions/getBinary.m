function [ f ] = getBinary( filePath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filePath);
f = {};
f.fileVersion = fread(fid, 1, 'int','l');
strlen = fread(fid, 1,'uchar','l');
fileName = fread(fid, strlen(1),'uchar','l');
f.fileName = char(fileName');
strlen = fread(fid, 1,'uchar','l');
sampleFreqString = fread(fid, strlen(1),'uchar','l');
sampleFreqString = char(sampleFreqString');
sampleFreq = fread(fid, 1,'double','l');
f.SampleFrequency = sampleFreq(1);
strlen = fread(fid, 1,'uchar','l');
userName = fread(fid, strlen(1),'uchar','l');
f.userName = char(userName');
strlen = fread(fid, 1,'uchar','l');
machineName = fread(fid, strlen(1),'uchar','l');
f.machineName = char(machineName');
%keyboard;


f.date = fread(fid, 1,'uint64','l');

calibrationParametersLabels = struct('beadRadius', [], 'temperature', [], 'viscosity', [], 'diffCoeffTheo', [], 'dragCoeffTheo', []);
fields = fieldnames(calibrationParametersLabels);

for ii = 1:numel(fields)
    calibrationParametersLabels.(fields{ii}) = fread(fid, 1,'double','l');
end

f.calibrationParametersLabels.beadRadius = calibrationParametersLabels.beadRadius*1E9;
f.calibrationParametersLabels.viscosity = calibrationParametersLabels.viscosity*1E-6;
f.calibrationParametersLabels.temperature = calibrationParametersLabels.temperature-273;

calibrationFitLabels = struct('alphaSinusoidal', [], 'beta', [], 'calibrationQuality', [], 'cornerFrequency', [], 'diffCoeffExp', [], 'dragCoeffExp', [], 'drivingAmplitude', [], 'drivingFrequency', [],...
    'forceFactor', [], 'peakHeight', [], 'trapStiffnessSinusoidal', [], 'stiffness', []);

fields = fieldnames(calibrationFitLabels);

for ii = 1:4
    for jj = 1:numel(fields)
        f.calibrationFitLabels.(fields{jj})(ii) = fread(fid, 1,'double','l');
    end
end

for ii = 1:4
    offsets(ii) = fread(fid, 1,'double','l');
end

f.measurementVariablesLabels = struct('time', [], 'index', [], 'PSD1VxDiff', [], 'PSD1VxSum', [], 'PSD1VyDiff', [], 'PSD1VySum', [], 'PSD2VxDiff', [], 'PSD2VxSum', [], 'PSD2VyDiff', [], 'PSD2VySum', [], 'MirrorX', [], 'MirrorY', [], 'Status', []);
fields = fieldnames(f.measurementVariablesLabels);
f.beadVariablesLabels = struct('time', [], 'frame', [], 'Bead1X', [], 'Bead1Y', [], 'Bead2X', [], 'Bead2Y', []);
fields2 = fieldnames(f.beadVariablesLabels);
checkByte = [];

ii = 1;
numIndx(ii) = ii;
checkByte(ii) = fread(fid, 1, 'int','l');
testRead = [];
while ~isempty(checkByte(ii))
    
    if checkByte(ii) == 17
        numIndx(end+1) = ii;
        
        if exist('startTimeFile')
            startTimeFile(end+1) = fread(fid, 1, 'double','l');
            startTimeInacc(end+1) = startTimeFile(end);
            startTimeFile(end) = startTimeFile(end-1) + 1/sampleFreq;
        else
            startTimeFile(1) = fread(fid, 1, 'double','l');
            startTimeInacc(1) = startTimeFile(1);
        end          
        
        if exist('nSamples')
            nSamples(end+1) = fread(fid, 1, 'int','l');            
        else
            nSamples(1) = fread(fid, 1, 'int','l');
        end
        
        startTimeFile(end:end+nSamples(end)-1) = startTimeFile(end):1/sampleFreq:startTimeFile(end)+(nSamples(end)-1)/sampleFreq;
        testRead = [testRead fread(fid, [12 nSamples(end) ], 'int','l')];
        
        
    else
        f.beadVariablesLabels.(fields2{1})(ii) = fread(fid, 1, 'double','l');
        beadFrametemp = fread(fid, 1,'uint64','l');        
        if ~isempty(beadFrametemp)
            f.beadVariablesLabels.(fields2{2})(ii) = beadFrametemp;
        end

        for jj = 3:numel(fields2)
            f.beadVariablesLabels.(fields2{jj})(ii) = fread(fid, 1,'double','l');
        end

    end
    ii = ii + 1;
    
    try
        checkByte(ii) = fread(fid, 1, 'int','l');
    catch exception
        break;
    end 
    
end 
startTimeFile = startTimeFile - startTimeFile(1);
if ~isempty(f.beadVariablesLabels.Bead1X)
    [iii, jjj] = find(f.beadVariablesLabels.Bead1X ~= 0);
    for jj = 1:numel(fields2)
        f.beadVariablesLabels.(fields2{jj}) = f.beadVariablesLabels.(fields2{jj})(jjj);
        f.beadData(:,jj) = f.beadVariablesLabels.(fields2{jj})';
    end 
    
    f.Data = [startTimeFile' testRead'];
else
    f.Data = [startTimeFile' testRead'];
end

fclose(fid);
numIndx = numIndx(2:end);
return

end

