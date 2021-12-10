clear all;
close all;
clc;
figure;

% open('20110212012T4-Time.fig');
%  h = gcf; %current figure handle
%  axesObjs = get(h, 'Children');  %axes handles
%       dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%       objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
%       time = get(dataObjs, 'XData');  %data from low-level grahics objects
% pLength = get(dataObjs, 'YData');
% close all;
% clear axesObjs dataObjs objTypes h
% figure;
% 
% open('20110212012T4.fig');
%  h = gcf; %current figure handle
%  axesObjs = get(h, 'Children');  %axes handles
%       dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%       objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
%       time = get(dataObjs, 'XData');  %data from low-level grahics objects
% pLength = get(dataObjs, 'YData');
% time = time{4};
% pLength = pLength{4};

open('GroELTest20060818028.fig');
 h = gcf; %current figure handle
 axesObjs = get(h, 'Children');  %axes handles
      dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
      objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
      time = get(dataObjs, 'XData');  %data from low-level grahics objects
pLength = get(dataObjs, 'YData');

close all;
clear axesObjs dataObjs objTypes h


figure;
plot(time,pLength);
force = pLength;
l = length(force);
a= [];

f1 = 0.5; % threshold for force jumps defining start and end points
f2 = 1; % threshold for step selection, minimal force difference (between min and max in interval [start end])
n = 25; % number of frames that correspond to half a second, upper limit for duration of unfolding event.

% First find steps where the force decreases and then encircle the
% unfolding event [startpoint endpoint]

i = 2;
while i <l-1
    if force(i)-force(i+1)> f1
       endpoint = i;
       while force(endpoint)-force(endpoint+1)>0
           endpoint = endpoint + 1;
           if endpoint-i >n || endpoint>l-1
               display('unfolding event too long, or end of trace reached')
               break
           end
       end
       startpoint = i;
       while force(startpoint-1)-force(startpoint)>0
           startpoint = startpoint - 1;
           if i - startpoint>n || startpoint <2
                display('unfolding event too long, or end of trace reached')
               break
           end
       end
       a = [a; startpoint endpoint];
       i = endpoint+1;
    else
       i = i+1;
    end
end

%check if points were found
if isempty(a)
    points = [];
    disp('no steps found')
    return
end

la = size(a);

% now try to figure out which steps are actual steps
% unfolding should consist of at least two frames, and the force step must
% be big enough, and force before unfolding above 6 pN
points = [];


for ii =1:la
    if abs(a(ii,2) - a(ii,1)) > 1 && abs(force(a(ii,1)) - force(a(ii,2))) > f2 && force(a(ii,2))>6
        points = [points; a(ii,:)];
    end
end

hold on;
for ii = 1:length(points)
    hold on;
    plot([time(points(ii,1)) time(points(ii,2))], [force(points(ii,1)) force(points(ii,2))], 'linewidth', 2, 'color', 'k')
end
hold off;