function [ btContainer ] = editLines( axesPos, figBord, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%btLinePropsCbStr = ['uisetlineprops(findall(' num2str(gca,99) ',''type'',''line''))'];
if exist('btContainer')
    for ii = 1:size(btContainer,2)
        delete(btContainer(ii));
        
    end
    clear btContainer;
end

lHand = findall(gca ,'type','line');
for ii = 1:size(lHand,1)


btLinePropsPos = [axesPos(1)+figBord(1)*axesPos(3)+0.003+(ii-1)*0.015,axesPos(2)+axesPos(4)-figBord(1)+0.013,0.01/0.7,0.01/0.5];
 
% Note: all the following code is just to have a specific cursor
% ^^^^ (HAND_CURSOR) when hovering over the button...
btLineprops = com.mathworks.mwswing.MJButton(num2str(ii));
btLineprops.setBorder([]);
btLineprops.setBackground(java.awt.Color.white);
btLineprops.setCursor(java.awt.Cursor(java.awt.Cursor.HAND_CURSOR));
btLineprops.setFlyOverAppearance(true);
btLineprops.setToolTipText('Modify properties of plot lines');
[dummy,btContainer(ii)] = javacomponent(btLineprops,[0 0 1 1],gcf); 
set(btLineprops, 'ActionPerformedCallback',['uisetlineprops(findall(' num2str(lHand(ii),99) ',''type'',''line''))']);
set(btContainer(ii), 'Units','Norm', 'Position',btLinePropsPos);
end

end

