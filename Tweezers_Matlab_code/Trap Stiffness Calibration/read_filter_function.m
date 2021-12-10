if isvalid(hfilter)
    val = get(hfilter,'Value');
else
    val = 1;
end

if val == 7,
    elec_filters = 1;
    want_filter = 0;
elseif val <= 6 & val > 1,
    Nfilter = val - 1;
    elec_filters = [];
    want_filter = 1;
    for i = 1 : Nfilter,
        eval(['f3dB' num2str(i) ' = str2num(get(hf3dB' num2str(i) ',''String''));'])
    elec_filters = [elec_filters '(1 ./ (1 + (f*1.e-3/' eval(['num2str(f3dB' num2str(i) ')']) ').^2)) .* '];
    end; %(for i = 1 : Nfilter)
    elec_filters = elec_filters(1:length(elec_filters)-3);
    
else
    if isvalid(hfilter)
        elec_filters = get(hfilter_formula,'String');
    else
       elec_filters    =   '(1 ./ (1 + (f/22000).^16))'; 
    end
    want_filter = 1;
    
end; %(if val <= 6 & val > 1)