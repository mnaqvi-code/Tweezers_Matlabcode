val = get(hdiode,'Value');

if val == 1
    want_alpha  =   0;
    want_diode  =   0;
elseif val == 2
    want_diode  =   1;
    want_alpha  =   0;
elseif val == 3
    want_diode  =   1;
    want_alpha  =   1;
end