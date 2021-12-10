%Choose the number of terms to include in the aliasing sum
val = get(alias_corr,'Value');

if (val == 1)
    n = 20;
end
if (val == 2)
    n = 30;
end
if (val == 3)
    n = 40;
end
if (val == 4)
    n = 10;
    close(figure(17));
end;