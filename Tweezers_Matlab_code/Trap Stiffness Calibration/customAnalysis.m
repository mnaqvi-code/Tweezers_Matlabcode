clear all;
close all;
clc;

extractData;
if exist('ans','var')
    data = ans;
    clear ans;
else
    break
end