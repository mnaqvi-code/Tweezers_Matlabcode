In 'fitsettings.m' change the lines

==

addpath('/Win2K-E/Mfiles');
addpath('/Win2K-E/Mfiles/Iva-perfect');

cd('/Win2K-E/Mfiles/Iva-perfect');

==

to

==

addpath('C:/directory-of-tweezercalib/');

cd('C:/directory-of-tweezercalib/');

==

assuming that all the .m-files of the program package are
in the directory  'C:/directory-of-tweezercalib/'


The instructions for executing the program can be found in
Section 5 (Application) of the original published manuscript,
Com. Phys. Commun. 159(2004)225.
