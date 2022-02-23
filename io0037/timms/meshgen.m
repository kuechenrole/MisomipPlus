meshOutPath='/work/ollie/orichter/MisomipPlus/io0037/fesommesh/1099.11/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0037/uadata/1099.11-Nodes8445-Ele16643-Tri3-kH1000-MismipPlus-io0037.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0037/fesommesh/meshgen.goodfile.1099.11';

cd meshgen;

disp('jigsaw2fesom');
jigsaw2fesomUa;

disp('makeTopoUa');
makeTopoUa;

%disp('ensure min overlap');
%ensure_minLayers_ovl;

disp('reorder');
reorder;

disp('makeFesom3d');
makeFesom3d;

pause(2);
disp('makeCavity');
makeCavity;

cd ..;

fid = fopen(goodfile_path,'w');
fid = fclose(fid);

exit;
