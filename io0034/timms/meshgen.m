meshOutPath='/work/ollie/orichter/MisomipPlus/io0034/fesommesh/1034.11/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0034/uadata/1034.11-Nodes8225-Ele16208-Tri3-kH1000-MismipPlus-io0034.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0034/fesommesh/meshgen.goodfile.1034.11';

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
