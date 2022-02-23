meshOutPath='/work/ollie/orichter/MisomipPlus/io0029/fesommesh/1081.08/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0029/uadata/1081.08-Nodes9778-Ele19300-Tri3-kH1000-MismipPlus-io0029.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0029/fesommesh/meshgen.goodfile.1081.08';

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
