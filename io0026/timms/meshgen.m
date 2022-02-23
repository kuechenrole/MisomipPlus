meshOutPath='/work/ollie/orichter/MisomipPlus/io0026/fesommesh/1001.03/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0026/uadata/1001.03-Nodes8028-Ele15810-Tri3-kH1000-MismipPlus-io0026.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0026/fesommesh/meshgen.goodfile.1001.03';

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
