meshOutPath='/work/ollie/orichter/MisomipPlus/io0008/fesommesh/1031.60/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0008/uadata/1031.60-Nodes8571-Ele16894-Tri3-kH1000-MismipPlus-io0008.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0008/fesommesh/meshgen.goodfile.1031.60';

cd meshgen;

disp('jigsaw2fesom');
jigsaw2fesomUa;

disp('makeTopoUa');
makeTopoUa;

disp('ensure min overlap');
ensure_minLayers_ovl;

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
