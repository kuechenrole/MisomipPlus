meshOutPath='/work/ollie/orichter/MisomipPlus/io0005/fesommesh/1028.00/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0005/uadata/0100000-Nodes8158-Ele16068-Tri3-kH1000-MismipPlus-io0001.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0005/fesommesh/meshgen.goodfile.1028.00';

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
