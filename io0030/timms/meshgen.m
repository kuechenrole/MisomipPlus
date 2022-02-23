meshOutPath='/work/ollie/orichter/MisomipPlus/io0030/fesommesh/1043.10/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0030/uadata/1043.10-Nodes8744-Ele17234-Tri3-kH1000-MismipPlus-io0030.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0030/fesommesh/meshgen.goodfile.1043.10';

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
