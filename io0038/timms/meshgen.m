meshOutPath='/work/ollie/orichter/MisomipPlus/io0038/fesommesh/1001.07/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0038/uadata/1001.07-Nodes8046-Ele15846-Tri3-kH1000-MismipPlus-io0038.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0038/fesommesh/meshgen.goodfile.1001.07';

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
