meshOutPath='/work/ollie/orichter/MisomipPlus/io0019/fesommesh/1035.10/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0019/uadata/1035.10-Nodes9013-Ele17776-Tri3-kH1000-MismipPlus-io0019.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0019/fesommesh/meshgen.goodfile.1035.10';

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
