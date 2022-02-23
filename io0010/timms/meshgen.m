meshOutPath='/work/ollie/orichter/MisomipPlus/io0010/fesommesh/1040.60/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0010/uadata/1040.60-Nodes9285-Ele18310-Tri3-kH1000-MismipPlus-io0010.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0010/fesommesh/meshgen.goodfile.1040.60';

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
