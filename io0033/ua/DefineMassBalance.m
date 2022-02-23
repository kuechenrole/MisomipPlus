function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

as=zeros(MUA.Nnodes,1)+0.3;

fesomMeltPath= '/work/ollie/orichter/MisomipPlus/io0033/fesomdata/io0033.1034.forcing.diag.nc';
fesomCoordPath= '/work/ollie/orichter/MisomipPlus/io0033/fesommesh/1034.11/nod2d.out';

rhofw = 1000;
rho_ice = 917;% already defined by default


wnetFes = ncread(fesomMeltPath,'wnet');
wnetFes = double(wnetFes);


fid=fopen(fesomCoordPath,'r');
n2d=fscanf(fid,'%g',1);
nodes=fscanf(fid, '%g', [4,n2d]);
fclose(fid);
xfes=transpose(nodes(2,:)*111000);
yfes=transpose(nodes(3,:)*111000);

xUa = MUA.coordinates(:,1);
yUa = MUA.coordinates(:,2);

interp = scatteredInterpolant(xfes,yfes,wnetFes,'linear','nearest');
wnetUa = interp(xUa,yUa);
wnetUa = wnetUa.*365.25*24*3600.*-1;
wnetUa = wnetUa.*(rhofw/rho_ice);

L=10e3 ;  % Smoothing length scale
[Blablub,wnetUa2]=HelmholtzEquation([],CtrlVar,MUA,1,L^2,wnetUa,0);

ab=wnetUa2.*(1-GF.node);

end
