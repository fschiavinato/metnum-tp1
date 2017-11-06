function show_normales(imagen, direcciones)

imagen = 'buda';
direcciones = '034';
%% 
filename_nx = strcat(imagen,'.',direcciones,'.x.normal');
filename_ny = strcat(imagen,'.',direcciones,'.y.normal');
filename_nz = strcat(imagen,'.',direcciones,'.z.normal');

x=dlmread(filename_nx); N=x;
y=dlmread(filename_ny); N(:,:,2)=y;
z=dlmread(filename_nz); N(:,:,3)=z;

[height,width,~] = size(N);

[X,Y] = meshgrid(1:width,1:height);
%%
p = 8; % paso / modificar para ver más o menos datos
%N = N(1:p:end, 1:p:end,:);
%X = X(1:p:end, 1:p:end,:);
%Y = Y(1:p:end, 1:p:end,:);
%quiver3(Z2, N(:,:,1),N(:,:,2),N(:,:,3))
h = figure;
quiver(X,-Y,N(:,:,1),N(:,:,2))

%% 
saveas(h,strcat(imagen,'.',direcciones,'.','normales'),'png');
%quiver(N(:,:,1),N(:,:,2))

end