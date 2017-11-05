function show_normales_caballo()
x=dlmread('budax'); N=x;
y=dlmread('buday'); N(:,:,2)=y;
z=dlmread('budaz'); N(:,:,3)=z;

[height,width,~] = size(N);

[X,Y] = meshgrid(1:width,1:height);


%%
%Z = zeros(size(N,1),size(N,2));
Z2 = dlmread('eg.txt');
p = 8; % paso / modificar para ver más o menos datos
%N = N(1:p:end, 1:p:end,:);
%X = X(1:p:end, 1:p:end,:);
%Y = Y(1:p:end, 1:p:end,:);
%quiver3(Z2, N(:,:,1),N(:,:,2),N(:,:,3))
quiver(X,Y,N(:,:,1),N(:,:,2))
%quiver(N(:,:,1),N(:,:,2))

end