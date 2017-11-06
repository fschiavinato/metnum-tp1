function show_depths()

imagen = 'buda';
direcciones = '034';
metodo = 'eg';
%% 
fname = strcat(imagen,'.',direcciones,'.',metodo,'.depth');
Z = dlmread(fname);

profundidad = figure
imagesc(Z)
colormap(gray)
%% 

saveas(profundidad, fname, 'png');
%fclose(fname);
end
