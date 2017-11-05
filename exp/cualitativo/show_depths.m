function show_depths()

%% 
fname = 'eg.txt';
Z = dlmread(fname);
%Z2 = read_depths(fname);

figure
imagesc(Z)
colormap(gray)

fclose(fp);
end
