%% Parameters
addpath('C:/Users/dorra/Documents/ETH/Thesis/Data/')


dem_file = 'M_SRTM3S_Projected.tif';

for i=1:length(lines)
   lines{i}.z_m = lines{i}.z_m * -1;
end

DEM = GRIDobj(dem_file);
DEM = fillsinks(DEM);

%%
figure(1)
clf
surf(DEM,'exaggerate',5)
alpha 0.5
hold on
for i=1:length(lines)
   plot3(lines{i}.x_m,lines{i}.y_m,lines{i}.z_m,'b')
   scatter3(lines{i}.x_m,lines{i}.y_m,lines{i}.z_m,'k')
end

title('Mean Posterior MHT Geometries (5x Exaggeration)')
xlabel('X [m]')
ylabel('Y [m]')