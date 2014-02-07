function post_proc(num_reg)

clc;
close force all;

%% Read in HDF5 data files
psi = cell(num_reg, 1);
dirs = cell(num_reg, 1);
wgts = cell(num_reg, 1);
for i = 1 : num_reg
    % Angular flux
    filename = sprintf('psi_reg%i', i-1);
    psi{i} = hdf5read(filename, 'dataset');
    % Quadrature directions
    filename = sprintf('dirs_reg%i', i-1);
    dirs{i} = hdf5read(filename, 'dataset');
    % Quadrature weights    
    filename = sprintf('wgts_reg%i', i-1);
    wgts{i} = hdf5read(filename, 'dataset');
end
% Region size    
reg_size = hdf5read('reg_size', 'dataset');
% Number of cells
num_cells = hdf5read('num_cells', 'dataset');
% Absorption cross section
abs_xs = hdf5read('abs_xs', 'dataset');
% External source
ext_source = hdf5read('ext_source', 'dataset');

%% Calculate balance parameter

% Incoming current
% Left boundary
inflow_left = 0;
j = 1;
for i = 1 : size(dirs{1},1)
    if dirs{1}(i, 1) > 0
        inflow_left = inflow_left+...
            psi{1}(1,1,j)*wgts{1}(i,1)*abs(dirs{1}(i,1))+...
            psi{1}(1,1,j+1)*wgts{1}(i,2)*abs(dirs{1}(i,2));
    end
    j = j + 2;
end
% Right boundary
inflow_right = 0;
j = 1;
for i = 1 : size(dirs{num_reg},1)
    if dirs{num_reg}(i, 1) < 0
        inflow_right = inflow_right+...
            psi{num_reg}(1,num_cells(num_reg)+1,j)*...
            wgts{num_reg}(i,1)*abs(dirs{num_reg}(i,1))+...
            psi{num_reg}(1,num_cells(num_reg)+1,j+1)*...
            wgts{num_reg}(i,2)*abs(dirs{num_reg}(i,2));
    end
    j = j + 2;
end

% Outgoing current
% Left boundary
outflow_left = 0;
j = 1;
for i = 1 : size(dirs{1},1)
    if dirs{1}(i, 1) < 0
        outflow_left = outflow_left+...
            psi{1}(1,1,j)*wgts{1}(i,1)*abs(dirs{1}(i,1))+...
            psi{1}(1,1,j+1)*wgts{1}(i,2)*abs(dirs{1}(i,2));
    end
    j = j + 2;
end
% Right boundary
outflow_right = 0;
j = 1;
for i = 1 : size(dirs{num_reg},1)
    if dirs{num_reg}(i, 1) > 0
        outflow_right = outflow_right+...
            psi{num_reg}(1,num_cells(num_reg)+1,j)*...
            wgts{num_reg}(i,1)*abs(dirs{num_reg}(i,1))+...
            psi{num_reg}(1,num_cells(num_reg)+1,j+1)*...
            wgts{num_reg}(i,2)*abs(dirs{num_reg}(i,2));
    end
    j = j + 2;
end

% Absorption rate and external source
abs_rate = 0;
ext_source_tot = 0;
for i = 1 : num_reg 
    cell_size = double(reg_size(i))/double(num_cells(i));
    counter = 1;
    for j = 1 : num_cells(i)
        left_phi = 0;
        right_phi = 0;
        r = 1;
        for k = 1 : size(dirs{i}, 1)           
            if dirs{i}(k,1) > 0
               left_phi = left_phi+psi{i}(2,j,r)*wgts{i}(k,1)+...
                   psi{i}(2,j,r+1)*wgts{i}(k,2);
               right_phi = right_phi+psi{i}(1,j+1,r)*wgts{i}(k,1) + ...
                   psi{i}(1,j+1,r+1)*wgts{i}(k,2);                
            else
               left_phi = left_phi+psi{i}(1,j,r)*wgts{i}(k,1)+...
                   psi{i}(1,j,r+1)*wgts{i}(k,2);
               right_phi = right_phi+psi{i}(2,j+1,r)*wgts{i}(k,1)+...
                   psi{i}(2,j+1,r+1)*wgts{i}(k,2);                
            end  
            ext_source_tot = ext_source_tot+(ext_source(i)*wgts{i}(k,1)+...
                ext_source(i)*wgts{i}(k,2))*cell_size;
            r = r + 2;
        end
        abs_rate = abs_rate + abs_xs(i)*((left_phi+right_phi)/2)*cell_size;
        phi_store{i}(counter) = (left_phi+right_phi)/2;
        counter = counter+1;
    end
end

% Balance parameter
bal = abs(inflow_right+inflow_left+ext_source_tot-outflow_right- ...
    outflow_left-abs_rate)/(inflow_right+inflow_left+ext_source_tot);
fprintf('The balance parameter is: %E\n', bal);

%% Plot psi and phi across each region
current_pos = 0;
for i = 1 : num_reg % Go through each region
    figure
    % Rearrange quadrature directions
    dirs_use = zeros(1, size(dirs{i},2)*2);
    counter = 1;
    for j = 1 : size(dirs{i},1)
        dirs_use(counter) = dirs{i}(j, 1);
        dirs_use(counter+1) = dirs{i}(j, 2);
        counter = counter+2;
    end
    % Edge locations
    edges = zeros(1, num_cells(i)+1);
    edges(1) = current_pos;
    delta = reg_size(i)/double(num_cells(i));  
    centers = zeros(1, num_cells(i));
    for j = 2 : num_cells(i)+1
        edges(j) = edges(j-1)+delta;
        centers(j-1) = edges(j-1)+delta/2;
    end
    % Plot each psi
    s = cell(1, size(psi{i}, 3));
    for j = 1 : size(psi{i}, 3)
      current_plot = plot(edges, psi{i}(1,:,j));
      set(current_plot, 'LineWidth', 1.25, 'MarkerSize', 10);
      s{j} = sprintf('\\mu=%.3f',dirs_use(j));
      hold all;                
    end
    % Set up plot options
    grid on
    legend(s);        
    set(gca, 'FontSize', 14)
    xlabel('cell edge position', 'FontSize', 18)
    y_name = sprintf('\\psi');
    ylabel(y_name, 'FontSize', 18)
    title_name = sprintf('Region %i Upstream \\psi', i);
    title(title_name, 'FontSize', 18);
    current_pos = current_pos+reg_size(i);
    % Plot phi
    figure
    current_plot = plot(centers, phi_store{i});
    set(current_plot, 'LineWidth', 1.25, 'MarkerSize', 5);
    grid on
    xlabel('cell center position', 'FontSize', 18)
    y_name = sprintf('\\phi');
    ylabel(y_name, 'FontSize', 18)
    title_name = sprintf('Region %i \\phi', i);
    title(title_name, 'FontSize', 18);
end
 
%% Plot psi at interfaces
for i = 1 : num_reg - 1
    % Rearrange quadrature directions for left region
    dirs_left = zeros(1, size(dirs{i},2)*2);
    counter = 1;
    for j = 1 : size(dirs{i},1)
        dirs_left(counter) = dirs{i}(j, 1);
        dirs_left(counter+1) = dirs{i}(j, 2);
        counter = counter+2;
    end
    % Gather left psi
    psi_left = zeros(1, size(psi{i}, 3));
    for j = 1 : size(psi{i}, 3)
        psi_left(j) = psi{i}(1, num_cells(i) + 1, j);
        
    end
    % Rearrange quadrature directions for right region
    dirs_right = zeros(1, size(dirs{i+1},2)*2);
    counter = 1;
    for j = 1 : size(dirs{i+1},1)
        dirs_right(counter) = dirs{i+1}(j, 1);
        dirs_right(counter+1) = dirs{i+1}(j, 2);
        counter = counter+2;
    end
    % Gather right psi
    psi_right = zeros(1, size(psi{i+1}, 3));
    for j = 1 : size(psi{i+1}, 3)
        psi_right(j) = psi{i+1}(1, 1, j);        
    end
    % Plot
    figure
    left_plot = plot(dirs_left, psi_left, 's');
    set(left_plot, 'LineWidth', 1.25, 'MarkerSize', 10);
    hold all;
    right_plot = plot(dirs_right, psi_right, 's');
    set(right_plot, 'LineWidth', 1.25, 'MarkerSize', 10);
    grid on
    x_name = sprintf('directional cosine \\mu');
    xlabel(x_name, 'FontSize', 18)
    y_name = sprintf('\\psi');
    ylabel(y_name, 'FontSize', 18)
    title_name = sprintf('Interface %i \\psi', i);
    title(title_name, 'FontSize', 18);
    legend_name_left = sprintf('left \\psi');
    legend_name_right = sprintf('right \\psi');
    legend(legend_name_left, legend_name_right);        

end

end