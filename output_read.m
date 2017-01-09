close all; clear; clc;

%% General info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('outputaid.txt', 'r');

% Cartesian topology
scrap = fscanf(fileID, '%s %s *'); % scrap strings
cart_top = fscanf(fileID, '%d', 3);
nprocs = prod(cart_top);

% topology reconstruction
procgrid = zeros(cart_top(1), cart_top(2), cart_top(3));
% procgrid(:) = 1:nprocs;     % ATTENZIONE: qui i processori vanno da 1 a nprocs

for k = 1 : cart_top(1)
    for j = 1 : cart_top(2)
        for i = 1 : cart_top(3)
            procgrid(i, j, k) = (i-1)*cart_top(1)*cart_top(2) + (j-1)*cart_top(2) + k;
        end
    end
end

% Arrays of shapes
% 1 ind. - process number
% 2 ind. - direction
% 3 ind. - variable (u, v, w, p)
shapes = zeros(nprocs, 3, 4);

% arrays of sizes
% 1 ind. - process number
% 2 ind. - variable (u, v, w, p)
sizes = zeros(nprocs, 4);

% Domain dimensions
scrap = fscanf(fileID, '%s %s *'); % scrap strings
L = fscanf(fileID, '%f', 3);

% shapes
for j = 1 : 4                          % for each component
    scrap = fscanf(fileID, '%s %s *'); % scrap strings
    for i = 1:nprocs                   % for each process
        shapes(i, :, j) = fscanf(fileID, '%i', 3);
    end
end

% sizes
for j = 1 : 4
    sizes(:, j) = shapes(:, 1, j).*shapes(:, 2, j).*shapes(:, 3, j);
end

% starts
starts = 0*sizes;
starts(1, :) = 1;
for i = 2 : nprocs
    starts(i, :) = sum(sizes(1:i-1, :), 1)+1;
end

fclose(fileID);


%% u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;

% reading
fid = fopen('u', 'r');
u = fread(fid, 'double');
fclose(fid);

uproc = cell(nprocs, 1);

% temporary 3D array to store the variable
for i = 1 : nprocs-1
    uproc{i} = reshape(u(starts(i, j) : starts(i+1, j)-1), shapes(i, :, j));
end
i = nprocs;
uproc{i} = reshape(u(starts(i, j) : end), shapes(i, :, j));

clear u;

U = data_assemble(uproc, procgrid);

clear uproc;


%% v %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 2;

% reading
fid = fopen('v', 'r');
v = fread(fid, 'double');
fclose(fid);

vproc = cell(nprocs, 1);

% temporary 3D array to store the variable
for i = 1 : nprocs-1
    vproc{i} = reshape(v(starts(i, j) : starts(i+1, j)-1), shapes(i, :, j));
end
i = nprocs;
vproc{i} = reshape(v(starts(i, j) : end), shapes(i, :, j));

clear v;

V = data_assemble(vproc, procgrid);

clear vproc;


%% w %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 3;

% reading
fid = fopen('w', 'r');
w = fread(fid, 'double');
fclose(fid);

wproc = cell(nprocs, 1);

% temporary 3D array to store the variable
for i = 1 : nprocs-1
    wproc{i} = reshape(w(starts(i, j) : starts(i+1, j)-1), shapes(i, :, j));
end
i = nprocs;
wproc{i} = reshape(w(starts(i, j) : end), shapes(i, :, j));

clear w;

W = data_assemble(wproc, procgrid);

clear wproc;


%% p %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 4;

% reading
fid = fopen('p', 'r');
p = fread(fid, 'double');
fclose(fid);

pproc = cell(nprocs, 1);

% temporary 3D array to store the variable
for i = 1 : nprocs-1
    pproc{i} = reshape(p(starts(i, j) : starts(i+1, j)-1), shapes(i, :, j));
end
i = nprocs;
pproc{i} = reshape(p(starts(i, j) : end), shapes(i, :, j));

clear p;

P = data_assemble(pproc, procgrid);

clear pproc;


%% Grid generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTENZIONE: per ora funziona solo con griglie equispaziate


%%%%%%%%%% u %%%%%%%%%%

xu = linspace(0, L(1), size(U, 1));

yu = zeros(1, size(U, 2));
dy = L(2)/(size(U, 2)-2);
yu(1) = 0;                  yu(end) = L(2);
yu(2:end-1) = linspace(dy/2, L(2)-dy/2, size(U, 2)-2);

zu = zeros(1, size(U, 3));
dz = L(3)/(size(U, 3)-2);
zu(1) = 0;                  zu(end) = L(3);
zu(2:end-1) = linspace(dz/2, L(3)-dz/2, size(U, 3)-2);


%%%%%%%%%% v %%%%%%%%%%

xv = zeros(1, size(V, 1));
dx = L(1)/(size(V, 1)-2);
xv(1) = 0;                  xv(end) = L(1);
xv(2:end-1) = linspace(dx/2, L(1)-dx/2, size(V, 1)-2);

yv = linspace(0, L(2), size(V, 2));

zv = zeros(1, size(V, 3));
dz = L(3)/(size(V, 3)-2);
zv(1) = 0;                  zv(end) = L(3);
zv(2:end-1) = linspace(dz/2, L(3)-dz/2, size(V, 3)-2);


%%%%%%%%%% w %%%%%%%%%%

xw = zeros(1, size(W, 1));
dx = L(1)/(size(W, 1)-2);
xw(1) = 0;                  xw(end) = L(1);
xw(2:end-1) = linspace(dx/2, L(1)-dx/2, size(W, 1)-2);

yw = zeros(1, size(W, 2));
dy = L(2)/(size(W, 2)-2);
yw(1) = 0;                  yw(end) = L(2);
yw(2:end-1) = linspace(dy/2, L(2)-dy/2, size(W, 2)-2);

zw = linspace(0, L(3), size(W, 3));





%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nz = floor(size(P, 3)/2);
Pg = P(:, :, nz);

figure
surf(Pg)


% nx = floor(size(U, 1)/2);
% nz = floor(size(U, 3)/2);
% 
% % [Yu, Xu] = meshgrid(yu, xu);
% 
% figure
% plot(U(nx, :, nz), yu, 'k'); hold on;
% % plot(U(nx, :, 16), yu, 'r');
% plot(U(nx, :, 8), yu, 'b');
% plot(U(nx, :, 4), yu, 'g');
% xlim([-1, 1]);
% xlabel('u');
% ylabel('y');
% 
% ny = floor(size(V, 1)/2);
% nz = floor(size(V, 3)/2);
% 
% figure
% plot(xv, V(:, ny, nz), 'k'); hold on;
% xlabel('x');
% ylabel('v');


% nz = floor(size(W, 3)/2);
% nz = 2
% 
% [Xw, Yw, Zw] = meshgrid(xw, yw, zw);
% 
% figure
% surf(Xw(:, :, nz), Yw(:, :, nz), W(:, :, nz)); shading interp;


% nx = floor(size(U, 1)/2);
% 
% [Xu, Yu, Zu] = meshgrid(xu, yu, zu);
% 
% 
% figure
% surf(Yu(nx, :, :), Zu(nx, :, :), U(nx, :, :));



















