% Create a 2D subduction setup with particles & temperature
function [] = setup2D (SubductionAngle)

SubductionAngle1 = SubductionAngle(1);
SubductionAngle2 = SubductionAngle(2);
% Tell the code where the LaMEM matlab routines are 
addpath /home/chris/software/LaMEM/matlab

LaMEM_Parallel_output   =    0;

RandomNoise             =   logical(0); % add random noise to particles?

LaMEM_input_file        =   'Subduction2D.dat';

%% Compute 3D grid, depending on whether we are on 1 or >1 processors
if ~LaMEM_Parallel_output 
    % In the other case, we create a setup for 1 processor and defined the
    % parameters here. 
    % Important: the resolution you use here should be identical to what
    % is specified in then *.dat file!
    %disp(['Creating setup for 1 processor using LaMEM file: ', LaMEM_input_file])

    % Read info from LaMEM input file & create 3D grid    
    [npart_x,npart_y,npart_z,Grid,X,Y,Z,W,L,H] =   LaMEM_ParseInputFile(LaMEM_input_file);
    
    nump_x              =   Grid.nel_x*npart_x;
    nump_y              =   Grid.nel_y*npart_y;
    nump_z              =   Grid.nel_z*npart_z;
    Parallel_partition  =   [];   % since we run this on one code
    
else
    % We perform a paralel simulation; or this a 'ProcessorPartitioning'
    % file shoule be created first by running LaMEM on the desired # of
    % processors as:
    %   mpiexec -n 2 ../../bin/opt/LaMEM -ParamFile Subduction2D_FreeSlip_MATLABParticles_Linear_DirectSolver.dat -mode save_grid
   % disp(['Creating setup in parallel using LaMEM file: ', LaMEM_input_file])
      
    % Define parallel partition file
    Parallel_partition                          =   'ProcessorPartitioning_4cpu_1.2.2.bin'
    
    % Load grid from parallel partitioning file
    [npart_x,npart_y,npart_z]  =   LaMEM_ParseInputFile(LaMEM_input_file);
    [X,Y,Z,xcoor,ycoor,zcoor,Xpart,Ypart,Zpart] =   FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition,RandomNoise);
    
    % Update other variables
    nump_x  = size(X,2);
    nump_y  = size(X,1);
    nump_z  = size(X,3);
    
    % Domain parameters
    W       =   xcoor(end)-xcoor(1);    % x-dir
    L       =   ycoor(end)-ycoor(1);    % y-dir
    H       =   zcoor(end)-zcoor(1);    % z-dir
end
%==========================================================================


%%
%==========================================================================
% SPECIFY PARAMETERS OF THE SLAB
%==========================================================================

dcrusts    = 25;         % depth subducting crust
dmliths = dcrusts+50;    % depth subducting mantle lithosphere

dcrusto    = 15;         % depth oceanic crust
dmlitho = dcrusto+40;    % depth oceanic mantle lithosphere  %% must not exceed other lithospheric depths

dcrust  = 20;        % depth continental crust middle of the plate
dmlith  = dcrust+50; % depth mantle lithosphere

m_x = max(X(:))/2;

l_x = 120;
p_x = m_x - l_x;
d_x = 2*l_x;

l_sl1 = 200;
l_sl2 = 200;

lslab1 =  200;        % length slab 1
lslab2 =  200;        % length slab 2
wslab1 =  100;        % width slab 1
wslab2 =  100;        % width slab 2

% extension of weak zone 1
d_weak = 12; % depth
h_weak = 50; % horizontal extension in y direction

% properties of subducting slab
ThermalAge_Myrs     =   50;     % Thermal age of the slab in Myrs
ThicknessCrust      =   10;

ThicknessSlab       =   400;    % Thickness of slab box; 

z_surface           =   0;      % initial free surface

T_mantle            =   1350;
T_surf              =   20;     


%==========================================================================
% DEFINE SLAB
%==========================================================================
%% mantle phase
Phase               =   ones(size(X));                 
Temp                =   T_mantle*ones(size(Phase));     

BoxSides            =   [min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(Z(:))  0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 1,'TempType','linear', 'topTemp',T_mantle,'botTemp', T_mantle); % Set slab to mantle lithosphere phase

%% oceanic plates around
BoxSides            =   [min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) -dmlitho  0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 3,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);  

BoxSides            =   [min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) -dcrusto  0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 5,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle,'thermalAge',40); % Set slab to mantle lithosphere phase
  

%% Add horizontal part of slab 
BoxSides            =   [p_x (p_x+d_x) min(Y(:)) max(Y(:)) -dmliths 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 3,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);                         % Set slab to mantle lithosphere phase

BoxSides            =   [p_x (p_x+d_x) min(Y(:)) max(Y(:)) -dcrusts 0]; % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 4,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle,'thermalAge',80);               % Add crust (will override the mantle lithosphere phase above)

% weak zones
BoxSides            =   [p_x (p_x+h_weak) min(Y(:)) max(Y(:)) -d_weak 0]; % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 6); % Add crust (will override the mantle lithosphere phase above)

BoxSides            =   [p_x+d_x-h_weak p_x+d_x min(Y(:)) max(Y(:)) -d_weak 0]; % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 6);               % Add crust (will override the mantle lithosphere phase above)


% %% slab 1
RotPt               =   [p_x,0,0];
BoxSides            =   [p_x-l_sl1 p_x min(Y(:)) max(Y(:)) -dmliths 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 3, 'RotationPoint',RotPt, 'DipAngle', -SubductionAngle1,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);          % mantle lithosphere

BoxSides            =   [p_x-l_sl1 p_x min(Y(:)) max(Y(:)) -dcrusts 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 4, 'RotationPoint',RotPt, 'DipAngle', -SubductionAngle1,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle,'thermalAge',80);       	% crust (will override the slab phase above

% weak zone
BoxSides            =   [p_x-l_sl1/2 p_x min(Y(:)) max(Y(:))  -d_weak 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 6, 'RotationPoint',RotPt, 'DipAngle', -SubductionAngle1);

% %% slab 2
RotPt               =   [p_x+d_x,0,0];
BoxSides            =   [p_x+d_x p_x+d_x+l_sl2 min(Y(:)) max(Y(:))  -dmliths 0];  % [Left Right Front Back Bottom Top] of the box% [Phase,Temp]        =   AddBox(Phase,Temp,Y,X,Z,BoxSides, 3, 'RotationPoint',RotPt, 'DipAngle', SubductionAngle2,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);          % mantle lithosphere 
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 3, 'RotationPoint',RotPt, 'DipAngle', SubductionAngle2,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle,'thermalAge',80);       	% crust (will override the slab phase above

BoxSides            =   [p_x+d_x p_x+d_x+l_sl2 min(Y(:)) max(Y(:))  -dcrusts 0];  % [Left Right Front Back Bottom Top] of the box% [Phase,Temp]        =   AddBox(Phase,Temp,Y,X,Z,BoxSides, 3, 'RotationPoint',RotPt, 'DipAngle', SubductionAngle2,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);          % mantle lithosphere 
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 4, 'RotationPoint',RotPt, 'DipAngle', SubductionAngle2,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle,'thermalAge',80);       	% crust (will override the slab phase above


% weak zone
BoxSides            =   [p_x+d_x p_x+d_x+l_sl1/2 min(Y(:)) max(Y(:))  -d_weak 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 6, 'RotationPoint',RotPt, 'DipAngle', SubductionAngle2);       	% crust (will override the slab phase above


%% inner continental plate
% BoxSides            =   [p_xc (p_xc+d_xc) p_yc (p_yc+d_yc) -dmlith 0];  % [Left Right Front Back Bottom Top] of the box
% [Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 3,'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);                         % Set slab to mantle lithosphere phase
% 
% BoxSides            =   [p_xc (p_xc+d_xc) p_yc (p_yc+d_yc) -dcrust 0]; % [Left Right Front Back Bottom Top] of the box
% [Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 2);               % Add crust (will override the mantle lithosphere phase above)


%% Set Mantle Lithosphere for mantle points that have temperatures < 1200 Celcius
%ind                 =    find(Temp<1200 & Phase==1);
%Phase(ind)          =    3;

%% Add sticky air
BoxSides            =   [min(X(:)) max(X(:))  min(Y(:)) max(Y(:)) 0 max(Z(:))];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 0,  'TempType','Constant','cstTemp',T_surf);                      % sticky air with constant temperature 


%==========================================================================
% PREPARE DATA FOR VISUALIZATION/OUTPUT (no need to change this)
%==========================================================================

% Prepare data for visualization/output
A = struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);

Phase    = permute(Phase,[2 1 3]);
Temp     = permute(Temp, [2 1 3]);

% Linear vectors containing coords
x        =  X(1,:,1);
y        =  Y(:,1,1);
z        =  Z(1,1,:);
X        =  permute(X,[2 1 3]);
Y        =  permute(Y,[2 1 3]);
Z        =  permute(Z,[2 1 3]);

A.W      =  W;
A.L      =  L;
A.H      =  H;
A.nump_x =  nump_x;
A.nump_y =  nump_y;
A.nump_z =  nump_z;
A.Phase  =  Phase;
A.Temp   =  Temp;
A.x      =  x(:);
A.y      =  y(:);
A.z      =  z(:);
A.Xpart  =  X;
A.Ypart  =  Y;
A.Zpart  =  Z;
A.npart_x=  npart_x;
A.npart_y=  npart_y;
A.npart_z=  npart_z;

% PARAVIEW VISUALIZATION
%FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option

% SAVE PARALLEL DATA (parallel)
FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition);

end


