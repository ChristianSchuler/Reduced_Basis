# Reduced Basis

With this application reduced basis solutions to Stokes flow with linear rheolgy solved by [LaMEM](https://bitbucket.org/bkaus/lamem/src/master/) can be created.
At the moment it is important, that the src files of [LaMEM](https://bitbucket.org/bkaus/lamem/src/master/) which are stored in the [LaMEM_src_files](LaMEM_src_files) are used.
It is feasible for 2D setups with linear rheology.

## Workflow

### LaMEM input file
- FSSA needs to be set to 0
- in the Petsc options sections the following option must be enabled:
   * -objects_dump
	* -dump_sol
	* -dump_precondition_matrixes
	* -dump_precondition_matrixes_prefix Ass
	* -dump_rhs
	* -dump_rho
	
### Create RB components

- All ingredients that are needed to access a reduced solution can be produced with the file [createRB.m](createRB.m).
- Steps:
   * specify LaMEM input file / filename of marker routine (e.g setup2D.m)  / number of cells in every direction / folder in which RB components should be stored
   * specify input parameter with appropriate ranges
   
   ```
   % parameter 2
   st   = 30;  % smallest parameter value
   en   = 40;  % largest parameter value
   n    = 14;  % parameter spacing
   par2 = linspace(st,en,n);
   ```
   * add parameter to list
    ```
    par = allcomb(par1,par2);
    ```
    * marker file must be adapted according to the number and type of input parameters (here for file **setup2D.m**). **par** is the list that was created by the     previous     step. In this example **par(1)** belongs to **par1** and is in **setup2D.m** treated as the variable **SubductionAngle1**
    ```
    function [] = setup2D (par)

    SubductionAngle1 = par(1);
    SubductionAngle2 = par(2);
    .
    .
    .
    end
    ```
   
   * specify tolerance for greedy algorithm
   ```
   tol  = 1e-6;    % tolerance for Greedy algorithm
   ```

### access RB solution
- the file [accessRB.m](accessRB.m) help to access a RB solution with the previously calculated components
- Steps:
  * specify LaMEM input file / filename of marker routine (e.g setup2D.m)  / number of cells in every direction / folder in which RB components were stored
  * specify if simulation is nondimensional or geo units were used.
     ```
    % units: geo or none
    unit = 'geo';
    ```
  * specify values for the parameters that were selceted in [createRB.m](createRB.m)
     ```
     par1 = 32.456;
     par2 = 36.342;
     ```
   * creates truth solution, RB solution (creaqted without DEIM) and DEIM solution
     ```
     % whole solutions
     u_t    % truth solution   
     u_RB   % RB solution      
     u_DEIM % DEIM solution
     
     % only velocity components
     v_t    % truth solution   
     v_RB   % RB solution      
     v_D    % DEIM solution
     
     ```
   * some plotting routines are added. Example:
     ```
     [V3d_t] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, v_t,'z','xz');  % z velocity is plotted in xz plane
     x     = linspace(0,coordx,nel_x);
     y     = linspace(0,coordz,nel_z+1);         % nel_z +1 velocity nodes in z direction
     [X,Y] = meshgrid(x,y);
     pcolor(X,Y,V3d_t(:,:,2).'); c = colorbar;   %   V3d_t(:,:,n) n decides about the number of the plane which should be plotted 
     ```
  ### 2D example: double subduction zone
