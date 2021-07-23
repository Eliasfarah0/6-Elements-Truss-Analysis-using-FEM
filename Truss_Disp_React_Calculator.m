clear
clc

% Connectivity Matrix
con_mat = [1 2 ; 1 3 ; 2 3 ; 2 4 ; 3 4 ; 3 5 ; 4 5];

% Number of Elements
num_ele = size(con_mat,1);

% Number of Nodes
num_nod = 5;

% Elements Degree of Freedom (DOF) 
ele_dof = [1 2 3 4 ; 1 2 5 6 ; 3 4 5 6 ; 3 4 7 8 ; 5 6 7 8 ; 5 6 9 10 ; 7 8 9 10];

% Assign zero values for the Displacement Matrix
disp_mat  = zeros(2*num_nod,1); 

% Assign zero values for the Reaction Forces Matrix
force_mat = zeros(2*num_nod,1);

% Assign zero values for the Assembled Stiffness Matrix
ass_stiff_mat = zeros(2*num_nod);

% User input the Length of the bars (Vertical & Horizental):
Length = input ('Enter a value for the length L of the bars (m): ');
fprintf ('\n')

% Nodes Coordinates
nod_coor = [0 0 ; 0 Length ; Length 0 ; Length Length ; 2*Length Length];

% User input the Area of the bars:
Area = input ('Enter a value for the cross-sectional area A of the bars (m^2): ');
fprintf ('\n')

% User input the Modulus of Elasticity of the bars:
Mod_Elasticity = input ('Enter a value for the modulus of elasticity E of the bars (N/m^2): ');
fprintf ('\n')

% User input Applied Loads Px = F3 & Py = F4 at node 2:
force_mat(3,1) = input ('Enter a value for the load Px applied on node 2 of the truss (N): ');
fprintf ('\n')
force_mat(4,1) = input ('Enter a value for the load Py applied on node 2 of the truss (N): ');
fprintf ('\n')

% Applied Loads F7 = 0 & F8 = 0 on node 4
force_mat(7,1) = 0;
force_mat(8,1) = 0;
force_mat(9,1) = 0; % F9_tangential = 0

% Boundary conditions u1 = v1 = u3 = v3 = 0
disp_mat (1,1)  = 0;
disp_mat (2,1)  = 0;
disp_mat (5,1)  = 0;
disp_mat (6,1)  = 0;
disp_mat (10,1)  = 0; % U_hat_10 = 0

% Elements Length Matrix
L(num_ele,1) = 0;

% Computation of the Assembled Stiffness Matrix
for u = 1:num_ele
    % Length of each element L = sqrt ( (X(node i+1) - X(node i))^2 + (Y(node i+1) - Y(node i))^2 )
    L(u , 1)      =  sqrt( (nod_coor(con_mat(u,2),1) - nod_coor(con_mat(u,1),1))^2  +  (nod_coor(con_mat(u,2),2) - nod_coor(con_mat(u,1),2))^2 );
    
    % cos(THETA)  = ( X(node i+1) - X(node i) ) / L
    cos           =  ( nod_coor(con_mat(u,2),1) - nod_coor(con_mat(u,1),1) )  /  L(u);
    
    % sin(THETA)  = ( Y(node i+1) - Y(node i) ) / L
    sin           =  ( nod_coor(con_mat(u,2),2) - nod_coor(con_mat(u,1),2) )  /  L(u);
    
    % Element Stiffness Matrix calculator
    ele_stiff_mat =  ( ((Area * Mod_Elasticity) / L(u))  *  [cos*cos cos*sin -cos*cos -cos*sin;cos*sin sin*sin -cos*sin -sin*sin;-cos*cos -cos*sin cos*cos cos*sin;-cos*sin -sin*sin cos*sin sin*sin]);
    
    % Extract the rows from ele_dof for each element and place them in a row matrix
    ele_dof_row   =  ele_dof(u , :);
    
    for i = 1:4
        for j = 1:4
            % Assemble the Stiffness Matrix of every element into the Assembled Stiffness Matrix
            ass_stiff_mat( ele_dof_row(1,i) , ele_dof_row(1,j) ) = ass_stiff_mat( ele_dof_row(1,i) , ele_dof_row(1,j) ) + ele_stiff_mat(i,j);
        end
    end
end

% DIRECT APPROACH - START
% Creating an Identity Matrix R
R = eye(2*num_nod);

% Adding to Node 5 (DOF 9 & 10) cos(theta = 45) = sin(theta = 45) = sqrt(2)/2
for i = 9:10
    for j = 9:10
    R(i,j) = R(i,j) + (sqrt(2)/2);
    end
end

% Adding the minus sign to Node 5 (DOF)
R(9,10) = -1 * R(9,10);

% Computing the Transformed Assembled Stiffness Matrix --> [K_hat] = [R_transpose] * [K] * [R] 
ass_stiff_mat = R.' * ass_stiff_mat * R;
% DIRECT APPROACH - END

% Known Forces on the nodes F3, F4, F7, F8 & F10
known_force_mat = [3 ; 4 ; 7 ; 8 ; 9];

% Assign zero values for the Reduced Displacement Matrix
red_disp_mat (size(known_force_mat) , 1) = 0;

% Assign zero values for the Reduced Reaction Forces Matrix
red_force_mat(size(known_force_mat) , 1) = 0;

% Assign zero values for the Reduced Assembled Stiffness Matrix
red_ass_stiff_mat(size(known_force_mat) , size(known_force_mat)) = 0;

% Filling the Reduced Displacement & Reaction Forces Matrix
for i = 1:size(known_force_mat , 1)
    red_disp_mat(i , 1)   =   disp_mat  (known_force_mat(i,1) , 1);
    red_force_mat(i , 1)  =   force_mat (known_force_mat(i,1) , 1);
end

% Filling the Reduced Assembled Stiffness Matrix
for i = 1:size(known_force_mat , 1)
    for j = 1:size(known_force_mat , 1)
        red_ass_stiff_mat(i , j) = ass_stiff_mat( known_force_mat(i , 1) , known_force_mat(j , 1) );
    end
end

% Solving the Reduced Linear System to Find the Displacements of each node ([K] * [U] = [F] --> [U] = [K] \ [F]) 
red_disp_mat = red_ass_stiff_mat \ red_force_mat;

% Replace the computed displacements in the original displacement matrix
for i = 1:size(known_force_mat  ,1)
    disp_mat(known_force_mat(i , 1) , 1) = red_disp_mat(i , 1); 
end

% Known Displacement on the nodes U1, U2, U5, U6 & U9
known_disp_mat = [1 ; 2 ; 5 ; 6 ; 10];

% Calculating the Unknown Reaction Forces on the remaining Nodes ([F] = [K] * [U])
for i = 1:size(known_disp_mat , 1)
    force_mat(known_disp_mat(i , 1) , 1) = ass_stiff_mat(known_disp_mat(i , 1) , :) * disp_mat; 
end

% Finding the Displacements U9 & U10 in the x-y plane using U_hat_9 & U_hat_10 
disp_mat (9,1)    =  disp_mat (9,1)  * (sqrt(2)/2) - disp_mat (10,1)  * (sqrt(2)/2);
disp_mat (10,1)   =  disp_mat (9,1)  * (sqrt(2)/2) + disp_mat (10,1)  * (sqrt(2)/2);

% Finding the Forces F9 & F10 in the x-y plane using F9_tangential & F10_normal
force_mat (9,1)   =  force_mat (9,1) * (sqrt(2)/2) - force_mat (10,1) * (sqrt(2)/2);
force_mat (10,1)  =  force_mat (9,1) * (sqrt(2)/2) + force_mat (10,1) * (sqrt(2)/2);

% Dipslay the Outputs
disp('The Assembled Stiffness Matrix is:')
disp(ass_stiff_mat)

disp('The Displacements at each Node are:')
for i = 1:2*num_nod
    X = ['U' , num2str(i) , ' = ' , num2str(disp_mat(i) * 10^3 , 5) , ' mm'];
    disp(X)
end
fprintf('\n')

disp('Reaction Forces on each Node are:')
        for i = 1:2*num_nod
    X = ['F' , num2str(i) , ' = ' , num2str(force_mat(i) * 10^-3 , 5) , ' kN'];
    disp(X)
        end
fprintf('\n')



% New Coordinates of the nodes after applying the Loads
new_nod_coor(num_nod , 2) = 0;
k = 0;
for i = 1:num_nod
    for j = 1:2
    k = k + 1;
    new_nod_coor(i,j) = nod_coor(i,j) + disp_mat(k,1);
    end
end

% Initial Truss plot with no Deformation (Blue) VS Final Truss plot with Deformation (Red)
for i = 1:num_ele
    x_old = [nod_coor(con_mat(i,1)     , 1) nod_coor(con_mat(i,2)     , 1)];
    y_old = [nod_coor(con_mat(i,1)     , 2) nod_coor(con_mat(i,2)     , 2)];
    x_new = [new_nod_coor(con_mat(i,1) , 1) new_nod_coor(con_mat(i,2) , 1)];
    y_new = [new_nod_coor(con_mat(i,1) , 2) new_nod_coor(con_mat(i,2) , 2)];
    plot(x_old , y_old , 'b' , x_new , y_new , 'r')
    title('Truss Initial (Blue) & Final (Red) Plot')
    xlabel ('x-axis (m)')
    ylabel ('y-axis (m)')
    hold on
    axis image
end