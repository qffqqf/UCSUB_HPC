function [ Nodes ] = calculate_nodenumbers( Direction, origin, coordinates)
%% Calc_Node_numbers
% This function automates assigning the nodes according to their location:
% top, bottom, left, right, topleft, topright, bottomleft & bottomright
% for a beam model unit cell.

% This function reads the coordinates matrix and iterates over the
% nodes as to find the appropiate group.
%
% In the second part of this function, the nodes to be visualised are
% assigned.

dispScatterPlot = false;
%% Coordinates transformation: The standard code is written in the assumption that
%there is wave propagation in the xz-plane (along positive x & z axis) so if your unit cell has been
%constructed in Patran for propagation in the xy- or yz-plane the
%code will transform the coordinates to an xz-orientation.  The code
%assumes the unit cell is always made in a right-handed system.

% The Direction.L1 direction = direction of Nodes.left
% The Direction.L2 direction = direction of Nodes.bottom
% The Direction.L3 direction = direction of Nodes.bottomleft

% User input can be: '+x','-x','+y','-y','+z','-z'

% origin is the local origin of the unit cell

%% coordinates
%column 1: number of node
%column 2: x-coordinates
%column 3: y-coordinates
%column 4: z-coordinates

% Check whether all Nodes are numbered in inceasing order.
if length(coordinates(:,1))~=coordinates(end,1)
    error('Max node number is larger then number of nodes in the model')
end

%% Translation of local origin of structure to global origin
% Nodes.origin = origin;
Origin_coord=coordinates(origin,[2:4]);

A = coordinates;
A(:,[2:4])=A(:,[2:4])-ones(coordinates(end,1),1)*Origin_coord;
A_temp = abs(A(:,[2,3,4]));

%A tolerance is used so that nodes that only differ by a slight margin are
%still included in the overlap
T = sortrows(roundn(A,-5),[2,3,4]);
Tol = abs(T(2,4)-T(1,4))/10;%10
% CLAUS: takes smalles step in z-direction and divide that by 10; this is assumed to be the tolerance

%We want to achieve following order Direction.L1 = Z, Direction.L2 = X, Direction.L3 = y
B_temp = A_temp;

if strcmp(Direction.L1,'+x')
    B_temp = A_temp(:,[2,3,1]);%x column should be z column
elseif strcmp(Direction.L1,'+y')  %x column should be y column
    B_temp = A_temp(:,[3,1,2]);%x column should be z column
elseif strcmp(Direction.L1,'+z')  %ok!
else
    error('Directions where not valid')
end


A_temp = B_temp;
A = [A(1:end,1) A_temp(1:end,:)];
%A = roundn(A,-5);   %to solve problems with rounding in Patran


%% Min and max search in coordinates to detect edges later on
Extremes.MinX = min(A(:,2));
Extremes.MaxX = max(A(:,2));
Extremes.MinY = min(A(:,3));
Extremes.MaxY = max(A(:,3));
Extremes.MinZ = min(A(:,4));
Extremes.MaxZ = max(A(:,4));


%% Filter out edges in y-direction of unit cell
% disp('Calculating corner edges...')
% Bottom Left: x is min, z is min , y decreases
Nodes.bottomleft = findEdge( A, Extremes, [2 4 -3] );

% Bottom Right: x is min, z is max , y decreases
Nodes.bottomright = findEdge( A, Extremes, [2 -4 -3] );

% Top Left: x is max, z is min pos, y decreases
Nodes.topleft = findEdge( A, Extremes, [-2 4 -3] );

% Top Right: x is max, z is max pos, y decreases
Nodes.topright = findEdge( A, Extremes, [-2 -4 -3] );


%% Filter out planes of bottom, top, left & right and delete double nodes already in corner edges
% disp('Finding planes...')
% Bottom plane: min x value
Nodes.bottom = findPlane (A, Nodes, 2);

% Top plane: max x value
Nodes.top = findPlane (A, Nodes, -2);

% Left plane: min z value
Nodes.left = findPlane (A, Nodes, 4);

% Right plane: max z value
Nodes.right = findPlane (A, Nodes, -4);



%% Check if there are facing nodes for nodes in planes and delete nodes that have none
% disp('Deleting non-facing nodes...')
% Bottom vs top
toKeepA = [];
toKeepB = [];

for i=1:length(Nodes.bottom)
    for j=1:length(Nodes.top)
        if abs(A(A(:,1)==Nodes.bottom(i),4) - A(A(:,1)==Nodes.top(j),4)) <= 2*Tol && abs(A(A(:,1)==Nodes.bottom(i),3) - A(A(:,1)==Nodes.top(j),3)) <= 2*Tol
            toKeepA = [toKeepA i];
            toKeepB = [toKeepB j];
        end
    end
end

Nodes.bottom = Nodes.bottom(:,toKeepA);
Nodes.top = Nodes.top(:,toKeepB);

% left vs right
toKeepA = [];
toKeepB = [];

for i=1:length(Nodes.left)
    for j=1:length(Nodes.right)
        if abs(A(A(:,1)==Nodes.left(i),2) - A(A(:,1)==Nodes.right(j),2)) <= 2*Tol && abs(A(A(:,1)==Nodes.left(i),3) - A(A(:,1)==Nodes.right(j),3)) <= 2*Tol
            toKeepA = [toKeepA i];
            toKeepB = [toKeepB j];
        end
    end
end
Nodes.left = Nodes.left(:,toKeepA);
Nodes.right = Nodes.right(:,toKeepB);

%% Get remaining, interior nodes
Nodes.interior = setdiff([1:length(coordinates(:,1))],[Nodes.left,Nodes.right,Nodes.bottom,Nodes.top,Nodes.bottomleft,Nodes.bottomright,Nodes.topleft,Nodes.topright]);

%% Construct matrices with coordinates of planes
% disp('Sorting node vectors...')
% Bottom
bottomCoord = zeros(length(Nodes.bottom),4);
for i=1:length(Nodes.bottom)
    bottomCoord(i,:) = A(A(:,1)==Nodes.bottom(i),:);
end

% Top
topCoord = zeros(length(Nodes.top),4);
for i=1:length(Nodes.top)
    topCoord(i,:) = A(A(:,1)==Nodes.top(i),:);
end

% Left
leftCoord = zeros(length(Nodes.left),4);
for i=1:length(Nodes.left)
    leftCoord(i,:) = A(A(:,1)==Nodes.left(i),:);
end

% Right
rightCoord = zeros(length(Nodes.right),4);
for i=1:length(Nodes.right)
    rightCoord(i,:) = A(A(:,1)==Nodes.right(i),:);
end

%% Sort the Nodes vectors of the planes again, to make sure they are in correct order

% Bottom
temp = sortrows(bottomCoord, [2 -3 4]);
Nodes.bottom = temp(:,1)';

% Top
temp = sortrows(topCoord, [2 -3 4]);
Nodes.top = temp(:,1)';

% Left
temp = sortrows(leftCoord, [4 -3 2]);
Nodes.left = temp(:,1)';

% Right
temp = sortrows(rightCoord, [4 -3 2]);
Nodes.right = temp(:,1)';


%% Check boundary conditions
% disp('Checking boundary conditions...')
% The boundary conditions will be checked for logical errors

% Check bottom and top

%A tolerance is used so that nodes that only differ by a slight margin are
%still included in the overlap
i = 1;
while i <= length(Nodes.bottom)
    node1 = Nodes.bottom(i);
    node2 = Nodes.top(i);
    
    % check if same y value or smaller than a tolerance
    % search for row with correct node number
    row1 = find(A(:,1)==node1);
    row2 = find(A(:,1)==node2);
    if abs(A(row1,3) - A(row2,3)) <= Tol
    else
        Diff = abs(A(row1,3) - A(row2,3))
        error('Error in boundary conditions for bottom vs top!');
    end
    
    % check if same z value or smaller than 0.5mm
    if abs(A(row1,4) - A(row2,4)) <= Tol
    else
        Diff = abs(A(row1,4) - A(row2,4))
        error('Error in boundary conditions for bottom vs top!');
    end
    
    % check if difference in x is length of the unit cell
    lengthX = min(A);
    lengthX = abs(lengthX(2));
    diffCheck = abs(A(row1,2) - A(row2,2));
    if diffCheck >= (lengthX - Tol)
    elseif diffCheck <= (lengthX + Tol)
    else
        Diff = abs(A(row1,2) - A(row2,2))
        error('Error in boundary conditions for bottom vs top!');
    end
    
    i = i + 1;
end


% Check bottom left and bottom right

i = 1;
while i <= length(Nodes.bottomleft)
    node1 = Nodes.bottomleft(i);
    node2 = Nodes.bottomright(i);
    row1 = find(A(:,1)==node1);
    row2 = find(A(:,1)==node2);
    
    % check if same y value or smaller than 0.5mm
    if abs(A(row1,3) - A(row2,3)) <= Tol
    else
        Diff = abs(A(row1,3) - A(row2,3))
        error('Error in boundary conditions for bottomLeft vs bottomRight!');
    end
    
    % check if same x value or smaller than 0.5mm
    if abs(A(row1,2) - A(row2,2)) <= Tol
    else
        Diff = abs(A(row1,2) - A(row2,2))
        error('Error in boundary conditions for bottomLeft vs bottomRight!');
    end
    
    % check if difference in z is length of the unit cell
    lengthZ = min(A);
    lengthZ = abs(lengthZ(4));
    diffCheck = abs(A(row1,4) - A(row2,4));
    if diffCheck >= (lengthZ - Tol)
    elseif diffCheck <= (lengthZ + Tol)
    else
        Diff = abs(A(row1,4) - A(row2,4))
        error('Error in boundary conditions for bottomLeft vs bottomRight!');
    end
    
    i = i + 1;
end



% Check left and right

i = 1;
while i <= length(Nodes.bottomleft)
    node1 = Nodes.left(i);
    node2 = Nodes.right(i);
    row1 = find(A(:,1)==node1);
    row2 = find(A(:,1)==node2);
    
    % check if same y value or smaller than 0.5mm
    if abs(A(row1,3) - A(row2,3)) <= Tol
    else
        Diff = abs(A(row1,3) - A(row2,3))
        error('Error in boundary conditions for left vs right!');
    end
    
    % check if same x value or smaller than 0.5mm
    if abs(A(row1,2) -A(row2,2)) <= Tol
    else
        Diff = abs(A(row1,2) -A(row2,2))
        error('Error in boundary conditions for left vs right!');
    end
    
    % check if difference in z is length of the unit cell
    lengthZ = min(A);
    lengthZ = abs(lengthZ(4));
    diffCheck = abs(A(row1,4) - A(row2,4));
    if diffCheck >= (lengthZ - Tol)
    elseif diffCheck <= (lengthZ + Tol)
    else
        Diff = abs(A(row1,4) - A(row2,4))
        error('Error in boundary conditions for left vs right!');
    end
    
    i = i + 1;
end


% Check bottom left and top left

i = 1;
while i <= length(Nodes.bottomleft)
    node1 = Nodes.bottomleft(i);
    node2 = Nodes.topleft(i);
    row1 = find(A(:,1)==node1);
    row2 = find(A(:,1)==node2);
    % check if same y value or smaller than 0.5mm
    if abs(A(row1,3) - A(row2,3)) <= Tol
    else
        Diff = abs(A(row1,3) - A(row2,3))
        error('Error in boundary conditions for bottomLeft vs topLeft!');
    end
    
    % check if same z value or smaller than 0.5mm
    if abs(A(row1,4) - A(row2,4)) <= Tol
    else
        Diff = abs(A(row1,4) - A(row2,4))
        error('Error in boundary conditions for bottomLeft vs topLeft!');
    end
    
    % check if difference in x is length of the unit cell
    lengthX = min(A);
    lengthX = abs(lengthX(2));
    diffCheck = abs(A(row1,2) - A(row2,2));
    if diffCheck >= (lengthX - Tol)
    elseif diffCheck <= (lengthX + Tol)
    else
        Diff = abs(A(row1,2) - A(row2,2))
        error('Error in boundary conditions for bottomLeft vs topLeft!');
    end
    
    i = i + 1;
end

% disp('...Calculation node numbers finished')

%% Make scatter plot for the boundary conditions

if dispScatterPlot == true
    %bottom
    x_bottom = [];
    y_bottom = [];
    z_bottom = [];
    for i=1:length(Nodes.bottom)
        x_bottom = [x_bottom; A(A(:,1)==Nodes.bottom(i),4)];
        y_bottom = [y_bottom; A(A(:,1)==Nodes.bottom(i),2)];
        z_bottom = [z_bottom; A(A(:,1)==Nodes.bottom(i),3)];
    end
    
    %top
    x_top = [];
    y_top = [];
    z_top = [];
    for i=1:length(Nodes.top)
        x_top = [x_top; A(A(:,1)==Nodes.top(i),4)];
        y_top = [y_top; A(A(:,1)==Nodes.top(i),2)];
        z_top = [z_top; A(A(:,1)==Nodes.top(i),3)];
    end
    
    %left
    x_left = [];
    y_left = [];
    z_left = [];
    for i=1:length(Nodes.left)
        x_left = [x_left; A(A(:,1)==Nodes.left(i),4)];
        y_left = [y_left; A(A(:,1)==Nodes.left(i),2)];
        z_left = [z_left; A(A(:,1)==Nodes.left(i),3)];
    end
    
    %right
    x_right = [];
    y_right = [];
    z_right = [];
    for i=1:length(Nodes.right)
        x_right = [x_right; A(A(:,1)==Nodes.right(i),4)];
        y_right = [y_right; A(A(:,1)==Nodes.right(i),2)];
        z_right = [z_right; A(A(:,1)==Nodes.right(i),3)];
    end
    
    %bottomright
    x_bottomright = [];
    y_bottomright = [];
    z_bottomright = [];
    for i=1:length(Nodes.bottomright)
        x_bottomright = [x_bottomright; A(A(:,1)==Nodes.bottomright(i),4)];
        y_bottomright = [y_bottomright; A(A(:,1)==Nodes.bottomright(i),2)];
        z_bottomright = [z_bottomright; A(A(:,1)==Nodes.bottomright(i),3)];
    end
    
    %bottomleft
    x_bottomleft = [];
    y_bottomleft = [];
    z_bottomleft = [];
    for i=1:length(Nodes.bottomleft)
        x_bottomleft = [x_bottomleft; A(A(:,1)==Nodes.bottomleft(i),4)];
        y_bottomleft = [y_bottomleft; A(A(:,1)==Nodes.bottomleft(i),2)];
        z_bottomleft = [z_bottomleft; A(A(:,1)==Nodes.bottomleft(i),3)];
    end
    %topright
    x_topright = [];
    y_topright = [];
    z_topright = [];
    for i=1:length(Nodes.topright)
        x_topright = [x_topright; A(A(:,1)==Nodes.topright(i),4)];
        y_topright = [y_topright; A(A(:,1)==Nodes.topright(i),2)];
        z_topright = [z_topright; A(A(:,1)==Nodes.topright(i),3)];
    end
    
    %topleft
    x_topleft = [];
    y_topleft = [];
    z_topleft = [];
    for i=1:length(Nodes.topleft)
        x_topleft = [x_topleft; A(A(:,1)==Nodes.topleft(i),4)];
        y_topleft = [y_topleft; A(A(:,1)==Nodes.topleft(i),2)];
        z_topleft = [z_topleft; A(A(:,1)==Nodes.topleft(i),3)];
    end
    
    figure('name','Scatter Nodes Boundary Condition')
    hold on
    scatter3(x_bottom,y_bottom,z_bottom,10,'r')
    scatter3(x_top,y_top,z_top,10,'b')
    scatter3(x_left,y_left,z_left,10,'k')
    scatter3(x_right,y_right,z_right,10,'m')
    scatter3(x_bottomright,y_bottomright,z_bottomright,10,'g')
    scatter3(x_bottomleft,y_bottomleft,z_bottomleft,10,'g')
    scatter3(x_topleft,y_topleft,z_topleft,10,'g')
    scatter3(x_topright,y_topright,z_topright,10,'g')
    legend('Bottom','Top','Left','Right','Corner edges')
    
    
    figure('name','Scatter Nodes Boundary Condition')
    hold on
    scatter3(coordinates(Nodes.bottom,2),coordinates(Nodes.bottom,3),coordinates(Nodes.bottom,4),10,'r')
    scatter3(coordinates(Nodes.top,2),coordinates(Nodes.top,3),coordinates(Nodes.top,4),10,'y')
    scatter3(coordinates(Nodes.left,2),coordinates(Nodes.left,3),coordinates(Nodes.left,4),10,'k')
    scatter3(coordinates(Nodes.right,2),coordinates(Nodes.right,3),coordinates(Nodes.right,4),10,'m')
    scatter3(coordinates(Nodes.bottomright,2),coordinates(Nodes.bottomright,3),coordinates(Nodes.bottomright,4),10,'g')
    scatter3(coordinates(Nodes.bottomleft,2),coordinates(Nodes.bottomleft,3),coordinates(Nodes.bottomleft,4),10,'g')
    scatter3(coordinates(Nodes.topleft,2),coordinates(Nodes.topleft,3),coordinates(Nodes.topleft,4),10,'g')
    scatter3(coordinates(Nodes.topright,2),coordinates(Nodes.topright,3),coordinates(Nodes.topright,4),10,'g')
    legend('Bottom','Top','Left','Right','Corner edges')
    scatter3(coordinates(:,2),coordinates(:,3),coordinates(:,4),10,'b')
    axis equal
    xlabel('L1')
    ylabel('L2')
    zlabel('L3')
    %%
end
end