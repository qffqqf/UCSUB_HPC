function mesh2vtk(filename, varargin)
% Required inputs:
%   element type 
%   nodes coordinates 
%   elements connectivity
% Optional inputs:
%   scalar field 
%   precision

%% Open and check
fid = fopen(filename, 'w'); 
if numel(varargin)<4 
    error('Not enough input arguments'); 
end
if any(strcmpi(varargin, 'precision')) 
    pres_id = find(strcmpi(varargin,'precision')); % find arg id
    precision = num2str(varargin{pres_id+1});
else
    precision = '10'; % default precision is 10^(-10)
end

%% Write headings
fprintf(fid, '# vtk DataFile Version 2.0\n'); 
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

%% Write nodes
nd_posi = varargin{2}; 
npoint = size(nd_posi,1);
fprintf(fid, ['POINTS ' num2str(npoint) ' float\n']);
spec = [repmat(['%0.', precision, 'f '], 1, 3), '\n'];
fprintf(fid, spec, nd_posi');

%% Write elements 
% Now only quadratic tetrahedral element is implemented
% Other types of elements to be added later...
ele_type = varargin{1};
connectivity = varargin{3};
nele = size(varargin{3},1);
switch ele_type
    case "C3D10"
        fprintf(fid,'\nCELLS %d %d\n', nele,11*nele);
        fprintf(fid,['10', repmat(' %d', 1, 10), '\n'], (connectivity-1)');
        fprintf(fid,'\nCELL_TYPES %d\n', nele);
        fprintf(fid,'%d\n', 24*ones(nele,1)); 
    case "C3D4"
        fprintf(fid,'\nCELLS %d %d\n', nele,5*nele);
        fprintf(fid,['4', repmat(' %d', 1, 4), '\n'], (connectivity-1)');
        fprintf(fid,'\nCELL_TYPES %d\n', nele);
        fprintf(fid,'%d\n', 10*ones(nele,1)); 
    case "C3D27" % 27 nodes -> 20 nodes
        fprintf(fid,'\nCELLS %d %d\n', nele,21*nele);
        fprintf(fid,['20', repmat(' %d', 1, 20), '\n'], (connectivity-1)');
        fprintf(fid,'\nCELL_TYPES %d\n', nele);
        fprintf(fid,'%d\n', 25*ones(nele,1)); % 25 is the ID for VTK_QUADRATIC_HEXAHEDRON
    case "S9" % 9 nodes -> 8 nodes
        fprintf(fid,'\nCELLS %d %d\n', nele,9*nele);
        fprintf(fid,['8', repmat(' %d', 1, 8), '\n'], (connectivity-1)');
        fprintf(fid,'\nCELL_TYPES %d\n', nele);
        fprintf(fid,'%d\n', 23*ones(nele,1)); 
end

%% Write fields (loop for different scalar fields...)
sc_id = find(strcmpi(varargin,'scalars'));
fprintf(fid, ['\nPOINT_DATA ' num2str(npoint)]);
if sc_id~=0
    for ii = 1:length(sc_id)
        title = varargin{sc_id(ii)+1};
        fprintf(fid, ['\nSCALARS ', title,' float\n']);
        fprintf(fid, 'LOOKUP_TABLE default\n');
        spec = ['%0.', precision, 'f\n'];
        fprintf(fid, spec, varargin{sc_id(ii) + 2});
    end
end
