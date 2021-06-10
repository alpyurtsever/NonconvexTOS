function [A, B, OPT, SOL] = qapread(fname)
% This file reads QAPLIB datasets. 

%% Preamble
% Optimal value and Solution will output empty if there is no output file.
OPT = [];
SOL = [];

% If the dataname includes the ".dat" filetype, remove it
if length(fname) >= 4
    if strcmp(fname(end-3:end),'.dat')
        fname(end-3:end) = [];
    end
end

%% Load data file
FID = fopen([fname,'.dat'], 'r');
if (FID == -1); error('File cannot be opened.'); end

TLINE = strtrim(fgetl(FID));
DIMENSION = sscanf(TLINE, '%d');

A = zeros([DIMENSION,DIMENSION]);
for INDEX = 1:DIMENSION
    ROW = [];
    while numel(ROW) < DIMENSION
        TLINE = strtrim(replace(fgetl(FID),',',' '));
        ROW = [ROW; sscanf(TLINE, '%f')]; %#ok
    end
    A(INDEX,:) = ROW';
end
 
B = zeros([DIMENSION,DIMENSION]);
for INDEX = 1:DIMENSION
    ROW = [];
    while numel(ROW) < DIMENSION
        TLINE = strtrim(replace(fgetl(FID),',',' '));
        ROW = [ROW; sscanf(TLINE, '%f')]; %#ok
    end
    B(INDEX,:) = ROW';
end

fclose(FID);

% Use sparse matrices if data matrices are sparse

if nnz(A)/numel(A) <= 0.25
    A = sparse(A);
end

if nnz(B)/numel(B) <= 0.25
    B = sparse(B);
end

%% Load solution file
if exist([fname,'.sln'],'file')
    FID = fopen([fname,'.sln'], 'r');
    if (FID == -1); error('File cannot be opened.'); end
    
    TLINE = strtrim(fgetl(FID));
    ROW = sscanf(TLINE, '%f');
    DIMENSION = round(ROW(1));
    OPT = ROW(2);
    
    ROW = [];
    while numel(ROW) < DIMENSION
        TLINE = strtrim(replace(fgetl(FID),',',' '));
        ROW = [ROW; sscanf(TLINE, '%f')]; %#ok
    end
    
    if any(ROW == 0)
        ROW = ROW + 1;
    end
    
    SOL = sparse((1:DIMENSION)',ROW,ones(DIMENSION,1),DIMENSION,DIMENSION);
    
    fclose(FID);
    
    % Check the convention used and fix potential data transpose issues
    if iprod(A'*SOL,SOL*B) ~= OPT
        if iprod(A*SOL,SOL*B) == OPT
            A = A';
        elseif iprod(B'*SOL,SOL*A) == OPT
            TMP = A;
            A = B;
            B = TMP;
        end
    end

end

end