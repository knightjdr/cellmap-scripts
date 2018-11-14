function write_matrix_file(fid, labels_rows, labels_cols, data, varargin)

% Inputs:
%
% * |labels_cols| = |m x n| cell array of labels, where |m| is the number of
% columns and |n| is the number of labels for each column.
%
% * |labels_row| = |p x q| cell array of labels, where |p| is the number of
% rows and |q| is the number of labels for each row.

%% Check inputs

PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};
end

PROGRESSTITLE = 'Writing the matrix file...';
if ~isempty(find(strcmpi('ProgressTitle', varargin)))
    PROGRESSTITLE = varargin{find(strcmpi('ProgressTitle', varargin))+1};
end

str = '\t%.3f';
if isinteger(data)
    str = '\t%d';
end

if ~iscell(labels_rows)
    warning('Row labels are numerical. Transforming them into strings...');
    labels_rows = cellfun(@num2str, num2cell(labels_rows), 'UniformOutput', false);
end

if ~iscell(labels_cols)
    warning('Column labels are numerical. Transforming them into strings...');
    labels_cols = cellfun(@num2str, num2cell(labels_cols), 'UniformOutput', false);
end

nrows = size(data,1);
ncols = size(data,2);

if size(labels_rows,1) ~= nrows
    error('The number of row labels does not match the number of rows in the data matrix. For a %d x %d data matrix, row labels are expected to be a %d x m cell array (with m > 0).', nrows, ncols, nrows);
end

if size(labels_cols,1) ~= ncols
    error('The number of column labels does not match the number of columns in the data matrix. For a %d x %d data matrix, column labels are expected to be a %d x p cell array (with p > 0).', nrows, ncols, ncols);
end


%%

if PROGRESSBAR
    w = waitbar(0, PROGRESSTITLE, 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(w, 'canceling', 0);
end

fprintf([PROGRESSTITLE '\n']);
fprintf(['|', blanks(100), '|\n']);
fprintf('|');


fprintf(fid,'ORF');

for n = 1 : size(labels_cols,2)
    if n > 1
        fprintf(fid, '\n');
    end
    for m = 1 : ncols
        if m == 1
            fprintf(fid, repmat('\t', 1, size(labels_rows,2)-1));
        end
        fprintf(fid, '\t%s',labels_cols{m,n});
    end
end

y = 0;

for p = 1 : nrows
    
    % Check for Cancel button press
    if PROGRESSBAR && getappdata(w,'canceling')
        delete(w);
        error('Canceled.');
    end
    
    % -- Progress bar ---
    if PROGRESSBAR
        waitbar(p/nrows, w);
    end
    
    x = fix(p * 100 / nrows);
    if x > y
        fprintf('*');
        y = x;
    end 
    % -- Progress bar ---
        
    for q = 1 : size(labels_rows,2)
        if q == 1
            fprintf(fid, '\n%s', labels_rows{p,q});
        else
            fprintf(fid, '\t%s', labels_rows{p,q});
        end
    end
    
    for m = 1 : ncols
        fprintf(fid, str, data(p,m));
    end
end

if PROGRESSBAR
    delete(w);
end

fprintf('|\n');
