function D = read_matrix_file(input_file, varargin)

%% Usage
%
% D = read_matrix_file(input_file, num_row_headers, num_col_headers)
% TREEVIEW FILE: num_row_headers = 2, num_col_headers = 2
% 

%% Check inputs

numrowheaders = 1;
numcolheaders = 1;

if nargin > 1 && isnumeric(varargin{1}) && isnumeric(varargin{2})
    numrowheaders = varargin{1};
    numcolheaders = varargin{2};
end

% Show GUI-style progress bar
PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};   
end

% Empty values
empty = '(none|empty|blank)';

%%

D.numrowheaders = numrowheaders;
D.numcolheaders = numcolheaders;

% D.data = dlmread(input_file,'\t', numrowheaders, numcolheaders);

fid = fopen(input_file,'r');
if fid == -1
    error(sprintf('Unable to open file %s\n', input_file));
    return;
end

% Get the number of lines, if possible
lineCount = -1;
if ~ispc
    [status, cmdout] = system(['wc -l ' input_file]);
    if status ~= 1
        scanCell = textscan(cmdout,'%u %s');
        lineCount = scanCell{1}; 
    end
end
lineCount = lineCount - numcolheaders;

% Get column headers
D.labels_col = {};
for i = 1 : numcolheaders
    
    line = '';
    
    % Skip lines if they start with '#'
    while isempty(line) || ~isempty(regexp(line, '^#','once'))
        line = fgetl(fid);
    end
    
    tmp = textscan(line,'%s','delimiter','\t');
    colheaders = tmp{1}(numrowheaders+1:end);
    if length(colheaders) > size(D.labels_col,1)
        tmp = D.labels_col;
        D.labels_col = cell(length(colheaders),size(D.labels_col,2));
        D.labels_col(1:size(tmp,1),1:size(tmp,2)) = tmp;
        D.labels_col(:,i) = colheaders;
    else
        D.labels_col(1:length(colheaders),i) = colheaders;
    end
end

% Get row headers and data
fmt = [repmat('%s ', 1, numrowheaders) repmat('%f ', 1, length(D.labels_col))];

D.labels_row = cell(1,numrowheaders);
D.data = [];

if lineCount > 0
    if PROGRESSBAR
        w = waitbar(0, 'Reading the file...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(w, 'canceling', 0);
    else
        fprintf('\nReading the file...\n');
        fprintf(['|', blanks(100), '|\n']);
        fprintf('|');
    end
end

r = 1; y = 0;
while 1
    line = fgetl(fid);
    if line == -1 
        break; 
    end
    
    if lineCount > 0
        
        % Check for Cancel button press
        if PROGRESSBAR && getappdata(w,'canceling')
            delete(w);
            error('Canceled.');
        end
        
        if PROGRESSBAR
            waitbar(r / lineCount, w);
        else
            x = round(r * 100 / lineCount);
            if x > y
                fprintf(repmat('*',1,x-y)); y = x;
            end
        end
        
    else
        
        if mod(r,100) == 0
            fprintf('Read %d lines\n', r);
        end
        
    end
    
    % Replace empty values with NaN
    line = regexprep(line, empty, 'NaN', 'ignorecase');
    
    t = textscan(line,fmt,'delimiter','\t');
    t(find(cellfun(@isempty,t))) = {NaN};
    
    D.labels_row(r,:) = cat(2,t{1:numrowheaders});
    D.data(r,:) = cell2mat(t(numrowheaders+1:end));
    
    r = r + 1;
    
end

if lineCount > 0
    if PROGRESSBAR
        delete(w);
    else
        fprintf('|\n');
    end
end