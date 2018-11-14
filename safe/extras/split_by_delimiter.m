function tokens = split_by_delimiter(delimiter, string)

if ~iscell(string) string = {string}; end

if strcmp(delimiter,'\t') delimiter = char(9); end

for i = 1 : length(string)

    remain = string{i};
    tokens{i,1} = {};
    
    j = 1;
    while true
        [str, remain] = strtok(remain, delimiter);
        if isempty(str) break; end
        tokens{i,j} = str;
        j = j + 1;
    end

end

