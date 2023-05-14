%MCC_GetContentIndex  Get content index from a query (pattern).
%
%   r = GetContentIndex(content, pattern)
%
%   INPUTS
%       content: text cell array
%       pattern: pattern string to find
%
%   OUTPUT
%       r : array of indices (found occurrences for pattern)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function r = MCC_GetContentIndex(content, pattern)
% Get content index
% NGR 2016 03 10
pat = lower(pattern);
r = [];
j = 1;
for i=1:length(content)
    s = lower(content{i});
    s = MCC_RemoveSeriesNumber(s);
    k = strfind(s, pat);
    if ~isempty(k)
        r(j) = i;
        return;
    end
end

if isempty(r)
    disp(' ');
    disp('No pattern found !');    
    MCC_DisplayContent(content);    
end