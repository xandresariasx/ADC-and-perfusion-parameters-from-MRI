%MCC_DisplayContent  Get content index from a query (pattern).
%
%   r = MCC_DisplayContent(content, pattern)
%
%   INPUTS
%       content: text cell array describing class content
%
%   OUTPUT
%       
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function MCC_DisplayContent(content)
% Display content
% NGR 2016 03 10
disp(' ');
disp('		Content ...');
for i=1:length(content)        
    disp(['				' num2str(i) ': ' content{i}]);    
end
disp(' ');
