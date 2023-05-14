% MCC_CloneObj  Clone an object.
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function new_obj = MCC_CloneObj(obj)

% Instantiate new object of the same class.
new_obj = feval(class(obj));

% Copy all non-private properties.
p = properties(obj);
for i = 1:length(p)
    new_obj.(p{i}) = obj.(p{i});
end
