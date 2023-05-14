% MCC_cTemplate  Class for template
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

classdef MCC_cTemplate < handle
    properties
        Cohort     % 4-D data (sequence of volumes) (x,y,z,k)
        Value
    end
    methods
        % Constructor
        function obj = MCC_cTemplate(val)
            if nargin > 0
                if ischar(val)
                    obj.Value = MCC_dicomreadvol(val);
                else
                    error('Value must be a string (DICOM path name)');
                end
            end
        end
                
        % Methods
        function r = roundOff(obj)
            r = round([obj.Value],2);
        end
        function r = multiplyBy(obj,n)
            r = [obj.Value] * n;
        end
        % Overlaod
        function r = plus(o1,o2)
            r = [o1.Value] + [o2.Value];
        end
    end
end