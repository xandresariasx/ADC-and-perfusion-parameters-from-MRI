% MCC_MutualInfomation  Normalized mutual information of two volumes.
%
%   MI = MCC_MutualInfomation(V2, V1, Nbins)
%
%       V1      : fixed volume
%       V2      : registered volume
%       Nbins   : number of bins
%       MI      : Mutual Information
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center


function [MI, H12] = MCC_MutualInfomation(V1, V2, Nbins)

% Normalization & Quantization
V1  = single(V1);
V2  = single(V2);
V1 =  uint8(round((Nbins-1)*(V1/max(V1(:)))));
V2 =  uint8(round((Nbins-1)*(V2/max(V2(:)))));
% Marginal histograms
H1 = zeros(1,Nbins);
H2 = H1;
% Joint histogram
H12 = zeros(Nbins);
% Compute 
N1 = length(V1(:));
N2 = length(V2(:));
for i1=1:N1
    n1 = V1(i1)+1;
    H1(n1) =  H1(n1)+1;
    i2 = i1;
    n2 =  V2(i2)+1;
    H2(n2) =  H2(n2)+1; % Negative MI bug fixed, NGR 2016 03 15
    H12(n1,n2) = H12(n1,n2)+1;
end
% Probability Density Function
H1=H1/N1; H2= H2/N2; H12=H12/(N1);
% Mutual Information
MI = 0;
for n1=1:Nbins
    for n2=1:Nbins
        if (H1(n1)>0) && (H2(n2)>0) && (H12(n1,n2)>0)       
            MI = MI + ...
                H12(n1,n2) * log2( ...
                H12(n1,n2)/(H1(n1)*H2(n2)) ...
                );
        end
    end
end
