% MCC_Histogram  Normalized histogram of a volume.
%
%   [ H1, N1 ] = MCC_Histogram(V1 , Nbins)
%
%       V1 : input volume
%       H1 : histogram
%       N1 : total number of bins
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function [ H1, X1 ] = MCC_Histogram(V1 , Nbins)

% Normalization & Quantization
V1 = single(V1);
Vmax = max(V1(:));
V1 =  uint8(round((Nbins-1)*(V1/Vmax)));
X1 = linspace(0,Vmax, Nbins);
% Marginal histograms
H1 = zeros(1,Nbins);
% Compute 
N1 = length(V1(:));
for i1=1:N1
    n1 = V1(i1)+1;
    H1(n1) =  H1(n1)+1;
end
% Probability Density Function
H1=H1/N1;
