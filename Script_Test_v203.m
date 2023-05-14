% Verify consistency between visit dates
% All visits must contain the number of sequences and types
% eg. : ADC, FLAIR, T1post, T1pre, T2

%-------------%
% PID: 686483 %
%-------------%

% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('Z:\Data_MRI\Institution\Study\Test_v203');
% Register
obj.Register('T2','elastic');

% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('Z:\Data_MRI\Institution\Study\Test_v203');
obj.Register_Date('T2','20110321');
