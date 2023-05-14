% Verify consistency between visit dates
% All visits must contain the number of sequences and types
% eg. : ADC, FLAIR, T1post, T1pre, T2

%-------------%
% PID: 686483 %
%-------------%

% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('U:\Data_MRI\Institution\Study\Test_v201_test');
% Register
obj.Register('T2','rigid');

% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('U:\Data_MRI\Institution\Study\Test_v201_test');
% Register
obj.Register('T2','elastic');

return;

% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('E:\Data\Moffitt\Visit_test\Brain_MRI\686483_elastic');
obj.Register_Date('T2','2015.12.17');

return;

