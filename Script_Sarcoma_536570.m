% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('C:\Temp\Institution\Study\536570\');
% Register
obj.Register_Date('AX_T1W','04_21_08');
%obj.Register('AX_T1W','rigid');

return;

% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('C:\Temp\Institution\Study\536570');
% Register
obj.Register('AX_T1W','elastic');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
close all; clear all; warning off;
% Create object
obj = MCC_cPatient('C:\Temp\Institution\Study\536570_2');
% Register
obj.Register('AX_T1W','elastic');





