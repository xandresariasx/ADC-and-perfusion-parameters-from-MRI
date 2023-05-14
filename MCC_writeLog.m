%MCC_writeLog Update log file with a message
%
%   function MCC_writeLog(fileName, ME, , version, info)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function MCC_writeLog(fileName, ME, info, bVerbose)

% Check input arguments
if (nargin==3)
    bVerbose = 1;
end
% Open
fid_log = fopen(fileName,'at+');

% Message
if bVerbose,
    disp('___________________________________________________________________________');
end
me = ...
    sprintf(...
    '\tClock: %s\n \t Ver.:\t%s\n \t Info:  %s\n \t   ID:\t%s\n \tMess.:\t%s\n \tStack:\n',...
    MCC_clock2str, ...
    MCC_GetVersionFormated, ...
    info, ...
    ME.identifier,...
    ME.message ...
    );
if bVerbose
    disp(me);
end
fprintf(fid_log, '%s', me);
% Stack
for i=1:length(ME.stack)
    st = ...
        sprintf(...
        '\t\t\t%s (line %s)\n', ...
        ME.stack(i).name, ...
        num2str(ME.stack(i).line) ...
        );
    if bVerbose,
        disp(st(1:end-1));
    end
    fprintf(fid_log, '%s', st);
end
if bVerbose
    disp('___________________________________________________________________________');
end
fprintf(fid_log, '___________________________________________________________________________\n');
fclose(fid_log);