% MCC_CurateDatabaseDICOM Curate a DICOM file database.
%
%   MCC_CurateDatabaseDICOM(scrFolder, dstFolder)
%
%   INPUTS
%       scrFolder: source folder
%       dstFolder: destination folder
%
%   All DICOM image files in the database are copied into a consistent
%   folder structure:
%
%   Dest. folder --- Patient ID 1 --- Study date 1 --- Series description 1
%                 |                |                |
%                 |                |                 - Series description 2
%                 |                |                |
%                 |                |                ...
%                 |                |
%                 |                 - Study date 2 ...
%                 |                |
%                 |                ...
%                 |
%                  - Patient ID 2  ...
%                 |
%                 ...
%
%   Typical rate: 8.5 gigabyte/hour
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function MCC_CurateDatabaseDICOM(scrFolder, dstFolder)

disp('Curate ...');
% Search for all folders and files
[subFolderList,~] = MCC_subdir(scrFolder);
% CSV file
sd = struct( ...
    'Clock','', ...
    'PatientID','',...
    'StudyDate','',...
    'SeriesDescription','',...
    'SeriesNumber','', ...
    'SeriesInstanceUID','', ...
    'TriggerTime', '' , ...     % NGR 2016 06 27
    'StudyInstanceUID', ''  ...	% NGR 2016 08 29
    );
if ~isdir(dstFolder)
    mkdir(dstFolder);
end
fileName_csv = [dstFolder '\Curation.csv'];
if exist(fileName_csv, 'file')
    fid_csv = fopen(fileName_csv,'at+');
    fprintf(fid_csv,'%s\n','');
else
    fid_csv = fopen(fileName_csv,'wt+');    
    fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s,%s\n', ...
        'Clock','Patient ID', 'Study Date', ...
        'Series Description', 'Series Number', ...
        'Trigger Time','EchoTime', 'StudyInstanceUID', ...
        'Source Folder');
end
% Log File
fileName_log = [dstFolder '\log.txt'];
% No subfolder case
if isempty(subFolderList)
   subFolderList{1}=scrFolder;
end
% Folder loop
for n=1:length(subFolderList)
    s = dir(subFolderList{n});
    iSD = 1;
    % File loop
    for i=3:length(s)
        if ~s(i).isdir
            CD =cd;
            scrFile = [subFolderList{n} '\' s(i).name];
            try
                md  = MCC_ReadMetadataDICOM(scrFile);
                bContinue = 1;
            catch
                bContinue = 0;
            end
            cd(CD);
            if bContinue
                try
                    N=length(md.SeriesDescription);
                catch
                    N=0; % Suspicious fragmentary file, might not be DICOM.
                    md.SeriesDescription = 'Unknown';N=1;
                end
                if N>0
                    % Imaging data only
                    try
%                         md.PixelSpacing;
%                         md.PatientID; md.StudyDate; md.SeriesDescription;
%                         md.SeriesInstanceUID; md.SeriesNumber;
                        bImage = 1;
                    catch ME
                        bImage = 0;
                        MCC_writeLog(fileName_log, ME, scrFile);
                    end
                    if bImage
                        % Update
                        sd(iSD).Clock = MCC_clock2str;
                        sd(iSD).PatientID = md.PatientID;
                        sd(iSD).StudyDate = md.StudyDate;
                        sd(iSD).SeriesDescription = md.SeriesDescription;
                        sd(iSD).SeriesNumber = md.SeriesNumber;
                        sd(iSD).SeriesInstanceUID = md.SeriesInstanceUID;                                                                        
                        try
                            sd(iSD).TriggerTime = md.TriggerTime; %NGR 2016 06 27                            
                        catch
                            md.TriggerTime = 0;                            
                            sd(iSD).TriggerTime = md.TriggerTime; %NGR 2016 06 27                            
                        end
                        try                            
                            sd(iSD).EchoTime = md.EchoTime;
                        catch                           
                            md.EchoTime = 0;                           
                            sd(iSD).EchoTime = md.EchoTime;
                        end
                        try
                            sd(iSD).AcquisitionTime = md.AcquisitionTime; %NGR 2016 08 19
                        catch
                            md.AcquisitionTime = 0;
                            sd(iSD).AcquisitionTime = md.AcquisitionTime;
                        end
                        try
                            sd(iSD).StudyInstanceUID = md.StudyInstanceUID;
                        catch                        
                            sd(iSD).StudyInstanceUID = ''; 
                            md.StudyInstanceUID = '';                            
                        end                        
                        
                        iSD = iSD + 1;
                        % Characters to avoid in folder name
                        ca = {...
                            '''','<','>',':','"','?', ...
                            '/','\','|','?','*'};                        
                        str = { ...                            
                            md.PatientID, ...
                            md.StudyDate, ...
                            [md.SeriesDescription  ...
                            ' [ S#' num2str(md.SeriesNumber) ...
                            '-' num2str(md.TriggerTime) ...
                            '-' num2str(md.EchoTime) ...
                            '-' num2str(md.AcquisitionTime) ... %NGR 2016 08 19
                            '-' num2str(md.StudyInstanceUID) ... %NGR 2016 08 19
                            ' ]']};
                        
                        % NGR 2016 09 14 , without acquisition time
                        str = { ...                            
                            md.PatientID, ...
                            md.StudyDate, ...
                            [md.SeriesDescription  ...
                            ' [ S#' num2str(md.SeriesNumber) ...
                            '-' num2str(md.TriggerTime) ...
                            '-' num2str(md.EchoTime) ...                            
                            '-' num2str(md.StudyInstanceUID) ... %NGR 2016 08 19
                            ' ]']};   
                        
                        for iCA = 1:length(ca)
                            for iSTR = 1:length(str)
                                str{iSTR} = ...
                                    strrep(str{iSTR}, ca(iCA), ' ');
                            end
                        end
                        fld = strjoin( ...
                            [dstFolder str{1} str{2} str{3}],'\');
                        % Create destination folder
                        if ~isdir(fld)
                            [SUCCESS,MESSAGE,~] = mkdir(fld);
                            if SUCCESS
                                % CSV File
                                fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s,%s%s\n', ...
                                    strrep(sd(end).Clock, ',', ' '), ...
                                    strrep(sd(end).PatientID, ',', ' '), ...
                                    strrep(sd(end).StudyDate, ',', ' '), ...
                                    strrep(sd(end).SeriesDescription, ',', ' '), ...
                                    num2str(sd(end).SeriesNumber), ...
                                    num2str(sd(end).TriggerTime), ...
                                    num2str(sd(end).EchoTime), ...
                                    num2str(sd(end).StudyInstanceUID), ...
                                    subFolderList{n} ...
                                    );                                
                            else
                                disp(MESSAGE);
                            end                            
                        end
                        % Copy file
                        dstFile = [fld '\' s(i).name];
                        disp([...
                            '               ' dstFile]);                        
                        if isempty(dir(dstFile))
                            try                                
                                copyfile(scrFile, dstFile);
                                bSuccess = 1;
                            catch ME
                                MCC_writeLog(fileName_log, ME, scrFile);
                                bSuccess = 0;
                            end
                            % Anonimize 
                            bSuccess = 0;
                            if bSuccess                                
                                try                                   
                                    dicomanon(scrFile, dstFile, 'keep',...
                                        {...
                                        'PatientID',...
                                        'PatientAge',...
                                        'PatientSex',...
                                        'PatientWeight',...                                        
                                        'StudyDescription',...
                                        'SeriesDescription'});
                                catch ME
                                    MCC_writeLog(fileName_log, ME, scrFile);                                   
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Curate duplicate slides

% Write CSV file
fclose(fid_csv);