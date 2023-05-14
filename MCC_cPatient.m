% MCC_cPatient  Class for patient (collection of study dates)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2018 © Moffitt Cancer Center

classdef MCC_cPatient < handle
    
    properties
        % Who
        desc='Radiomics'% Description (eg. mpMRI)
        name            % Name (Patient ID)
        content         % List of study dates
        % Variables
        StudyDate       % Study date sequence (5D) - array of 4D objects
        Machine         %
        % Location
        root            % Root directory
        InstituteID     % Institute ID
        DatabaseID      % Database ID
        CohortID        % Cohort ID
        PatientID       % Patient ID
        Address         % Address (Radiomics language)
    end
    
    properties (Access = private)
        bDisplay = 1;   % Display execution line and figures
    end
    
    methods
        
        %-----------------------------------------------------------------%
        % Constructor
        %-----------------------------------------------------------------%
        
        function obj = MCC_cPatient(root)
            % Code
            if nargin > 0
                if ischar(root)
                    % Path
                    obj.root = root;
                    % Location
                    cellTxt = strsplit(root,'\');
                    obj.name = cellTxt{end};
                    obj.InstituteID =  cellTxt{end-3};
                    obj.DatabaseID =  cellTxt{end-2};
                    obj.CohortID = cellTxt{end-1};
                    obj.PatientID = cellTxt{end};
                    obj.Address = lower([ ...
                        obj.InstituteID '.' ...
                        obj.DatabaseID '.' ...
                        obj.CohortID '.' ...
                        obj.PatientID  ...
                        ]);
                    % Load study dates
                    obj.Load;
                else
                    % Update
                    error('Value must be a string (DICOM name)');
                end
            end
        end
        
        %-----------------------------------------------------------------%
        % Processing
        %-----------------------------------------------------------------%
        
        function Process(obj, VolRef)
            % Process data
            % NGR 2016 08 23
            
            if ~exist('var','VolRef')
                VolRef = 1;
            end
            for i=1:obj.GetNbObjects
                obj.StudyDate(i).Process(VolRef);
            end
        end
        
        function Resample(obj, iVolRef)
            % Resample all study dates
            % NGR 2016 03 04
            
            for i=1:obj.GetNbObjects
                if nargin==2 && isnumeric(iVolRef)
                    obj.StudyDate(i).Resample(iVolRef);
                else
                    obj.StudyDate(i).Resample;
                end
            end
        end
        
        function Register(obj, VolRef, options, options2, options3)
            % Register all study dates
            % NGR 2016 03 08, NGR 2016 09 12
            
            if ~exist('options2')       % 2018 06 27 
                options2 = [];
            end 
            
            
            for i=1:obj.GetNbObjects
                % Display
%                 obj.StudyDate(i).Display(...          % 1/2/2019
%                     MCC_eView.Axial, MCC_eType.Original);
%                 obj.StudyDate(i).Display(...
%                     MCC_eView.Coronal, MCC_eType.Original);
%                 obj.StudyDate(i).Display(...
%                     MCC_eView.Sagittal, MCC_eType.Original);
                % Register
                %obj.StudyDate(i).Register(VolRef, options);
                %try
                if nargin==3
                    obj.StudyDate(i).Register(VolRef{i}, options, options2); % 2018 06 07 
                else
                    obj.StudyDate(i).Register(VolRef{i}, options, options2, options3);
                end
%                 catch
%                     nop=1;
%                 end
                % Display
%                 obj.StudyDate(i).Display(...              % 1/2/2019
%                     MCC_eView.Axial, MCC_eType.Registered);
                % obj.StudyDate(i).Display(MCC_eView.RGB);
                % Clear memory
                obj.StudyDate(i).ClearVolSeqData;
                close all;
            end
        end
        
        function Register_orgi(obj, VolRef, options)
            % Register all study dates
            % NGR 2016 03 08, NGR 2016 09 12
            
            for i=1:obj.GetNbObjects
                % Display
                obj.StudyDate(i).Display(...
                    MCC_eView.Axial, MCC_eType.Original);
                obj.StudyDate(i).Display(...
                    MCC_eView.Coronal, MCC_eType.Original);
                obj.StudyDate(i).Display(...
                    MCC_eView.Sagittal, MCC_eType.Original);
                % Register
                obj.StudyDate(i).Register(VolRef, options);                 
                % Display
                obj.StudyDate(i).Display(...
                    MCC_eView.Axial, MCC_eType.Registered);
                % obj.StudyDate(i).Display(MCC_eView.RGB);
                % Clear memory
                obj.StudyDate(i).ClearVolSeqData;
                close all;
            end
        end

        
        function Flipup(obj, eType)
            % Flip volumes in up / down direction
            % NGR 2018 041 19 v205
            
            for i=1:obj.GetNbObjects            
                obj.StudyDate(i).Flipud(eType);
                obj.StudyDate(i).Write(eType);
            end
        end
        
       function Register_Date_orgi(obj, VolRef, options)
            % Flip volumes in up / down direction
            % NGR 2018 041 19 v205
            
            %obj.Display_Registered_Date(options);return;
            
            % Load data
            CD = cd;
            Nd = obj.GetNbObjects;
            eT_Reg = MCC_eType.Registered;
            % Identify reference visit
            iDref = obj.Identify_Baseline_Date(options);
            % Patient registration loop
            for iD=1:Nd
                if obj.StudyDate(iD).IsVolSeqEmpty(eT_Reg)
                    obj.StudyDate(iD).Load;
                    obj.StudyDate(iD).Load(eT_Reg);
                    if obj.StudyDate(iD).IsVolSeqEmpty(eT_Reg)
                        % Register if no registered data
                        obj.StudyDate(iD).Register(VolRef);
                        close all;
                    end
                end
                iVolRef(iD) = obj.StudyDate(iD).iVolRef_Registered; % NGR 2018 04 17 v206
                subFolderList{iD} = ...
                    obj.StudyDate(iD).vsRegistered(iVolRef(iD)).root;
            end
            % Registration
            scrFolder=[obj.root '\' obj.StudyDate(iDref).StudyDateID];
            x = MCC_c4D(scrFolder, subFolderList);
            for iD=1:Nd                
                if iD == 1
                    x.vsOriginal = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                else
                    x.vsOriginal(iD) = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                end
            end
            x.Register_orgi(iDref, 'date');
            % Manage results
            switch iDref
                case 1
                    iaSel = [1:Nd];
                case Nd
                    iaSel = [Nd 1:(Nd-1)];
                otherwise
                    iaSel = [iDref 1:(iDref-1) (iDref+1):Nd];
            end
            res_mm_o = x.vsRegistered(iDref).res_mm; % NGR 2017 09 30 v310
            for iD = iaSel
                % Resample
                obj.StudyDate(iD).vsRegistered_Date = ...
                    obj.StudyDate(iD).CloneVolSeq(....
                    obj.StudyDate(iD).vsRegistered);
                obj.StudyDate(iD).Resample_Type(MCC_eType.Registered_Date, res_mm_o)
                obj.StudyDate(iD).eTypeCur = MCC_eType.Registered_Date;
                % Transform
                R = x.Volref3d(x.vsRegistered(iDref));
                obj.StudyDate(iD).ApplyTransform(...
                    MCC_eType.Registered_Date, x.tform{iD}, R);
                obj.StudyDate(iD).iVolRef_Registered_Date =  ...
                    obj.StudyDate(iD).iVolRef_Registered;
                obj.StudyDate(iD).eTypeCur = MCC_eType.Registered_Date;
                % ---------------------------------
                % Display comparison for all visits
                % ---------------------------------
                y = MCC_CloneObj(obj.StudyDate(iD));
                % Data
                y.vsOriginal(iVolRef(iD)).data  =  ...
                    obj.StudyDate(iDref).vsOriginal(iVolRef(iD)).data;
                y.vsOriginal = ...
                    obj.StudyDate(iD).CloneVolSeq(....
                    obj.StudyDate(iD).vsRegistered);
                y.vsResampled = y.Resample(iVolRef(iD));
                % Compare
                y.eTypeCur = MCC_eType.Registered_Date;
                y.Compare;
                close all;
                % ---------------------------
                % Display baseline comparison
                % ---------------------------
                y = MCC_CloneObj(obj.StudyDate(iD));
                y.vsOriginal = [];
                y.vsResampled = [];
                y.vsRegistered = [];
                y.vsRegistered_Date = [];
                y.content = [];
                % Data
                y.vsOriginal  =  MCC_CloneObj( ...
                    obj.StudyDate(iDref).vsOriginal(iVolRef(iD)));
                y.vsOriginal(2)  =  MCC_CloneObj( ...
                    obj.StudyDate(iD).vsOriginal(iVolRef(iD)));
                y.vsResampled = y.Resample(1);
                y.vsRegistered  =  MCC_CloneObj( ...
                    obj.StudyDate(iDref).vsRegistered(iVolRef(iD)));
                y.vsRegistered(2)  =  MCC_CloneObj( ...
                    obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                y.vsRegistered_Date  =  MCC_CloneObj( ...
                    obj.StudyDate(iDref).vsRegistered_Date(iVolRef(iD)));
                y.vsRegistered_Date(2)  =  MCC_CloneObj( ...
                    obj.StudyDate(iD).vsRegistered_Date(iVolRef(iD)));
                % Name
                y.vsOriginal(1).ImagingModeID = [ ...
                    obj.StudyDate(iD).vsOriginal(iVolRef(iD)).ImagingModeID ...
                    ' baseline'];
                % Compare
                y.eTypeCur = MCC_eType.Registered_Date;
                y.iVolRef_Registered = 1; y.iVolRef_Registered_Date = 1; % NGR 2018 04 11 v203
                try
                    y.Compare;
                catch ME
                    disp(ME);
                    MCC_writeLog([CD '\log.txt'], ...
                        ME, '');
                    cd(CD);
                end
                close all;
                % Display                                           12/3 no slow down pc with figs
%                 obj.StudyDate(iD).Display(...
%                     MCC_eView.All, MCC_eType.Original);
%                 obj.StudyDate(iD).Display(...
%                     MCC_eView.All, MCC_eType.Registered);
                try
                    obj.StudyDate(iD).iVolRef_Original = obj.StudyDate(iD).iVolRef_Registered;
                    obj.StudyDate(iD).iVolRef_Resampled = obj.StudyDate(iD).iVolRef_Registered;
                    obj.StudyDate(iD).Resample;
                    obj.StudyDate(iD).eTypeCur = MCC_eType.Registered;
                    obj.StudyDate(iD).Compare;
                catch ME
                    disp(ME);
                    MCC_writeLog([CD '\log.txt'], ...
                        ME, '');
                    cd(CD);
                end
                close all;
                obj.StudyDate(iD).Display(...
                    MCC_eView.All, MCC_eType.Registered_Date);
                close all;
                % DICOM output
                obj.StudyDate(iD).eTypeCur = MCC_eType.Registered_Date;
                obj.StudyDate(iD).Write;
            end
        end

        
        function Register_Date(obj, VolRef, options, options2)
            % Flip volumes in up / down direction
            % NGR 2018 041 19 v205
            
            %obj.Display_Registered_Date(options);return;
            
            % Load data
            CD = cd;
            Nd = obj.GetNbObjects;
            eT_Reg = MCC_eType.Registered;
            % Identify reference visit
            iDref = obj.Identify_Baseline_Date(options);
            % Patient registration loop
            for iD=1:Nd
                if obj.StudyDate(iD).IsVolSeqEmpty(eT_Reg)
                    obj.StudyDate(iD).Load;
                    obj.StudyDate(iD).Load(eT_Reg);
                    if obj.StudyDate(iD).IsVolSeqEmpty(eT_Reg)
                        % Register if no registered data
                        %obj.StudyDate(iD).Register(VolRef);
                        obj.StudyDate(iD).Register(VolRef(iD),'elastic');
                        close all;
                    end
                end
                %%% 7/19/2018
                Ref=cell(size(obj.StudyDate(iD).content));
                Ref(:)=VolRef(iD);
                aux=cellfun(@isequal, obj.StudyDate(iD).content,Ref);
                obj.StudyDate(iD).iVolRef_Registered=find(aux);
                %%%%%%
                try
                    iVolRef(iD) = obj.StudyDate(iD).iVolRef_Registered; % NGR 2018 04 17 v206
                catch
                    nop=1;
                end
                subFolderList{iD} = ...
                    obj.StudyDate(iD).vsRegistered(iVolRef(iD)).root;
            end
            % Registration
            scrFolder=[obj.root '\' obj.StudyDate(iDref).StudyDateID];
            x = MCC_c4D(scrFolder, subFolderList);
            for iD=1:Nd                
                if iD == 1
                    x.vsOriginal = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                else
                    x.vsOriginal(iD) = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                end
            end
            x.Register(iDref, 'date',options2);
            % Manage results
            switch iDref
                case 1
                    iaSel = [1:Nd];
                case Nd
                    iaSel = [Nd 1:(Nd-1)];
                otherwise
                    iaSel = [iDref 1:(iDref-1) (iDref+1):Nd];
            end
            res_mm_o = x.vsRegistered(iDref).res_mm; % NGR 2017 09 30 v310
            for iD = iaSel
                for J=1:length(obj.StudyDate(iD).vsRegistered)                      %10_29_2018 update metadata in obj
                    obj.StudyDate(iD).vsRegistered(J).metadata=x.vsRegistered(iD).metadata;
                end
                % Resample
                obj.StudyDate(iD).vsRegistered_Date = ...
                    obj.StudyDate(iD).CloneVolSeq(....
                    obj.StudyDate(iD).vsRegistered);
                obj.StudyDate(iD).Resample_Type(MCC_eType.Registered_Date, res_mm_o)
                obj.StudyDate(iD).eTypeCur = MCC_eType.Registered_Date;
                % Transform
                R = x.Volref3d(x.vsRegistered(iDref));
                try
                    obj.StudyDate(iD).ApplyTransform(...
                        MCC_eType.Registered_Date, x.tform{iD}, R);
                catch
                    obj.StudyDate(iD).ApplyTransform(...
                        MCC_eType.Registered_Date, x.tform{iD}{1}, R);
                end
                obj.StudyDate(iD).iVolRef_Registered_Date =  ...
                    obj.StudyDate(iD).iVolRef_Registered;
                obj.StudyDate(iD).eTypeCur = MCC_eType.Registered_Date;
                % ---------------------------------
                % Display comparison for all visits
                % ---------------------------------
                try
                    y = MCC_CloneObj(obj.StudyDate(iD));
                    % Data                
                    y.vsOriginal(iVolRef(iD)).data  =  ...
                        obj.StudyDate(iDref).vsOriginal(iVolRef(iD)).data;
                    y.vsOriginal = ...
                        obj.StudyDate(iD).CloneVolSeq(....
                        obj.StudyDate(iD).vsRegistered);
                    y.vsResampled = y.Resample(iVolRef(iD));
                    % Compare
                    y.eTypeCur = MCC_eType.Registered_Date;
                    y.Compare;
                    close all;
                    % ---------------------------
                    % Display baseline comparison
                    % ---------------------------
                    y = MCC_CloneObj(obj.StudyDate(iD));
                    y.vsOriginal = [];
                    y.vsResampled = [];
                    y.vsRegistered = [];
                    y.vsRegistered_Date = [];
                    y.content = [];
                    % Data
                    y.vsOriginal  =  MCC_CloneObj( ...
                        obj.StudyDate(iDref).vsOriginal(iVolRef(iD)));
                    y.vsOriginal(2)  =  MCC_CloneObj( ...
                        obj.StudyDate(iD).vsOriginal(iVolRef(iD)));
                    y.vsResampled = y.Resample(1);
                    y.vsRegistered  =  MCC_CloneObj( ...
                        obj.StudyDate(iDref).vsRegistered(iVolRef(iD)));
                    y.vsRegistered(2)  =  MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                    y.vsRegistered_Date  =  MCC_CloneObj( ...
                        obj.StudyDate(iDref).vsRegistered_Date(iVolRef(iD)));
                    y.vsRegistered_Date(2)  =  MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered_Date(iVolRef(iD)));
                    % Name
                    y.vsOriginal(1).ImagingModeID = [ ...
                        obj.StudyDate(iD).vsOriginal(iVolRef(iD)).ImagingModeID ...
                        ' baseline'];
                    % Compare
                    y.eTypeCur = MCC_eType.Registered_Date;
                    y.iVolRef_Registered = 1; y.iVolRef_Registered_Date = 1; % NGR 2018 04 11 v203
                    try
                        y.Compare;
                    catch ME
                        disp(ME);
                        MCC_writeLog([CD '\log.txt'], ...
                            ME, '');
                        cd(CD);
                    end
                    close all;
                end
                % Display
%                 obj.StudyDate(iD).Display(...
%                     MCC_eView.All, MCC_eType.Original);
%                 obj.StudyDate(iD).Display(...
%                     MCC_eView.All, MCC_eType.Registered);
                try
                    obj.StudyDate(iD).iVolRef_Original = obj.StudyDate(iD).iVolRef_Registered;
                    obj.StudyDate(iD).iVolRef_Resampled = obj.StudyDate(iD).iVolRef_Registered;
                    obj.StudyDate(iD).Resample;
                    obj.StudyDate(iD).eTypeCur = MCC_eType.Registered;
                    obj.StudyDate(iD).Compare;
                catch ME
                    disp(ME);
                    MCC_writeLog([CD '\log.txt'], ...
                        ME, '');
                    cd(CD);
                end
                close all;
                obj.StudyDate(iD).Display(...
                    MCC_eView.All, MCC_eType.Registered_Date);
                close all;
                % DICOM output
                obj.StudyDate(iD).eTypeCur = MCC_eType.Registered_Date;
                obj.StudyDate(iD).Write;
            end
        end
        
        function Register_Date_Local(obj, VolRef, options, options2,VOI)
            CD = cd;
            Nd = obj.GetNbObjects;
            eT_Reg = MCC_eType.Registered;
            % Identify reference visit
            iDref = obj.Identify_Baseline_Date(options);
            % Patient registration loop
            for iD=1:Nd
                if obj.StudyDate(iD).IsVolSeqEmpty(eT_Reg)
                    obj.StudyDate(iD).Load;
                    obj.StudyDate(iD).Load(eT_Reg);
                    if obj.StudyDate(iD).IsVolSeqEmpty(eT_Reg)
                        % Register if no registered data
                        %obj.StudyDate(iD).Register(VolRef);
                        obj.StudyDate(iD).Register(VolRef(iD),'elastic');
                        close all;
                    end
                end
                %%% 7/19/2018
                Ref=cell(size(obj.StudyDate(iD).content));
                Ref(:)=VolRef(iD);
                aux=cellfun(@isequal, obj.StudyDate(iD).content,Ref);
                obj.StudyDate(iD).iVolRef_Registered=find(aux);
                %%%%%%
                try
                    iVolRef(iD) = obj.StudyDate(iD).iVolRef_Registered; % NGR 2018 04 17 v206
                catch
                    nop=1;
                end
                subFolderList{iD} = ...
                    obj.StudyDate(iD).vsRegistered(iVolRef(iD)).root;
            end
            % Loca Registration VOI
            scrFolder=[obj.root '\' obj.StudyDate(iDref).StudyDateID];
            x = MCC_c4D(scrFolder, subFolderList);
            for iD=1:Nd                
                if iD == 1
                    x.vsOriginal = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                     x.vsRegistered = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                else
                    x.vsOriginal(iD) = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                    x.vsRegistered(iD) = MCC_CloneObj( ...
                        obj.StudyDate(iD).vsRegistered(iVolRef(iD)));
                end
            end
            x.LoadData;
            x.LoadWorkspace;
            for I=1:length(VOI)
                x.VOI=VOI;
                x.VOI_mm(I,:)=VOI(I).xyzDxDyDz_mm;                
%                 x.Register({VolRef{iDref}},[],[],['voi' num2str(I)]);
                 % Parameters
                eT = MCC_eType.Parametric;
                eV = MCC_eView.Axial;
                iROI = I;
                eval(['eT  = MCC_eType.Registered_VOI' num2str(iROI) ';']);
                x.tType = options2;    % 7/6/2020
                        %{'none','rigid','similarity','affine'};
                ext = x.GetVolSeqLabel(eT);  
                x.Resample(VolRef{iDref});
                vsI = x.GetVolSeq(MCC_eType.Registered);  
                iVolRef = x.iVolRef_Resampled;
                Vf = x.vsResampled(iVolRef);
                for iVol = 1:length(x.vsResampled);
                    Vm = vsI(iVol);
                    Vt = x.ForceFOV(Vf, Vm);
                    vsI(iVol).data = Vt.data;
                end
                 % Crop
                VOIx = x.VOI_mm(iROI,:);
                VOIx = x.InteresctFOVandVOI(vsI, VOIx);
                vsC = x.CropVolSeqVOI(vsI, VOIx);
                 % Register
                vsO = x.RegisterVolSeq(vsC, iVolRef);                
                 % Crop and apply transformation to other images
                Vaux=MCC_c3D;
                iVaux=[];
                for II=1:Nd     
                    for iVol = 1:obj.StudyDate(II).GetNbObjects 
                        if iVol~=obj.StudyDate(II).iVolRef_Registered
                            Vaux(end+1)=obj.StudyDate(II).vsRegistered(iVol);
                            iVaux(end+1)=iVol;
                        end        
                    end
                end
                Vaux(1)=[];
                VauxC=x.CropVolSeqVOI(Vaux, VOIx);                
                Ro = x.Volref3d(vsC(iVolRef));
                Vt=cell(0);
                III=1;
                for II=1:Nd
                    for iVol = 1:obj.StudyDate(II).GetNbObjects-1 
                        V = VauxC(III).data;
                        R = x.Volref3d(vsC(II));
                        Vt{III} = imwarp(V ,R, x.tform{II}, 'bicubic','OutputView',Ro);
                        III=III+1;
                    end    
                end
                % Fill VOI and save image in original obj
                III=1;
                for II=1:Nd
                    if iROI==1
                        eval(['obj.StudyDate(II).vsRegistered_VOI' num2str(iROI)...
                            '= obj.StudyDate(II).CloneVolSeq(obj.StudyDate(II).vsRegistered);']);
                    else
                        eval(['obj.StudyDate(II).vsRegistered_VOI' num2str(iROI)...
                            '= obj.StudyDate(II).CloneVolSeq(obj.StudyDate(II).vsRegistered_VOI'...
                            num2str(iROI-1) ');']);  % 2/22
                    end
                    for iVol = 1:obj.StudyDate(II).GetNbObjects 
                        if iVol==obj.StudyDate(II).iVolRef_Registered
                            Vi=vsO(II).data;
                            iVol2=iVol;
                        else
                            Vi = Vt{III}; 
                            iVol2=iVaux(III);
                            III=III+1;
                        end
                        eval(['obj.StudyDate(II).vsRegistered_VOI' num2str(iROI)...
                            '(iVol2).FillVOI(Vi, VOIx);']);                        
                    end
                end
                % Fill gaps to visually appeal 
                for II=1:Nd
                    for III=1:numel(obj.StudyDate(II).vsRegistered)                        
                        eval(['Vaux2=obj.StudyDate(II).vsRegistered_VOI' num2str(iROI) '(III).data;'])
                        Vaux2(Vaux2==0)=obj.StudyDate(II).vsRegistered(III).data(Vaux2==0);
                        eval(['obj.StudyDate(II).vsRegistered_VOI' num2str(iROI) '(III).data=Vaux2;'])
                    end                   
                end   
                % Write
                for II=1:Nd
                    obj.StudyDate(II).Write(eT); 
                end
                
                R=obj.StudyDate(1).vsRegistered_VOI1(1).GetCropVec(VOIx);
                Reg_ROI=zeros(size(obj.StudyDate(1).vsRegistered_VOI1(1).data));
                Reg_ROI(R{1},R{2},R{3}) = 1;
                save([obj.root '\roi' num2str(iROI) '.mat'],'Reg_ROI')                 
            end
        end
        
        function Display_Registered_Date( obj , options )
            % Display Registered Dates
            % NGR 2017 09 30 v310
            
            iDref = obj.Identify_Baseline_Date(options);
            Nd = obj.GetNbObjects;
            % Load
            for iD=1:Nd
                % Load study
                obj.StudyDate(iD).LoadData;
                % Orginal
                obj.StudyDate(iD).Display(...
                    MCC_eView.All, MCC_eType.Original);
                % Registered
                obj.StudyDate(iD).Display(...
                    MCC_eView.All, MCC_eType.Registered);
                obj.StudyDate(iD).eTypeCur = MCC_eType.Registered;
                obj.StudyDate(iD).Compare;
                % Registered Dates
                obj.StudyDate(iD).Display(...
                    MCC_eView.All, MCC_eType.Registered_Date);
                % Close figures to save memory
                close all;
            end
        end
        
        function iDref = Identify_Baseline_Date(obj, options)
            % Identify reference visit
            % NGR 2017 09 30 v310
            
            Nd = obj.GetNbObjects;
            iDref = 1; % Base line (first visit by default)
            for iD=1:Nd
                txt = strsplit(obj.StudyDate(iD).root,'\');
                if strcmp(lower(txt{end}),lower(options))
                    iDref = iD;
                end
            end
        end
        
        function Segment(obj, VolRef)
            % Segment all
            % NGR 2016 04 27
            
            for i=1:obj.GetNbObjects
                % Display
                obj.StudyDate(i).Display(...
                    MCC_eView.Axial, MCC_eType.Original);
                % Segment
                if nargin==2
                    obj.StudyDate(i).Segment(VolRef);
                else
                    obj.StudyDate(i).Segment;
                end
                % Display
                obj.StudyDate(i).Display4Q(...
                    MCC_eType.Segmented, [0 0 0], 0, 1);
                % Clear memory
                obj.StudyDate(i).ClearVolSeqData;
                close all;
            end
        end
        
        function Quantify(obj, options)
            % Quantify all
            % NGR 2016 09 12
            
            for i=1:obj.GetNbObjects
                try
                    obj.StudyDate(i).Quantify(options);
                    obj.StudyDate(i).ClearVolSeq;
                    close all;
                catch ME
                    disp(ME);
                    CD = cd;
                    MCC_writeLog([CD '\log.txt'], ...
                        ME, ' ');
                    cd(CD);
                end
            end
        end
        
        function r = GetVOI(obj)
            % Get VOI set
            % NGR 2016 06 24
            
            k = 1;
            for i=1:obj.GetNbObjects
                VOI = obj.StudyDate(i).GetVOI;
                if ~isempty(VOI)
                    for j = 1:length(VOI)
                        r(k) = VOI(j);
                        k = k + 1;
                    end
                end
            end
        end
        
        function Localize(obj, Organ)
            % Localize organ
            % NGR 2016 06 30
            
            switch nargin
                case 1
                    Organ = 'Kidney';
            end
            % Localize
            for i=1 : obj.GetNbObjects
                obj.StudyDate(i).Localize(Organ);
                % Clear memory
                obj.StudyDate(i).ClearVolSeqData;
                close all;
            end
        end
        
        function Classify(obj, options)
            % Classify all
            % NGR 2016 11 03
            
            if ~exist('options')
                options = [];
            end
            for i=1:obj.GetNbObjects
                obj.StudyDate(i).Classify(options);
                obj.StudyDate(i).ClearVolSeq;
                close all;
            end
        end
        
        %-----------------------------------------------------------------%
        % Study dates
        %-----------------------------------------------------------------%
        
        function Load(obj)
            % Load all study dates
            % NGR 2016 03 04
            
            CD = cd;
            subFolderList = MCC_SubFolderList(obj.root);
            j = 1;
            if ~isempty(subFolderList)
                for i=1:length(subFolderList)
                    %                     try
                    obj.LoadOne(subFolderList{i}, j);
                    obj.content{j} = subFolderList{i};
                    j = j + 1;
                    %                     catch ME
                    %                         disp(ME);
                    %                         MCC_writeLog([CD '\log.txt'], ...
                    %                             ME, subFolderList{i});
                    %                         cd(CD);
                    %                     end
                end
            end
        end
        
        function LoadOne(obj, StudyDateID, index)
            % Load a study date
            % NGR 2016 03 04
            
            scrFolder = [obj.root '\' StudyDateID];
            if isempty(obj.StudyDate)
                obj.StudyDate = ...
                    MCC_c4D(scrFolder);
                obj.content{1} = StudyDateID;
            else
                obj.StudyDate(index) = ...
                    MCC_c4D(scrFolder);
                obj.content{index} = StudyDateID;
            end
        end
        
        %-----------------------------------------------------------------%
        % Export
        %-----------------------------------------------------------------%
        
        function ExportDICOMtags(obj)
            % Export DICOM tags
            % NGR 2017 08 02 v285
            
            for i=1:obj.GetNbObjects
                obj.StudyDate(i).ExportDICOMtags;
            end
        end
        
        %-----------------------------------------------------------------%
        % Display
        %-----------------------------------------------------------------%
        
        function Compare(obj)
            % Display registration for visual comparison
            % NGR 2016 03 10
            
            for i=1:length(obj.StudyDate)
                obj.StudyDate(i).Compare;
            end
        end
        
        function Display(obj, eView, eType)
            % Display central sections
            % NGR 2016 03 10
            
            switch nargin
                case 1
                    for i=1:obj.GetNbObjects
                        obj.StudyDate(i).Display(...
                            MCC_eView.Axial, MCC_eType.Original);
                        obj.StudyDate(i).Display(...
                            MCC_eView.Coronal, MCC_eType.Original);
                        obj.StudyDate(i).Display(...
                            MCC_eView.Sagittal, MCC_eType.Original);
                        obj.StudyDate(i).Display(...
                            MCC_eView.Planar, MCC_eType.Original);
                    end
                case 2
                    for i=1:obj.GetNbObjects
                        obj.StudyDate(i).Display(eView);
                        pause
                    end
                case 3
                    for i=1:obj.GetNbObjects
                        obj.StudyDate(i).Display(eView, eType);
                        %pause
                    end
            end
            close all;
        end
        
        %-----------------------------------------------------------------%
        % Language
        %-----------------------------------------------------------------%
        
        function r = GetContentIndex(obj, pattern)
            % Get content index from a query
            % NGR 2016 03 10
            
            r = MCC_GetContentIndex(obj.content, pattern);
        end
        
        function DisplayContent(obj)
            % Display content
            % NGR 2016 03 10
            
            DisplayContent(obj.content);
        end
        
        function r = GetNbObjects(obj)
            % Get number of objects
            % NGR 2016 03 23
            
            r = length(obj.StudyDate);
        end
        
    end
    
    methods (Access = private)
        
    end
    
end