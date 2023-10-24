% MCC_c4D  Class for 4D processing
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2018 © Moffitt Cancer Center

classdef MCC_c4D < handle
    
    properties
        % Who
        desc='Radiomics'% Description (eg. mpMRI)
        name            % Name (study date)
        content;        % Imaging modes (sub-folder list)
        contentSeg;     % Segmentation results
        time            % Acquisition time (hhmmss.ffffff)
        % Variables
        vsOriginal      % Volume sequence (4D) - array of 3D objects
        vsResampled     % Volume sequence (4D) - resampled
        vsRegistered    % Volume sequence (4D) - registered
        vsRegistered_Elastic  % Volume sequence (4D) - elastic registered
        vsRegistered_VOI1     % Volume sequence (4D) - locally registered
        vsRegistered_VOI2     % Volume sequence (4D) - locally registered
        vsRegistered_VOI3     % Volume sequence (4D) - locally registered
        vsRegistered_VOI4     % Volume sequence (4D) - locally registered
        vsRegistered_VOI5     % Volume sequence (4D) - locally registered   8/9/2020
        vsRegistered_VOI6     % Volume sequence (4D) - locally registered   8/9/2020
        vsRegistered_VOI7     % Volume sequence (4D) - locally registered   8/9/2020
        vsRegistered_VOI8     % Volume sequence (4D) - locally registered   8/9/2020
        vsRegistered_Date    % 4D registered along study Date, NGR 2018 03 21
        vsReady         % Volume sequence (4D) - ready (study date regist.)
        vsSegmented     % Volume sequence (4D) - segmented
        vsLocalized     % Volume sequence (4D) - localized (organ level)
        vsParametric    % Volume sequence (4D) - parametric volumes
        vsTemporary     % Volume sequence (4D) - temporary
        v_RGB           % RGB volume (3D)
        iVolRef_Original          % Reference volume index for original
        iVolRef_Resampled         % Reference volume index for resampling
        iVolRef_Registered        % Reference volume index for registration
        iVolRef_Registered_Elastic
        iVolRef_Registered_VOI1   % Reference volume index for registration
        iVolRef_Registered_VOI2   % Reference volume index for registration
        iVolRef_Registered_VOI3   % Reference volume index for registration
        iVolRef_Registered_VOI4   % Reference volume index for registration
        iVolRef_Registered_VOI5   % Reference volume index for registration  8/9/2020
        iVolRef_Registered_VOI6   % Reference volume index for registration  8/9/2020
        iVolRef_Registered_VOI7   % Reference volume index for registration  8/9/2020
        iVolRef_Registered_VOI8   % Reference volume index for registration  8/9/2020
        iVolRef_Ready=1;          % Reference volume index for registration
        iVolRef_Segmented         % Reference volume index for segmentation
        iVolRef_Localized         % Reference volume index for registration (local)
        iVolRef_Registered_Date   % NGR 2018 03 21
        eTypeCur        % Current type (original, resampled, registered...)
        eViewCur        % Current view (axial, coronal, sagittal)
        tExe=0;         % Execution time for each operation
        rMI             % Mutual Information (registration)
        MIpc            % Improvement % of normalized MI (registration)
        tform           % Geometric transformation (array)
        tType = {'none','rigid','similarity','affine'}; % Reg. scheme
        net             % Machine neural network
        MaxAxis         % Maximum axis
        % VOI
        hROI            % ROI handlers
        hROIl           % ROI handlers (localization)
        VOI_mm          % Volume Of Interest [x y z dx dy dz]
        VOIsl_mm        % Localized VOI
        VOI             % VOI class (structure)
        pSeeds          % Seeds points for region growing
        % Display
        viewpoint       % Viewpoint (planar/volume rendering)
        vPlaneCur       % Vector of orthogonal plane
        vPlaneCur_mm    % Vector of orthogonal plane in [mm]
        vShiftPlaneCur = [0 0 0];  % Vector of orthogonal plane
        maskReg;        % Registration mask (processing area)
        % Organ
        eOrganCur = MCC_eOrgan.Kidney; %Any; % Organ of interest
        % Location
        root            % Root directory
        InstituteID     % Institute ID
        DatabaseID      % Database ID
        CohortID        % Cohort ID
        PatientID       % Patient ID
        StudyDateID     % Study date ID
        Address         % Adress (Radiomics language)
        % Version
        version = MCC_GetVersionFormated;
        signature=[MCC_GetArtistSignature ' - ' MCC_GetInstituteSignature];
    end
    
    properties (Access = private)
        res_fact        % Resize factors (for resampling)
        bContiguous=[]; % Contiguous slice flag
    end
    
    events
        StateChange
    end
    
    methods (Access = private)
        % Event handler
        % NGR 2016 05 17
        function evntCb(~,~,evnt,varargin)
            disp(['Number of inputs: ',num2str(nargin)])
            disp(evnt.EventName)
            disp(varargin{:})
        end
    end
    
    methods
        
        %-----------------------------------------------------------------%
        % Events
        %-----------------------------------------------------------------%
        
        function triggerEvnt(obj)
            % Event trigger
            % NGR 2016 05 17
            notify(obj,'StateChange')
        end
        
        %-----------------------------------------------------------------%
        % Constructor
        %-----------------------------------------------------------------%
        
        function obj = MCC_c4D(root, subfolder_list) % NGR 2018 03 21
            % Code
            if nargin > 0
                if ~isempty(root)
                    % Event
                    t = datestr(now);
                    addlistener(obj, 'StateChange', ...
                        @(src,evnt)obj.evntCb(src,evnt,t));
                    % Root
                    obj.root = root;
                    obj.content = MCC_SubFolderList(root);
                    if exist('subfolder_list','var')
                        % Registration along study date, NGR 2018 03 21
                        obj.content = subfolder_list;
                    else
                        obj.content = MCC_SubFolderList(root);
                    end
                    % Location
                    cellTxt = strsplit(root,'\');
                    obj.name = cellTxt{end};
                    if numel(cellTxt)>4   % 9/29/2020
                        obj.InstituteID =  cellTxt{end-4};
                        obj.DatabaseID =  cellTxt{end-3};
                        obj.CohortID = cellTxt{end-2};
                    else
                        obj.InstituteID =  '';
                        obj.DatabaseID =  '';
                        obj.CohortID = '';
                    end
                    obj.PatientID = cellTxt{end-1};
                    obj.StudyDateID = cellTxt{end};
                    obj.Address = lower([ ...
                        obj.InstituteID '.' ...
                        obj.DatabaseID '.' ...
                        obj.CohortID '.' ...
                        obj.PatientID '.' ...
                        obj.StudyDateID ...
                        ]);
                    % Organ (FMX), NGR 2016 09 09 v195                    
                    bFound = strfind(lower(obj.StudyDateID), 'tubes');
                    if ~isempty(bFound)
                        obj.eOrganCur = MCC_eOrgan.Tubes_FMX;
                    end                    
                    bFound = strfind(lower(obj.StudyDateID), 'brain');
                    if ~isempty(bFound)
                        obj.eOrganCur = MCC_eOrgan.Brain_FMX;
                    end
                    bFound = strfind(lower(obj.StudyDateID), 'body');
                    if ~isempty(bFound)
                        obj.eOrganCur = MCC_eOrgan.Body_FMX;
                    end
                    bFound = strfind(lower(obj.StudyDateID), 'thoraco');
                    if ~isempty(bFound)
                        obj.eOrganCur = MCC_eOrgan.Body_FMX;
                    end
                    bFound = strfind(lower(obj.StudyDateID), 'abdomen');
                    if ~isempty(bFound)
                        obj.eOrganCur = MCC_eOrgan.Body_FMX;
                    end
                else
                    obj.content = subfolder_list;
                end
                % Update
                obj.eTypeCur = MCC_eType.Original;
                obj.eViewCur = MCC_eType.All;
            end
        end
        
        function ExludeCertainContent(obj, toEx)
            % Exlude certain content from processing
            % NGR 2016 07 07
            k = 1;
            bEx = 0;
            for i=1:length(obj.content)
                txt = lower(obj.content{i});
                for j=1:length(toEx)
                    if ~isempty(strfind(txt,lower(toEx{j})))
                        bEx = 1;
                        break;
                    end
                end
                if ~bEx
                    x{k} = obj.content{i};
                    k = k + 1;
                end
                bEx = 0;
            end
            obj.content = x;
        end
        
        %-----------------------------------------------------------------%
        % Processing
        %-----------------------------------------------------------------%
                
        function Quantify(obj, options)            
            % Process FMX study
            % NGR 2016 08 29
            
            if strcmp(options,'fmx')                
                % FMX pilot study   
                obj.LoadData;                            
                obj.Quantify_FMX;                                
                obj.SaveWorkspace;
            else
                disp('''fmx'' only available');
            end
            
        end
        
        function vsO = Resample(obj, VolRef)
            % Resampling (uniform spatial resolution across modalities)
            % NGR 2016 05 04
            obj.ClearVolSeq(MCC_eType.Resampled);
            % Check data
            obj.LoadData(MCC_eType.Original);
            % Input arguments
            switch nargin
                case 1
                    if isempty(obj.iVolRef_Resampled)
                        iVolRef = obj.GetLowestResolutionIndex;
                    else
                        iVolRef = obj.iVolRef_Resampled;
                    end
                case 2
                    if ~isnumeric(VolRef)
                        iVolRef = obj.GetContentIndex(VolRef);
                        if isempty(iVolRef)
                            iVolRef = obj.GetLowestResolutionIndex;
                        end
                    else
                        iVolRef = VolRef;
                    end
            end
            obj.iVolRef_Resampled = iVolRef;
            % Input
            vsI = obj.CloneVolSeq(obj.vsOriginal);
            % Resize
            res_mm_o  = vsI(iVolRef).res_mm;
            vsO = obj.ResampleVolSeq(vsI, res_mm_o);
            % Update
            vsO(obj.iVolRef_Resampled).name = ...
                [vsO(iVolRef).name ' [Reference]'];
            obj.vsResampled = obj.CloneVolSeq(vsO);
            obj.eTypeCur = MCC_eType.Resampled;
        end
        
        function Resample_Type(obj, eType, res_mm_o)
            % Resampling (uniform spatial resolution across modalities)
            % NGR 2018 03 21
            
            vsI = obj.GetVolSeq(eType);
            vsO = obj.ResampleVolSeq(vsI, res_mm_o);
            txt = ['obj.vs' obj.GetVolSeqLabel(eType) ' = obj.CloneVolSeq(vsO);'];
            eval(txt);
        end

        function vsO = Register_orgi(obj, VolRef, options)
            % Registration (volume sequence)
            % NGR 2016 02 09, NGR 2016 06 22
            
            if ~exist('options')
                options = [];
            end     

            if ~isempty(options) % NGR 2018 03 20
                % Register whole volume (no VOI)
                obj.ClearVolSeq(MCC_eType.Registered);
                % Check data
                obj.LoadData(MCC_eType.Original);
                % Input arguments
                switch nargin
                    case 1
                        if isempty(obj.iVolRef_Registered)
                            iVolRef = obj.GetLowestResolutionIndex;
                        else
                            iVolRef = obj.iVolRef_Registered;
                        end
                    case 2
                        if ~isnumeric(VolRef)
                            iVolRef = obj.GetContentIndex(VolRef);
                            if isempty(iVolRef)
                                iVolRef = obj.GetLowestResolutionIndex;
                            end
                        else
                            iVolRef = VolRef;
                        end
                    case {3, 4}
                        if ~isnumeric(VolRef)
                            iVolRef = obj.GetContentIndex(VolRef);
                            if isempty(iVolRef)
                                iVolRef = obj.GetLowestResolutionIndex;
                            end
                        else
                            iVolRef = VolRef;
                        end
                        % NGR 2018 03 20
                        switch options
                            case 'none'                 % AA 1/26/2021
                                obj.tType = {'none'};
                                eT = MCC_eType.Registered; 
                            case 'rigid'
                                obj.tType = {'none','rigid'};
                                eT = MCC_eType.Registered;     % NGR 2018 04 09
                            case 'elastic'
                                obj.tType = {'none','rigid','similarity','affine'};
                                eT = MCC_eType.Registered;     % NGR 2018 04 09
                            case 'date'
                                obj.tType = {'none','rigid','similarity','affine'};
                                eT = MCC_eType.Registered_Date;     % NGR 2018 03 21
                        end
                                
                end
                obj.iVolRef_Registered = iVolRef;
                % Resample data (if already not)
                if obj.IsVolSeqEmpty(MCC_eType.Resampled)
                    obj.Resample(iVolRef);
                end
                obj.eViewCur = MCC_eView.All;
                % Registration
                vsI = obj.CloneVolSeq(obj.vsResampled);
                %obj.tType = {'none','rigid'}; % Registration scheme                
                vsO = obj.RegisterVolSeq(vsI, iVolRef);
                obj.vsRegistered = obj.CloneVolSeq(vsO);
                % Upate, NGR 2018 03 21
                switch eT
                    case MCC_eType.Registered
                        obj.eTypeCur = eT;
                        obj.Compare;
                        obj.Write;
                    case MCC_eType.Registered_Date 
                        obj.eTypeCur = eT;
                        obj.vsRegistered_Date = ...
                            obj.CloneVolSeq(obj.vsRegistered);
                end
            else
               % Local registration 
               obj.Register_Local(VolRef, options);
            end                       
        end

        
        function vsO = Register(obj, VolRef, options, options2, options3)
            % Registration (volume sequence)
            % NGR 2016 02 09, NGR 2016 06 22
            
            if ~exist('options')
                options = [];
            end     
            if ~exist('options2')
                options2 = [];
            end 
            if ~isempty(options) % NGR 2018 03 20
                % Register whole volume (no VOI)
                obj.ClearVolSeq(MCC_eType.Registered);
                % Check data
                obj.LoadData(MCC_eType.Original);
                % Input arguments
                switch nargin
                    case 1
                        if isempty(obj.iVolRef_Registered)
                            iVolRef = obj.GetLowestResolutionIndex;
                        else
                            iVolRef = obj.iVolRef_Registered;
                        end
                    case 2
                        if ~isnumeric(VolRef)
                            iVolRef = obj.GetContentIndex(VolRef);
                            if isempty(iVolRef)
                                iVolRef = obj.GetLowestResolutionIndex;
                            end
                        else
                            iVolRef = VolRef;
                        end
                    case {3, 4}
                        if ~isnumeric(VolRef)
                            iVolRef = obj.GetContentIndex(VolRef);
                            if isempty(iVolRef)
                                iVolRef = obj.GetLowestResolutionIndex;
                            end
                        else
                            iVolRef = VolRef;
                        end
                        % NGR 2018 03 20
                        eT = MCC_eType.Registered;  %9/3/2020
                        switch options
                            case 'none'                 % AA 1/26/2021
                                obj.tType = {'none'};
                                eT = MCC_eType.Registered;
                            case 'rigid'
                                obj.tType = {'none','rigid'};
                                eT = MCC_eType.Registered;     % NGR 2018 04 09
                            case 'elastic'
                                obj.tType = {'none','rigid','similarity','affine'};
                                eT = MCC_eType.Registered;     % NGR 2018 04 09
                            case 'date'
                                obj.tType = options2;
                                eT = MCC_eType.Registered_Date;     % NGR 2018 03 21
                        end
                                
                end
                obj.iVolRef_Registered = iVolRef;
                % Resample data (if already not)
                if obj.IsVolSeqEmpty(MCC_eType.Resampled)
                    obj.Resample(iVolRef);
                end
                obj.eViewCur = MCC_eView.All;
                % Registration
                vsI = obj.CloneVolSeq(obj.vsResampled);
                %obj.tType = {'none','rigid'}; % Registration scheme                
                vsO = obj.RegisterVolSeq(vsI, iVolRef);
                obj.vsRegistered = obj.CloneVolSeq(vsO);
                % Upate, NGR 2018 03 21
                switch eT
                    case MCC_eType.Registered
                        obj.eTypeCur = eT;
                        %obj.Compare;           %1/2/2019
                        obj.Write;
                    case MCC_eType.Registered_Date 
                        obj.eTypeCur = eT;
                        obj.vsRegistered_Date = ...
                            obj.CloneVolSeq(obj.vsRegistered);
                end
            else
               % Local registration 
               obj.Register_Local(VolRef, options3);
            end                       
        end
        
        function Register_Local(obj, VolRef, options)  
            % Local registration (selectable VOI)
            % NGR 2016 09 12
            
            % Parameters
            eT = MCC_eType.Parametric;
            eV = MCC_eView.Axial;
                                
            % VOI sorting
            
            iROI = str2num(options(end));
            
            eval(['eT  = MCC_eType.Registered_VOI' num2str(iROI) ';']);
            obj.tType = {'none','rigid','similarity','affine'};
            ext = obj.GetVolSeqLabel(eT);
            % FMX
            switch obj.eOrganCur
                case MCC_eOrgan.Tubes_FMX
                    switch iROI
                        case 1  % Tubes
                            obj.tType = {'none','rigid'};
                        case 2  % Tubes
                            obj.tType = {'none','rigid'};
                    end                
                case MCC_eOrgan.Brain_FMX
                    switch iROI
                        case 1  % Tissue
                            obj.tType = {'none','rigid'};
                        case 2  % Tubes
                            obj.tType = {'none','rigid'};
                    end
                case MCC_eOrgan.Body_FMX
                    switch iROI
                        case 1  % Tissue
                            obj.tType = {'none','rigid','similarity','affine'};
                        case 2  % Tubes
                            obj.tType = {'none','rigid'};
                    end
            end            
            % Check data
            obj.LoadData;
            obj.LoadWorkspace;
            if isempty(obj.VOI_mm)
                return;
            end
            % Select data
            obj.Resample(VolRef);
            if obj.IsVolSeqEmpty(MCC_eType.Registered)               
                vsI = obj.GetVolSeq(MCC_eType.Resampled);
            else
                if iROI==1
                    vsI = obj.GetVolSeq(MCC_eType.Registered);  
                else
                    eval(['vsI = obj.GetVolSeq(MCC_eType.Registered_VOI'...
                        num2str(iROI-1) ');'])
                    %vsI = obj.GetVolSeq(MCC_eType.Registered_VOI1);
                end
                    
            end            
%             if iROI == 2 % Tube registration   %2/22
%                 vsI = obj.GetVolSeq(MCC_eType.Resampled);
%             end
            % Force same FOV     
            iVolRef = obj.iVolRef_Resampled;
            Vf = obj.vsResampled(iVolRef);
            for iVol = 1:length(obj.vsResampled);
                Vm = vsI(iVol);
                Vt = obj.ForceFOV(Vf, Vm);
                vsI(iVol).data = Vt.data;
            end
            % Crop
            VOIx = obj.VOI_mm(iROI,:);
            VOIx = obj.InteresctFOVandVOI(vsI, VOIx);
            vsC = obj.CropVolSeqVOI(vsI, VOIx);
            % Register
            vsO = obj.RegisterVolSeq(vsC, iVolRef);
            vsC = obj.CloneVolSeq(vsO);     
            if iROI==1
                eval(['obj.vs' ext '= obj.CloneVolSeq(obj.vsRegistered);']);
            else
                eval(['obj.vs' ext '= obj.CloneVolSeq(obj.vsRegistered_VOI'...
                    num2str(iROI-1) ');']);  % 2/22
            end
            for iVol = 1:obj.GetNbObjects %2:obj.GetNbObjects 2/22
                Vi = vsC(iVol).data;                 
                eval(['obj.vs' ext '(iVol).FillVOI(Vi, VOIx);']);                
            end
            % Mask            
            if 0
                v = MCC_CloneObj(obj.vsResampled(1));
                v.data = zeros(size(v.data));
                v.FillVOI(obj.maskReg, VOIx);
                obj.maskReg = uint16(v.data);
            end
            % Upate              
            eval(['obj.iVolRef_' ext ' = iVolRef;']);
            obj.eTypeCur  = eT;
            %obj.Compare;
            for I=1:obj.GetNbObjects          % 11/6/2020 fill gaps
                if iROI==1
                    etOld='Registered'; 
                else
                    etOld=['Registered_VOI' num2str(iROI-1)]; 
                end
                etS=['Registered_VOI' num2str(iROI)];
                eval(['obj.vs' etS '(I).data(obj.vs' etS '(I).data==0)=obj.vs' etOld '(I).data(obj.vs' etS '(I).data==0);'])
            end
            obj.Write;  
            % Save roi registration
            R=obj.vsRegistered_VOI1(1).GetCropVec(VOIx);
            Reg_ROI=zeros(size(obj.vsRegistered_VOI1(1).data));
            Reg_ROI(R{1},R{2},R{3}) = 1;
            save([obj.root '\Processed\' ext '\roi.mat'],'Reg_ROI') 
        end
                        
        %-----------------------------------------------------------------%
        % Quantification
        %-----------------------------------------------------------------%
                
        function Quantify_FMX(obj)
            % Compute T2*, R2*
            % NGR 2016 08 23

            % Parameters
            eT = MCC_eType.Parametric;
            if strfind(lower(obj.StudyDateID), 'coronal') % NGR 2016 10 20
                eV = MCC_eView.Coronal;
            else
                eV = MCC_eView.Axial;
            end            
            if isempty(obj.vsRegistered_VOI1)
                vs    = obj.GetVolSeq(MCC_eType.Registered);
            else
                vs = obj.GetVolSeq(MCC_eType.Registered_VOI1);
            end
            % Correct image contour size, NGR 2017 01 11
            if 0
                s1 = size(obj.vsOriginal(1).dataRGB{1}(:,:,1));
                s2 = size(obj.vsOriginal(obj.iVolRef_Registered).data);                
                x = ((s1(1)/2)-s2(1)/2):((s1(1)/2)+s2(1)/2-1);
                y = ((s2(2)/2)-s1(2)/2):((s2(2)/2)+s1(2)/2-1);
                Nc = length(obj.vsOriginal(1).dataRGB);
                obj.vsOriginal(1).data = zeros(s2);
                for iC = 1:Nc,
                    Ic = obj.vsOriginal(1).dataRGB{iC};
                    x0 = uint8(zeros([s2(1:2) 3]));
                    x0(:,y,1) = Ic(x,:,1);
                    x0(:,y,2) = Ic(x,:,2);
                    x0(:,y,3) = Ic(x,:,3);
                    obj.vsOriginal(1).dataRGB{iC} = x0;
                    % RGB countour to mask
                    NN = size(x0);
                    roi = zeros([NN(1) NN(2)]);
                    for l = 1:NN(1)
                        for c= 1:NN(2)
                            if (x0(l,c,1) == x0(l,c,2)) && (x0(l,c,2) == x0(l,c,3))
                                roi(l,c) = 0;
                            else
                                roi(l,c) = 1;
                            end
                        end
                    end
                    mask = imfill(roi,'holes');
                    if sum(mask(:)>0)
                        iZ = obj.vsOriginal(1).metadata{iC}.InstanceNumber;
                        obj.vsOriginal(1).data(:,:,iZ) = mask;
                    end
                end
            end
            % Quantify tubes
            switch obj.eOrganCur
                case MCC_eOrgan.Tubes_FMX                    
                    Nseeds = 9; 
                    if eV == MCC_eView.Coronal % Plasma tubes
                        Nseeds = 2;     
                    end
                case MCC_eOrgan.Brain_FMX
                    Nseeds = 2;
                case MCC_eOrgan.Body_FMX
                    Nseeds = 8; % (4:MM002-001-0022-LAK)
                otherwise
                    Nseeds = 4;
            end
            
            % Tubes analysis
            bTubes =1;
            if bTubes
                CD = cd;
%                try
                switch eV
                    case MCC_eView.Axial
                        out = obj.Quantify_FMX_Tubes_Seeds(Nseeds);
                    case MCC_eView.Coronal
                        out = obj.Quantify_FMX_Tubes_Seeds_CoronalView(Nseeds);
                        eV = MCC_eView.Coronal;
                end
                Npix_tubes = out(:,1);
                avROI_GOF_tubes = out(:,2);
                avROI_T2s_tubes = out(:,3);
                avR2s_tubes = out(:,4);
%                 catch ME
%                     MCC_writeLog([CD '\log.txt'], ...
%                         ME, '', 1);
%                     cd(CD);
%                     bTubes = 0;
%             end
            end
            % Paramatric volumes                             
            obj.ComputeParametricVolumes_FMX(1);
            vS0  = obj.vsParametric(1).data;
            vT2s = obj.vsParametric(2).data;
            vR2s = obj.vsParametric(3).data;
            vGOF = obj.vsParametric(4).data;
            % Parametric volumes (for display)
            vS0_disp  = obj.vsParametric(1).data; 
            vT2s_disp = obj.vsParametric(2).data; 
            vR2s_disp = obj.vsParametric(3).data; 
            vGOF_disp = obj.vsParametric(4).data;
            % Quantification
            Nroi = length(obj.vsOriginal(1).dataRGB);
            count = 0;  
            Npix      = [];
            avROI_T2s = [];
            avROI_GOF = [];
            iSlice    = [];
            iROIcount = 0;
            for iROI = 1:Nroi
                iPlane = obj.vsOriginal(1).metadata{iROI}.InstanceNumber;
                mask = obj.vsOriginal(1).data(:,:,iPlane);
                % Muiltiple masks
                CC = bwconncomp(mask);
                if CC.NumObjects>1
                    for iCC = 1:CC.NumObjects
                        tmp = zeros(size(mask));
                        tmp(CC.PixelIdxList{iCC}) = 1;
                        masks(:,:,iCC) =  tmp;
                    end
                    iROIcount = iROIcount + CC.NumObjects;
                else
                    masks = mask;
                    iROIcount = iROIcount + 1;
                end
                for iCC = 1:CC.NumObjects
                    % Countour
                    se = strel('disk',1); % Morpho structuring element
                    contour = masks(:,:,iCC) - imerode(masks(:,:,iCC),se);
                    % Find connected components in binary image
                    index = iROIcount - CC.NumObjects + iCC;
                    % Masking
                    I_T2s = vT2s(:,:,iPlane);
                    Iroi_T2s = I_T2s(find(masks(:,:,iCC)));                    
                    I_GOF = vGOF(:,:,iPlane);
                    Iroi_GOF = I_GOF(find(masks(:,:,iCC)));
                    % Statistics
                    Imask = masks(:,:,iCC);
                    Npix(index)=sum(Imask(:));                    
                    % Average
                    [val,count] = obj.MeanMaskNaN(Iroi_T2s);
                    Npix_mask(index) = count;
                    if ~isempty(val)
                        avROI_T2s(index) = obj.MeanMaskNaN(Iroi_T2s);                        
                    else
                        avROI_T2s(index) = 0;                        
                    end
                    val = obj.MeanMaskNaN(Iroi_GOF);
                    obj.MeanMaskNaN(Iroi_GOF);                    
                    if ~isempty(val)
                        avROI_GOF(index) = obj.MeanMaskNaN(Iroi_GOF);                    
                    else
                        avROI_GOF(index) = 0;
                    end    
                    iSlice(index) = iPlane;
                    
                    %%%%%%%%%%%
                    % Display %
                    %%%%%%%%%%%                                    
                    
                    % Dimensions
                    R = obj.Volref3d(vs(2));
                    RI = imref2d(size(I_T2s));
                    RI.XWorldLimits = R.XWorldLimits;
                    RI.YWorldLimits = R.YWorldLimits;
                    % Figure
                    fig = figure(...
                        'units','normalized',...
                        'outerposition',[0 0 1 1],...
                        'Color',[0 0 0], ...
                        'WindowScrollWheelFcn', @MCC_WindowScrollWheelFcn, ...
                        'WindowButtonDownFcn', @MCC_WindowButtonDownFcn, ...
                        'WindowButtonUpFcn', @MCC_WindowButtonUpFcn ...
                        );
                    % Annotation [x y dX dY]
                    %
                    % ClinicalTrials.gov ID: NCT01770353
                    %                        Merrimack Pharmaceuticals
                    %
                    % ClinicalTrials.gov ID: NCT02367196
                    %                        Celgene Corporation
                    annotation(fig,'textbox',...
                        [0 0.88 0.17 0.12],...
                        'String',...
                        { ...
                        ['ClinicalTrials.gov ID: NCT01770353'],...
                        ['PID  : ' strrep(obj.PatientID,'_', ' ')],...
                        ['Study: ' strrep(obj.StudyDateID, '_', ' ')], ...
                        ['Type : ' strrep(eT.char,'_',' ')], ...
                        ['View : ' eV.char], ...
                        ['Comment: FMX MRI' ] ...
                        },...
                        'FitBoxToText','off',...
                        'BackgroundColor',[0 0 0], ...
                        'Color', [238, 130, 238] / 255, ... % Violet (web)
                        'FontSize',14 ...
                        );
                    annotation(fig,'textbox',...
                        [0.842940301120441 0.0332839174765179 0.525210083564018 0.0303522262078559],...
                        'String', ...
                        { ...
                        ['Radiomics language (' MCC_GetVersionFormated ')'],...
                        MCC_GetInstituteSignature,...
                        }, ...
                        'FitBoxToText','off', ...
                        'Color',[0 1 1], ...
                        'FontSize',16  ...
                        );
                    % ROI
                    subplot(221);
                    Id = obj.vsOriginal(1).dataRGB{iROI};
                    try 
                        Id = Id(1:RI.ImageSize(1),1:RI.ImageSize(2),:);                                        
                        imshow(Id,RI);
                        title(['Slice # = ' num2str(iPlane)]);
                        xlabel('[mm]');ylabel('[mm]');
                        cb = colorbar; cb.Color = [0 0 0];
                        obj.SetAxis;
                        obj.DrawROIs(MCC_eView.Axial);
                    catch
                        imshow(Id);
                        xlabel('[pix]');ylabel('[pix]');
                        cb = colorbar; cb.Color = [0 0 0];
                        obj.SetAxis;                        
                    end
                    % T2*
                    subplot(222);
                    Id = vT2s_disp(:,:,iPlane);
                    Id(find(contour))=100;
                    Id = obj.RecplaceNaNbyZero(Id);
                    imshow(Id, RI, [0 100]);
                    title(['T2* map (ms) : ROI_ ' num2str(index) ' = ' sprintf('%1.2f',avROI_T2s(index))]);
                    xlabel('[mm]');ylabel('[mm]');
                    cb = colorbar; cb.Color = [1 1 1];
                    obj.SetAxis;                    
                    % R^2
                    subplot(223);
                    Id = vGOF_disp(:,:,iPlane);
                    Id(find(contour))=1;
                    Id = obj.RecplaceNaNbyZero(Id);
                    imshow(Id, RI, [0 1]);
                    title(['R^2 : ROI_ ' num2str(index) ' = ' sprintf('%1.2f',avROI_GOF(index))]);
                    xlabel('[mm]');ylabel('[mm]');
                    cb = colorbar; cb.Color = [1 1 1];
                    obj.SetAxis;                    
                    % R2*
                    subplot(224);
                    Id = vR2s_disp(:,:,iPlane);
                    Id(find(contour))=300;
                    Id = obj.RecplaceNaNbyZero(Id);
                    imshow(Id, RI, [0 300]); 
                    title(['R2* (1/s) : ROI_ ' num2str(index) ' = ' sprintf('%1.2f',1e3/avROI_T2s(index))]);
                    xlabel('[mm]');ylabel('[mm]');
                    cb = colorbar; cb.Color = [1 1 1];
                    obj.SetAxis;                    
                    % Colormap
                    mapJ = jet(256);
                    mapJ(1,:) = [0 0 0];
                    colormap(mapJ);
                    % Save
                    obj.SaveFigure(gcf,[],eT,MCC_eView.Axial,0,index);
                end
            end
            % Save CSV file
            fid_csv = obj.CreateQuantifyCSV;            
            switch eV
                case MCC_eView.Axial
                    Apix = ... % pixel area
                        obj.vsResampled(2).res_mm(1)*obj.vsResampled(2).res_mm(2);
                case MCC_eView.Coronal
                    Apix = ... % pixel area
                        obj.vsResampled(2).res_mm(1)*obj.vsResampled(2).res_mm(2);                    
            end
            Nroi = length(avROI_GOF);
            for i=1:Nroi
                fprintf(fid_csv, ...
                    '%d,%d,%d,%1.2f,%d,%1.0f,%1.2f,%1.2f,%1.2f\n', ...
                    i, ...
                    iSlice(i), ...
                    Npix(i), ...
                    Npix(i)*Apix, ...
                    Npix_mask(i), ...
                    100*Npix_mask(i)/Npix(i), ...
                    avROI_GOF(i), ...
                    avROI_T2s(i), ...
                    1e3/avROI_T2s(i) ...                    
                    );
            end
            if bTubes
                Ntubes = length(avROI_GOF_tubes);
                for i=1:length(avROI_GOF_tubes) %  % Tubes
                    fprintf(fid_csv, ...
                        ',,%d,,,,%1.2f,%1.2f,%1.2f, Tube %d\n', ...
                        Npix_tubes(i), ...
                        avROI_GOF_tubes(i), ...
                        avROI_T2s_tubes(i), ...
                        avR2s_tubes(i), ...
                        i ...
                        );
                end
            else
                Ntubes = 0;
            end
            % DICOM tags for quality control (QC), NGR 2016 09 16
            N = 34 - (Nroi + Ntubes + 8);
            for l = 1:N
                fprintf(fid_csv, '%s\n', ' ');
            end
            obj.SaveDICOMtags_FMX(fid_csv);
            % Close
            fclose(fid_csv);
        end
        

        function SaveDICOMtags_FMX(obj,fid_csv)
            % Save DICOM tags for QC
            % NGR 2016 09 16
            
             bNew = 0;
            if ~exist('fid_csv')
                fileN_csv = [ ...
                    obj.root '\Processed\' ...
                    obj.MakeFileName(MCC_eType.DICOM,MCC_eView.All,0) ...
                    '.csv' ];
                fid_csv = fopen(fileN_csv,'wt+');
                bNew = 1;
            end            
            N = length(obj.vsOriginal);
            fprintf(fid_csv, '%s',...
                'Group Element,Description,');
            for i = 2:N
                fprintf(fid_csv, 'TE = %1.1f ms Scan,', obj.vsOriginal(i).metadata{1}.EchoTime);
            end
            
            fprintf(fid_csv, '\n(0018 0087),MagneticFieldStrength,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.MagneticFieldStrength);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end            
            fprintf(fid_csv, '\n\n(0018 0020),ScanningSequence,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%s,', obj.vsOriginal(i).metadata{1}.ScanningSequence);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 0021),SequenceVariant,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%s,', obj.vsOriginal(i).metadata{1}.SequenceVariant);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 0022),ScanOptions,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%s,', obj.vsOriginal(i).metadata{1}.ScanOptions);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 0024),SequenceName,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%s,', obj.vsOriginal(i).metadata{1}.SequenceName);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end            
            fprintf(fid_csv, '\n\n(0018 0050),SliceThickness,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.SliceThickness);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 0088),SpacingBetweenSlices,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.SpacingBetweenSlices);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end            
            fprintf(fid_csv, '\n\n(0018 0050),RepetitionTime,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.SliceThickness);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 0081),EchoTime,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.EchoTime);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 1314),FlipAngle,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.FlipAngle);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end            
            fprintf(fid_csv, '\n\n(0018 1312),InPlanePhaseEncodingDirection,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%s,', obj.vsOriginal(i).metadata{1}.InPlanePhaseEncodingDirection);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 0089),NumberOfPhaseEncodingSteps,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.NumberOfPhaseEncodingSteps);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 0095),PixelBandwidth,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%f,', obj.vsOriginal(i).metadata{1}.PixelBandwidth);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0018 1310),AcquisitionMatrix,');
            for i = 2:N
                try
                    x = obj.vsOriginal(i).metadata{1}.AcquisitionMatrix;
                    fprintf(fid_csv, '%d\\%d\\%d\\%d,', x(1),x(2),x(3),x(4));
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end            
            fprintf(fid_csv, '\n\n(0028 0010),Rows,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%d,', obj.vsOriginal(i).metadata{1}.Rows);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0028 0011),Columns,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%d,', obj.vsOriginal(i).metadata{1}.Columns);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0028 0100),BitsAllocated,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%d,', obj.vsOriginal(i).metadata{1}.BitsAllocated);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0028 0101),BitsStored,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%d,', obj.vsOriginal(i).metadata{1}.BitsStored);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end
            fprintf(fid_csv, '\n(0028 0102),HighBit,');
            for i = 2:N
                try
                    fprintf(fid_csv, '%d,', obj.vsOriginal(i).metadata{1}.HighBit);
                catch
                    fprintf(fid_csv, '%s,', ' ');
                end
            end            
            % Close
            if bNew
                fclose(fid_csv);
            end
        end
        
        function out = Quantify_FMX_Tubes_Seeds(obj, Nseeds)
            % Compute T2* in tubes - Region growing solution
            % NGR 2016 09 09
            
            % Parameters
            switch obj.eOrganCur
                case MCC_eOrgan.Tubes_FMX
                    pcTh = 0.05;  
                case MCC_eOrgan.Brain_FMX
                    pcTh = 0.05;                    
                case MCC_eOrgan.Body_FMX
                    pcTh = 0.08;  %pcTh = 0.12; 
            end
            % Extract VOI            
            iVOI = 2;
            if isempty(obj.vsRegistered_VOI2)
                Vroi = obj.ExtractVolVOI(...
                    MCC_eType.Registered_VOI2, ...
                    2, obj.VOI(iVOI).xyzDxDyDz_mm);                
            else
                Vroi = obj.ExtractVolVOI(...
                    MCC_eType.Registered, ...
                    2, obj.VOI(iVOI).xyzDxDyDz_mm);
            end
            Vroi = double(Vroi);
            szM = round(size(Vroi)/2);
            % Fitler out noise
            th=pcTh*max(Vroi(:));
            Vb = Vroi>th;                        
            % Erosion
            Vbe = zeros(size(Vb));
            sz = size(Vbe);
            switch obj.eOrganCur
                case MCC_eOrgan.Brain_FMX
                    Nse = 5;
                case MCC_eOrgan.Body_FMX
                    Nse = 3;
                otherwise
                    Nse = 3;
            end
            se = strel('disk', Nse);
            for z=1:sz(3)
                if z>=1 && z<=size(Vb,3)
                    Vbe(:,:,z) = imerode(Vb(:,:,z),se);
                end
            end            
            Vbe = double(Vbe);            
            % Capture seed points   
            if isempty(obj.pSeeds)
                figure;
                I = Vbe(:,:,szM(3));
                himage = imagesc(I); colormap gray; axis image;
                title(['Select ' num2str(Nseeds) ' seed points']);
                p = ginput(Nseeds);
                obj.pSeeds = [];       
                bCaptured = 1;
            else
                bCaptured = 0;
            end
            % Region growing
            if ~isempty(obj.pSeeds)
                Nseeds = size(obj.pSeeds,1)
            end            
            Vtubes = zeros(size(Vroi));
            for iSeed = 1:Nseeds
                if bCaptured
                    initPos(1) = round(axes2pix(size(I, 2), get(himage, 'XData'), p(iSeed,2)));
                    initPos(2) = round(axes2pix(size(I, 1), get(himage, 'YData'), p(iSeed,1)));
                    initPos(3) = szM(3);
                    obj.pSeeds(iSeed,:) = initPos;
                else
                    initPos = obj.pSeeds(iSeed,:);
                end               
                thresVal = 0.5; 
                [poly, mask] =  MCC_regionGrowing(....
                    squeeze(Vbe), initPos , thresVal, ...
                    Inf, [], false, false);              
                idx = find(mask);
                Vtubes(idx) = iSeed;
            end
            % Vertical erosion (top to bottom)             
            for iT=1:max(Vtubes(:))
                count = 0;
                for z=1:sz(3)
                    tmp = (Vtubes(:,:,z)==iT);                    
                    if sum(tmp(:)) && count<=2
                        idx = find(tmp);
                        I = Vtubes(:,:,z);
                        I(idx) = 0;
                        Vtubes(:,:,z) = I;
                        count = count + 1;
                    end
                end
            end
            % Vertical erosion (bottom to top)            
            for iT=1:max(Vtubes(:))
                count = 0;
                for z=sz(3):-1:1
                    tmp = (Vtubes(:,:,z)==iT);                    
                    if sum(tmp(:)) && count<=2
                        idx = find(tmp);
                        I = Vtubes(:,:,z);
                        I(idx) = 0;
                        Vtubes(:,:,z) = I;
                        count = count + 1;
                    end
                end
            end            
            % Parametric volumes            
            Vroi_param = obj.ComputeParametricVolumes_FMX(iVOI);
            vS0  = Vroi_param(:,:,:,1); 
            vT2s = Vroi_param(:,:,:,2);
            vR2s = Vroi_param(:,:,:,3);
            vGOF = Vroi_param(:,:,:,4);            
            % Display
            I = obj.RecplaceNaNbyZero(vT2s(:,:,szM(3)));
            T = Vtubes(:,:,szM(3))>=1;
            figure;            
            subplot(311);
            imagesc(Vroi(:,:,szM(3)));colormap hot;axis image;colorbar
            title('Shortest TE');
            subplot(312);
            imagesc(Vtubes(:,:,szM(3)));colormap hot;axis image;colorbar
            title('Tubes (detection)')
            subplot(313);            
            imagesc(I.*T);colormap hot;axis image;colorbar
            title('T2* (ms) in tubes');
            % Save figure
            eT = MCC_eType.Parametric;
            eV = MCC_eView.Axial;            
            fmt = 'png';
            fld = [obj.root '\Processed'];
            fn = [ ...
                fld '\' ...
                obj.MakeFileName(eT, eV) ...
                'Tubes.' fmt];            
            capt = getframe(gcf);
            imwrite(...
                capt.cdata, fn , fmt, 'Comment',MCC_GetArtistSignature);  
            % Update parametric volumes
            for iVol=1:length(obj.vsParametric)                
                Vroi = obj.ExtractVolVOI(...
                    MCC_eType.Parametric, ...
                    iVol, obj.VOI(iVOI).xyzDxDyDz_mm);                
                obj.vsParametric(iVol).FillVOI(...
                    Vroi.*(Vtubes>=1), obj.VOI_mm(iVOI,:));
            end            
            % Average
            Nvoi = max(Vtubes(:));            
            for iROI = 1:Nvoi                
                Iroi_T2s = vT2s(find(Vtubes(:)==iROI));
                Iroi_GOF = vGOF(find(Vtubes(:)==iROI));       
                Npix(iROI) = (sum(Vtubes(:)==iROI));                
                avROI_T2s(iROI) = obj.MeanMaskNaN(Iroi_T2s);
                avROI_GOF(iROI) = obj.MeanMaskNaN(Iroi_GOF); 
                avR2s(iROI) = 1e3/avROI_T2s(iROI);                
            end            
            out(:,1) = Npix;
            out(:,2) = avROI_GOF;
            out(:,3) = avROI_T2s;
            out(:,4) = avR2s;            
        end
        
        function out = Quantify_FMX_Tubes_Seeds_CoronalView(obj, Nseeds)
            % Compute T2* in tubes - Region growing solution
            % MM002-002-0034-NAK\20160426_Body_0h_preFMX
            % NGR 2016 09 09
            
            % Parameters
            switch obj.eOrganCur
                case MCC_eOrgan.Tubes_FMX
                    pcTh = 0.05; % pcTh = 0.02;
                case MCC_eOrgan.Brain_FMX
                    pcTh = 0.05;
                case MCC_eOrgan.Body_FMX
                    pcTh = 0.08; % pcTh = 0.05;
            end            
            % Extract VOI            
            iVOI = 2;
            if isempty(obj.vsRegistered_VOI2)
                Vroi = obj.ExtractVolVOI(...
                    MCC_eType.Registered_VOI2, ...
                    2, obj.VOI(iVOI).xyzDxDyDz_mm);                
            else
                Vroi = obj.ExtractVolVOI(...
                    MCC_eType.Registered, ...
                    2, obj.VOI(iVOI).xyzDxDyDz_mm);
            end
            % Change volume orientation
            Vroi = obj.SwitchView(Vroi, ...
                MCC_eView.Coronal, MCC_eView.Axial);            
            sz =size(Vroi);
            szM = round(sz/2);
            % Fitler out noise
            th=pcTh*max(Vroi(:));
            Vb = Vroi>th;                        
            % Erosion
            Vbe = zeros(size(Vb));
            Nse = 0;
            if Nse
                se = strel('disk', Nse);
                for z=1:sz(3)
                    if z>=1 && z<=size(Vb,3)
                        Vbe(:,:,z) = imerode(Vb(:,:,z),se);
                    end
                end
            else
                Vbe = Vb;
            end
            
            % MM002-002-0042-L-C\20160727_Body_16-24h_postFMX_Coronal
            % Vbe(:,end,:) = 0; %@NGR            
            % Vbe(:,(end-4:end),:) = 0; %@NGR
            
            % Capture seed points            
            if isempty(obj.pSeeds)
                figure;                
                I = Vbe(:,:,szM(3));
                subplot(211);
                imagesc(Vroi(:,:,szM(3))');                
                subplot(212);
                himage = imagesc(I); colormap gray;
                title(['Select ' num2str(Nseeds) ' seed points']);
                p = ginput(Nseeds);
                obj.pSeeds = [];       
                bCaptured = 1;
            else
                bCaptured = 0;
            end
            % Region growing
            Vtubes = zeros(size(Vroi));
            for iSeed = 1:Nseeds
                if bCaptured
                    initPos(1) = round(axes2pix(size(I, 2), get(himage, 'XData'), p(iSeed,2)));
                    initPos(2) = round(axes2pix(size(I, 1), get(himage, 'YData'), p(iSeed,1)));
                    initPos(3) = szM(3);
                    obj.pSeeds(iSeed,:) = initPos;
                else
                    initPos = obj.pSeeds(iSeed,:);
                end               
                thresVal = 0.5; 
                [poly, mask] =  MCC_regionGrowing(....
                    squeeze(Vbe), initPos , thresVal, ...
                    Inf, [], false, false);              
                idx = find(mask);
                Vtubes(idx) = iSeed;
            end            
            % Vertical erosion (top to bottom)             
            for iT=1:max(Vtubes(:))
                count = 0;
                for z=1:sz(3)
                    tmp = (Vtubes(:,:,z)==iT);                    
                    if sum(tmp(:)) && count<=2
                        idx = find(tmp);
                        I = Vtubes(:,:,z);
                        I(idx) = 0;
                        Vtubes(:,:,z) = I;
                        count = count + 1;
                    end
                end
            end
            % Vertical erosion (bottom to top)            
            for iT=1:max(Vtubes(:))
                count = 0;
                for z=sz(3):-1:1
                    tmp = (Vtubes(:,:,z)==iT);                    
                    if sum(tmp(:)) && count<=2
                        idx = find(tmp);
                        I = Vtubes(:,:,z);
                        I(idx) = 0;
                        Vtubes(:,:,z) = I;
                        count = count + 1;
                    end
                end
            end            
            % Parametric volumes            
            Vroi_param = obj.ComputeParametricVolumes_FMX(iVOI);
            % Rotate volume
            for iP = 1:size(Vroi_param,4)
                Vroi_param_(:,:,:,iP) = obj.SwitchView(...
                    Vroi_param(:,:,:,iP), ...
                    MCC_eView.Coronal, MCC_eView.Axial);
                Vroi_param_(:,:,:,iP) = obj.SwitchView(...
                    Vroi_param(:,:,:,iP), ...
                    MCC_eView.Coronal, MCC_eView.Axial);
            end            
            vS0  = Vroi_param_(:,:,:,1); 
            vT2s = Vroi_param_(:,:,:,2);
            vR2s = Vroi_param_(:,:,:,3);
            vGOF = Vroi_param_(:,:,:,4);
            % Display            
            I = obj.RecplaceNaNbyZero(vT2s(:,:,szM(3)));
            T = Vtubes(:,:,szM(3))>=1;
            figure;            
            subplot(311);
            imagesc(Vroi(:,:,szM(3))');colormap hot;axis image;colorbar
            title('Shortest TE');
            subplot(312);
            imagesc(Vtubes(:,:,szM(3))');colormap hot;axis image;colorbar
            title('Tubes (detection)')
            subplot(313);            
            imagesc((I.*T)');colormap hot;axis image;colorbar
            title('T2* (ms) in tubes');
            % Save figure
            eT = MCC_eType.Parametric;
            eV = MCC_eView.Axial;            
            fmt = 'png';
            fld = [obj.root '\Processed'];
            fn = [ ...
                fld '\' ...
                obj.MakeFileName(eT, eV) ...
                'Tubes.' fmt];            
            capt = getframe(gcf);
            imwrite(...
                capt.cdata, fn , fmt, 'Comment',MCC_GetArtistSignature);          
            % Average
            Nvoi = max(Vtubes(:));            
            for iROI = 1:Nvoi                
                Iroi_T2s = vT2s(find(Vtubes(:)==iROI));
                Iroi_GOF = vGOF(find(Vtubes(:)==iROI));       
                Npix(iROI) = (sum(Vtubes(:)==iROI));                
                avROI_T2s(iROI) = obj.MeanMaskNaN(Iroi_T2s);
                avROI_GOF(iROI) = obj.MeanMaskNaN(Iroi_GOF); 
                avR2s(iROI) = 1e3/avROI_T2s(iROI);                
            end            
            out(:,1) = Npix;
            out(:,2) = avROI_GOF;
            out(:,3) = avROI_T2s;
            out(:,4) = avR2s;    
            % Change tubes orientation
            Vtubes = obj.SwitchView(...
                    Vtubes, ...
                    MCC_eView.Axial, MCC_eView.Coronal);
            % Update parametric volumes
            for iVol=1:length(obj.vsParametric)
                Vroi = obj.ExtractVolVOI(...
                    MCC_eType.Parametric, ...
                    iVol, obj.VOI(iVOI).xyzDxDyDz_mm);                
                obj.vsParametric(iVol).FillVOI(...
                    Vroi.*(Vtubes>=1), obj.VOI_mm(iVOI,:));
            end
        end
        
        function Vo = SwitchView(obj, Vi, eVs, eVd)
            % Rotate volume according to view 
            % NGR 2016 09 16
                        
            % Coronal -> Axial
            if eVs == MCC_eView.Coronal && eVd == MCC_eView.Axial
                sz = size(Vi);
                Vo = zeros([ sz(2) sz(3) sz(1)]);
                for iZ = 1:size(Vi,1)
                    Vo(:,:,iZ) = squeeze(Vi(iZ,:,:));
                end                
            end
            % Axial -> Coronal
            if eVs == MCC_eView.Axial && eVd == MCC_eView.Coronal
                sz = size(Vi);
                %Vo = zeros([ sz(3) sz(2) sz(1)]);
                for iZ = 1:size(Vi,3)
                    Vo(iZ,:,:) = squeeze(Vi(:,:,iZ));
                end                
            end
        end
        
        function Vroi_param = ComputeParametricVolumes_FMX(obj, iVOI)
            % Compute parametric volumes
            % NGR 2016 09 14
            
            % Parameters
            pcTEmask = 0.05;
            thGOF = 0.2;
            % Extract shortest TE
            eT = eval(['MCC_eType.Registered_VOI' num2str(iVOI)]);
            Vroi = obj.ExtractVolVOI(...
                eT, 2, obj.VOI(iVOI).xyzDxDyDz_mm);
            sz = size(Vroi);
            Vroi = double(Vroi);
            % Extract volume
            Nvol  = obj.GetNbObjects;
            seqVroi = zeros([size(Vroi) Nvol-1]);
            for iVol=2:Nvol
%                 if iVol ==7
%                     iVol
%                 end
                Vvoi = obj.ExtractVolVOI(...
                eT, iVol, obj.VOI(iVOI).xyzDxDyDz_mm);
                seqVroi(:,:,:,iVol-1) = Vvoi(:,:,:);
                %seqVroi(1:34,:,:,iVol-1) = Vvoi(1:34,:,:);
            end
            % Parametric volumes initialization
            Vroi_param = zeros([sz 1]);
            TEmax = max(obj.vsRegistered(2).data(:)); % shortest TE            
            % Echo time vector
            TE = zeros(1,Nvol-1);
            for iVol=2:Nvol  % Note: 1 is countour volume index
                TE(iVol-1) = ...
                    double(obj.vsOriginal(iVol).metadata{1}.EchoTime);
            end
            % Compute parametric volumes
            for i=1:sz(1)
                for j=1:sz(2)
                    for k=1:sz(3)
                        Y = squeeze(seqVroi(i,j,k,:))';
                        % Linear regression in the log domain
                        x = TE; y = log(Y);
                        N = length(x);
                        covXY = (1/N) * sum((x-mean(x)).*(y-mean(y)));
                        varX  = (1/N) * sum((x-mean(x)).^2);
                        % Parameters
                        a = covXY / varX;
                        b = mean(y) - a*mean(x);
                        y_ = a*x + b;
                        % FMX parameters
                        S0  = exp(b);
                        T2s = -1/a;
                        S   = S0 * exp(-TE/T2s);
                        % Goodness of fit (R-squared)
                        SSR = sum((y_ - mean(y)).^2);
                        SST = sum((y  - mean(y)).^2);
                        GOF = SSR/SST;
                        if isnan(GOF)
                            GOF = 0;
                        else
                            if GOF<0
                                GOF = 0;
                            end
                        end
                        % Masking
                        thMask = pcTEmask*TEmax;
                        if (...
                                ~isnan(T2s) && ...
                                T2s>0       && ...
                                GOF>thGOF   && ...
                                Vroi(i,j,k)>thMask ...
                                )
                            Vroi_param(i,j,k,1) = S0;
                            Vroi_param(i,j,k,2) = T2s; 
                            Vroi_param(i,j,k,3) = 1/(T2s*1e-3);
                        else
                            Vroi_param(i,j,k,1) = NaN;
                            Vroi_param(i,j,k,2) = NaN;
                            Vroi_param(i,j,k,3) = NaN;
                        end
                        if Vroi(i,j,k)>thMask
                            Vroi_param(i,j,k,4) = GOF;
                        end
                    end
                end
            end
            % Paramatric volumes            
            names = {'S0', 'T2s','vR2s','GOF'};
            sz = size(obj.vsRegistered(2).data);  
            if isempty(obj.vsParametric)                
                obj.vsParametric = MCC_CloneObj(obj.vsRegistered(2));                
                bInstanciate = 1;
            else
                bInstanciate = 0;
            end
            for iVol=1:length(names)
                if bInstanciate
                    obj.vsParametric(iVol) = MCC_CloneObj(obj.vsRegistered(2));
                    obj.vsParametric(iVol).data = zeros(sz);
                    obj.vsParametric(iVol).name = names{iVol};
                    obj.vsParametric(iVol).eType = MCC_eType.Parametric;
                    obj.vsParametric(iVol).ImagingModeID = [];
                end
                obj.vsParametric(iVol).FillVOI(...
                    Vroi_param(:,:,:,iVol), obj.VOI_mm(iVOI,:));
            end            
        end                
        
        function fid_csv = CreateQuantifyCSV(obj)
            % Create CSV file for FMX quantification
            % NGr 2016 09 12, 2016 10 04
            fileN_csv = [ ...
                obj.root '\Processed\' ...
                obj.MakeFileName(MCC_eType.Quantified,MCC_eView.All,0) ...
                '.csv' ];
            fid_csv = fopen(fileN_csv,'wt+');
            fprintf(fid_csv,...
                [ ...
                'Version:,%s\nClock:,%s\n\n'...
                'ClinicalTrials.gov ID:,NCT01770353\n' ...
                'PID:,%s\nExam:,%s\n\n' ...
                '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'], ...
                MCC_GetVersionFormated, ...
                MCC_clock2str, ...
                strrep(obj.PatientID, ',', ' '), ...
                strrep(obj.StudyDateID, ',', ' '), ...
                'ROI #', ...
                'Slice #', ...
                '# of pixels', ...
                'Area (mm^2)', ...
                '# of mask pixels', ...
                '% of mask pixels', ...
                'R^2', ...
                'T2* (ms)', ...
                'R2* (1/s)', ...
                'Comments' ...
                );
        end

        function Io = RecplaceNaNbyZero(obj, Ii)
            % Recplace NaN by Zero (image)
            % NGR 2016 09 02
            Io = Ii;
            for n=1:length(Ii(:))
                if isnan(Ii(n))
                    Io(n) = 0;
                end
            end
        end
        
        function [out, count]= MeanMaskNaN(obj, in)
            % Average of the elements (excluding NaN)
            % NGR 2016 09 07
            in = in(:);
            data = [];
            count = 0;
            for id = 1:length(in)
                if ~isnan(in(id))
                    count = count + 1;
                    data(count) = in(id);                                        
                end
            end
            if ~isempty(data)
                out = mean(data);
            else
                out = [];
                count = 0;
            end
        end
        
        %-----------------------------------------------------------------%
        % Files
        %-----------------------------------------------------------------%
        
        function UpdateImagingModes(obj)
            % Update imaging modes (work in progress)
            % 2016 03 09
            for i=1:length(obj.content)
                fld = [obj.root '\' obj.content{i}];
                lst = dir(fld);
                CD = cd;
                if (length(lst)>2)
                    try
                        file = [fld '\' lst(3).name];
                        md = MCC_ReadMetadataDICOM(file);
                        obj.content{i} = md.SeriesDescription;
                    catch ME
                        MCC_writeLog([CD '\log.txt'], ...
                            ME, file);
                        cd(CD);
                    end
                end
            end
        end
        
        function Load(obj, eType, bFirstSlice)
            % Load DICOM data from files
            % NGR 2016 02 17
            if nargin<=2
                bFirstSlice = 0;
            end
            iP = 1;
            % Input arguments
            if (nargin<=1)
                eType = MCC_eType.Original;
            end
            % Contiguous slices
            obj.bContiguous{obj.GetNbObjects} = [];
            % Source folder
            fld = obj.GetSubDirPath(1,eType);
            if ~isdir(fld), return; end
            obj.DisplayExeLine(iP);
            MCC_Disp(fld,5);
            ext = obj.GetVolSeqLabel(eType);
            eval(['obj.vs' ext ...
                '=MCC_c3D(fld, bFirstSlice);']);
            eval( ['obj.vs' ext '(1).eType = eType;']);
            % Volume sequence loop
            for i=2:obj.GetNbObjects(eType)
                fld = obj.GetSubDirPath(i,eType);
                if ~isdir(fld), return; end
                MCC_Disp(fld,5);
                eval( ...
                    ['obj.vs' ext '(' num2str(i) ...
                    ')=MCC_c3D(fld, bFirstSlice);']);
                if eval(['obj.vs' ext '(' num2str(i) ').bLoadSuccess'])
                    eval(['obj.vs' ext '(' num2str(i) ').eType = eType;']);
                else
                    % Flag failure
                    obj.content{i} = [];
                end
            end
            % Keep success only
            if ...
                    ~bFirstSlice && ...
                    eType~=MCC_eType.Segmented && ...
                    ~isempty(obj.vsOriginal(1).name)
                j = 1;
                for i=1:obj.GetNbObjects
                    if ~isempty(obj.content{i})
                        x{j} = obj.content{i};
                        y(j) = eval(...
                            ['MCC_CloneObj(obj.vs' ext ...
                            '(' num2str(i) '));']);
                        j = j + 1;
                    end
                end
                obj.content = x;
                eval(['obj.vs' ext ' = [];']);
                eval(['obj.vs' ext ' = obj.CloneVolSeq(y);']);
            end
            % Mask management (radio therapy data structure)
            obj.CreateMasks;
            % Acquisition time vector
            if eType == MCC_eType.Original
                for i=1:obj.GetNbObjects
                    try
                        str = obj.vsOriginal(i).metadata{1}.AcquisitionTime;
                        obj.time(i) = ...
                            str2num(str(1:2))*3600 + ...
                            str2num(str(3:4))*60 + ...
                            str2num(str(5:end));
                    catch
                        obj.time(i) = i;
                    end
                end
                obj.time = obj.time - min(obj.time);
                [obj.time, ind] = sort(obj.time);
                obj.content = obj.content(ind);
                % Re-arrange data
                vs = obj.GetVolSeq(eType);
                for i=1:obj.GetNbObjects
                    obj.vsOriginal(i) = MCC_CloneObj(vs(ind(i)));
                end
            end
            % Update
            switch eType
                case MCC_eType.Original
                    obj.name = obj.vsOriginal(1).metadata{1}.StudyDate;
                    obj.PatientID = obj.vsOriginal(1).metadata{1}.PatientID;
                case MCC_eType.Resampled
                    obj.iVolRef_Resampled = obj.GetVolRefIndex(eType);
                case MCC_eType.Registered
                    obj.iVolRef_Registered = obj.GetVolRefIndex(eType);
            end
            obj.eTypeCur = eType;
            obj.DisplayExeTime(iP);
        end
        
        function LoadData(obj, eT)
            % Load data form database
            % NGR 2016 04 15
            if nargin == 1
                eT = [ ...
                    MCC_eType.Original ...
                    MCC_eType.Registered ...
                    MCC_eType.Registered_Date ... % NGR 2018 03 21
                    MCC_eType.Registered_Elastic ...                    
                    MCC_eType.Registered_VOI1 ...
                    MCC_eType.Registered_VOI2 ...
                    MCC_eType.Localized ...
                    MCC_eType.Segmented ...
                    ];
            end
            for i = 1:length(eT)
                if ...
                        obj.IsVolSeqDataEmpty(eT(i)) && ...
                        isdir(obj.GetSubDirPath(0, eT(i)))
                    obj.Load(eT(i));
                    if eT(i) == MCC_eType.Registered
                        if ~isempty(obj.iVolRef_Registered)
                            obj.Resample(obj.iVolRef_Registered);
                        else
                            obj.Resample;
                        end
                    end
                    obj.LoadWorkspace;
                    obj.eTypeCur = eT(i);
                end
            end
        end
        
        function Write(obj, eType)
            % Write DICOM files (in subdirectories, eg. Registered)
            % NGR 2016 02 08
            if ~isempty(obj.root)
                iP = 7;
                % Display
                obj.DisplayExeLine(iP);
                if (nargin<=1)
                    eType = obj.eTypeCur;
                end
                vs = obj.GetVolSeq(eType);
                % Write in database
                for i=1:obj.GetNbObjects(eType)
                    fld = obj.MakeSubDir(i,eType);
                    MCC_Disp(fld,5);
                    try
                        % Data
                        vs(i).Write(fld);
                    catch ME
                        MCC_writeLog('log.txt', ME, ...
                            obj.GetSubDirPath(i,eType));
                    end
                end
                % Mask
                if ~isempty(obj.maskReg)
                    vs(1).data = obj.maskReg;
                    vs(1).metadata{1}.SeriesDescription = 'maskReg';
                    fld = [obj.root '\Processed\Mask'];
                    % Cretate or remove content
                    if ~isdir(fld)
                        mkdir(fld);
                    else
                        CD = cd;
                        cd(fld);
                        delete(['*.*']);
                        cd(CD);
                    end
                    vs(1).Write(fld);
                end
                % Update
                obj.DisplayExeTime(iP);
            end
        end
        
        function fld = MakeSubDir(obj, iVol, eType, name)
            % Make a subdirectory in the database
            % NGR 2016 03 15
            if iVol
                switch nargin
                    case 1
                        fld = obj.GetSubDirPath(0, obj.eTypeCur);
                    case 2
                        fld = obj.GetSubDirPath(iVol, obj.eTypeCur);
                    case 3
                        fld = obj.GetSubDirPath(iVol, eType);
                    case 4
                        fld = obj.GetSubDirPath(iVol, eType, name);
                end
                % Cretate or remove content
                if ~isdir(fld)
                    mkdir(fld);
                else
                    CD = cd;
                    cd(fld);
                    delete(['*.dcm']);
                    delete(['*.jpg']);
                    cd(CD);
                end
            else
                fld = obj.root;
            end
        end
        
        function fld = GetSubDirPath(obj, iVol, eType, name)
            % Get subdirectory full path in the database
            % NGR 2016 03 31
            if iVol
                switch eType
                    case MCC_eType.Original
                        fld= [...
                            obj.root '\' ...
                            obj.content{iVol} ...
                            ];
                    case MCC_eType.Segmented
                        if obj.IsVolSeqEmpty(eType);
                            % Read mode
                            sd = MCC_subdir([ ...
                                obj.root '\Processed\' ...
                                char(eType) ...
                                ]);
                            for i=1:length(sd)
                                cellTxt = strsplit(sd{i},'\');
                                obj.contentSeg{i} = cellTxt{end};
                            end
                        end
                        s = obj.contentSeg{iVol};
                        fld= [...
                            obj.root '\Processed\' ...
                            char(eType) '\' ...
                            s ...
                            ];
                    otherwise
                        s=obj.vsOriginal(iVol).ImagingModeID;
                        fld= [...
                            obj.root '\Processed\' ...
                            char(eType) '\' ...
                            s ...
                            ];
                end
            else
                if eType == MCC_eType.Original
                    fld = obj.root;
                else
                    fld = [obj.root '\Processed\' char(eType)];
                end
            end
            if nargin==4 && ~isempty(name)
                fld = [fld '\' name];
            end
        end
        
        function WriteCSV(obj, fld, MI, method)
            % Write results in a CSV file
            % NGR 2016 03 16
            fileN_csv = [MCC_GetDateFileFormat 'Registration.csv'];
            if exist(fileN_csv, 'file')
                fid_csv = fopen(fileN_csv,'at+');
            else
                fid_csv = fopen(fileN_csv,'wt+');
                fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s\n', ...
                    'Clock','Name', ...
                    'MI (ini.)', 'MI (final)', 'Relative (%)', ...
                    'Method index','Method name');
            end
            % Write
            fprintf(fid_csv, '%s,%s,%1.2f,%1.2f,%1.2f,%d,%s\n', ...
                strrep(MCC_clock2str, ',', ' '), ...
                fld, ...
                MI(1), MI(2), MI(3), MI(4), ...
                method ...
                );
            % Close
            fclose(fid_csv);
        end
        
        function Save(obj)
            % Write results as DICOM files
            % NGR 2016 02 08
            iP = 6;
            obj.DisplayExeLine(iP);
            save([obj.root '\MCC_c4D']);
            obj.DisplayExeTime(iP);
        end
        
        function CreateMasks(obj)
            % Create masks
            % NGR 2016 04 25
            
            % Are DICOM masks already created?
            for iVol=1:obj.GetNbObjects
                if ~isempty(strfind(obj.content{iVol},'Mask'))
                    return;
                end
            end
            % Generate DICOM mask
            bVolSeqReload = 0;
            for iVol=1:obj.GetNbObjects
                scrFld = [obj.root '\' obj.content{iVol} '\RTSS'];
                if isdir(scrFld)
                    MCC_Disp('Create masks ...');
                    % Mask loop
                    d = dir(scrFld); d([d.isdir])= []; N = length(d);
                    for i=1:N
                        fname = d(i).name;
                        if ~strcmp(fname,'.') && ~strcmp(fname,'..')
                            % Load DICOM headers
                            rtssfile = [scrFld '\' fname];
                            rtssheader = MCC_ReadMetadataDICOM(rtssfile);
                            path = [obj.root '\' obj.content{iVol}];
                            imageheaders = ...
                                loadDicomImageInfo(...
                                path, rtssheader.StudyInstanceUID);
                            % Read contour sequences
                            contours = ...
                                ReadRTstructures(rtssheader, imageheaders);
                            % Save segmentations
                            for i=1:size(contours,2)
                                name = ['Mask ' contours(i).ROIName];
                                dstFld = fullfile(obj.root, name);
                                MCC_Disp(dstFld,5);
                                if ~exist(dstFld, 'dir')
                                    mkdir(dstFld);
                                end
                                v = MCC_CloneObj(obj.vsOriginal(iVol));
                                v.root = dstFld;
                                Vmask = 255*...
                                    single(contours(i).Segmentation) / ...
                                    max(contours(i).Segmentation(:));
                                v.data = uint16(Vmask);
                                v.name = name;
                                v.Write(dstFld);
                                bVolSeqReload = 1;
                            end
                        end
                    end
                end
            end
            % Reload data
            if bVolSeqReload
                i1 = obj.iVolRef_Original;
                i2 = obj.iVolRef_Resampled;
                i3 = obj.iVolRef_Registered;
                obj.content = MCC_SubFolderList(obj.root);
                obj.ClearVolSeq;
                obj.Load(MCC_eType.Original);
                obj.iVolRef_Original = i1;
                obj.iVolRef_Resampled = i2;
                obj.iVolRef_Registered = i3;
            end
        end
        
        %-----------------------------------------------------------------%
        % Display
        %-----------------------------------------------------------------%
        
        function DisplayExeLine(obj, iMode)
            % Display execution line
            % NGR 2016 02 17
            tic;
            disp(' ');
            switch iMode
                case 1
                    MCC_Disp('Load DICOM ...');
                case 2
                    MCC_Disp('Resample ...');
                case 3
                    MCC_Disp('Register ...');
                case 4
                    MCC_Disp('Segment ...');
                case 5
                    MCC_Disp('Habitat ...');
                case 6
                    MCC_Disp('Save ...');
                case 7
                    MCC_Disp('Write DICOM ...');
                case 8
                    MCC_Disp('Write JPG ...');
                otherwise
                    disp('? ...');
            end
        end
        
        function DisplayExeTime(obj, iMode)
            % Display execution time
            % NGR 2016 02 17
            obj.tExe(iMode) = toc;
            MCC_Disp([datestr(obj.tExe(iMode)/86400, 'HH:MM:SS.FFF')]);
            disp(' ');
        end
        
        function Display(obj, eView, eType)
            % Display volume sequence (central slice)
            % NGR 2016 02 05
            obj.LoadData;
            switch nargin
                case 1
                    eView  = MCC_eView.Planar;
                    eType = obj.eTypeCur;
                case 2
                    eType = obj.eTypeCur;
            end
            obj.eViewCur = eView;
            vs = obj.GetVolSeq(eType);
            Nvol = obj.GetNbVol(eType);
            N = Nvol;
            
            %eType = MCC_eType.Registered;
            
            % Special cases
            switch eView
                case MCC_eView.RGB
                    obj.DisplayRGB(0,0);
                    return
                case MCC_eView.FourQ
                    obj.Figure(eType, MCC_eView.VOIs);
                    if isempty(obj.vShiftPlaneCur)
                        obj.vShiftPlaneCur = [0 0 0];
                    end
                    obj.Display4Q(eType, obj.vShiftPlaneCur, 0);
                    return;
            end
            % Figure
            obj.Figure(eType, eView);
            colormap gray;
            % Subplot matrix
            txtSubp = '1,1';
            if (N>=1) && (N<=2)
                txtSubp = '1,2,';
            end
            if (N>2) && (N<=4)
                txtSubp = '2,2,';
            end
            if (N>4) && (N<=6)
                txtSubp = '3,2,';
            end
            if (N>6) && (N<=9)
                txtSubp = '3,3,';
            end
            if (N>9) && (N<=12)
                txtSubp = '4,3,';
            end
            if (N>12)
                txtSubp = '4,4,';
            end
            % Display
            if Nvol>=16
                Nvol = 16;
            end
            for i=1:Nvol
                v = vs(i);
                if N>=1
                    ax(i) = eval(['subplot(' txtSubp num2str(i) ');']);
                end
                switch eView
                    case MCC_eView.Planar
                        % Volume slice plot
                        s = size(v.data);
                        r = v.res_mm;
                        j=1;lx=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                        j=2;ly=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                        j=3;lz=linspace(s(j)*r(j)/2, -s(j)*r(j)/2, s(j));
                        [x,y,z] = meshgrid(ly,lx,lz);
                        slice(x,y,z,single(v.data),0,0,0);
                        shading flat;axis image;
                        % Axis
                        obj.SetAxis;
                    case MCC_eView.Volumic
                        Nc = 8;
                        % 3D rendering
                        obj.Figure(eType, eView);
                        %V = v.data;
                        th = 1;
                        D = vs(Nc-1).data;
                        s = size(v.data);
                        r = v.res_mm;
                        j=1;lx=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                        j=2;ly=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                        j=3;lz=linspace(s(j)*r(j)/2, -s(j)*r(j)/2, s(j));
                        [x, y, z] = meshgrid(ly,lx,lz);
                        %p = patch(isosurface(x,y,z, V, 0.5), 'FaceColor', 'red', 'EdgeColor', 'none');
                        se_sz = 10;
                        % Kidney
                        D = vs(Nc-1).data==Nc;
                        D = single(D);
                        for iS = 1:size(D,3)
                            se = strel('ball',se_sz,se_sz);
                            D(:,:,iS) = imerode(D(:,:,iS),se);
                            D(:,:,iS) = imdilate(D(:,:,iS),se);
                        end
                        D_ = D;
                        for iS = 1:size(D,3)
                            D(:,:,iS) = D_(:,:,size(D,3)-iS+1);
                        end
                        p1 = patch(isosurface(D),'FaceColor','red',...
                            'EdgeColor','red');
                        % Liver
                        D = vs(Nc-1).data==Nc-1;
                        D = single(D);
                        for iS = 1:size(D,3)
                            se = strel('ball',se_sz,se_sz);
                            D(:,:,iS) = imerode(D(:,:,iS),se);
                            D(:,:,iS) = imdilate(D(:,:,iS),se);
                        end
                        D_ = D;
                        for iS = 1:size(D,3)
                            D(:,:,iS) = D_(:,:,size(D,3)-iS+1);
                        end
                        p2 = patch(isosurface(D),'FaceColor','green',...
                            'EdgeColor','none');
                        % Liver
                        D = vs(Nc-1).data==Nc-1;
                        D = single(D);
                        for iS = 1:size(D,3)
                            se = strel('ball',se_sz,se_sz);
                            D(:,:,iS) = imerode(D(:,:,iS),se);
                            D(:,:,iS) = imdilate(D(:,:,iS),se);
                        end
                        D_ = D;
                        for iS = 1:size(D,3)
                            D(:,:,iS) = D_(:,:,size(D,3)-iS+1);
                        end
                        p3 = patch(isosurface(D),'FaceColor','blue',...
                            'EdgeColor','none');
                        axis tight
                        daspect([1,1,obj.vsRegistered(1).res_mm(1)/obj.vsRegistered(1).res_mm(3)])
                        camlight
                        lighting gouraud
                        isonormals(D,p1);hold on;
                        isonormals(D,p2);
                        isonormals(D,p3);
                        % Axis
                        obj.SetAxis;
                        axis off;
                        % Video
                        % Prepare the new file.
                        vidObj = VideoWriter('Cluster_3D.mp4');
                        open(vidObj);
                        for af=-180:180
                            view([-45 af 45]);
                            drawnow;
                            %writeVideo(vidObj,getframe(gcf));
                            %X(:,:,af) = frame2im(F);
                        end
                        % Close the file.
                        close(vidObj);
                        return
                        
                    otherwise
                        % Image
                        iPlane = obj.GetCentralOrthoPlaneIndex(eType, eView, i, 0);
                        [I,RI]=obj.GetOrthoPlane(eType,eView,i,iPlane);
                        % Colormap
                        if eType == MCC_eType.Segmented
                            V = obj.vsSegmented(i).data;
                            Nc = max(V(:)); % nb of clusters
                            map = jet(256);
                            map(1,:) = [0 0 0];
                            I = round(256*I/Nc);
                            I(find(I==min(I(:))))=1;
                            Irgb = map(I(:),:);
                            Irgb = reshape(Irgb,[size(I) 3]);
                            imshow(Irgb,RI);
                        else
                            imshow(I,RI,[min(I(:)) 0.75*max(I(:))]);
                        end
                        % Axis
                        obj.SetAxis;
                        axis(obj.GetMaxAxis(eType, eView));
                        if eType == MCC_eType.Ready
                            %axis([-100 100 -100 100]);
                        end
                end
                obj.SetLabels(eView);
                % Detect reference
                iVolRef = obj.GetVolRefIndex(eType);
                txtRef = '';
                if ~isempty(iVolRef)
                    if (i == iVolRef)
                        if eType ~= MCC_eType.Original
                            txtRef = ' [Reference]';
                        end
                    end
                end
                % Title
                title([ ...
                    strrep(MCC_RemoveSeriesNumber(v.ImagingModeID),'_', ...
                    ' ') '  (' ...
                    num2str(v.res_mm(1),2) ' x ' ...
                    num2str(v.res_mm(2),2) ' x ' ...
                    num2str(v.res_mm(3),2) ' = ' ...
                    num2str(v.GetVoxelSize_mm3,2) ' mm^3 ' ...
                    ') ' txtRef]);
            end
            % Save figure
            obj.SaveFigure(gcf, [] ,eType, eView, 0, 0);
            %linkprop([ax(1), ax(2) ax(3)],'View'); % WORK IN DEBUGING MODE
        end
                
        function Display4Q_AllSlices_(obj, eT)
            % Display registration for visual comparison
            % NGR 2016 05 09
            vs = obj.GetVolSeq(eT);
            N  = vs(1).GetNbSlices;
            R = floor(N/2)-1;
            for j=-R:R
                obj.Figure(eT, MCC_eView.All);
                obj.Display4Q(eType, [0 0 j], 1);
                close;
            end
        end
        
        function Display4Q(obj, ...
                eT, ...
                vShiftPlane, ...
                bAllSlices, bSaveFig)
            % Display RBG image (3 views)
            % NGR 2016 05 09, 2016 07 19
            
            if nargin<=4
                bSaveFig = 0;
            end
            % Data
            obj.LoadData(eT);
            if obj.IsVolSeqDataEmpty(eT)
                disp(' ');
                disp('No processed data !');
                disp(' ');
                return;
            end
            eV = MCC_eView.Axial;
            obj.vShiftPlaneCur = vShiftPlane;
            obj.vsTemporary = obj.GetVolSeq(eT);
            
            % Volume Of Interest (VOI)
            %iC  = obj.GetContentIndex('T2');
            iC  = obj.iVolRef_Registered;
            obj.InitializeVOI();
            for i =1:length(obj.VOI)
                if isempty(obj.VOI_mm)
                    s = size(obj.vsTemporary(iC).data);
                    r_mm = obj.vsTemporary(iC).res_mm;
                    s_mm = s .* r_mm;
                    o_mm = -[s_mm(2) s_mm(1) s_mm(3)]/2 + 10;
                    obj.VOI(i).xyzDxDyDz_mm =[ ...
                        (o_mm(1)+(i-1)*50) ...
                        o_mm(2) + 35 ...
                        o_mm(3) ...
                        25 25 25];
                end
            end
            % Display
            if isempty(obj.MaxAxis)
                ax = obj.GetMaxAxis(MCC_eType.Original, eV);
                obj.MaxAxis= ax;
            else
                ax = obj.MaxAxis;
            end
            for  k =1:3
                switch k
                    case 1
                        eV = MCC_eView.Axial;
                        subplot('221');
                    case 2
                        eV = MCC_eView.Coronal;
                        subplot('223');
                    case 3
                        eV = MCC_eView.Sagittal;
                        subplot('222');
                end
                % Plane index
                iPlane = obj.GetCentralOrthoPlaneIndex(...
                    MCC_eType.Temporary, eV, iC);
                switch k
                    case 1
                        iP = iPlane+vShiftPlane(3);
                        iPsaved = iP;
                    case 2
                        iP = iPlane+vShiftPlane(1);
                        iPsaved = iP;
                    case 3
                        iP = iPlane+vShiftPlane(2);
                        iPsaved = iP;
                end
                % Image
                [I,RI]=obj.GetOrthoPlane(...
                    MCC_eType.Temporary, eV, iC, iP);
                % Display
                f = 0.75;
                imshow(I,RI,[min(I(:)) f*max(I(:))]);
                % Axis
                obj.SetAxis;
                axis(obj.GetMaxAxis(eT, eV));
                axis(ax);
                % Title and line
                res_mm= obj.vsTemporary(1).res_mm;
                switch k
                    case 1
                        %  Axial
                        title(char(eV));
                        axColor = [1 1 1];
                        hlColor = [0 1 0];
                        vlColor = [1 1 0];
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(1);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(3) = iP;
                        % VOI
                        obj.DrawROIs(eV);
                    case 2
                        %  Coronal
                        title(char(eV));
                        axColor = [0 1 0];
                        hlColor = [1 1 1];
                        vlColor = [1 1 0];
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(1) = iP;
                        % ROI
                        obj.DrawROIs(eV);
                    case 3
                        %  Sagital
                        title(char(eV));
                        axColor = [1 1 0];
                        hlColor = [1 1 1];
                        vlColor = [0 1 0];
                        XXres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(2) = iP;
                        % ROI
                        obj.DrawROIs(eV);
                end
                % Axis
                obj.SetLabels(eV);
                obj.SetAxis(axColor);
                axis(ax);
            end
            % Planar
            bSaveFig = 1;
            if 0
                subplot('224');
                V = obj.vsTemporary(iC).data;
                s = size(V);
                r = obj.vsTemporary.res_mm;
                j=1;lx=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                j=2;ly=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                j=3;lz=linspace(s(j)*r(j)/2, -s(j)*r(j)/2, s(j));
                [x,y,z] = meshgrid(ly,lx,lz);
                vP(1) = (obj.vPlaneCur(1) - round(size(V,1)/2))*r(1);
                vP(2) = (obj.vPlaneCur(2) - round(size(V,2)/2))*r(2);
                vP(3) = (obj.vPlaneCur(3) - round(size(V,3)/2))*r(3);
                if eT == MCC_eType.Segmented
                    V(find(V==0)) = 1;
                    V(find(I==min(V(:)))) = 1;
                    slice(x,y,z,single(V), vP(2),vP(1),vP(3));
                    view([45 -45 45]);
                    colormap(map);
                else
                    slice(x,y,z,single(V), vP(2),vP(1),-vP(3));
                    %view([45 -45 45])
                    colormap gray;
                end
                shading flat;axis image;
                obj.SetAxis;
                xlabel('x [mm]');ylabel('y [mm]');zlabel('z [mm]');
                colorbar;
                % Draw box
                obj.DrawVOI;
                % Title
                vs = obj.vsTemporary; txtRef = '';
                title([ ...
                    strrep(...
                    MCC_RemoveSeriesNumber(vs(iC).ImagingModeID),'_', ...
                    ' ') '  (' ...
                    num2str(vs(iC).res_mm(1),2) ' x ' ...
                    num2str(vs(iC).res_mm(2),2) ' x ' ...
                    num2str(vs(iC).res_mm(3),2) ' = ' ...
                    num2str(vs(iC).GetVoxelSize_mm3,2) ' mm^3 ' ...
                    ') ' txtRef]);
            end
            % Save figure
            if bSaveFig
                if bAllSlices
                    obj.SaveFigure(gcf, [] , eT, MCC_eView.VOIs, 0, iPsaved);
                else
                    obj.SaveFigure(gcf, [] , eT, MCC_eView.VOIs, 0);
                end
            end
        end
        
        function Display4Q_(obj, ...
                eT, ...
                vShiftPlane, ...
                bAllSlices, bSaveFig)
            % Display RBG image (3 views)
            % NGR 2016 05 09
            
            eT = MCC_eType.Original; % @NGR 2016 06 02
            
            if nargin<=4
                bSaveFig = 0;
            end
            
            % Organ specificities
            switch obj.eOrganCur                
                case MCC_eOrgan.Prostate
                    eT = MCC_eType.Resampled; % @NGR 2016 06 02
                    if obj.IsVolSeqDataEmpty(eT)
                        obj.Resample('t2');
                    end
                otherwise
            end
            % Data
            obj.LoadData;
            if obj.IsVolSeqDataEmpty(eT)
                disp(' ');
                disp('No processed data !');
                disp(' ');
                return;
            end
            eV = MCC_eView.Axial;
            obj.vShiftPlaneCur = vShiftPlane;            
            obj.vsTemporary = obj.GetVolSeq(eT);
            % Number of clusters
            if eT == MCC_eType.Segmented
                iC = length(obj.vsTemporary);
            else
                iC = 1;
            end
            % Volume Of Interest (VOI)
            if isempty(obj.VOI)
                s = size(obj.vsTemporary(iC).data);
                r_mm = obj.vsTemporary(iC).res_mm;
                s_mm = s .* r_mm;
                o_mm = -[s_mm(2) s_mm(1) s_mm(3)]/2 + 5;
                obj.VOI = [o_mm 50 50 50];
            end
            % Display
            if isempty(obj.MaxAxis)
                ax = obj.GetMaxAxis(MCC_eType.Original, eV);
                obj.MaxAxis= ax;
            else
                ax = obj.MaxAxis;
            end
            for  k =1:3
                switch k
                    case 1
                        eV = MCC_eView.Axial;
                        subplot('221');
                    case 2
                        eV = MCC_eView.Coronal;
                        subplot('223');
                    case 3
                        eV = MCC_eView.Sagittal;
                        subplot('222');
                end
                % Edge outlining
                if 0
                    obj.vsTemporary = obj.Edge(obj.vsTemporary);
                    obj.vsTemporary = obj.MultiplyByVol(...
                        obj.GetVolSeq(eT), ...
                        uint16(not(obj.vsTemporary(1).data)));
                    obj.vsTemporary(1).data(find(obj.vsTemporary(1).data==0))=1;
                end
                % Plane index
                iPlane = obj.GetCentralOrthoPlaneIndex(...
                    MCC_eType.Temporary, eV, 1);
                switch k
                    case 1
                        iP = iPlane+vShiftPlane(3);
                        iPsaved = iP;
                    case 2
                        iP = iPlane+vShiftPlane(1);
                        iPsaved = iP;
                    case 3
                        iP = iPlane+vShiftPlane(2);
                        iPsaved = iP;
                end
                % Image
                [I,RI]=obj.GetOrthoPlane(...
                    MCC_eType.Temporary, eV, 1, iP);
                % Colormap
                if eT == MCC_eType.Segmented
                    % Erosion of clusters
                    if 0
                        obj.vsTemporary = obj.ErodeClusters(obj.vsTemporary);
                    end
                    % Slice data
                    V = obj.vsTemporary(iC).data;
                    [I,RI]=obj.GetOrthoPlane(...
                        MCC_eType.Temporary, eV, 1, iP);
                    % Color mapping
                    Nc = max(V(:)); % nb of clusters
                    map = obj.Colormap('jet',Nc);
                    I(find(I==0)) = 1;
                    I(find(I==min(I(:))))=1;
                    % RGB
                    Irgb = map(round(I(:)),:);
                    Irgb = reshape(Irgb,[size(I) 3]);
                    % Display
                    imshow(Irgb, RI);
                else
                    %I = histeq(I);
                    imshow(I,RI,[min(I(:)) 0.5*max(I(:))]);
                end
                % Axis
                obj.SetAxis;
                axis(obj.GetMaxAxis(eT, eV));
                axis(ax);
                % Title and line
                res_mm= obj.vsTemporary(1).res_mm;
                % Inference settings
                VOI = obj.VOI;
                crop_mm(1,:) = [VOI(1) VOI(1) + VOI(4)];
                crop_mm(2,:) = [VOI(2) VOI(2) + VOI(5)];
                crop_mm(3,:) = [VOI(3) VOI(3) + VOI(6)];
                scoreTh = 0.70;
                switch k
                    case 1
                        %  Axial
                        title(char(eV));
                        axColor = [1 1 1];
                        hlColor = [0 1 0];
                        vlColor = [1 1 0];
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(1);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(3) = iP;
                        % ROI
                        ROI_mm = [ ...
                            obj.VOI(1) obj.VOI(2) ...
                            obj.VOI(4) obj.VOI(5) ...
                            ];
                        if isempty(obj.hROI)
                            obj.hROI = imrect(gca, ROI_mm);
                        else
                            obj.hROI = imrect(gca, ROI_mm);
                        end
                        fcn = makeConstrainToRectFcn(...
                            'imrect',get(gca,'XLim'),get(gca,'YLim'));
                        setPositionConstraintFcn(obj.hROI(k),fcn);
                        % Inference
                        if ~isempty(obj.net)
                            Iroi = obj.vsTemporary(iC).CropImage(...
                                MCC_eView.Axial, crop_mm, iP);
                            iClass = k;
                            res = obj.Infere(Iroi);
                            scores = squeeze(gather(res(end).x));
                            [scoreClassBest, iClassBest] = max(scores);
                            disp('___________________________________________');
                            sprintf( ...
                                '  %s (%1.2f %%) %s (%1.2f %%) %s (%1.2f %%)',...
                                obj.net.meta.classes.name{1},...
                                100*scores(1), ...
                                obj.net.meta.classes.name{2},...
                                100*scores(2), ...
                                obj.net.meta.classes.name{3},...
                                100*scores(3) ...
                                )
                            if iClassBest == iClass && scoreClassBest > scoreTh
                                txtBox = sprintf('  %s (%1.2f %%)',...
                                    obj.net.meta.classes.name{iClassBest},...
                                    100*scoreClassBest)
                                obj.hROI(k).setColor([1 0 0]);
                            else
                                txtBox =  ' ';
                            end
                            text(...
                                ROI_mm(1)+ROI_mm(3),...
                                ROI_mm(2), ...
                                crop_mm(3,2),txtBox,'Color','red');
                        end
                    case 2
                        %  Coronal
                        title(char(eV));
                        axColor = [0 1 0];
                        hlColor = [1 1 1];
                        vlColor = [1 1 0];
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(1) = iP;
                        % VOI
                        ROI_mm = [ ...
                            obj.VOI(1) obj.VOI(3) ...
                            obj.VOI(4) obj.VOI(6) ...
                            ];
                        obj.hROI(k) = imrect(gca, ROI_mm);
                        fcn = makeConstrainToRectFcn(...
                            'imrect',get(gca,'XLim'),get(gca,'YLim'));
                        setPositionConstraintFcn(obj.hROI(k),fcn);
                        % Inference
                        if ~isempty(obj.net)
                            Iroi = obj.vsTemporary(iC).CropImage(...
                                MCC_eView.Coronal, crop_mm, iP);
                            iClass = k; % k
                            res = obj.Infere(Iroi);
                            scores = squeeze(gather(res(end).x));
                            [scoreClassBest, iClassBest] = max(scores);
                            sprintf( ...
                                '  %s (%1.2f %%) %s (%1.2f %%) %s (%1.2f %%)',...
                                obj.net.meta.classes.name{1},...
                                100*scores(1), ...
                                obj.net.meta.classes.name{2},...
                                100*scores(2), ...
                                obj.net.meta.classes.name{3},...
                                100*scores(3) ...
                                )
                            if iClassBest == iClass && scoreClassBest > scoreTh
                                txtBox = sprintf('  %s (%1.2f %%)',...
                                    obj.net.meta.classes.name{iClassBest},...
                                    100*scoreClassBest)
                                obj.hROI(k).setColor([1 0 0]);
                            else
                                txtBox =  ' ';
                            end
                            text(...
                                ROI_mm(1)+ROI_mm(3),...
                                ROI_mm(2), ...
                                crop_mm(3,2),txtBox,'Color','red');
                        end
                    case 3
                        %  Sagital
                        title(char(eV));
                        axColor = [1 1 0];
                        hlColor = [1 1 1];
                        vlColor = [0 1 0];
                        XXres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(2) = iP;
                        % VOI
                        ROI_mm = [ ...
                            obj.VOI(2) obj.VOI(3) ...
                            obj.VOI(5) obj.VOI(6) ...
                            ];
                        obj.hROI(k) = imrect(gca, ROI_mm);
                        fcn = makeConstrainToRectFcn(...
                            'imrect',get(gca,'XLim'),get(gca,'YLim'));
                        setPositionConstraintFcn(obj.hROI(k),fcn);
                        % Inference
                        if ~isempty(obj.net)
                            Iroi = obj.vsTemporary(iC).CropImage(...
                                MCC_eView.Sagittal, crop_mm, iP);
                            iClass = k; % k
                            res = obj.Infere(Iroi);
                            scores = squeeze(gather(res(end).x));
                            [scoreClassBest, iClassBest] = max(scores);
                            sprintf( ...
                                '  %s (%1.2f %%) %s (%1.2f %%) %s (%1.2f %%)',...
                                obj.net.meta.classes.name{1},...
                                100*scores(1), ...
                                obj.net.meta.classes.name{2},...
                                100*scores(2), ...
                                obj.net.meta.classes.name{3},...
                                100*scores(3) ...
                                )
                            if iClassBest == iClass && scoreClassBest > scoreTh
                                txtBox = sprintf('  %s (%1.2f %%)',...
                                    obj.net.meta.classes.name{iClassBest},...
                                    100*scoreClassBest)
                                obj.hROI(k).setColor([1 0 0]);
                            else
                                txtBox =  ' ';
                            end
                            text(...
                                ROI_mm(1)+ROI_mm(3),...
                                ROI_mm(2), ...
                                crop_mm(3,2),txtBox,'Color','red');
                        end
                end
                obj.SetLabels(eV);
                % Axis
                obj.SetAxis(axColor);
                axis(ax);
            end
            % Planar
            if 0
                subplot('224');
                V = obj.vsTemporary(iC).data;
                s = size(V);
                r = obj.vsTemporary.res_mm;
                j=1;lx=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                j=2;ly=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                j=3;lz=linspace(s(j)*r(j)/2, -s(j)*r(j)/2, s(j));
                [x,y,z] = meshgrid(ly,lx,lz);
                vP(1) = (obj.vPlaneCur(1) - round(size(V,1)/2))*r(1);
                vP(2) = (obj.vPlaneCur(2) - round(size(V,2)/2))*r(2);
                vP(3) = (obj.vPlaneCur(3) - round(size(V,3)/2))*r(3);
                if eT == MCC_eType.Segmented
                    V(find(V==0)) = 1;
                    V(find(I==min(V(:)))) = 1;
                    slice(x,y,z,single(V), vP(2),vP(1),vP(3));
                    view([45 -45 45]);
                    colormap(map);
                else
                    slice(x,y,z,single(V), vP(2),vP(1),-vP(3));
                    view([45 -45 45])
                    colormap gray;
                end
                shading flat;axis image;
                obj.SetAxis;
                xlabel('x [mm]');ylabel('y [mm]');zlabel('z [mm]');
                colorbar;
                if ~isempty(obj.viewpoint)
                    %view(obj.viewpoint);
                end
                % Draw box
                VOImm = obj.VOI;
                VOImm(3) = -VOImm(3);
                VOImm(6) = -VOImm(6);
                obj.Draw3Dbox(VOImm);
            end
            % Volume
            if 0
                V = obj.vsTemporary.data;
                Nc = max(V(:))
                D = (V==Nc);
                D = single(D);
                p1 = patch(isosurface(D),'FaceColor','none',...
                    'EdgeColor','red');
                daspect(...
                    [1,1,obj.vsRegistered(1).res_mm(1)/obj.vsRegistered(1).res_mm(3)]);
                camlight
                lighting gouraud
                isonormals(D,p1);hold on;
            end
            % Title
            vs = obj.vsTemporary; txtRef = '';
            title([ ...
                strrep(MCC_RemoveSeriesNumber(vs(iC).ImagingModeID),'_', ...
                ' ') '  (' ...
                num2str(vs(iC).res_mm(1),2) ' x ' ...
                num2str(vs(iC).res_mm(2),2) ' x ' ...
                num2str(vs(iC).res_mm(3),2) ' = ' ...
                num2str(vs(iC).GetVoxelSize_mm3,2) ' mm^3 ' ...
                ') ' txtRef]);
            % Save figure
            if bSaveFig
                if bAllSlices
                    obj.SaveFigure(gcf, [] ,eT, MCC_eView.All, 0, iPsaved);
                else
                    obj.SaveFigure(gcf, [] ,eT, MCC_eView.All, 0, 0);
                end
            end
            obj.SaveFigure(gcf, [] ,eT, MCC_eView.All, 0, 0);
        end
        
        function map = Colormap(obj, mapType, N)
            % RGB mixing model
            % NGR 2016 05 20
            switch mapType
                case 'rgb'
                    map = zeros(1,3);
                    map(1,1:3) = [ 0 0 0 ];     % Black
                    map(2,1:3) = [ 0 0 1 ];     % Bleu
                    map(3,1:3) = [ 0 1 1 ];     % Cyan
                    map(4,1:3) = [ 0 1 0 ];     % Green
                    map(5,1:3) = [ 1 1 0 ];     % Yellow
                    map(6,1:3) = [ 1 0 0 ];     % Red
                    map(7,1:3) = [ 1 0 1 ];     % Magenta
                    map(8,1:3) = [ 1 1 1 ];     % White
                case 'jet'
                    %map = jet(N-1);
                    %map = [0 0 0 ; map];
                    
                    N = N-1;
                    mapJ = jet(256);
                    idx = round((((1:N)/N)*256));
                    map = mapJ(idx,:);
                    map = [0 0 0 ; map];
                    
                case 'hot'
                    map = hot(N);
                otherwise
                    map = gray(N);
            end
        end
        
        function Display4Q_Localization(obj, ...
                eT, ...
                vShiftPlane, ...
                bAllSlices, bSaveFig)
            % Display RBG image (3 views)
            % NGR 2016 05 09
            
            if nargin<=4
                bSaveFig = 0;
            end
            
            obj.LoadData;
            if obj.IsVolSeqDataEmpty(eT)
                disp(' ');
                disp('No processed data !');
                disp(' ');
                return;
            end
            eV = MCC_eView.Axial;
            obj.vShiftPlaneCur = vShiftPlane;
            % Data
            obj.vsTemporary = obj.GetVolSeq(eT);
            
            % Number of clusters
            if eT == MCC_eType.Segmented
                iC = length(obj.vsTemporary);
            else
                iC = 1;
            end
            
            % Volume Of Interest (VOI)
            if isempty(obj.VOI)
                s = size(obj.vsTemporary(iC).data);
                r_mm = obj.vsTemporary(iC).res_mm;
                s_mm = s .* r_mm;
                o_mm = -[s_mm(2) s_mm(1) s_mm(3)]/2 + 5;
                obj.VOI = [o_mm 50 50 50];
            end
            % Display
            if isempty(obj.MaxAxis)
                ax = obj.GetMaxAxis(MCC_eType.Original, eV);
                obj.MaxAxis= ax;
            else
                ax = obj.MaxAxis;
            end
            for  k =1:3
                switch k
                    case 1
                        eV = MCC_eView.Axial;
                        subplot('221');
                    case 2
                        eV = MCC_eView.Coronal;
                        subplot('223');
                    case 3
                        eV = MCC_eView.Sagittal;
                        subplot('222');
                end
                % Edge outlining
                if 0
                    obj.vsTemporary = obj.Edge(obj.vsTemporary);
                    obj.vsTemporary = obj.MultiplyByVol(...
                        obj.GetVolSeq(eT), ...
                        uint16(not(obj.vsTemporary(1).data)));
                    obj.vsTemporary(1).data(find(obj.vsTemporary(1).data==0))=1;
                end
                % Plane index
                iPlane = obj.GetCentralOrthoPlaneIndex(...
                    MCC_eType.Temporary, eV, 1);
                switch k
                    case 1
                        iP = iPlane+vShiftPlane(3);
                        iPsaved = iP;
                    case 2
                        iP = iPlane+vShiftPlane(1);
                        iPsaved = iP;
                    case 3
                        iP = iPlane+vShiftPlane(2);
                        iPsaved = iP;
                end
                % Image
                [I,RI]=obj.GetOrthoPlane(...
                    MCC_eType.Temporary, eV, 1, iP);
                % Display
                imshow(I,RI,[min(I(:)) 0.5*max(I(:))]);
                % Axis
                obj.SetAxis;
                %axis(obj.GetMaxAxis(eT, eV));
                axis(ax);
                % Title and line
                res_mm= obj.vsTemporary(1).res_mm;
                % Inference settings
                VOI = obj.VOI;
                crop_mm(1,:) = [VOI(1) VOI(1) + VOI(4)];
                crop_mm(2,:) = [VOI(2) VOI(2) + VOI(5)];
                crop_mm(3,:) = [VOI(3) VOI(3) + VOI(6)];
                scoreTh = 0.70;
                switch k
                    case 1
                        %  Axial
                        title(char(eV));
                        axColor = [1 1 1];
                        hlColor = [0 1 0];
                        vlColor = [1 1 0];
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(1);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(3) = iP;
                        % ROI
                        ROI_mm = [ ...
                            obj.VOI(1) obj.VOI(2) ...
                            obj.VOI(4) obj.VOI(5) ...
                            ];
                        obj.Draw2Dbox(ROI_mm,[1 1 1]);
                        % VOI localization
                        [~, vid]  = sort(obj.VOIsl_mm.score);
                        [~, iMax] = max(obj.VOIsl_mm.score);
                        for iVOI = vid
                            ROIl_mm = [ ...
                                obj.VOIsl_mm.data(iVOI,1) obj.VOIsl_mm.data(iVOI,2) ...
                                obj.VOIsl_mm.data(iVOI,4) obj.VOIsl_mm.data(iVOI,5) ...
                                ];
                            c = [obj.VOIsl_mm.score(iVOI)/max(obj.VOIsl_mm.score) 0 1];
                            c = [0 0 1];
                            obj.Draw2Dbox(ROIl_mm,c);
                        end
                        iVOI = iMax;
                        ROIl_mm = [ ...
                            obj.VOIsl_mm.data(iVOI,1) obj.VOIsl_mm.data(iVOI,2) ...
                            obj.VOIsl_mm.data(iVOI,4) obj.VOIsl_mm.data(iVOI,5) ...
                            ];
                        obj.Draw2Dbox(ROIl_mm,[1 0 1]);
                    case 2
                        %  Coronal
                        title(char(eV));
                        axColor = [0 1 0];
                        hlColor = [1 1 1];
                        vlColor = [1 1 0];
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(1) = iP;
                        % VOI
                        ROI_mm = [ ...
                            obj.VOI(1) obj.VOI(3) ...
                            obj.VOI(4) obj.VOI(6) ...
                            ];
                        obj.Draw2Dbox(ROI_mm,[1 1 1]);
                        % VOI localization
                        [~, vid] = sort(obj.VOIsl_mm.score);
                        for iVOI = vid
                            
                            ROIl_mm = [ ...
                                obj.VOIsl_mm.data(iVOI,1) obj.VOIsl_mm.data(iVOI,3) ...
                                obj.VOIsl_mm.data(iVOI,4) obj.VOIsl_mm.data(iVOI,6) ...
                                ];
                            c = [obj.VOIsl_mm.score(iVOI)/max(obj.VOIsl_mm.score) 0 1];
                            c = [0 0 1];
                            obj.Draw2Dbox(ROIl_mm,c);
                        end
                        iVOI = iMax;
                        ROIl_mm = [ ...
                            obj.VOIsl_mm.data(iVOI,1) obj.VOIsl_mm.data(iVOI,3) ...
                            obj.VOIsl_mm.data(iVOI,4) obj.VOIsl_mm.data(iVOI,6) ...
                            ];
                        obj.Draw2Dbox(ROIl_mm,[1 0 1]);
                    case 3
                        %  Sagital
                        title(char(eV));
                        axColor = [1 1 0];
                        hlColor = [1 1 1];
                        vlColor = [0 1 0];
                        XXres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                        obj.vPlaneCur(2) = iP;
                        % VOI
                        ROI_mm = [ ...
                            obj.VOI(2) obj.VOI(3) ...
                            obj.VOI(5) obj.VOI(6) ...
                            ];
                        obj.Draw2Dbox(ROI_mm,[1 1 1]);
                        % VOI localization
                        [~, vid] = sort(obj.VOIsl_mm.score);
                        for iVOI = vid
                            ROIl_mm = [ ...
                                obj.VOIsl_mm.data(iVOI,2) obj.VOIsl_mm.data(iVOI,3) ...
                                obj.VOIsl_mm.data(iVOI,5) obj.VOIsl_mm.data(iVOI,6) ...
                                ];
                            c = [obj.VOIsl_mm.score(iVOI)/max(obj.VOIsl_mm.score) 0 1];
                            c = [0 0 1];
                            obj.Draw2Dbox(ROIl_mm,c);
                        end
                        iVOI = iMax;
                        ROIl_mm = [ ...
                            obj.VOIsl_mm.data(iVOI,2) obj.VOIsl_mm.data(iVOI,3) ...
                            obj.VOIsl_mm.data(iVOI,5) obj.VOIsl_mm.data(iVOI,6) ...
                            ];
                        obj.Draw2Dbox(ROIl_mm,[1 0 1]);
                end
                obj.SetLabels(eV);
                % Axis
                obj.SetAxis(axColor);
                axis(ax);
            end
            % Planar
            if 1
                subplot('224');
                V = obj.vsTemporary(iC).data;
                s = size(V);
                r = obj.vsTemporary.res_mm;
                j=1;lx=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                j=2;ly=linspace(-s(j)*r(j)/2, s(j)*r(j)/2, s(j));
                j=3;lz=linspace(s(j)*r(j)/2, -s(j)*r(j)/2, s(j));
                [x,y,z] = meshgrid(ly,lx,lz);
                vP(1) = (obj.vPlaneCur(1) - round(size(V,1)/2))*r(1);
                vP(2) = (obj.vPlaneCur(2) - round(size(V,2)/2))*r(2);
                vP(3) = (obj.vPlaneCur(3) - round(size(V,3)/2))*r(3);
                if eT == MCC_eType.Segmented
                    V(V==0) = 1;
                    V(find(I==min(V(:)))) = 1;
                    slice(x,y,z,single(V), vP(2),vP(1),vP(3));
                    view([45 -45 45]);
                    colormap(map);
                else
                    slice(x,y,z,single(V), vP(2),vP(1),-vP(3));
                    view([45 -45 45])
                    colormap gray;
                end
                shading flat;axis image;
                obj.SetAxis;
                axis([-250 250 -250 250 -250 250])
                xlabel('x [mm]');ylabel('y [mm]');zlabel('z [mm]');
                colorbar;
                if ~isempty(obj.viewpoint)
                    %view(obj.viewpoint);
                end
                % Draw box
                VOImm = obj.VOI;
                VOImm(3) = -VOImm(3);
                VOImm(6) = -VOImm(6);
                obj.Draw3Dbox(VOImm, [1 1 1]);
                % Draw box (localization)
                [~, vid] = sort(obj.VOIsl_mm.score);
                for iVOI = vid,
                    VOIlmm = obj.VOIsl_mm.data(iVOI,:);
                    VOIlmm(3) = -VOIlmm(3);
                    VOIlmm(6) = -VOIlmm(6);
                    c = [obj.VOIsl_mm.score(iVOI)/max(obj.VOIsl_mm.score) 0 1];
                    c = [0 0 1];
                    obj.Draw3Dbox(VOIlmm,c);
                    if iVOI == iMax
                        obj.Draw3Dbox(VOIlmm,[1 0 1]);
                    end
                end
                
            end
            % Title
            vs = obj.vsTemporary; txtRef = '';
            title([ ...
                strrep(MCC_RemoveSeriesNumber(vs(iC).ImagingModeID),'_', ...
                ' ') '  (' ...
                num2str(vs(iC).res_mm(1),2) ' x ' ...
                num2str(vs(iC).res_mm(2),2) ' x ' ...
                num2str(vs(iC).res_mm(3),2) ' = ' ...
                num2str(vs(iC).GetVoxelSize_mm3,2) ' mm^3 ' ...
                ') ' txtRef]);
            % Save figure
            bSaveFig = 1;
            if bSaveFig
                if bAllSlices
                    obj.SaveFigure(gcf, [] ,eT, MCC_eView.All, 0, iPsaved);
                else
                    obj.SaveFigure(gcf, [] ,eT, MCC_eView.All, 0, 0);
                end
            end
        end
        
        function Compare(obj)
            % Display registration for visual comparison
            % NGR 2016 03 16
            
            obj.LoadData;
            % Volume loop
            for iVol = 1:obj.GetNbVol
                if iVol~=obj.iVolRef_Registered
                    ext = obj.GetVolSeqLabel(obj.eTypeCur);
                    eval(['N = obj.vs' ext '(iVol).GetNbSlices;']);                    
                    R = floor(N/2)-1;
                    if 1
                        % Central sclice
                        obj.ComparePlane(iVol,obj.vShiftPlaneCur,0);
                        %close;
                    else
                        % All sclices
                        for j=-R:R
                            obj.ComparePlane(iVol,j,1);
                            close;
                        end
                    end
                end
            end
        end
        
        function ComparePlane(obj, iVol, vShiftPlane, bAllSlices)
            % Display registration for visual comparison
            % NGR 2016 03 10
            bEdge = 0;
            thSat1 = inf;
            thSat2 = inf;
            eV = MCC_eView.Axial;
            % Display
            j = 1;
            obj.Figure(obj.eTypeCur, MCC_eView.All);
            ax = obj.GetMaxAxis(MCC_eType.Original, MCC_eView.Axial);
            res_mm = obj.vsResampled(1).res_mm;
            for  k =1:3
                switch k
                    case 1
                        eV = MCC_eView.Axial;
                        iSI = 3;
                    case 2
                        eV = MCC_eView.Coronal;
                        iSI = 1;
                    case 3
                        eV = MCC_eView.Sagittal;
                        iSI = 2;
                end
                %--------------%
                % Unregistered %
                %--------------%
                subplot(['23' num2str(j)]); j = j + 1;
                eTun = MCC_eType.Resampled;
                bSliceLoc = 0;
                % ROI, NGR 2016 09 13   
                switch obj.eTypeCur
                    case MCC_eType.Registered_VOI1
                        r = obj.vsRegistered_VOI1(2).res_mm;
                        iVOI = 1;
                        vShiftPlane(2) = round((obj.VOI_mm(iVOI,1) + obj.VOI_mm(iVOI,4)/2)/r(1));
                        vShiftPlane(1) = round((obj.VOI_mm(iVOI,2) + obj.VOI_mm(iVOI,5)/2)/r(2));
                        vShiftPlane(3) = round((obj.VOI_mm(iVOI,3) + obj.VOI_mm(iVOI,6)/2)/r(3));
                    case MCC_eType.Registered_VOI2
                        r = obj.vsRegistered_VOI2(2).res_mm;
                        iVOI = 2;
                        vShiftPlane(2) = round((obj.VOI_mm(iVOI,1) + obj.VOI_mm(iVOI,4)/2)/r(1));
                        vShiftPlane(1) = round((obj.VOI_mm(iVOI,2) + obj.VOI_mm(iVOI,5)/2)/r(2));
                        vShiftPlane(3) = round((obj.VOI_mm(iVOI,3) + obj.VOI_mm(iVOI,6)/2)/r(3));
                end
                % Fixed
                iPlane = obj.GetCentralOrthoPlaneIndex( ....
                    eTun, eV, obj.iVolRef_Registered, bSliceLoc);
                iP = iPlane+vShiftPlane(iSI);
                [I1, RI1] = obj.GetOrthoPlane(eTun,eV,obj.iVolRef_Registered,iP);
                % Moving            
                iPlane = obj.GetCentralOrthoPlaneIndex(...
                    eTun, eV, iVol, bSliceLoc);
                iP = iPlane+vShiftPlane(iSI);
                if (iP<1)
                    [I2, RI2] = obj.GetOrthoPlane(...
                        eTun, eV, iVol, iPlane, bSliceLoc);
                    I2 = zeros(size(I2));
                else
                    [I2, RI2] = obj.GetOrthoPlane(...
                        eTun, eV, iVol, iP, bSliceLoc);
                end
                % Edges
                if bEdge
                    I1 = edge(I1,'Canny');I2 = edge(I2,'Canny');
                end
                % Saturation
                if thSat1
                    I1(find(I1>thSat1)) = thSat1;
                end
                if thSat2
                    I2(find(I2>thSat2)) = thSat2;
                end
                % Show pair
                imshowpair(...
                    I1, RI1, I2, RI2,...
                    'Scaling','independent',...
                    'ColorChannels','red-cyan');
                % Cursor location
                YYres = [vShiftPlane(1) vShiftPlane(1)]*obj.vsResampled(iVol).res_mm(3);
                % Title and line
                switch k
                    case 1                        
                        title(char(eV)); %title('Axial');
                        axColor = [1 1 1];
                        hlColor = [0 1 0];
                        vlColor = [1 1 0];
                        % Lines
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(1);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);                  
                    case 2
                        title(char(eV));
                        axColor = [0 1 0];
                        hlColor = [1 1 1];
                        vlColor = [1 1 0];
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                    case 3
                        title(char(eV));
                        axColor = [1 1 0];
                        hlColor = [1 1 1];
                        vlColor = [0 1 0];
                        XXres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                end
                % ROI
                switch obj.eTypeCur
                    case MCC_eType.Registered_VOI1
                        obj.DrawROI(eV,1);
                    case MCC_eType.Registered_VOI2
                        obj.DrawROI(eV,2);
                end
                obj.SetLabels(eV);
                % Axis
                obj.SetAxis(axColor);
                axis(ax);                
                %------------%
                % Registered %
                %------------%
                subplot(['23' num2str(j+2)]);
                eT = obj.eTypeCur;
                switch eT
                    case MCC_eType.Registered
                        iVolRef = obj.iVolRef_Registered;
                    case MCC_eType.Registered_Elastic
                        iVolRef = obj.iVolRef_Registered_Elastic;                        
                    case MCC_eType.Localized                        
                        iVolRef = obj.iVolRef_Resampled;
                    case MCC_eType.Registered_VOI1
                        iVolRef = obj.iVolRef_Registered_VOI1;
                    case MCC_eType.Registered_VOI2
                        iVolRef = obj.iVolRef_Registered_VOI2; 
                    case MCC_eType.Registered_VOI3
                        iVolRef = obj.iVolRef_Registered_VOI3; 
                    case MCC_eType.Registered_VOI4
                        iVolRef = obj.iVolRef_Registered_VOI4;
                    case MCC_eType.Registered_VOI5
                        iVolRef = obj.iVolRef_Registered_VOI5;
                    case MCC_eType.Registered_VOI6
                        iVolRef = obj.iVolRef_Registered_VOI6;
                    case MCC_eType.Registered_VOI7
                        iVolRef = obj.iVolRef_Registered_VOI7;
                    case MCC_eType.Registered_VOI8
                        iVolRef = obj.iVolRef_Registered_VOI8;
                    case MCC_eType.Registered_Date
                        iVolRef = obj.iVolRef_Registered_Date;  % NGR 2018 03 21                        
                end
                iPlane = obj.GetCentralOrthoPlaneIndex(...
                    eT, eV, iVolRef);
                iP = iPlane + vShiftPlane(iSI);
                % Fixed                                
                [I1, RI1] = obj.GetOrthoPlane(eT, eV, iVolRef, iP);
                % Moving
                [I2, RI2] = obj.GetOrthoPlane(eT, eV, iVol, iP);
                % Edges
                if bEdge
                    I1 = edge(I1,'Canny');I2 = edge(I2,'Canny');
                end
                % Saturation
                if thSat1
                    I1(find(I1>thSat1)) = thSat1;
                end
                if thSat2
                    I2(find(I2>thSat2)) = thSat2;
                end
                % Show pair
                imshowpair(...
                    I1, RI1, I2, RI2,...
                    'Scaling','independent',...
                    'ColorChannels','red-cyan');
                xlabel('x [mm]');ylabel('y [mm]');
                % Cursor location
                YYreg = [vShiftPlane(1) vShiftPlane(1)]*obj.vsResampled(iVol).res_mm(3);
                % Title and line
                switch k
                    case 1
                        str1 = strrep(MCC_RemoveSeriesNumber(obj.vsOriginal(iVolRef).ImagingModeID),'_',' ');
                        str2 = strrep(MCC_RemoveSeriesNumber(obj.vsOriginal(iVol).ImagingModeID),'_',' ');
                        title([str1 ' [fixed] - '  str2 ' [moving]' ]);
                        % Lines
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(1);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                    case 2
                        XXres = [vShiftPlane(2) vShiftPlane(2)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                    case 3
                        XXres = [vShiftPlane(1) vShiftPlane(1)] ...
                            * res_mm(2);
                        YYres = [vShiftPlane(3) vShiftPlane(3)] ...
                            * res_mm(3);
                        line(XXres,[ax(3) ax(4)],'Color',vlColor);
                        line([ax(1) ax(2)],YYres,'Color',hlColor);
                end
                % ROI
                switch obj.eTypeCur
                    case MCC_eType.Registered_VOI1
                        obj.DrawROI(eV,1);
                    case MCC_eType.Registered_VOI2
                        obj.DrawROI(eV,2);
                end
                obj.SetLabels(eV);
                % Axis
                obj.SetAxis(axColor);
                axis(ax);
            end
            % Save figure
            iPlane = obj.GetCentralOrthoPlaneIndex(...
                eT, MCC_eView.Axial, iVolRef);
            iP = iPlane+vShiftPlane(1);
            if bAllSlices
                obj.SaveFigure(gcf, 'Slices', eT, MCC_eView.All, iVol, iP);
            else
                obj.SaveFigure(gcf, [], eT, MCC_eView.All, iVol, 0);
            end
        end
        
        function fig = Figure(obj, eType, eView)
            % Create figure
            % NGR 2016 03 16
            fig = figure(...
                'units','normalized',...
                'outerposition',[0 0 1 1],...
                'Color',[0 0 0], ...
                'WindowScrollWheelFcn', @MCC_WindowScrollWheelFcn, ...
                'WindowButtonDownFcn', @MCC_WindowButtonDownFcn, ...
                'WindowButtonUpFcn', @MCC_WindowButtonUpFcn ...
                );
            % Set global varaiable (for slice navigation)
            global eViewAxis
            eViewAxis = MCC_eView.Axial;
            MCC_setGlobalx(obj);
            % Annotation
            annotation(fig,'textbox',...
                [0 0.88 0.17 0.12],...
                'String',...
                { ...
                ['           Co.: ' strrep(obj.CohortID,'_', ' ')],...
                ['           PID: ' strrep(obj.PatientID,'_', ' ')],...
                ['Study date: ' obj.name], ...
                ['         Type: ' strrep(eType.char,'_',' ')], ...
                ['         View: ' eView.char] ...
                },...
                'FitBoxToText','off',...
                'BackgroundColor',[0 0 0], ...
                'Color', [238, 130, 238] / 255, ... % Violet (web)
                'FontSize',14 ...
                );
            annotation(fig,'textbox',...
                [0.842940301120441 0.0332839174765179 0.525210083564018 0.0303522262078559],...
                'String', ...
                { ...
                ['Radiomics language (' MCC_GetVersionFormated ')'],...
                MCC_GetInstituteSignature,...
                }, ...
                'FitBoxToText','off', ...
                'Color',[0 1 1], ...
                'FontSize',16  ...
                );
        end
        
        function SetAxis(obj, RGB)
            % Set axis (white)
            % NGR 2016 03 15
            if nargin<=1
                RGB = [ 1 1 1];
            end
            h = gca;
            set(h,...
                'XColor', RGB, ...
                'YColor', RGB, ...
                'ZColor', RGB, ...
                'LineWidth', 1, ...
                'FontSize', 14, ...
                'Color', [ 0 0 0] ...
                );
            h.Title.Color = RGB;
            grid on;
        end
        
        function SetLabels(obj, eView)
            % Set labels
            % NGR 2016 04 06
            switch eView
                case MCC_eView.Planar
                    xlabel('x [mm]');ylabel('y [mm]');
                    zlabel('z [mm]');
                case MCC_eView.Axial
                    xlabel('x [mm]');ylabel('y [mm]');
                case MCC_eView.Coronal
                    xlabel('x [mm]');ylabel('z [mm]');
                case MCC_eView.Sagittal
                    xlabel('y [mm]');ylabel('z [mm]');
            end
        end
        
        function SaveFigure(obj, fig, subpath, eType, eView, iVol, iPlane)
            % Capture figure and save as a JPG
            % NGR 2016 03 15
            if ~isempty(obj.root)
                capt = getframe(fig);
                fld = [obj.root '\Processed'];
                if ~isempty(subpath)
                    fld = [fld '\' subpath];
                end                
                if ~isdir(fld)
                    mkdir(fld);
                else
                    % delete([fld '\*.*']); %@
                end
                fmt = 'png';                
                if exist('iPlane')
                    str = obj.MakeFileName(eType, eView, iVol , iPlane);
                else
                    str = obj.MakeFileName(eType, eView, iVol);
                end                
                fn = [ ...
                    fld '\' ...
                    str ...
                    '.' fmt];
                imwrite(...
                    capt.cdata, fn , fmt, 'Comment',MCC_GetArtistSignature);
            end
        end
        
        function r = MakeFileName(obj, eType, eView, iVol, iPlane)
            % Make customized file name
            % NGR 2016 03 16
            if exist('iVol')
                if iVol
                    ImMod = ...
                        MCC_RemoveSeriesNumber(obj.vsOriginal(iVol).ImagingModeID);
                else
                    ImMod = 'All';
                end
            else
                ImMod = '';
            end
            if exist('iPlane')
                if iPlane<10
                    sPlane = ['_00' num2str(iPlane)];
                end
                if iPlane>=10 && iPlane<100
                    sPlane = ['_0' num2str(iPlane)];
                end
                if iPlane>=100 && iPlane<1000
                    sPlane = ['_' num2str(iPlane)];
                end
            else
                sPlane = '';
            end
            r{1} = [ ...
                MCC_GetDateFileFormat() '_' ...
                'v' num2str(MCC_GetVersion) '_' ...
                eType.char '_' ...
                eView.char '_' ...
                obj.PatientID '_' ...
                obj.StudyDateID '_' ...
                ImMod, ...
                sPlane ... %@NGR 2016 08 24
                ];
            r = strjoin(r,'');
        end
        
        function DisplayMontage(obj, iVol, eType)
            % Display multiple image frames as rectangular montage
            % NGR 2016 02 15
            if (nargin<=1)
                eType = obj.eTypeCur;
                iVol = obj.GetVolRefIndex(eType);
            end
            if (nargin==2)
                eType = obj.obj.eTypeCur;
            end
            figure('units','normalized','outerposition',[0 0 1 1]);
            vs = obj.GetVolSeq(eType);
            Nz = vs(iVol).GetNbSlices;
            s = [ ...
                size(vs(iVol).data,1) ...
                size(vs(iVol).data,2) ...
                1 ...
                size(vs(iVol).data,3)];
            V = zeros(s,'uint16');
            for i=1:Nz
                V(:,:,1,i) = vs(iVol).data(:,:,i);
            end
            montage(V, 'DisplayRange', [0 max(V(:))]);
            title(vs(iVol).content);
        end
        
        %-----------------------------------------------------------------%
        % Workspace
        %-----------------------------------------------------------------%
        
        function LoadWorkspace(obj)
            % Load workspace
            % NGR 2016 05 19
            fld = [obj.root '\Processed'];
            if isdir(fld)
                fid = fopen([fld '\Workspace.txt'],'r');
                if fid~=-1
                    formatSpec = '%c';
                    A = fscanf(fid,formatSpec);
                    B = strsplit(A,'\n');
                    for i=1:length(B)
                        eval(B{i});
                    end
                    fclose(fid);
                end
            end
            % Machine
            if exist('cnn.mat', 'file')
                obj.net = load('cnn.mat');
                obj.net = vl_simplenn_move(obj.net, 'gpu');
            end
            % Reset shift
            obj.vShiftPlaneCur = [0 0 0];
            % VOIs
            % Backward compatibility, NGR 2016 09 09 v195
            if isempty(obj.VOI_mm) &&  ~isempty(obj.VOI)
                for i =1:length(obj.VOI)
                    obj.VOI_mm(i,:)= obj.VOI(i).xyzDxDyDz_mm;
                end
            end        
        end
        
        function SaveWorkspace(obj)
            % Save workspace
            % NGR 2016 05 19
            fld = [obj.root '\Processed'];
            if ~isdir(fld)
                mkdir(fld);
            end
            fid = fopen([fld '\Workspace.txt'],'wt+');
            % Organ
            str=['obj.eOrganCur = MCC_eOrgan.' char(obj.eOrganCur) ';'];
            fprintf(fid,'%s\n',str);
            % VOIs            
            for i=1:length(obj.VOI)
                str = [...
                    'obj.VOI(' num2str(i) ').label = ''' ...
                    obj.VOI(i).label ''';'];
                fprintf(fid,'%s\n',str);
                str = [ ...
                    'obj.VOI(' num2str(i) ').xyzDxDyDz_mm = ' ...
                    obj.Mat2str(obj.VOI(i).xyzDxDyDz_mm) ];
                fprintf(fid,'%s\n',str);
                str = [...
                    'obj.VOI(' num2str(i) ').index = ' ...
                    num2str(obj.VOI(i).index) ';'];
                fprintf(fid,'%s\n',str);
                str = [...
                    'obj.VOI(' num2str(i) ').RGB = ' ...
                    obj.Mat2str(obj.VOI(i).RGB) ';'];
                fprintf(fid,'%s\n',str);
            end
            % Current plane
            if 0
                str = [...
                    'obj.vPlaneCur = [' num2str(obj.vPlaneCur) '];'];
                fprintf(fid,'%s\n',str);
                % Shift plane
                str = [...
                    'obj.vShiftPlaneCur = [' num2str(obj.vShiftPlaneCur) '];'];
                fprintf(fid,'%s\n',str);
            end
            % Seeds points
            % NGR 2016 09 09 v195
            str = [...
                'obj.pSeeds = ' ...
                obj.Mat2str(obj.pSeeds) ';'];
            fprintf(fid,'%s\n',str);
            % Close
            fclose(fid);
        end
        
        function str = Mat2str(obj, m)
            % Matrix to string conversion
            % NGR 2016 07 22
            N = size(m,1);
            str = [];
            for i = 1:N
                if i<N
                    str = [ str num2str(m(i,:)) '; ' ];
                else
                    str = [ str num2str(m(i,:)) ];
                end
            end
            str=['[ ' str ' ];'];
        end
        
        %-----------------------------------------------------------------%
        % VOI & ROI
        %-----------------------------------------------------------------%
        
        function InitializeVOI(obj)
            % Instanciate VOI object
            % NGR 2016 07 19
            if isempty(obj.VOI)
                % Color
                RGB(1,:) = [1 0 0];
                RGB(2,:) = [0 1 0];
                RGB(3,:) = [0 1 1];
                RGB(4,:) = [1 1 0];
                RGB(5,:) = [1 1 1];
                % Organ specifities
                switch obj.eOrganCur
                    case MCC_eOrgan.Kidney
                        labels = {...
                            'Suspicious','Left kidney', ...
                            'Right kidney','Paravertebral muscle'};
                    case MCC_eOrgan.Prostate
                        labels = {'Prostate'};
                    case MCC_eOrgan.Brain_FMX
                        labels = {'Head', 'Tubes'};
                    case MCC_eOrgan.Body_FMX
                        labels = {'Body', 'Tubes'};
                    otherwise
                        labels = {'VOI 1', 'VOI 2'};
                end
                % Number of VOIs
                N = length(labels);
                % Structure
                obj.VOI = struct('xyzDxDyDz_mm',[],'label',[],'RGB',[]);
                for i =1:N
                    if ~isempty(obj.VOI_mm)
                        obj.VOI(i).xyzDxDyDz_mm = obj.VOI_mm(i,:);
                    end
                    obj.VOI(i).index = i;
                    obj.VOI(i).label = labels{i};
                    obj.VOI(i).RGB = RGB(i,:);
                end
            end
        end
        
        function DrawROIs(obj, eView)
            % Draw editable ROI
            % NGR 2016 07 19
            if ~isempty(obj.VOI)                
                for iROI=1:length(obj.VOI)                    
                    obj.DrawROI(eView, iROI);
                end
            end
        end
        
        function DrawROI(obj, eView, iROI)
            % Draw editable ROI
            % NGR 2016 07 19
            if ~isempty(obj.VOI)
                switch eView
                    case MCC_eView.Axial
                        idx = [1 2 4 5];
                        k = 1;
                    case MCC_eView.Coronal
                        idx = [1 3 4 6];
                        k = 2;
                    case MCC_eView.Sagittal
                        idx = [2 3 5 6];
                        k = 3;
                end
                ROI_mm = obj.VOI(iROI).xyzDxDyDz_mm(idx);
                if  isempty(obj.hROI)
                    obj.hROI.h(k) = imrect(gca, ROI_mm);
                else
                    obj.hROI(iROI).h(k) = imrect(gca, ROI_mm);
                end
                fcn = makeConstrainToRectFcn(...
                    'imrect',get(gca,'XLim'),get(gca,'YLim'));
                %setPositionConstraintFcn(obj.hROI(iROI).h(k),fcn);
                obj.hROI(iROI).h(k).setColor(obj.VOI(iROI).RGB);
                if k==1
                    text(...
                        ROI_mm(1),...
                        ROI_mm(2)-25, ...
                        obj.VOI(iROI).label,'Color',obj.VOI(iROI).RGB);
                end
            end
        end
        
        function DrawVOI(obj)
            % Draw VOI
            % NGR 2016 07 19
            if ~isempty(obj.VOI)
                for i=1:length(obj.VOI)
                    VOImm = obj.VOI(i).xyzDxDyDz_mm(i,:);
                    VOImm(3) = -VOImm(3);
                    VOImm(6) = -VOImm(6);
                    obj.Draw3Dbox(VOImm, obj.VOI(i).RGB);
                end
            end
        end
        
        function Draw3Dbox(obj, VOI_mm, color)
            % Draw 3D box [mm]
            % NGR 2016 05 19
            switch nargin
                case 2
                    color = [1 1 1];
            end
            c = VOI_mm(1:3) + VOI_mm(4:6)/2;
            d = VOI_mm(4:6)/2;
            %  -
            %
            %
            line(...
                [c(1)-d(1) c(1)+d(1)], ...
                [c(2)-d(2) c(2)-d(2)], ...
                [c(3)-d(3) c(3)-d(3)], ...
                'Color',color,'LineWidth',2);
            % |
            %  -
            line(...
                [c(1)-d(1) c(1)-d(1)], ...
                [c(2)-d(2) c(2)+d(2)], ...
                [c(3)-d(3) c(3)-d(3)], ...
                'Color',color,'LineWidth',2);
            %  -
            % |
            %  -
            line(...
                [c(1)-d(1) c(1)+d(1)], ...
                [c(2)+d(2) c(2)+d(2)], ...
                [c(3)-d(3) c(3)-d(3)], ...
                'Color',color,'LineWidth',2);
            %  -
            % | |
            %  -
            line(...
                [c(1)+d(1) c(1)+d(1)], ...
                [c(2)-d(2) c(2)+d(2)], ...
                [c(3)-d(3) c(3)-d(3)], ...
                'Color',color,'LineWidth',2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -
            %
            %
            line(...
                [c(1)-d(1) c(1)+d(1)], ...
                [c(2)-d(2) c(2)-d(2)], ...
                [c(3)+d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
            %
            % |
            %  -
            line(...
                [c(1)-d(1) c(1)-d(1)], ...
                [c(2)-d(2) c(2)+d(2)], ...
                [c(3)+d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
            %  -
            % |
            %  -
            line(...
                [c(1)-d(1) c(1)+d(1)], ...
                [c(2)+d(2) c(2)+d(2)], ...
                [c(3)+d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
            %  -
            % | |
            %  -
            line(...
                [c(1)+d(1) c(1)+d(1)], ...
                [c(2)-d(2) c(2)+d(2)], ...
                [c(3)+d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % |
            line(...
                [c(1)-d(1) c(1)-d(1)], ...
                [c(2)-d(2) c(2)-d(2)], ...
                [c(3)-d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
            % ||
            line(...
                [c(1)+d(1) c(1)+d(1)], ...
                [c(2)-d(2) c(2)-d(2)], ...
                [c(3)-d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
            % |||
            line(...
                [c(1)-d(1) c(1)-d(1)], ...
                [c(2)+d(2) c(2)+d(2)], ...
                [c(3)-d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
            % ||||
            line(...
                [c(1)+d(1) c(1)+d(1)], ...
                [c(2)+d(2) c(2)+d(2)], ...
                [c(3)-d(3) c(3)+d(3)], ...
                'Color',color,'LineWidth',2);
        end
        
        function Draw2Dbox(obj, ROI_mm, color)
            % Draw 2D box [mm] [x y dx dy]
            % NGR 2016 07 07
            switch nargin
                case 2
                    color = [1 1 1];
            end
            c = ROI_mm(1:2) + ROI_mm(3:4)/2;
            d = ROI_mm(3:4)/2;
            %  -
            %
            %
            line(...
                [c(1)-d(1) c(1)+d(1)], ...
                [c(2)-d(2) c(2)-d(2)], ...
                'Color',color,'LineWidth',2);
            % |
            %  -
            line(...
                [c(1)-d(1) c(1)-d(1)], ...
                [c(2)-d(2) c(2)+d(2)], ...
                'Color',color,'LineWidth',2);
            %  -
            % |
            %  -
            line(...
                [c(1)-d(1) c(1)+d(1)], ...
                [c(2)+d(2) c(2)+d(2)], ...
                'Color',color,'LineWidth',2);
            %  -
            % | |
            %  -
            line(...
                [c(1)+d(1) c(1)+d(1)], ...
                [c(2)-d(2) c(2)+d(2)], ...
                'Color',color,'LineWidth',2);
        end
        
        %-----------------------------------------------------------------%
        % Volume
        %-----------------------------------------------------------------%
        
        function ClearVolSeq(obj, eType)
            % Clear a volume sequence from memory
            % NGR 2016 03 04
            if (nargin<=1)
                % All
                obj.vsOriginal = []; obj.vsResampled = [];
                obj.vsRegistered = []; obj.vsSegmented = [];
            else
                ext = obj.GetVolSeqLabel(eType);
                eval(['obj.vs' ext ' = [];']);
            end
        end
        
        function r = IsVolSeqEmpty(obj, eType)
            % Is a volume sequence empty ?
            % NGR 2016 03 04
            r = eval(['isempty(obj.vs' obj.GetVolSeqLabel(eType) ');']);
        end
        
        function ClearVolSeqData(obj, eType)
            % Clear a volume sequence (voxels) from memory
            % NGR 2016 03 04
            if (nargin<=1)
                % All
                for i = 1:length(obj.vsOriginal)
                    obj.vsOriginal(i).data = [];
                    obj.vsResampled(i).data = [];
                    obj.vsRegistered(i).data = [];
                end
            else
                ext = obj.GetVolSeqLabel(eType);
                N = obj.GetNbVol(eType);
                for i = 1:N
                    eval(['obj.vs' ext '(' num2str(i) ').data = [];']);
                end
            end
        end
        
        function r = IsVolSeqDataEmpty(obj, eType)
            % Is a volume sequence data (voxels) empty ?
            % NGR 2016 03 04
            r = 0;
            if ~obj.IsVolSeqEmpty(eType)
                ext = obj.GetVolSeqLabel(eType);
                N = eval(['length(obj.vs' ext ')']);
                for i = 1:N
                    v =  eval(['obj.vs' ext '(' num2str(i) ')']);
                    if v.IsVolEmpty
                        r = 1;
                        return;
                    end
                end
            else
                r = 1;
            end
        end
        
        function r = GetVolSeqLabel(obj, eType)
            % Get volume label extension
            % NGR 2016 03 04            
            r = char(eType);
        end
        
        function r = GetVolRefIndex(obj, eType)
            % Get index of reference volume (resampling & registration)
            % NGR 2016 02 16
            N = eval(['length(obj.vs' obj.GetVolSeqLabel(eType) ');']);
            for i=1:N
                str = eval([ ...
                    'obj.vs' obj.GetVolSeqLabel(eType) ...
                    '(' num2str(i) ').metadata{1}.SeriesDescription;']);
                k = strfind(str,'[Reference]');
                if ~isempty(k)
                    r = i;
                    eval(['obj.iVolRef_' obj.GetVolSeqLabel(eType) '=i;']);
                    return;
                end
            end
            r = [];
        end
        
        function r = GetMaxAxis(obj, eType, eView)
            % Get the maximum axis settings [mmm]
            % NGR 2016 03 16
            dimMax = [0 0];
            for iVol = 1:obj.GetNbVol(eType)
                [~, RI] = obj.GetCentralOrthoPlane(eType, eView, iVol);
                dimMax(1) = max(dimMax(1),max(RI.XWorldLimits));
                dimMax(2) = max(dimMax(2),max(RI.YWorldLimits));
            end
            r = [-dimMax(1) dimMax(1) -dimMax(2) dimMax(2)];
        end
        
        function r = GetHighestRes(obj)
            % Get highest spatial resolution
            % NGR 2016 02 04
            res_mm = [inf inf inf];
            sPx = inf;sPy = inf;sPz = inf;
            for i=1:obj.GetNbObjects
                for j=1:3
                    res_mm(j) = ...
                        min(res_mm(j),obj.vsOriginal(i).res_mm(j));
                end
                sPx  = min([sPx obj.vsOriginal(i).res_mm(1)]);
                sPy  = min([sPy obj.vsOriginal(i).res_mm(2)]);
                sPz  = min([sPz obj.vsOriginal(i).res_mm(3)]);
            end
            % Update
            res_mm(1) = sPx;
            res_mm(2) = sPy;
            res_mm(3) = sPz;
            % Return
            r = res_mm;
        end
        
        function r = GetHighestResolutionIndex(obj)
            % Get highest spatial resolution
            % NGR 2016 02 04
            for i=1:obj.GetNbObjects
                res(i) = obj.vsOriginal(i).GetVoxelSize_mm3;
            end
            [~, r] = min(res);
        end
        
        function r = GetLowestResolutionIndex(obj)
            % Get highest spatial resolution
            % NGR 2016 05 06
            for i=1:obj.GetNbObjects
                res(i) = obj.vsOriginal(i).GetVoxelSize_mm3;
            end
            [~, r] = max(res);
        end
        
        function [r, map] = GetVolIndexed(obj, eType)
            % Convert RGB volume to indexed volume
            % NGR 2016 02 11
            vRGB = obj.GetVolRBG(eType);
            s = size(vRGB);
            Ir = squeeze(vRGB(:,:,:,1)); I(:,:,1) = Ir(:,:);
            Ig = squeeze(vRGB(:,:,:,2)); I(:,:,2) = Ig(:,:);
            Ib = squeeze(vRGB(:,:,:,3)); I(:,:,3) = Ib(:,:);
            [X , map] = rgb2ind(I,256);
            r = reshape(X(:),s(1:3));
        end
        
        function r = GetVolRBG(obj, eType)
            % Get RGB volume (1st three channels, normalization by mean)
            % NGR 2016 02 11
            vs = obj.GetVolSeq(eType);
            r = zeros([size(vs(1).data) 3],'single');
            for i=1:3
                r(:,:,:,i) = ...
                    single(vs(i).data) / single(mean(vs(i).data(:))) / 3;
            end
        end
        
        function r = GetCentralOrthoPlaneIndex(obj, ...
                eType, eView, iVol, bSliceLoc)
            % Get central orthogonal plane index
            % NGR 2016 02 09
            if nargin<5
                bSliceLoc = 0;
            end
            vs = obj.GetVolSeq(eType);
            switch  eView
                case MCC_eView.Axial
                    if bSliceLoc
                        if ...
                                eType == MCC_eType.Original || ...
                                eType == MCC_eType.Resampled
                            deltaZ = obj.GetDeltaZ(eType, iVol);
                            iShift = round(deltaZ/vs(iVol).res_mm(3));
                        else
                            iShift = 0;
                        end
                        Nz = size(vs(iVol).data,3);
                        r = round(Nz/2)+iShift;
                        % Out of range
                        if r<1 || r>Nz
                            r = [];
                        end
                    else
                        r = round(size(vs(iVol).data,3)/2);
                    end
                case MCC_eView.Coronal
                    r = round(size(vs(iVol).data,1)/2);
                case MCC_eView.Sagittal
                    r = round(size(vs(iVol).data,2)/2);
                otherwise
                    r = round(size(vs(iVol).data,3)/2);
            end
        end
        
        function [I, RI] = GetCentralOrthoPlane(...
                obj, eType, eView, iVol, bSliceLoc)
            % Get central orthogonal plane
            % NGR 2016 02 09
            if nargin == 4
                bSliceLoc = 0;
            end
            vs = obj.GetVolSeq(eType);
            iPlane = obj.GetCentralOrthoPlaneIndex( ...
                eType, eView, iVol, bSliceLoc);
            [I, RI] = obj.GetOrthoPlane(...
                eType, eView, iVol, iPlane, bSliceLoc);
        end
        
        function [I, RI] = GetOrthoPlane( ...
                obj, eType, eView, iVol, iPlane, bSliceLoc)
            % Get an orthogonal plane
            % NGR 2016 02 09
            if nargin == 5
                bSliceLoc = 0;
            end
            if isempty(iVol)
                iVol = 1;
            end
            vs = obj.GetVolSeq(eType);
            if bSliceLoc
                R = obj.Volref3d(vs(iVol), iVol);
            else
                R = obj.Volref3d(vs(iVol));
            end
            switch eView
                case MCC_eView.Axial
                    % Axial
                    if isempty(iPlane)
                        I = zeros(size((vs(iVol).data(:,:,1))));
                    else
                        Np = size(vs(iVol).data,3);
                        if iPlane > Np || iPlane < 1
                            % Out of bounds
                            % NGR 2016 09 13                            
                            I = squeeze(vs(iVol).data(:,:,1));
                            I = zeros(size(I));
                        else
                            % In bounds
                            I = squeeze(vs(iVol).data(:,:,iPlane));
                        end
                    end
                    limX_mm = R.XWorldLimits;
                    limY_mm = R.YWorldLimits;
                case MCC_eView.Coronal
                    % Coronal
                    if isempty(iPlane)
                        I = zeros(size((vs(iVol).data(1,:,:))))';
                    else
                        Np = size(vs(iVol).data,1);
                        if iPlane > Np || iPlane < 1
                            % Out of bounds                                                                                  
                            I = squeeze(vs(iVol).data(1,:,:))';
                            I = zeros(size(I));
                        else
                            % In bounds
                            I = squeeze(vs(iVol).data(iPlane,:,:))';
                        end
                    end
                    limX_mm = R.XWorldLimits;
                    limY_mm = R.ZWorldLimits;
                case MCC_eView.Sagittal
                    % Sagital
                    if isempty(iPlane)
                        I = zeros(size((vs(iVol).data(:,1,:))))';
                    else
                        Np = size(vs(iVol).data,2);
                        if iPlane > Np || iPlane < 1
                            % Out of bounds
                            iPlane = Np;
                            I = squeeze(vs(iVol).data(:,1,:))';
                            I = zeros(size(I));
                        else
                            % In bounds
                            I = squeeze(vs(iVol).data(:,iPlane,:))';
                        end                     
                    end
                    limX_mm = R.YWorldLimits;
                    limY_mm = R.ZWorldLimits;
                    % Default (axial)
                otherwise
                    if isempty(iPlane)
                        I = zeros(size((vs(iVol).data(:,:,1))));
                    else
                        I = squeeze(vs(iVol).data(:,:,iPlane));
                    end
                    limX_mm = R.XWorldLimits;
                    limY_mm = R.YWorldLimits;
            end
            RI = imref2d(size(I));
            RI.XWorldLimits = limX_mm;
            RI.YWorldLimits = limY_mm;
        end
        
        function r = AdjustImaView(obj, I)
            % Adjust image view for display (cornoal and sagittal case)
            % NGR 2016 02 18
            r = imrotate(I,90);
        end
        
        function vs = GetVolSeq(obj, eType)
            % Get a volume sequence (orginal, resampled, registered ...)
            % NGR 2016 02 05
            r = eval(['obj.vs' obj.GetVolSeqLabel(eType) ';']);
            vs = obj.CloneVolSeq(r);
        end
        
        function SetVolSeq(obj, eType, vs)
            % Set a volume sequence (orginal, resampled, registered ...)
            % NGR 2016 05 05
            eval([...
                'obj.vs' obj.GetVolSeqLabel(eType) ...
                ' = obj.CloneVolSeq(vs);' ...
                ]);
        end
        
        function vs = CloneVolSeq(obj, vsI)
            % Clone a volume sequence
            % NGR 2016 02 19
            vs = MCC_CloneObj(vsI);
            for i=1:length(vsI)
                vs(i) = MCC_CloneObj(vsI(i));
            end
        end
        
        function r = GetNbVol(obj, eType)
            % Get number of volumes
            % NGR 2016 03 01 04
            if (nargin<=1), eType = MCC_eType.Original; end;
            r = eval(['length(obj.vs' obj.GetVolSeqLabel(eType) ');']);
        end
        
        function r = GetVolSeqSize(obj, eType)
            % Get volume sequence size
            % NGR 206 02 09
            Nvol = obj.GetNbVol;
            r = zeros([Nvol 3]);
            for i = 1:obj.Nvol
                switch  eType
                    case MCC_eType.Original
                        r(i,:) = size(obj.vsOriginal(1).data);
                    case MCC_eType.Resampled
                        r(i,:) = size(obj.vsResampled(1).data);
                    case MCC_eType.Registered
                        r(i,:) = size(obj.vsRegistered(1).data);
                    case MCC_eType.Segmented
                        r(i,:) = size(obj.vsSegmented(1).data);
                    otherwise
                        r = obj.vsOriginal;
                end
            end
        end
        
        function ApplyTransform(obj, eType, tform, Ro)
            % Apply a geometrical transform
            % NGR 2018 03 21
            
            vs = obj.GetVolSeq(eType);
            ext = obj.GetVolSeqLabel(eType);
            for i=1:obj.GetNbObjects
                V = vs(i).data;
                R = obj.Volref3d(vs(i));
                if exist('Ro','var')
                    Vt = imwarp(V ,R, tform, 'bicubic','OutputView',Ro);
                else
                    Vt = imwarp(V ,R, tform, 'bicubic','OutputView',R);
                end
                eval(['obj.vs' ext '(' num2str(i) ').data = Vt;']);
                eval(['obj.vs' ext '(' num2str(i) ').eType = eType;']);
            end
        end
        
        function CreateReadyVolSeq(obj, tform)
            % Create ready sequence (registered study dates)
            % NGR 2016 04 12
            eT = MCC_eType.Ready;
            obj.vsReady = obj.CloneVolSeq(obj.vsRegistered);
            obj.ApplyTransform(eT, tform);
            obj.Write(eT);
            obj.eTypeCur = eT;
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
            MCC_DisplayContent(obj.content);
        end
        
        function r = GetNbObjects(obj, eType)
            % Get number of objects
            % NGR 2016 03 23
            if nargin == 2
                switch eType
                    case MCC_eType.Segmented
                        r = length(obj.contentSeg);
                    otherwise
                        r = length(obj.content);
                end
            else
                r = length(obj.content);
            end
        end
        
        function r = GetVOI(obj, eType, VolRef, eOrgan)
            % Get VOI
            % NGR 2016 06 24
            switch nargin
                case 1
                    eType  = MCC_eType.Registered;
                    VolRef = 1;
                    eOrgan = MCC_eOrgan.Kidney;
                case 2
                    VolRef = 1;
                    eOrgan = MCC_eOrgan.Kidney;
                case 3
                    eOrgan = MCC_eOrgan.Kidney;
            end
            % Get
            obj.LoadWorkspace;
            obj.InitializeVOI;
            r = obj.VOI;
        end
                
    end
    
    methods (Access = public)   % NGR 2018 03 21
        
        %-----------------------------------------------------------------%
        % 3D processing
        %-----------------------------------------------------------------%
        
        function [tF, MI, Vr] = Volregtform(obj, ...
                Vf, Rf, Vm, Rm, tType, tFI, Nbits)
            % Registration (volume)
            % NGR 2016 02 23
            
            bDispOpt = 0;       % Verbose optimization mode
            Nbins = 2^Nbits;    % # of histogram bins
            Niter = 1000;   % 1000       
            Npl   = 3;   % 3    
            
            %Niter = 1; Npl   = 1; % @@@
            
            % All dimensions of the fixed and moving images should be
            % greater than or equal to 16 for 'PyramidLevels' 3.
            sMin = min([size(Vf) , size(Vm)]);
            if sMin<16
                Npl   = 1;
            end
            
            if (sum(Vf(:)==0)+sum(Vf(:)==1))==length(Vf(:)) % Binary ?
                % Edge-based registration
                % Optimizer
                [opt,met] = imregconfig('monomodal');
                opt.MaximumIterations = Niter;
                opt.MaximumStepLength = opt.MaximumStepLength/2;
                % Registration
                tF = imregtform( ...
                    Vm,Rm, ...
                    Vf,Rf,...
                    tType, opt, met, ...
                    'InitialTransformation', tFI,...
                    'PyramidLevels', Npl, ...
                    'DisplayOptimization', bDispOpt);
                Vr = imwarp(Vm,Rm, tF,'bicubic','OutputView',Rf);
                MI = MCC_MutualInfomation(Vf, Vr, Nbins);
            else
                % Intensity-based and gradient-based registration
                iniRadiusF = 3;
                % Optimizer
                [opt,met] = imregconfig('multimodal');
                opt.InitialRadius = opt.InitialRadius / iniRadiusF;
                met.NumberOfHistogramBins = Nbins;
                opt.MaximumIterations = Niter;
                % Registration
                tF = imregtform( ...
                    Vm,Rm, ...
                    Vf,Rf,...
                    tType, opt, met, ...
                    'InitialTransformation', tFI,...
                    'PyramidLevels', Npl, ...
                    'DisplayOptimization', bDispOpt);
                Vr = imwarp(Vm,Rm, tF,'bicubic','OutputView',Rf);
                MI = MCC_MutualInfomation(Vf, Vr, Nbins);
            end
        end
        
         function [tF, MI, Vr] = VolregtformTest(obj, ...
                Vf, Rf, Vm, Rm, tType, tFI, Nbits,i,j,k)
            % Registration (volume)
            % NGR 2016 02 23
            
            if j==2 && k==3
%                 save('Part1Volregtform.mat')  
                nop=1;
            end
            
            
            bDispOpt = 0;       % Verbose optimization mode
            Nbins = 2^Nbits;    % # of histogram bins
            Niter = 1000;   % 1000       
            Npl   = 3;   % 3    
            
%             Niter = 1; Npl   = 1; % @@@
            
            % All dimensions of the fixed and moving images should be
            % greater than or equal to 16 for 'PyramidLevels' 3.
            sMin = min([size(Vf) , size(Vm)]);
            if sMin<16
                Npl   = 1;
            end
            
            if (sum(Vf(:)==0)+sum(Vf(:)==1))==length(Vf(:)) % Binary ?
                % Edge-based registration
                % Optimizer
                [opt,met] = imregconfig('monomodal');
                opt.MaximumIterations = Niter;
                opt.MaximumStepLength = opt.MaximumStepLength/2;
                
                if j==2 && k==3
%                     save('Part2Op1Volregtform.mat')  
                    nop=1;
                end
                % Registration
                if ~strcmp(version('-release'),'2018b')
                    tF = imregtform( ...
                        Vm,Rm, ...
                        Vf,Rf,...
                        tType, opt, met, ...
                        'InitialTransformation', tFI,...
                        'PyramidLevels', Npl, ...
                        'DisplayOptimization', bDispOpt);
                else
                    filein=[num2str(i) '_' num2str(j) '_' num2str(k) '.mat'];
                    save(filein,'Vm','Rm', ...
                            'Vf','Rf',...
                            'tType', 'opt', 'met', ...
                             'tFI',...
                             'Npl', ...
                             'bDispOpt');
                    aux=strsplit(filein,'.');
                    system(['imregtform2 ' filein ' ' aux{1} 'out.mat']);
                    load([aux{1} 'out.mat']);
                    delete([aux{1} 'out.mat']);
                    delete(filein);
                end
                Vr = imwarp(Vm,Rm, tF,'bicubic','OutputView',Rf);
                MI = MCC_MutualInfomation(Vf, Vr, Nbins);
                if j==2 && k==3
%                     save('Part3Op1Volregtform.mat')  
                    nop=1;
                end
            else
                % Intensity-based and gradient-based registration
                iniRadiusF = 3;
                % Optimizer
                [opt,met] = imregconfig('multimodal');
                opt.InitialRadius = opt.InitialRadius / iniRadiusF;
                met.NumberOfHistogramBins = Nbins;
                opt.MaximumIterations = Niter;
                if j==2 && k==3
%                     save('Part2Op2Volregtform.mat')  
                    nop=1;
                end
                % Registration
                if ~strcmp(version('-release'),'2018b')
                    tF = imregtform( ...
                        Vm,Rm, ...
                        Vf,Rf,...
                        tType, opt, met, ...
                        'InitialTransformation', tFI,...
                        'PyramidLevels', Npl, ...
                        'DisplayOptimization', bDispOpt);
                else
                
                    filein=[num2str(i) '_' num2str(j) '_' num2str(k) '.mat'];
                    save(filein,'Vm','Rm', ...
                            'Vf','Rf',...
                            'tType', 'opt', 'met', ...
                             'tFI',...
                             'Npl', ...
                             'bDispOpt');
                    aux=strsplit(filein,'.');
                    system(['imregtform2 ' filein ' ' aux{1} 'out.mat']);
                    load([aux{1} 'out.mat']);
                    delete([aux{1} 'out.mat']);
                    delete(filein)  ;      
                end
                
                Vr = imwarp(Vm,Rm, tF,'bicubic','OutputView',Rf);
                MI = MCC_MutualInfomation(Vf, Vr, Nbins);
                if j==2 && k==3
%                     save('Part3Op2Volregtform.mat')  
                    nop=1;
                end
            end
            
        end
        
        function R = Volref3d(obj, c3D, iVol)
            % Reference volume to world coordinates
            % NGR 2016 02 23
            if nargin==2,
                iVol = 0;
            end
            V  = c3D.data;
            r = c3D.res_mm;
            R = imref3d(size(V),r(2),r(1),r(3));
            % Center
            lm = R.XWorldLimits(2)-R.XWorldLimits(1);
            R.XWorldLimits = [-lm/2 lm/2];
            lm = R.YWorldLimits(2)-R.YWorldLimits(1);
            R.YWorldLimits = [-lm/2 lm/2];
            lm = R.ZWorldLimits(2)-R.ZWorldLimits(1);
            R.ZWorldLimits = [-lm/2 lm/2];
            if iVol
                deltaZ = obj.GetDeltaZ(c3D.eType, iVol);
                R.ZWorldLimits = R.ZWorldLimits + deltaZ;
            end
        end
        
        function r = GetDeltaZ(obj, eType, iVol)
            % Get delta Z based on slice location
            % NGR 2016 03 23
            if ...
                    eType == MCC_eType.Original || ...
                    eType == MCC_eType.Resampled
                % Resampled ?
                if isempty(obj.iVolRef_Resampled)
                    obj.iVolRef_Resampled = 1;
                end
                if isempty(obj.iVolRef_Registered)
                    obj.iVolRef_Registered = 1;
                end
                % Compute delta Z [mm]
                z1 = obj.vsOriginal(obj.iVolRef_Registered).GetSliceLocation(1);
                z2 = obj.vsOriginal(iVol).GetSliceLocation(1);
                deltaZ = z1 - z2;
                r = deltaZ;
                % Slice contiguousity ?
                b1 = obj.vsOriginal(obj.iVolRef_Registered).CheckSliceContiguousity;
                b2 = obj.vsOriginal(iVol).CheckSliceContiguousity;
                if ~(b1 && b2)
                    % Warning
                    if isempty(obj.bContiguous{iVol})
                        me.identifier = 'MATLAB:UndefinedFunction';
                        me.message = 'Warning: none contiguous slices';
                        me.cause = {};
                        me.stack.file = 'MCC_c4D.m';
                        me.stack.name = 'GetDeltaZ';
                        me.stack.line = 0;
                        MCC_writeLog(...
                            [cd '\log.txt'], me, obj.vsOriginal(iVol).root, 0);
                        r = 0;
                        obj.bContiguous{iVol} = 0;
                    end
                else
                    obj.bContiguous{iVol} = 1;
                end
            else
                r = 0;
            end
        end
        
        %-----------------------------------------------------------------%
        % 4D processing (flipud)
        %-----------------------------------------------------------------%
        
        function Flipud(obj, eT)
            % Flip volumes in up / down direction
            % NGR 2018 041 19 v205
            
            obj.Load;
            obj.Load(eT);
            ext = obj.GetVolSeqLabel(eT);
            for i = 1:obj.GetNbObjects
                cmd = ['obj.vs' ext '(' num2str(i) ').Flipud;'];
                eval(cmd);
            end
        end

        %-----------------------------------------------------------------%
        % 4D processing (resampling)
        %-----------------------------------------------------------------%
        
        function vsO = ResampleVolSeq(obj, vsI, res_mm_o)
            % Resampling of a volume sequence (res_mm_o: output resolution)
            % NGR 2016 06 21
            iP = 2;
            % Display
            obj.DisplayExeLine(iP);
            % Instantiation
            vsO = obj.CloneVolSeq(vsI);
            % Volume loop
            for iVol=1:length(vsI)
                if ~strcmp(vsO(iVol).name,'Contours') % NGR 2016 10 21
                    vsO(iVol).Resample(res_mm_o);
                end
            end
            % Execution time
            obj.DisplayExeTime(iP);
        end
        
        %-----------------------------------------------------------------%
        % 4D processing (registration)
        %-----------------------------------------------------------------%
        
        function vsO = RegisterVolSeq(obj, vsI, iVolRef)
            % Registration of a volume sequence (iVolRef: reference indexs)
            % NGR 2016 06 22
            
            % FMX, NGR 2016 11 22
            % No registration for tubes
            if obj.eOrganCur == MCC_eOrgan.Tubes_FMX                
                vsO = obj.CloneVolSeq(vsI);
                return;
            end
            
            iP = 3;
            % Display
            obj.DisplayExeLine(iP);
            % Settings
            Nbits = 5;          % quantization (5 bits -> 32 levels)
            Nbins = 2^Nbits;    % number of histogram bins
            % Instantiation
            vsO = obj.CloneVolSeq(vsI);
            % Mask
%             masks = ones([size(vsO(1).data) length(vsO)]);
            % Instensity-based
            vs(1,:) = obj.CloneVolSeq(vsI);
%             vs(1,:) = obj.NormalizeMean(vs(1,:));
            vs(1,:) = obj.NormalizeMeanV2(vs(1,:));            
            r = -inf;
            for i=1:length(vs(1,:)), r = max(r,max(vs(1,i).data(:))); end
            vs(1,:) = obj.Multiply(vs(1,:), (Nbins-1)/r );
%             vs(1,:) = obj.Quantize(vs(1,:));
            vs(1,:) = obj.QuantizeV2(vs(1,:));
            % Gradient-based
            vs(2,:) = obj.CloneVolSeq(vsI);
            vs(2,:) = obj.AbsGradient(vs(2,:));
%             vs(2,:) = obj.NormalizeMean(vs(2,:));
            vs(2,:) = obj.NormalizeMeanV2(vs(2,:));
            r = -inf;
            for i=1:length(vs(2,:)), r = max(r,max(vs(2,i).data(:))); end
            vs(2,:) = obj.Multiply(vs(2,:), (Nbins-1)/r );
%             vs(2,:) = obj.Quantize(vs(2,:));
            vs(2,:) = obj.QuantizeV2(vs(2,:));
            % Edge-based
            vs(3,:) = obj.Edge(obj.CloneVolSeq(vsI));
%             vs(3,:) = obj.Quantize(vs(3,:));
            vs(3,:) = obj.QuantizeV2(vs(3,:));
            % Fixed volume (reference)
            Vf{1} = vs(1,iVolRef).data;
            Vf{2} = vs(2,iVolRef).data;
            Vf{3} = vs(3,iVolRef).data;
            Rf    = obj.Volref3d(vs(1,iVolRef));
            % Transforms
            tType = obj.tType;
            % Mutual information
            Nv = length(vsI);
            MI = zeros(Nv,length(tType),size(vs,1));
            pc = zeros(Nv,length(tType),size(vs,1));
            for i=1:Nv
                for k=1:length(tType)
                    tF(i,k,1:size(vs,1))=affine3d(eye(4));
                end
            end
            % Perfect alignement (reference with herself)
            Vf0 = vsI(iVolRef).data;
            MIref =  MCC_MutualInfomation(Vf0, Vf0, Nbins);
            % Volume loop (i)
            x=vsO;
            tforms=[];
            aux1=cell(Nv,1);
            aux2=cell(Nv,1);
            parfor i=1:Nv % FMX (Skip contour), NGR 2016 09 05   If multicore parfor             
                 %J(i)=batch(@register_eachVolFun,2,{obj,vsI,vs,i,tType,Rf,Vf0,Nbins,iVolRef,MI,tF,Vf,Nbits,x,MIref});
                    [aux1{i},aux2{i}]=register_eachVolFun(obj,vsI,vs,i,tType,Rf,Vf0,Nbins,iVolRef,MI,tF,Vf,Nbits,x,MIref);
%                     vsO(i)=aux1;
%                     tforms{i}=aux2;
            end                    %8/3
            for i=1:Nv
%                 wait(J(i));
%                 diary(J(i))
%                 outputs=fetchOutputs(J(i));
%                 vsO(i)=outputs{1};
%                 obj.tform{i}=outputs(2);
                vsO(i)=aux1{i};
                obj.tform{i}=aux2{i};
%                 delete(J(i))
            end

            % Intersected mask
            masks = ones([size(vsO(1).data) length(vsO)]);
            sz = size(masks);
            mask = ones(sz(1:3));
            for iVol=1:obj.GetNbObjects
                mask = mask .* masks(:,:,:,iVol);
            end
            %obj.maskReg = mask; % @NGR
            obj.maskReg =[];
            % Update
            obj.rMI = MI;
            obj.MIpc = pc;
            % Execution time
            obj.DisplayExeTime(iP);
        end
        
        
        function vsOi=register_eachVol(obj,vsI,vs,i,tType,Rf,Vf0,Nbins,iVolRef,MI,tF,Vf,Nbits,vsO)

        %save(['IterationS1_' num2str(i) '.mat'])     
        % Storage (original data)
        Vm0 = vsI(i).data;
        % Moving volumes
        Vm{1} = vs(1,i).data;
        Vm{2} = vs(2,i).data;
        Vm{3} = vs(3,i).data;
        Rm = obj.Volref3d(vs(1,i));
        % Is centered volume a better initialization ?
        Vr0 = imwarp(...
            Vm0 ,Rm, affine3d(eye(4)),...
            'bicubic','OutputView',Rf);
        MIc(i) = MCC_MutualInfomation(Vf0, Vr0, Nbins);
        Rm0 = obj.Volref3d(vs(1,i), iVolRef);
        Vr0 = imwarp(...
            Vm0 ,Rm0, affine3d(eye(4)), ...
            'bicubic','OutputView',Rf);
        MI0(i) = MCC_MutualInfomation(Vf0, Vr0, Nbins);
        if MI0(i)>MIc(i)
            Rm = Rm0;
        end
        %save(['IterationS2_' num2str(i) '.mat'])  
        % Reference volume
        if (i==iVolRef) || strcmp(vsI(i).name,'Contours') % NGR 2016 10 21
            MutInfo =  MCC_MutualInfomation(Vf0, Vf0, Nbins);
            for iType = 1:length(tType)
                for iData = 1:length(Vm)
                    MI(i,iType,iData) = MutInfo;
                end
            end
            % Update
            obj.tform{i} = ...
                affine3d(eye(4));
            if (i==iVolRef) % NGR 2016 10 21
                vsO(i).name = ...
                    [vsO(i).name ' [Reference]'];
            end
            vsO(i).UpdateMetadata;
            vsO(i).eType = MCC_eType.Registered;
        else
            % Verbose
            tic;
            scrF = vsI(i).root;
            %%%% 7/19
            MCC_Disp(sprintf('Fixed image:'), 5);
            MCC_Disp(vsI(iVolRef).root,5);
            MCC_Disp(sprintf('Moving image:'), 5);
            %%%%%%%%
            MCC_Disp(scrF,5);
            MCC_Disp(sprintf('Type        Inte.       Grad.       Edg.        Inte.       Grad.       Edg.'), 5);
            MCC_Disp(sprintf('----        ----------------------------        ----------------------------'), 5);
            % Type loop (j)
            for j=1:length(tType)
                MCC_Disp(tType{j},5);
                % Best previous transform for initialization
                MI_ = squeeze(MI(i,:,:));
                tF_ = squeeze(tF(i,:,:));
                [~, ind] = max(MI_(:));
                tFini = tF(ind);
                %save(['IterationS3_' num2str(i) '.mat'])  
                % Data loop
                for k = 1:length(Vm)
                    if j==1
                        % Non transformation case
                        tF(i,j,k) = tFini;
                    else
                        try
                            % Transformation search
                            if isempty(...
                                    strfind(obj.content{i},'Mask'))
                                % Similarity & Affine transformations limited to
                                % the smallest FOV
                                if ...
                                        j>=3 && ...
                                        sum(size(Vf{k})-size(Vm{k}))~=0 ...
                                        % Dimensions
                                    Rx(1)=max(Rm.XWorldLimits(1),Rf.XWorldLimits(1));
                                    Rx(2)=min(Rm.XWorldLimits(2),Rf.XWorldLimits(2));
                                    Ry(1)=max(Rm.YWorldLimits(1),Rf.YWorldLimits(1));
                                    Ry(2)=min(Rm.YWorldLimits(2),Rf.YWorldLimits(2));
                                    Rz(1)=max(Rm.ZWorldLimits(1),Rf.ZWorldLimits(1));
                                    Rz(2)=min(Rm.ZWorldLimits(2),Rf.ZWorldLimits(2));
                                    crop_mm = [Rx ; Ry ; Rz];
                                    % Fixed volume (crop)
                                    Vfc = MCC_CloneObj(vsI(iVolRef));
                                    Vfc.data = Vf{k};
                                    Vfc.Crop(crop_mm);
                                    Rfc = obj.Volref3d(Vfc);
                                    % Moving volume (crop)
                                    Vmc = MCC_CloneObj(vsI(1));
                                    Vmc.data = Vm{k};
                                    Vmc.Crop(crop_mm);
                                    Rmc = obj.Volref3d(Vmc);
                                    % Estimate transformation
                                    tF(i,j,k) = ...
                                        obj.Volregtform( ...
                                        Vfc.data, Rfc, ...
                                        Vmc.data, Rmc, ...
                                        tType{j}, ...
                                        tFini, Nbits);
                                else
                                    % Estimate transformation
                                    tF(i,j,k) = ...
                                        obj.Volregtform( ...
                                        Vf{k}, Rf, ...
                                        Vm{k}, Rm, ...
                                        tType{j}, ...
                                        tFini,Nbits);
                                    %save(['IterationS4_' num2str(i) '.mat'])  
                                end
                            else
                                tF(i,j,k) = obj.tform{1};
                            end
                        catch ME
                            MCC_writeLog([cd '\log.txt'],ME, ...
                                ['(' tType{j} ' registration)' ...
                                vsI(i).root]);
                        end
                    end
                    % Register
                    Vr{k} = imwarp( ...
                        Vm0 ,Rm, ...
                        tF(i,j,k),'bicubic','OutputView',Rf);
                    % Mutual information (MI)
                    MI(i,j,k) = MCC_MutualInfomation(...
                        Vf0, Vr{k}, Nbins);
                    % Relative MI
                    pc(i,j,k) = ...
                        100*(MI(i,j,k)-MI(i,1,1))/MI(i,1,1);
                    %save(['IterationS5_' num2str(i) '.mat'])  
                end
                % Echo
                c = { ....
                    sprintf('\t\t\t%1.2f', MI(i,j,1)), ...
                    sprintf('\t\t%1.2f',   MI(i,j,2)), ...
                    sprintf('\t\t%1.2f',   MI(i,j,3)), ...
                    sprintf('\t\t%1.1f%%', pc(i,j,1)), ...
                    sprintf('\t\t%1.1f%%', pc(i,j,2)), ...
                    sprintf('\t\t%1.1f%%', pc(i,j,3)) ...
                    };
                MCC_Disp(strjoin(c,''),5);
            end
            % Select the best and update
            MI_ = squeeze(MI(i,:,:));
            pc_ = squeeze(pc(i,:,:));
            tF_ = squeeze(tF(i,:,:));
            [~, ind]=max(MI_(:));
            %save(['IterationS6_' num2str(i) '.mat'])  
            % Apply final traformation
            if isempty(strfind(obj.content{i},'Mask'))
                obj.tform{i} = tF_(ind);
                % Organ specifities
                switch obj.eOrganCur
                    case MCC_eOrgan.Kidney
                    case MCC_eOrgan.Prostate
                        %obj.tform{i} = affine3d(eye(4));
                    otherwise
                end
                % Mask
                Vrm = imwarp( ...
                    single(Vm0),Rm, obj.tform{i}, ...
                    'bicubic','OutputView',Rf, 'FillValues', -1);
                try
                    masks(:,:,:,i) = (Vrm>=0);
                end
                % Trans
                vsO(i).data = uint16( ...
                    imwarp( ...
                    Vm0,Rm, obj.tform{i}, ...
                    'bicubic','OutputView',Rf));
                %save(['IterationS7_' num2str(i) '.mat'])  
            else
                % CT mask
                R = obj.Volref3d(vsI(i));
                V = single(vsI(i).data);
                V = 255 * V / max(V(:));
                % Transformation
                Vmask_CT = imwarp( ...
                    V, R, obj.tform{1}, ...
                    'bicubic','OutputView',R);
                % Croping and resampling
                Vmask_MR = imwarp( ...
                    Vmask_CT, R, affine3d(eye(4)), ...
                    'bicubic','OutputView',Rf);
                vsO(i).data = uint16(Vmask_MR);
            end
            % Selection
            selT = {...
                'none' ,'rigid', 'simi.', 'affi.', ...
                'none' ,'rigid', 'simi.', 'affi.',...
                'none' ,'rigid', 'simi.'  'affi.'...
                };
            selD = {...
                'inte.','inte.', 'inte.', 'inte.',...
                'grad.','grad.', 'grad.', 'grad.',...
                'edge' ,'edge' , 'edge'   'edge.'...
                };
            MCC_Disp( [ ...
                selT{ind} ' (' selD{ind} ')', ...
                sprintf(' = %1.2f (%1.1f%%)', ...
                MI_(ind), pc_(ind))] , 26);
            disp(' ');
            % Update
            vsO(i).UpdateMetadata;
            vsO(i).eType = MCC_eType.Registered;
            %save(['IterationS8_' num2str(i) '.mat'])  
            % CSV file
            indF=strfind(obj.root, obj.InstituteID);
            fileN_csv=[obj.root '\Registration.csv'];
            if exist(fileN_csv, 'file')
                fid_csv = fopen(fileN_csv,'at+');
            else
                fid_csv = fopen(fileN_csv,'wt+');
                try
                fprintf(fid_csv,...
                    '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                    'Version', ...
                    'Year', ...
                    'Month', ...
                    'Day', ...
                    'Hour', ...
                    'Min', ...
                    'Sec.', ...
                    'Institute', ...
                    'Database', ...
                    'Cohort', ...
                    'Patient', ...
                    'Study date', ...
                    'Imaging Mode', ...
                    'Imaging Mode {r}', ...
                    'Vx [mm3] {i}', ...
                    'Vx [mm3] {f}', ...
                    'MI {r}', ...
                    'MI {i0}', ...
                    'MI {i}', ...
                    'MI {f}', ...
                    'MI {(f-i0)/i0)} [%]', ...
                    'MI {(f-i)/i)} [%]', ...
                    'MI {i0/r} [%]', ...
                    'MI {f/r} [%]', ...
                    'Method', ...
                    'Transform', ...
                    'Data', ...
                    'Time [s]', ...
                    'Comments', ...
                    'Path' ...
                    );
                catch
                    nop=1;
                end
            end
            clk = strsplit(MCC_clock2str,' ');
            try                  % 8/2
            fprintf(fid_csv, ...
                '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%s,%s,%1.2f,%s,%s\n', ...
                obj.version, ...
                clk{2}, clk{3}, clk{4}, clk{5}, clk{6}, clk{7}, ...
                obj.InstituteID, ...
                obj.DatabaseID, ...
                obj.CohortID, ...
                obj.PatientID, ...
                obj.StudyDateID, ...
                strrep(vsI(i).name, ',', ' '), ...
                strrep(vsI(iVolRef).name, ',', ' '), ...
                vsI(i).GetVoxelSize_mm3, ...
                vsO(i).GetVoxelSize_mm3, ...
                MIref, ...
                MI0(i), ...
                MI_(1), ...
                MI_(ind), ...
                100*(MI_(ind)-MI0(i))/MI0(i), ...
                pc_(ind), ...
                100*MI0(i)/MIref, ...
                100*MI_(ind)/MIref, ...
                ind, ...
                selT{ind}, ...
                selD{ind}, ...
                toc, ...
                ' ', ...
                scrF ...
                );
                fclose(fid_csv);
            end
        end
        %save(['IterationS9_' num2str(i) '.mat'])  
        vsOi=vsO(i);
%         if exist(['IterationPar_' num2str(i) '.mat'])
%             save(['IterationParD2_' num2str(i) '.mat'])
%         else
%             save(['IterationPar_' num2str(i) '.mat'])
%         end
        
        end
        
        %-----------------------------------------------------------------%
        % 4D processing (misc.)
        %-----------------------------------------------------------------%
        
       
        function Vo = ForceFOV(obj, Vf, Vm)
            % Force the same FOV as the reference (Vf)
            % NGR 2016 06 29
            Rf = obj.Volref3d(Vf);
            Rm = obj.Volref3d(Vm);
            x = Vm.data;
            data = imwarp(x ,Rm, ...
                affine3d(eye(4)),'bicubic','OutputView',Rf);
            Vo = MCC_CloneObj(Vm);
            Vo.data = data;
        end
        
        function r = Adjust(obj, vsI)
            % Intensity adjustement (volume sequence)
            % NGR 2016 02 24
            
            % 1% of data is saturated at low and high intensities of vsI
            r = obj.CloneVolSeq(vsI);
            for i=1:length(r);
                for j=1:size(r(i).data,3)
                    r(i).data(:,:,j) =  imadjust(vsI(i).data(:,:,j));
                end
            end
        end
        
        function r = Edge(obj, vsI, method)
            % Edge detector (volume sequence)
            % Method:
            %        'Canny', 'Sobel'
            % NGR 2016 02 98
            if nargin<3
                method = 'Canny';
            end
            r = obj.CloneVolSeq(vsI);
            for i=1:length(r);
                for j=1:size(r(i).data,3)
                    r(i).data(:,:,j) =  edge(double(vsI(i).data(:,:,j)),method);
                end
            end
        end
        
        function r = NormalizeMax(obj, vsI)
            % Intensity normalization (volume sequence)
            % NGR 2016 02 18
            r = obj.CloneVolSeq(vsI);
            % Normalize volume per volume
            for i=1:length(r);
                m = single(max(vsI(i).data(:)));
                r(i).data = ...
                    single(vsI(i).data) / m;
            end
        end
        
        function r = NormalizeMean(obj, vsI)
            % Intensity normalization (volume sequence)
            % NGR 2016 02 18
            r = obj.CloneVolSeq(vsI);
            % Normalize volume per volume
            for i=1:length(r);
                m = single(mean(vsI(i).data(:)));
                r(i).data = ...
                    single(vsI(i).data) / m;
            end
        end
        
        function r = NormalizeMeanV2(obj, vsI)
            % Intensity normalization (volume sequence)
            % NGR 2016 02 18
            r = obj.CloneVolSeq(vsI);
            % Normalize volume per volume
            for i=1:length(r);
                m = double(mean(vsI(i).data(:)));
                r(i).data = ...
                    double(vsI(i).data) / m;
            end
        end
        
        function r = Multiply(obj, vsI, val)
            % Scalar multiplication (volume sequence)
            % NGR 2016 02 22
            r = obj.CloneVolSeq(vsI);
            % Normalize volume per volume
            for i=1:length(r);
                r(i).data = vsI(i).data * val;
            end
        end
        
        function r = MultiplyByVol(obj, vsI, Vm)
            % Volume sequence multiplied by a volume
            % NGR 2016 05 11
            r = obj.CloneVolSeq(vsI);
            % Normalize volume per volume
            for i=1:length(r);
                r(i).data = vsI(i).data .* Vm;
            end
        end
        
        function r = MultiplyByVol_GPU(obj, vsI, Vm)
            % Volume sequence multiplied by a volume (GPU version)
            % NGR 2016 05 27
            r = obj.CloneVolSeq(vsI);
            V1 = gpuArray(Vm);
            % Normalize volume per volume
            for i=1:length(r);
                V2 = gpuArray(vsI(i).data);
                V3 = V1 .* V2;
                r(i).data = gather(V3);
            end
        end
        
        function r = Sqrt(obj, vsI)
            % Sqrt of volume sequence (GPU)
            % NGR 2016 05 27
            r = obj.CloneVolSeq(vsI);
            % Normalize volume per volume
            for i=1:length(r);
                % Create array on GPU
                r(i).data = sqrt(single(vsI(i).data));
            end
        end
        
        function r = Sqrt_GPU(obj, vsI)
            % Sqrt of volume sequence (GPU)
            % NGR 2016 05 27
            r = obj.CloneVolSeq(vsI);
            % Normalize volume per volume
            for i=1:length(r);
                % Create array on GPU
                Vgpu = gpuArray(single(vsI(i).data));
                Vgpu = sqrt(Vgpu);
                r(i).data = gather(Vgpu);
            end
        end
        
        function r = Quantize(obj, vsI)
            % Quantization (volume sequence)
            % NGR 2016 02 22
            r = obj.CloneVolSeq(vsI);
            for i=1:length(r)
                r(i).data = uint16(round(vsI(i).data));
            end
        end
        
        function r = QuantizeV2(obj, vsI)
            % Quantization (volume sequence)
            % NGR 2016 02 22
            r = obj.CloneVolSeq(vsI);
            for i=1:length(r)
%                 r(i).data = uint32(round(100000*vsI(i).data));
                 r(i).data = uint32(round(vsI(i).data));
            end
        end
        
        function r = AbsGradient(obj, vsI)
            % Numerical gradient (volume sequence)
            % NGR 2016 02 17
            r = obj.CloneVolSeq(vsI);
            for i=1:length(r)
                [Gx, Gy, Gz] = gradient(single(vsI(i).data));
                r(i).data = sqrt(Gx.^2+Gy.^2+Gz.^2);
            end
        end
        
        function r = Abs(obj, vsI)
            % Lastgest voxel intensity (volume sequence)
            % NGR 2016 02 17
            r = obj.CloneVolSeq(vsI);
            for i=1:length(vsI)
                r(i).data = abs(r(i).data);
            end
        end
        
        function r = Max(obj, eType)
            % Lastgest voxel intensity (volume sequence)
            % NGR 2016 02 17
            r = -inf;
            vs = obj.GetVolSeq(eType);
            for i=1:length(vs)
                r = max(r,max(vs(i).data(:)));
            end
        end
        
        function r = Min(obj, eType)
            % Smallest voxel intensity (volume sequence)
            % NGR 2016 02 17
            r = +inf;
            vs = obj.GetVolSeq(eType);
            for i=1:length(vs)
                r = min(r,min(vs(i).data(:)));
            end
        end
        
        function r = Not(obj, vsI)
            % Logical not
            % NGR 2016 05 11
            r = obj.CloneVolSeq(vsI);
            for i=1:length(r);
                for j=1:size(r(i).data,3)
                    r(i).data(:,:,j) =  not(vsI(i).data(:,:,j));
                end
            end
        end
        
        function r = GetVolSeqMinFOV(obj, vsI)
            % Get smallest Field of View (FOV)
            % NGR 2016 05 05
            for iVol = 1:length(vsI)
                R = obj.Volref3d(vsI(iVol));
                % Lower limit
                R1(iVol,1) = R.XWorldLimits(1);
                R1(iVol,2) = R.YWorldLimits(1);
                R1(iVol,3) = R.ZWorldLimits(1);
                % Upper limit
                R2(iVol,1) = R.XWorldLimits(2);
                R2(iVol,2) = R.YWorldLimits(2);
                R2(iVol,3) = R.ZWorldLimits(2);
            end
            R1c = max(R1);
            R2c = min(R2);
            r = obj.ConvertCrop2VOI([R1c ; R2c]');
        end
        
        function r = GetMinFOV(obj, eType)
            % Get smallest Field of View (FOV)
            % NGR 2016 05 05
            vsI = obj.GetVolSeq(eType);
            r = obj.GetVolSeqMinFOV(vsI);
        end
        
        function vsO = CropVolSeq(obj, vsI, VOI_mm)
            % 4D crop
            % NGR 2016 05 05
            vsO = obj.CloneVolSeq(vsI);
            % Crop
            for iVol = 1:length(vsI)
                vsO(iVol).Crop(VOI_mm);
            end
        end
        
        function Vo = CropVol(obj, eT, iVol, VOI_mm)
            % 3D crop
            % NGR 2016 06 02
            vs = obj.GetVolSeq(eT);
            Vo = vs(iVol).Crop(VOI_mm);
        end
        
        function vsO = CropVolSeqVOI(obj, vsI, VOI_mm)
            % 4D crop VOI [x y z dx dy dz]
            % NGR 2016 06 22
            vsO = obj.CloneVolSeq(vsI);
            % Smallest commum VOI
            VOI = obj.InteresctFOVandVOI(vsI, VOI_mm);
            % Crop
            for iVol = 1:length(vsI)
                vsO(iVol).Crop(VOI);
            end
        end
        
        function r = InteresctFOVandVOI(obj, vsI, VOI_mm)
            % Itersection between FOV and VOI
            % NGR 2016 06 29
            FOV_mm = obj.GetVolSeqMinFOV(vsI);
            FOV = obj.ConvertVOI2Crop(FOV_mm);
            crop_mm = obj.ConvertVOI2Crop(VOI_mm);
            c = zeros(3,2);
            for i = 1:3
                c(i,1) = max(FOV(i,1),crop_mm(i,1));
                c(i,2) = min(FOV(i,2),crop_mm(i,2));
            end
            r = obj.ConvertCrop2VOI(c);
        end
        
        function Vout = ExtractVolVOI(obj, eT, iVol, VOI_mm)
            % 3D crop VOI [x y z dx dy dz]
            % NGR 2016 07 07
            vs = obj.GetVolSeq(eT);
            Vout = vs(iVol).ExtractVOI(VOI_mm);
        end
        
        function CropVolVOI(obj, eT, iVol, VOI_mm)
            % 3D crop VOI [x y z dx dy dz]
            % NGR 2016 06 22
            ext = obj.GetVolSeqLabel(eT);
            eval(['obj.vs' ext ...
                '(iVol).CropVOI(VOI_mm); ']);
        end
        
        function r = ConvertCrop2VOI(obj, crop_mm)
            % Covert VOI (Volume Of Interest) [x y z dx dy dz]
            % NGR 2016 06 29
            VOI = zeros(1,6);
            for i = 1:3
                VOI(i) = crop_mm(i,1);
                VOI(i+3) = crop_mm(i,2) - crop_mm(i,1);
            end
            r = VOI;
        end
        
        function r = ConvertVOI2Crop(obj, VOI_mm)
            % Covert VOI (Volume Of Interest) [x y z dx dy dz]
            % NGR 2016 06 23
            VOI = VOI_mm;
            crop_mm(1,:) = [VOI(1) VOI(1) + VOI(4)];
            crop_mm(2,:) = [VOI(2) VOI(2) + VOI(5)];
            crop_mm(3,:) = [VOI(3) VOI(3) + VOI(6)];
            r = crop_mm;
        end
        
        %-----------------------------------------------------------------%
        % Misc
        %-----------------------------------------------------------------%
        
        function Disp(obj, txt, Ntab)
            % Display text
            % NGR 2016 03 10
            if nargin<=2
                Ntab = 5;
            end
            for i=1:Ntab
                txt = ['    ' txt];
            end
            disp(txt);
        end
        
    end
end