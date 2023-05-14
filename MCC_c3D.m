% MCC_c3D Class for 3D processing
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center
% 5/30/2018


classdef MCC_c3D < handle
    
    properties
        % Who
        desc='Radiomics'% Description (eg. mpMRI)
        name            % Series description
        % Variable
        zFlipud = 0;    % Flip up/down in z direction, NGR 2018 04 11 v203
        data            % Volume (3D)
        dataRGB         % RGB data (countours)
        eType           % Volume type (Original, Resampled, ...)
        res_mm          % Spatial resolution in mm (array of 3 elements)
        z               % Slice location in mm
        TriggerTime     % Get "trigger time" (DICOM tag) vector
        filename        % DICOM file name
        metadata        % DICOM metadata
        bLoadSuccess    % DICOM read success
        modality        % Imaging modality (CT, MR, ...)
        % Location
        root            % Root directory
        InstituteID     % Institute ID
        DatabaseID      % Database ID
        CohortID        % Cohort ID
        PatientID       % Patient ID
        StudyDateID     % Study date ID
        ImagingModeID   % Imaging mode
        Address         % Adress (Radiomics language)
        % Version
        version = MCC_GetVersionFormated;
        signature=[MCC_GetArtistSignature ' - ' MCC_GetInstituteSignature];
    end
    
    methods
        
        %-----------------------------------------------------------------%
        % Constructor
        %-----------------------------------------------------------------%
        
        function obj = MCC_c3D(root, bFirstSlice)
            if nargin > 0
                if ischar(root)
                    % Root
                    obj.root = root;
                    cellTxt = strsplit(root,'\');
                    if ~isempty(strfind(root,'\Processed\'))
                        % Location
                        obj.name = cellTxt{end};
                        obj.InstituteID =  cellTxt{end-5};
                        obj.DatabaseID =  cellTxt{end-4};
                        obj.CohortID = cellTxt{end-3};
                        obj.PatientID = cellTxt{end-2};
                        obj.StudyDateID = cellTxt{end-1};
                        obj.ImagingModeID = cellTxt{end};
                        obj.Address = lower([ ...
                            obj.InstituteID '.' ...
                            obj.DatabaseID '.' ...
                            obj.CohortID '.' ...
                            obj.PatientID '.' ...
                            obj.StudyDateID '.' ...
                            obj.ImagingModeID
                            ]);
                    end
                    obj.ImagingModeID = cellTxt{end};
                    % Sort DICOM file based on instance number
                    d = dir(root);
                    d([d.isdir])= [];
                    N = length(d);
                    iFail = 0;
                    j = 1;
                    for i=1:N
                        fname = d(i).name;
                        if ~strcmp(fname,'.') && ~strcmp(fname,'..')
                            filename = [root '\' fname];
                            CD = cd;
                            try
                                info =  MCC_ReadMetadataDICOM(filename);
                                if bFirstSlice,
                                    obj.metadata{1} = info;
                                    obj.bLoadSuccess = 1;
                                    return;
                                end
                                pos = info.InstanceNumber;
                                sd(j) = ...
                                    struct('imagename',d(i).name, ...
                                    'instance',pos);
                                j = j+ 1;
                            catch
                                iFail = iFail + 1;
                            end
                            cd(CD);
                        end
                    end
                    if ~exist('sd')
                        obj.bLoadSuccess = 0;
                        return;
                    end
                    [~, order] = sort([sd(:).instance],'ascend');
                    fileNames = sd(order).';
                    % DICOM header
                    obj.root = root;
                    obj.metadata{1} = MCC_ReadMetadataDICOM(...
                        [obj.root '\' fileNames(1).imagename]);
                    if exist('obj.metadata{1}.SeriesDescription','var')==0  % 12/15/2021 In case SeriesDescription doesn't exist
                        obj.metadata{1}.SeriesDescription='';
                    end
                    obj.name = obj.metadata{1}.SeriesDescription;
                    obj.modality = obj.metadata{1}.Modality;
                    % Volume dimensions
                    dim(1) = obj.metadata{1}.Rows;
                    dim(2) = obj.metadata{1}.Columns;
                    dim(3) = length(fileNames);
                    % Axial resolution [mm]
                    try
                        obj.res_mm(2) = obj.metadata{1}.PixelSpacing(2); %% @@@
                        obj.res_mm(1) = obj.metadata{1}.PixelSpacing(1); %% @@@
                    catch
                        obj.metadata{1}.PixelSpacing(1) = 1;%0.9375;
                        obj.metadata{1}.PixelSpacing(2) = 1;%0.9375;
                        obj.res_mm(1) = obj.metadata{1}.PixelSpacing(1);
                        obj.res_mm(2) = obj.metadata{1}.PixelSpacing(2);
                    end
                    % Memory allocation
                    obj.data = zeros(dim,'uint32');         % 8/2/1018
%                     obj.data = zeros(dim,'uint16');
                    % Volume read (slices)
                    j = 1;
                    for i=1:dim(3)
                        obj.filename = fileNames(i).imagename;
                        try
                            obj.metadata{i} = ...
                                MCC_ReadMetadataDICOM( ...
                                [obj.root '\' fileNames(i).imagename]);
                            if exist('obj.metadata{i}.SeriesDescription','var')==0  % 12/15/2021 In case SeriesDescription doesn't exist
                                obj.metadata{i}.SeriesDescription='';
                            end
                            % mint Lesion managment, NGR 2016 08 22
                            if strfind(obj.metadata{i}.SeriesDescription,'mint Lesion') 
                                % RGB countour to mask
                                x0 = dicomread( ...
                                    [obj.root '\' fileNames(i).imagename]);
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
                                    obj.data(:,:,sd(i).instance) = mask;
                                    % RGB data
                                    obj.dataRGB{j} = dicomread(...
                                        [ ...
                                        obj.root '\' ...
                                        fileNames(i).imagename ...
                                        ]);
                                    j = j + 1 ;
                                end
                                if 0
                                    figure;
                                    subplot(221);
                                    imagesc(x0);axis image; axis off;
                                    subplot(223)
                                    imagesc(roi);axis image; axis off;
                                    subplot(224)
                                    imagesc(mask);axis image; axis off;
                                    colormap gray;
                                end
                            else
                                obj.data(:,:,i) = dicomread( ...
                                    [obj.root '\' fileNames(i).imagename]);
                            end
                        catch ME
                            MCC_writeLog([cd '\log.txt'],ME, obj.root);
                            obj.bLoadSuccess = 0;
                            return;
                        end
                    end
                    % Coronal/Sagittal resolution
                    obj.res_mm(3) = obj.metadata{1}.SliceThickness;
                    % 9/3/2021 One case SliceThickness was different to SpacingBetweenSlices
                    if isfield(obj.metadata{1},'SpacingBetweenSlices')  % 2/9/2022 check if defined
                        if obj.metadata{1}.SpacingBetweenSlices~=obj.metadata{1}.SliceThickness
                            obj.res_mm(3) = abs( ...
                                obj.metadata{2}.SliceLocation - ...
                                obj.metadata{1}.SliceLocation   ...
                                );
                        end
                    end
%                     try      % 2/5/2021 This was causing a small
%                     numerical issue
%                         z = abs( ...
%                             obj.metadata{2}.SliceLocation - ...
%                             obj.metadata{1}.SliceLocation   ...
%                             );
%                         if z
%                             obj.res_mm(3) = z;
%                         else
%                             obj.res_mm(3) = obj.metadata{1}.SliceThickness;
%                         end
%                     catch
%                         try
%                             obj.res_mm(3) = obj.metadata{1}.SliceThickness;
%                         catch
%                             obj.metadata{1}.SliceThickness = 4.9963;
%                             obj.res_mm(3) = obj.metadata{1}.SliceThickness;
%                         end
%                     end
                    % Update
                    obj.eType = MCC_eType.Original;
                    obj.PatientID = obj.metadata{1}.PatientID;
                    obj.StudyDateID = obj.metadata{1}.StudyDate;
                    obj.z = obj.GetSliceLocationVector;
                    slope = obj.z(end) - obj.z(1);
                    
                    % NGR 2018 03 20
                    if obj.zFlipud,
                        % Slope
                        for i=1:size(obj.data,3),
                            n = size(obj.data,3)-i+1;
                            data(:,:,n) = obj.data(:,:,i);
                            metadata{n} = obj.metadata{i};
                            z(n) = obj.z(i);
                        end
                        obj.data = data;
                        obj.z = z;
                        for i=1:size(obj.data,3),
                            obj.metadata{i} = metadata{i};
                        end
                    else
                        if slope<0     % slope>0, I changed it because it was flipping the images
                            %
                            for i=1:size(obj.data,3),
                                n = size(obj.data,3)-i+1;
                                data(:,:,n) = obj.data(:,:,i);
                                metadata{n} = obj.metadata{i};
                                z(n) = obj.z(i);
                            end
                            obj.data = data;
                            obj.z = z;
                            for i=1:size(obj.data,3),
                                obj.metadata{i} = metadata{i};
                            end
                        end
                    end
                    
                    obj.name = obj.ImagingModeID;
                    obj.bLoadSuccess = 1;
                    %  Curation
                    %obj.CurateDuplicateSlices; % @NGR
                    %x = obj.GetTriggerTimeVector
                else
                    error('Value must be a string (DICOM path name)');
                end
            end
        end
        
        %-----------------------------------------------------------------%
        % Flip volume in up/down direction
        %-----------------------------------------------------------------%
        
        function Flipud(obj)
            % Flip volume in up/down direction
            % NGR 2018 04 16

            for i=1:size(obj.data,3),
                n = size(obj.data,3)-i+1;
                data(:,:,n) = obj.data(:,:,i);
                metadata{n} = obj.metadata{i};
                z(n) = obj.z(i);
            end
            obj.data = data;
            obj.z = z;
            for i=1:size(obj.data,3),
                obj.metadata{i} = metadata{i};
            end
        end
        
        %-----------------------------------------------------------------%
        % Resampling
        %-----------------------------------------------------------------%
        
        function CurateDuplicateSlices(obj)
            % Curate duplicate slices (Prostate, WaterFat)
            % NGR 2016 07 04
            for i=1:length(obj.z)
                for j=1:length(obj.z)
                    dx(j) = obj.z(i)-obj.z(j);
                end
                idx(i,:) = find((dx==0));
            end
            if size(idx,2) > 1
                c = find(idx(:,1)==1);
                l = c(2)-1;
                splt = strsplit(obj.root,'\');
                for i = 1:size(idx,2)
                    R = idx(1:l,i);
                    % Slice #
                    for j = 1:length(R)
                        fld = fullfile(...
                            splt{1:end-1},[splt{end} '_' num2str(i)]);
                        if ~isdir(fld)
                            mkdir(fld);
                        end
                        scrFile = char(obj.metadata{R(j)}.Filename);
                        spltFN = strsplit(scrFile,'\'); filename = spltFN(end);
                        dstFile = char(fullfile(fld, filename));
                        copyfile(scrFile, dstFile);
                    end
                end
                % Remove source folder
                delete([obj.root '\*.*'])
                rmdir(obj.root);
            end
        end
        
        %-----------------------------------------------------------------%
        % Resampling
        %-----------------------------------------------------------------%
        
        function Vo = Resample(obj, res_mm_o)
            % Resampling a volume (res_mm_o: output resolution)
            % NGR 2016 06 22
            Vi = obj.data;
            res_mm_i = obj.res_mm;
            if (sum(res_mm_i - res_mm_o))
                % Echo
                MCC_Disp(obj.root,5);
                % Input
                Si = size(Vi);
                lX = (Si(2)-1)*res_mm_i(2);
                lY = (Si(1)-1)*res_mm_i(1);
                lZ = (Si(3)-1)*res_mm_i(3);
                vXi = linspace(0,lX, Si(2));
                vYi = linspace(0,lY, Si(1));
                vZi = linspace(0,lZ, Si(3));
                [Xi,Yi,Zi] = meshgrid(vXi,vYi,vZi);
                % Output
                So = obj.GetOutSize(res_mm_o);
                vXo = linspace(0, lX, So(2));
                vYo = linspace(0, lY, So(1));
                vZo = linspace(0, lZ, So(3));
                [Xo,Yo,Zo] = meshgrid(vXo,vYo,vZo);
                % 3-D spline interpolation (resampling)
                Vo = interp3(Xi,Yi,Zi,single(Vi),Xo,Yo,Zo,'spline');
                Vo = uint32(Vo);
                % Update
                obj.data = Vo;
                obj.res_mm = res_mm_o;
                obj.UpdateMetadata;
                obj.eType = MCC_eType.Resampled;
            end
        end
        
        function r = GetOutSize(obj, res_mm_o)
            % Get out volume size for resampling
            % NGR 2016 06 22
            Vi = obj.data;
            % Compute
            R = obj.GetResizeFactor(res_mm_o);
            r = zeros([1 3]);
            for j=1:3,
                r(j) = round(size(Vi,j)*R(j));
            end
        end
        
        function r = GetResizeFactor(obj, res_mm_o)
            % Get resize factor for resampling
            % NGR 2016 06 22
            res_mm_i = obj.res_mm;
            % Compute
            r = zeros([1 3]);
            for j=1:3,
                r(j) = res_mm_i(j) / res_mm_o(j);
            end
        end
        
        %-----------------------------------------------------------------%
        % Crop
        %-----------------------------------------------------------------%
        
        function Crop(obj, VOI_mm)
            % Crop VOI (Volume Of Interest) [x y z dx dy dz]
            % NGR 2016 06 22
            R = obj.GetCropVec(VOI_mm);
            obj.data = obj.data(R{1},R{2},R{3});
            obj.UpdateMetadata;
        end
        
        function Vo = ExtractVOI(obj, VOI_mm)
            % Extract VOI (Volume Of Interest) [x y z dx dy dz]
            % NGR 2016 07 08
            R = obj.GetCropVec(VOI_mm);
            Vo = obj.data(R{1},R{2},R{3});
        end
        
        function R = GetCropVec(obj, VOI_mm)
            % Get cropping box coordinates
            % VOI: [x y z dx dy dz]
            % NGR  2016 06 23
            c = VOI_mm([2 1 3]);    % coordinate
            d = VOI_mm([5 4 6]);    % dimension
            for i=1:3
                % Lower limit
                r(1) = round(c(i)/obj.res_mm(i) + size(obj.data,i)/2);
                if r(1)<1
                    r(1) = 1;
                end
                if r(1)>size(obj.data,i)
                    r(1) = size(obj.data,i);
                end
                % Upper limit
                %r(2) = round(d(i)/obj.res_mm(i)) + r(1);
                r(2) = round(d(i)/round(obj.res_mm(i),5)) + r(1);  %9/30/2020
                
                if r(2)<1
                    r(2) = 1;
                end
                if r(2)>size(obj.data,i)
                    r(2) = size(obj.data,i);
                end
                R{i} = r(1):r(2);
            end
        end
        
        function Io = CropImage(obj, eType, VOI_mm, iPlane)
            % Crop image
            % NGR 2016 04 15
            R = obj.GetCropVec(VOI_mm);
            switch eType
                case MCC_eView.Axial
                    R{3} = iPlane;
                case MCC_eView.Sagittal
                    R{2} = iPlane;
                case MCC_eView.Coronal
                    R{1} = iPlane;
            end
            Io = squeeze(obj.data(R{1},R{2},R{3}));
        end
        
        %-----------------------------------------------------------------%
        % Crop
        %-----------------------------------------------------------------%
        
        function FillVOI(obj, Vi, VOI_mm)
            % Fill a VOI with an input volume
            % 2016 06 23
            R = obj.GetCropVec(VOI_mm);
            obj.data(R{1},R{2},R{3}) = Vi;
        end
        
        %-----------------------------------------------------------------%
        % Misc.
        %-----------------------------------------------------------------%
        
        function r = IsVolEmpty(obj)
            % Is volume data empty ?
            % NGR 2016 03 04
            r = 0;
            if isempty(obj.data)
                r = 1;
            end
        end
        
        function r = GetVoxelSize_mm3(obj)
            % Get voxel size in mm^3
            r = obj.res_mm(1) * obj.res_mm(2) * obj.res_mm(3);
        end
        
        function r = GetNbSlices(obj)
            % Get number of slices
            r = size(obj.data,3);
        end
        
        function r = GetSliceLocation(obj, iSlice)
            % Get slice location
            % NGR 2016 03 24
            try
                r = obj.metadata{iSlice}.SliceLocation;
            catch
                r = -(iSlice-1)*obj.res_mm(3);
            end
        end
        
        function r = CheckSliceContiguousity(obj)
            % Check unifromity of slice location
            % NGR 2016 03 24
            sl = obj.GetSliceLocationVector;
            val = std(diff(sl));
            if abs(val) < 1e-3
                r = 1;
            else
                r = 0;
            end
        end
        
        function r = GetSliceLocationVector(obj)
            % Get slice location vector
            % NGR 2016 03 30
            N = obj.GetNbSlices;
            r = zeros(1,N);
            for i=1:N
                r(i) = obj.GetSliceLocation(i);
            end
            obj.z = r;
        end
        
        function r = GetTriggerTimeVector(obj)
            % Get trigger time vetor vector
            % NGR 2016 07 27
            N = obj.GetNbSlices;
            %r = zeros(1,N);
            for i=1:1
                r(i) = obj.metadata{i}.TriggerTime;
            end
            obj.TriggerTime = r;
        end
        
        function r = GetSize(obj)
            % Get number of slices
            r = size(obj.data);
        end
        
        function UpdateMetadata(obj)
            % Update DICOM metadata
            % NGR 2016 02 08
            Nz = obj.GetNbSlices;
            slope = obj.z(end) - obj.z(1);
            for i = 1:Nz
                if i>length(obj.metadata)
                    obj.metadata{i} = obj.metadata{1};
                end
                % Description
                obj.metadata{i}.SeriesDescription = obj.name;
                % Dimensions
                obj.metadata{i}.Rows    = size(obj.data,1);
                obj.metadata{i}.Height  = obj.metadata{i}.Rows;
                obj.metadata{i}.Columns = size(obj.data,2);
                obj.metadata{i}.Width   = obj.metadata{i}.Columns;
                % Spatial resolution
                obj.metadata{i}.PixelSpacing(1) = obj.res_mm(1);
                obj.metadata{i}.PixelSpacing(2) = obj.res_mm(2);
                obj.metadata{i}.SliceThickness  = obj.res_mm(3);
                % Slice location
                if slope>0
                    obj.metadata{i}.SliceLocation = ...
                        obj.z(1)+(i-1)*obj.res_mm(3);
                else
                    obj.metadata{i}.SliceLocation = ...
                        obj.z(1)-(i-1)*obj.res_mm(3);
                end
            end
        end
        
        %-----------------------------------------------------------------%
        % Files
        %-----------------------------------------------------------------%
        
        function obj = Write(obj, pathIn)
            % Write images as DICOM files
            % NGR 2016 02 08
            
            % Base file name
            bfn = [ ...
                MCC_GetDateFileFormat '_' ...
                'v' num2str(MCC_GetVersion) '_' ...
                char(obj.eType) '_' obj.PatientID '_' obj.StudyDateID '_s'];

            % Slice loop
            for i =1:obj.GetNbSlices
                
                % DICOM
                if i<10
                    s = ['00' num2str(i)];
                end
                if i>=10 && i<100
                    s = ['0' num2str(i)];
                end
                if i>=100 && i<1000
                    s = num2str(i);
                end
                dstF = [pathIn '\' bfn s '.dcm'];
                dstF = [pathIn '\' s '.dcm']; % Prostate, @NGR 2016 06 29
                obj.metadata{i}.InstanceNumber = i;
                %I = uint16(obj.data(:,:,i));
                I = uint32(obj.data(:,:,i));    % 8/2
                
                %obj.metadata{i}.PixelSpacing = [ 4.99938086 4.99938086 ];
                %obj.metadata{i}.SliceThickness =  1.6406;
                %obj.metadata{i}.SliceLocation = 0;
                
                try
%                     dicomwrite(I, dstF, obj.metadata{i});
                    obj.metadata{i}.SpacingBetweenSlices=obj.metadata{i}.SliceThickness;
                    dicomwrite(I, dstF, obj.metadata{i}, 'MultiframeSingleFile',false,'CreateMode','copy');  % 7/5/2018
                catch
                    try
                        dicomwrite(I, dstF, obj.metadata{i}, 'MultiframeSingleFile',false,'CreateMode','Create');  % 7/7/2020
                    catch
                        try    %12/6/2020
                            dicomwrite(I, dstF, ...
                                'Format', obj.metadata{i}.Format, ...
                                'Width', obj.metadata{i}.Width, ...
                                'Height', obj.metadata{i}.Height, ...
                                'BitDepth', 16, ...
                                'StudyDate', obj.metadata{i}.StudyDate, ...
                                'AcquisitionTime', obj.metadata{i}.AcquisitionTime, ...
                                'Modality', obj.metadata{i}.Modality, ... % No PixelSpacing if removed !
                                'Manufacturer', ['Radiomics Language (' obj.version ')'], ...
                                'SeriesDescription', obj.metadata{i}.SeriesDescription,...
                                'PatientID', obj.metadata{i}.PatientID, ...
                                'SliceThickness', obj.metadata{i}.SliceThickness, ...
                                'StudyInstanceUID', obj.metadata{i}.StudyInstanceUID, ...
                                'SeriesInstanceUID', obj.metadata{i}.SeriesInstanceUID, ...
                                'InstanceNumber', obj.metadata{i}.InstanceNumber, ...
                                'SeriesNumber', obj.metadata{i}.SeriesNumber, ...
                                'FrameOfReferenceUID', obj.metadata{i}.FrameOfReferenceUID, ...
                                'SliceLocation', obj.metadata{i}.SliceLocation, ...
                                'Rows', obj.metadata{i}.Rows, ...
                                'Columns', obj.metadata{i}.Columns, ...
                                'PixelSpacing', obj.metadata{i}.PixelSpacing, ...
                                'RepetitionTime', obj.metadata{i}.RepetitionTime, ...
                                'FlipAngle', obj.metadata{i}.FlipAngle, ...
                                'TriggerTime', obj.metadata{i}.TriggerTime, ...
                                'AcquisitionTime', obj.metadata{i}.AcquisitionTime ...
                                );
                        catch
                            dicomwrite(I, dstF, ...
                                'Format', obj.metadata{i}.Format, ...
                                'Width', obj.metadata{i}.Width, ...
                                'Height', obj.metadata{i}.Height, ...
                                'BitDepth', 16, ...
                                'StudyDate', obj.metadata{i}.StudyDate, ...
                                'AcquisitionTime', obj.metadata{i}.AcquisitionTime, ...
                                'Modality', obj.metadata{i}.Modality, ... % No PixelSpacing if removed !
                                'Manufacturer', ['Radiomics Language (' obj.version ')'], ...
                                'SeriesDescription', obj.metadata{i}.SeriesDescription,...
                                'PatientID', obj.metadata{i}.PatientID, ...
                                'SliceThickness', obj.metadata{i}.SliceThickness, ...
                                'StudyInstanceUID', obj.metadata{i}.StudyInstanceUID, ...
                                'SeriesInstanceUID', obj.metadata{i}.SeriesInstanceUID, ...
                                'InstanceNumber', obj.metadata{i}.InstanceNumber, ...
                                'SeriesNumber', obj.metadata{i}.SeriesNumber, ...
                                'FrameOfReferenceUID', obj.metadata{i}.FrameOfReferenceUID, ...
                                'SliceLocation', obj.metadata{i}.SliceLocation, ...
                                'Rows', obj.metadata{i}.Rows, ...
                                'Columns', obj.metadata{i}.Columns, ...
                                'PixelSpacing', obj.metadata{i}.PixelSpacing ...
                                );
                        end
                    end
                end
            end
        end
        
        function obj = WriteJPG(obj, pathIn)
            % Write images as DICOM files
            % NGR 2016 02 08
            
            % Base file name
            bfn = [char(obj.eType) '_' obj.PatientID '_' obj.StudyDateID '_'];
            % Normalization
            v = single(obj.data);
            v = 255*v/max(v(:));
            % JPG
            for i =1:obj.GetNbSlices
                dstF = ...
                    [MCC_GetDateFileFormat '_' ...
                    pathIn '\' bfn num2str(i) '.jpg'];
                imwrite(uint8(v(:,:,i)),dstF,'jpg','Comment',dstF);
            end
        end
        
        function obj = ConvertMasks(obj)
            % Convert PGM images (mask format) into DICOM ones
            % NGR 2016 04 18
            d = dir(obj.root);
            d(~[d.isdir])= [];
            N = length(d);
            x = strsplit(obj.root,'\');
            root = strjoin(x(1:end-1),'\');
            MCC_Disp('Convert masks ...');
            for i=1:N
                fname = d(i).name;
                if ~strcmp(fname,'.') && ~strcmp(fname,'..')
                    fld = [obj.root '\' fname];
                    MCC_Disp(fld,5);
                    for i = 1:obj.GetNbSlices
                        I = imread([fld '\' num2str(i) '.pgm']);
                        obj.data(:,:,i) = 255*(I>0);
                    end
                    dstFld = [root '\' fname];
                    if ~isdir(dstFld)
                        mkdir(dstFld);
                    end
                    MCC_Disp(dstFld,5);
                    obj.Write(dstFld);
                end
            end
            MCC_Disp(' ');
        end
        
        %-----------------------------------------------------------------%
        % Overload
        %-----------------------------------------------------------------%
        
        % Test
        function r = plus(o1,o2)
            r = [o1.data] + [o2.data];
        end
        
    end
end