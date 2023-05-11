function Main(PatientName, Folder, Options)
% Main file to process one patient
    
[DoDWI,CurateDWI,BiasCorrectDWI,AnnotateDWI,...
    GlobalRegisterDWI,LocalRegisterDWI,roiFactorDW,GetParDWI,GlobalDateRegisterDWI,...
    LocalDateRegisterDWI,roiFactorLDW,UseSegmentationForLocalRegDWI,DoT1W,...
    CurateT1W,BiasCorrectT1,AnnotateT1W,GlobalRegisterT1W,SpetialRegT1,LocalRegisterT1W,...
    roiFactorT1,IncludeArteryLR,GetT1,CorrectFA,WeightT1,AirTh,ScaleT1Ims,RefTiss,gammaT1,gammaM0,...
    solver,CorrectBiasBetweenSlices,AssumeT1s,GetConc,ContrastBrand,ScaleMo,AssumeT1sC,AssumeT1Artery,GetPar,...
    GlobalDateRegisterT1,LocalDateRegisterT1,roiFactorLT1,UseSegmentationForLocalRegT1,...
    DoT1WDWIReg,roiFactorLT1DW,RotateViewDW]=LoadMainOptions(Options);


FolderDWI='(^(?=.*diff)(?!.*(Apparent|Exponential|PosDisp)).*)|(Axial(.*)DWI)|(Ax(.*)DWI)';
Folder_T1=['(PRE T1_vibe)|((AX 3D)(((?!(SPGR|DCE)).)*$))|'...
    '(3d Non CE RSSG)|(^DCE((?!(WITH CONTRAST)).)*$)|(Ax 3D SPGR T1)|(TI WI(.*)Dynamic)|(^(?=.*fl3d_ce-1-fa)(?!.*(DYNAMIC|PosDisp)).*)'];

Folder_DCE=['(T1_vibe)|(AX 3D DCE 30 \+C)|'...
    '(3d PRE\/POST  RSSG  30fa)|(^DCE 30 WITH CONTRAST)|(Ax 3D SPGR T1 30)|(TI WI 30 DCE \+C)|(^(?=.*DYNAMIC_fl3d_ce-1-fa)(?!.*(PosDisp)).*)'];

if exist([Folder PatientName])~=0
    Dates=AdjustDirVariable(dir([Folder PatientName])); 
    Dates=Dates([Dates(:).isdir]);
    Dates={Dates.name};
end

WriteFolder='C:\Users\Andres\Desktop\Moffitt\Images\'; 

if exist('Dates','var')==0
    [~,~,RAW]=xlsread('C:\Users\Andres\Desktop\Moffitt\MATLAB\Project3\TableCurationV5.xls');
    aux=RAW(2:end,1);
    aux(cell2mat(cellfun(@(x) isnumeric(x),aux,'UniformOutput',false)))={' '};  
    Ind0=find(contains(aux,PatientName));
    IndL=find(~contains(aux(Ind0+1:end),' '),1)+Ind0-1;
    if ~isempty(Ind0) & isempty(IndL)
        IndL=numel(aux);
    end    
    Dates={RAW{Ind0+1:IndL+1,6}};   
    Dates=cellfun(@num2str, Dates,'UniformOutput',false);
    Dates(cell2mat(cellfun(@(x) strcmp(x,'NaN'),Dates,'UniformOutput',false)))=[];
end

Path=genpath('C:\Users\Andres\Desktop\Moffitt\MATLAB\');
addpath(Path);

[~,PCNAME]=system('hostname');
PCNAME=PCNAME(1:end-1);

close all force
warning off

% Read files
if exist([WriteFolder PatientName '\FilesInfo.mat'],'file')==0 | (CurateDWI & DoDWI) |...
        (CurateT1W & DoT1W)
    I=1;
    FilesFolderNames=[]; Infos=[]; SeriesDescription=[];
    disp('Reading files')
    for Date=Dates
        %[FilesFolderNames{I},Infos{I},SeriesDescription{I}]=ReadPatientFilesP3([Folder PatientName '_' Date{1} '\**\*']) ; 
        [FilesFolderNames{I},Infos{I},SeriesDescription{I}]=ReadPatientFilesP3([Folder PatientName '\' Date{1} '\**\*']); 
        I=I+1;
    end
    mkdir([WriteFolder PatientName])
    save([WriteFolder PatientName '\FilesInfo.mat'],'FilesFolderNames','Infos','SeriesDescription'); 
else
    load([WriteFolder PatientName '\FilesInfo.mat']);
end


if DoDWI
    % Curate DW-MRI    
    if CurateDWI
        try, rmdir([WriteFolder PatientName '\DWI\'],'s');  end 
        
        I=1;
        for Date=Dates
            CurateDataDWIP3V3(WriteFolder, Infos{I},SeriesDescription{I},PatientName,Date{1},FolderDWI,BiasCorrectDWI);
            [Templates{I},TemplateDate{I},LocalRegTemplate{I}]=GetDW_RegistrationTemplatesP3([WriteFolder PatientName '\DWI\' Date{1}]);
            I=I+1;
        end
        if exist([WriteFolder PatientName '\DWI\'])~=0
            save([WriteFolder PatientName '\DWI\Templates.mat'], 'LocalRegTemplate', 'TemplateDate', 'Templates');
        end
    end       
    % Global Register DW-MRI  
    if GlobalRegisterDWI   
        if ~strcmp(PCNAME,'PHYTM7M2') 
%             delete(gcp('nocreate'))
%             parpool(5);
        end
        disp('Global Reg. ADC')
        I=1;
        if exist([WriteFolder PatientName '\DWI\\Templates.mat'])~=0
            load([WriteFolder PatientName '\DWI\\Templates.mat']);
            for Date=Dates  
                if ~isempty(TemplateDate{I})
                    obj = MCC_cPatient([WriteFolder PatientName '\DWI\' Date{1}]);
                    obj.Register_Date(Templates{I},TemplateDate{I}, {'none','rigid','similarity','affine'});
                    clear obj
                    OrganizeRegisteredDWIP3([WriteFolder PatientName '\DWI\' Date{1}],WriteFolder) 
                end
                I=I+1;
            end
        end
    end
    % Copy masks from network
    if UseNetworkMasks
        disp('Copying ADC masks from Network')
        for Date=Dates
            if exist(['M:\dept\IRAT_Research\Andres Arias\PEGPH20 2020 Results\Results\' PatientName '\DWI\Registered\'...
                    Date{1} '\Global\Masks.mat'])~=0 
                copyfile(['M:\dept\IRAT_Research\Andres Arias\PEGPH20 2020 Results\Results\' PatientName '\DWI\Registered\'...
                    Date{1} '\Global\Masks.mat'],  [WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\Masks.mat'])           
%                 if InvertMasks
%                     InvertMarsk([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\Masks.mat'])
%                 end
            end
        end
    end 
    % Annotate DWI
    if AnnotateDWI
        disp('Annotation ADC')
        if exist([WriteFolder PatientName '\DWI\\Templates.mat'])~=0
            load([WriteFolder PatientName '\DWI\Templates.mat']);
            I=1;
            for Date=Dates 
                if ~isempty(TemplateDate{I})
                    f=0;
                    while f==0
                        try                                      
                            %MakeAnnotationsTemplateP3([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\'],I,LocalRegTemplate{I});
                            MakeAnnotationsTemplateV2P3([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\'],I,LocalRegTemplate{I});
                            f=1;
                        catch
                            disp('Error try again')
                            pause
                        end
                    end
                end
                I=I+1;
            end
        end
    end
    % Local registration 
    if LocalRegisterDWI
        if ~strcmp(PCNAME,'PHYTM7M2') 
%             delete(gcp('nocreate'))
%             parpool(5);
        end    
        disp('Local Reg. ADC')
        I=1;
        if exist([WriteFolder PatientName '\DWI\Templates.mat'])~=0
            load([WriteFolder PatientName '\DWI\Templates.mat']);
            for Date=Dates   
                if exist([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\'])~=0
%                     VisualizeLocalRegistrationDWIP3([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\'],...
%                         LocalRegTemplate{I},roiFactorDW); % Before all images were registered to template, now only b0s
                    VisualizeLocalRegistrationDWIP3V2([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\'],...
                        Templates{I},LocalRegTemplate{I},roiFactorDW)
                end
                I=I+1;
            end
            fid=fopen([WriteFolder PatientName '\DWI\Registered\roiFactorDW.txt'],'w');
            fprintf(fid, num2str(roiFactorDW));
            fclose(fid);
        end
    end
    % Obtain ADC maps
    if GetParDWI  
        disp('Computing ADC maps')
        I=1;
        for Date=Dates
            ADCMAPP3V2([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Local\'])            
            I=I+1;
        end
    end
    % Global Date registration
    if GlobalDateRegisterDWI 
        disp('Global Date registration ADC')
        try, rmdir([WriteFolder PatientName '\DWI\DateRegistered\'],'s'); end
        load([WriteFolder PatientName '\DWI\Templates.mat']);
        I=1;
        Ind=logical(zeros(size(Dates)));
        for Date=Dates
            if ~isempty(LocalRegTemplate{I}) &...
                    exist([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Local\ADC.mat'])
                ADC=load([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Local\ADC.mat']);
                if ~isempty(ADC.ADCav)
                    MoveFilesForDateRegistrationDW([WriteFolder PatientName],Date{1},LocalRegTemplate{I})
                    Ind(I)=1;
                end
            end
            I=I+1;
        end
        obj = MCC_cPatient([WriteFolder PatientName '\DWI\DateRegistered\']);
        obj.Register_Date(LocalRegTemplate(Ind),obj.content{1},...
            {'none','rigid','similarity','affine'});
        MetaCorrectionDW([WriteFolder PatientName],Dates(Ind))
    end
    % Local Date registration
    if LocalDateRegisterDWI 
        disp('Local Date registration ADC')
        try, rmdir([WriteFolder PatientName '\DWI\DateRegisteredLocal\'],'s'); end
        load([WriteFolder PatientName '\DWI\Templates.mat']);
        I=1;
        Inds=logical(zeros(size(Dates)));
        for Date=Dates
            if ~isempty(LocalRegTemplate{I}) &...
                    exist([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Local\ADC.mat'])
                ADC=load([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Local\ADC.mat']);
                if ~isempty(ADC.ADCav)
                    MoveFilesForLocalDateRegistrationDW([WriteFolder PatientName],Date{1},LocalRegTemplate{I})
                    Inds(I)=1;
                end
            end
            I=I+1;
        end
        Ind=find(Inds,1);
        VOI=GenerateVOIforLocalDateRegistration([WriteFolder PatientName],...
            [WriteFolder PatientName '\DWI\Registered\' Dates{Ind} '\Global\'],...
            [WriteFolder PatientName '\DWI\DateRegisteredLocal\' Dates{Ind}], LocalRegTemplate{Ind},roiFactorLDW);
        obj = MCC_cPatient([WriteFolder PatientName '\DWI\DateRegisteredLocal\']);
        if UseSegmentationForLocalRegDWI
            obj.Register_Date_Local(cellstr(repmat('TumorSegmentationDW',[sum(Inds),1]))',Dates{Ind},...
                {'none','rigid','similarity','affine'},VOI);
        else
            obj.Register_Date_Local(LocalRegTemplate(Inds),Dates{Ind},...
                {'none','rigid','similarity','affine'},VOI);
        end
        fid=fopen([WriteFolder PatientName '\DWI\DateRegisteredLocal\roiFactorLDW.txt'],'w');
        fprintf(fid, [num2str(UseSegmentationForLocalRegDWI) ' ' num2str(roiFactorLDW)]);
        fclose(fid);
    end
end

if DoT1W  
    % Curate T1
    if CurateT1W
        try, rmdir([WriteFolder PatientName '\T1W\'],'s');  end
        I=1;
        for Date=Dates
            CurateT1WP3V2(WriteFolder,Infos{I},SeriesDescription{I},PatientName,Date{1},Folder_T1,BiasCorrectT1); 
            CurateDCEP3V2(Infos{I},SeriesDescription{I},PatientName,Date{1},Folder_DCE,WriteFolder,BiasCorrectT1);
            I=I+1;
        end
    end
    % Copy masks from network
    if UseNetworkMasks
        disp('Copying T1 masks from network')
        for Date=Dates            
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0 &&...
                    exist(['M:\dept\IRAT_Research\Andres Arias\PEGPH20 2020 Results\Results\' PatientName '\T1W\' Date{1} '\Masks.mat'])
                copyfile(['M:\dept\IRAT_Research\Andres Arias\PEGPH20 2020 Results\Results\' PatientName '\T1W\' Date{1} '\Masks.mat'],...
                    [WriteFolder PatientName '\T1W\' Date{1} '\Masks.mat'])
%                 if InvertMasks
%                     InvertMarsk([WriteFolder PatientName '\T1W\' Date{1} '\Masks.mat'])
%                 end
            end
        end
    end 
    % Annotate T1
    if AnnotateT1W 
        disp('Annotate T1 masks')
        Ind=cellfun(@(x) exist([WriteFolder PatientName '\T1W\' x '\'])~=0, Dates);
        I=1;
        for Date=Dates(Ind)  
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0  %DCEav\'])~=0  %9/14/2020
                f=0;
                while f==0
                    try
%                         MakeAnnotationsTemplateP3([WriteFolder PatientName '\T1W\' Date{1} '\'],I);
                        MakeAnnotationsTemplateV2P3([WriteFolder PatientName '\T1W\' Date{1} '\'],I);
                        f=1;
                    catch 
                        disp(' ');disp('Error try again')
                        pause
                    end
                end
            end
            I=I+1;
        end
    end
    % Global Registration T1
    if GlobalRegisterT1W
        disp('Global Reg. T1')
        if exist([WriteFolder PatientName '\T1W\'])~=0
            I=1;
            if ~exist([WriteFolder PatientName '\T1W\ImageTemplates.mat'])
                for Date=Dates
                    if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0
                        ImageTemplate{I}=GetT1RegistrationTemplateP3([WriteFolder PatientName '\T1W\' Date{1} '\']);  
                    end
                    I=I+1;
                end
                ImageTemplate(cellfun(@isempty,ImageTemplate))=[];
                save([WriteFolder PatientName '\T1W\ImageTemplates.mat'], 'ImageTemplate');
            else
                load([WriteFolder PatientName '\T1W\ImageTemplates.mat']);
            end
            if ~strcmp(PCNAME,'PHYTM7M2') 
%                 delete(gcp('nocreate'))
%                 parpool(4); %6
            end
            try, rmdir([WriteFolder PatientName '\T1W\Registered\'],'s');  end
            try, rmdir([WriteFolder PatientName '\T1W\DateRegistered\'],'s');  end
            try, rmdir([WriteFolder PatientName '\T1W\DateRegisteredLocal\'],'s');  end            
            obj = MCC_cPatient([WriteFolder PatientName '\T1W\']);
            for I=1:numel(obj.StudyDate)        % No register average DCE
                obj.StudyDate(I).content(strcmp(obj.StudyDate(I).content,'DCEav'))=[];
            end          
            obj.Register(ImageTemplate,'elastic');           
        end
    end
    % Special Global registration for T1w images
    if SpetialRegT1
        disp('Spetial registration T1w')
        SpetialGlobalRegistrationV3(WriteFolder,Dates,PatientName)
        fid=fopen([WriteFolder PatientName '\T1W\Registered\SpetialRegT1.txt'],'w');
        fprintf(fid, 'True');
        fclose(fid);
    end
    % Local registration 
    if LocalRegisterT1W
        if ~strcmp(PCNAME,'PHYTM7M2') 
            %delete(gcp('nocreate'))
            %parpool(5);
        end
        disp('Local registration T1w')
        load([WriteFolder PatientName '\T1W\ImageTemplates.mat']);
        I=1;
        for Date=Dates    
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0
                VisualizeLocalRegistrationT1P3([WriteFolder PatientName '\T1W\' Date{1} '\Processed\Registered'],...
                    ImageTemplate{I},roiFactorT1, IncludeArteryLR);
                 I=I+1;
            end           
        end
        fid=fopen([WriteFolder PatientName '\T1W\Registered\roiFactorT1.txt'],'w');
        fprintf(fid, ['Artery: ' num2str(IncludeArteryLR) char(10) num2str(roiFactorT1)]);
        fclose(fid);
    end
    % Obtain T1 maps
    if GetT1
        if ~strcmp(PCNAME,'PHYTM7M2') 
%             delete(gcp('nocreate'))
%             parpool(10);
        end
        disp('Computing T1 maps')
        for Date=Dates
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0
                T1mapP3([WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\'],[WriteFolder PatientName '\T1W\' Date{1} '\'],...
                    RefTiss, ScaleT1Ims, gammaT1, gammaM0, CorrectBiasBetweenSlices,AssumeT1s,AirTh,solver,Date{1},Dates,...
                    WriteFolder,PatientName,CorrectFA,WeightT1); 
            end
        end
        fid=fopen([WriteFolder PatientName '\T1W\Registered\T1mapppingParameters.txt'],'w');
        fprintf(fid, ['CorrectFA: ' num2str(CorrectFA) char(10) 'Scale: ' num2str(ScaleT1Ims) char(10)...
            'RefTiss: ' RefTiss char(10) 'Gammas: ' num2str(gammaT1) ', ' num2str(gammaM0) char(10)...
            'BiasCorrection: ' num2str(CorrectBiasBetweenSlices) char(10) 'AssumeT1: ' num2str(AssumeT1s) char(10)...
            'AirTh: ' num2str(AirTh) char(10) 'Solver: ' solver]);
        fclose(fid);
    end
    % Get Concentration maps
    if GetConc
        disp('Computing C(t)')
        load([WriteFolder PatientName '\T1W\ImageTemplates.mat']);
        I=1;
        for Date=Dates      
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\DCE_t=1\'])~=0   % 9/15/2020
                [r,ContrastBrand]=Get_rContrast([WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\'],ContrastBrand);
                [Cons_NW_NL,Conds_NW_NL,Times]=GetConcentrationsP3([WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\'],...
                                    [WriteFolder PatientName '\T1W\' Date{1} '\'],r,ImageTemplate{I}, ScaleMo, AssumeT1sC, AssumeT1Artery);                 
                save([WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\Concentrations.mat'], 'Cons_NW_NL', 'Conds_NW_NL',...
                        'Times','r');
                I=I+1;
            end
        end
        save([WriteFolder PatientName '\T1W\Registered\AssumeT1sC.mat'],'AssumeT1sC')
        fid=fopen([WriteFolder PatientName '\T1W\Registered\ConcentrationParameteres.txt'],'w');
        fprintf(fid,['ContrastBrand: ' ContrastBrand char(10) 'r: ' num2str(r) char(10) 'AssumedT1: ' num2str(AssumeT1sC)...
            char(10) 'AssumedT1 Artery: ' num2str(AssumeT1Artery)]);
        fclose(fid);
    end
    % Get DCE parametric Maps
    if GetPar
       if ~strcmp(PCNAME,'PHYTM7M2') 
            delete(gcp('nocreate'))
%             parpool(5);
       end
        disp('AUC maps')
        for Date=Dates 
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\DCE_t=1\'])~=0
                AUC_NW=Get90sAUCP3([WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\'],...
                    [WriteFolder PatientName '\T1W\' Date{1} '\']);   
                save([WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\90sAUCs.mat'],...
                    'AUC_NW');
                close all force
            end
        end
        disp('Tofts maps')
        for Date=Dates  
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\DCE_t=1\'])~=0
                Folder2=[WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\'];
                [Vp_NW_NL,Ve_NW_NL,Kt_NW_NL]=GetKinetcisParametersP3(WriteFolder,Folder2,...
                    [WriteFolder PatientName '\T1W\' Date{1} '\']);   
                save([WriteFolder PatientName '\T1W\Registered\' Date{1} '\Local\Perfussion_Parameters_Maps.mat'],...
                    'Vp_NW_NL','Ve_NW_NL','Kt_NW_NL');
                %close all force
            end
        end
    end
    % Global Date registration
    if GlobalDateRegisterT1
        disp('Global date registration')
        try,rmdir([WriteFolder PatientName '\T1W\DateRegistered\'],'s');end
        load([WriteFolder PatientName '\T1W\ImageTemplates.mat']);
        Ind=cellfun(@(x) exist([WriteFolder PatientName '\T1W\' x '\'])~=0, Dates);       
        I=1;
        for Date=Dates
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0
                %MoveFilesForDateRegistrationT1([WriteFolder PatientName],Date{1},ImageTemplate{I})   
                MoveFilesForDateRegistrationT1V2([WriteFolder PatientName],Date{1},ImageTemplate{I},Dates(Ind),ImageTemplate) 
                I=I+1;
            end
        end
        obj = MCC_cPatient([WriteFolder PatientName '\T1W\DateRegistered\']);
        obj.Register_Date(ImageTemplate,obj.content{1},...
            {'none','rigid','similarity','affine'});
        MetaCorrectionT1([WriteFolder PatientName],obj.content)
    end
    % Local Date registration
    if LocalDateRegisterT1 
        disp('Local date registration')
        try,rmdir([WriteFolder PatientName '\T1W\DateRegisteredLocal\'],'s');end
        load([WriteFolder PatientName '\T1W\ImageTemplates.mat']);
        I=1;
        for Date=Dates
            if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0
                MoveFilesForLocalDateRegistrationT1([WriteFolder PatientName],Date{1},ImageTemplate{I})
                if I==1
                    DateI=Date{1};
                end                    
                I=I+1;
            end
        end
        VOI=GenerateVOIforLocalDateRegistration([WriteFolder PatientName],...
            [WriteFolder PatientName '\T1W\' DateI '\'],...
            [WriteFolder PatientName '\T1W\DateRegisteredLocal\' DateI], ImageTemplate{1},roiFactorLT1);
        obj = MCC_cPatient([WriteFolder PatientName '\T1W\DateRegisteredLocal\']);
        if UseSegmentationForLocalRegT1
            obj.Register_Date_Local(cellstr(repmat('TumorSegmentationT1',[numel(ImageTemplate),1]))',DateI,...
                {'none','rigid','similarity','affine'},VOI);
        else
            obj.Register_Date_Local(ImageTemplate,DateI,...
                {'none','rigid','similarity','affine'},VOI);
        end
        fid=fopen([WriteFolder PatientName '\T1W\DateRegisteredLocal\RegistrationParameters.txt'],'w');
        fprintf(fid, ['ROI Factor: ' num2str(roiFactorLT1) char(10) 'Use Seg.: ' num2str(UseSegmentationForLocalRegT1)]);
        fclose(fid);
    end
end
% T1-DW registration
if DoT1WDWIReg
    disp('T1-ADC registration')
    try,rmdir([WriteFolder PatientName '\T1-DW Registered\'],'s');end
    load([WriteFolder PatientName '\T1W\ImageTemplates.mat']);
    load([WriteFolder PatientName '\DWI\Templates.mat']);
    I=1;
    Inds=logical(zeros(size(Dates)));
    for Date=Dates
        if ~isempty(LocalRegTemplate{I}) 
            ADC=load([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Local\ADC.mat']);
            if ~isempty(ADC.ADCav)
                Inds(I)=1;
            end
        end
        I=I+1;
    end
    Ind=find(Inds,1);
    DatesDCE=cell(0);
    for Date=Dates
        if exist([WriteFolder PatientName '\T1W\' Date{1} '\'])~=0
            DatesDCE{end+1}= Date{1};    
        end
    end
    IndDCE=find(contains(DatesDCE,Dates{Ind}),1);
    MoveFilesForDWT1RegistrationV2([WriteFolder PatientName],Dates(Inds),...
        LocalRegTemplate{Ind},ImageTemplate{IndDCE},RotateViewDW)
    obj = MCC_cPatient([WriteFolder PatientName '\T1-DW Registered\']);
    obj.Register_Date({LocalRegTemplate{Ind} ImageTemplate{IndDCE}},'T1W',...
       {'none'}); %{'none','rigid','similarity','affine'});
    MoveFilesForDWT1LocalRegistration([WriteFolder PatientName],Dates(Inds),...
        LocalRegTemplate{Ind},ImageTemplate{IndDCE})
    VOI=GenerateVOIforLocalDateRegistration([WriteFolder PatientName],...
        [WriteFolder PatientName '\T1W\' Dates{Ind} '\'],...
        [WriteFolder PatientName '\T1W\DateRegisteredLocal\' Dates{Ind}], ImageTemplate{IndDCE},roiFactorLT1DW);   
    obj = MCC_cPatient([WriteFolder PatientName '\T1-DW Registered\']);
    obj.Register_Date_Local({'TumorSegmentationDW' 'TumorSegmentationT1'},'T1W',...
        {'none','rigid','similarity','affine'},VOI);
    
    fid=fopen([WriteFolder PatientName '\T1-DW Registered\RegistrationParameters.txt'],'w');
    fprintf(fid, ['ROI Factor: ' num2str(roiFactorLT1DW)]);
    fclose(fid);
end








function SlicePos=CheckAndLoadFiles(File)
if exist(File,'file')~=0
    load(File);
else
    SlicePos=[0 0];
    disp('No SlicePositions/Ks file');
end




