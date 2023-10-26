function Main(PatientName, Folder, WriteFolder, Options)
% Main file to process one patient
    
[DoDWI,CurateDWI,AnnotateDWI,...
    GlobalRegisterDWI,LocalRegisterDWI,roiFactorDW,GetParDWI,GlobalDateRegisterDWI,...
    LocalDateRegisterDWI,roiFactorLDW,UseSegmentationForLocalRegDWI,DoT1W,...
    CurateT1W,AnnotateT1W,GlobalRegisterT1W,SpetialRegT1,LocalRegisterT1W,...
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

close all force
warning off

% Read files
if exist([WriteFolder PatientName filesep 'FilesInfo.mat'],'file')==0 | (CurateDWI & DoDWI) |...
        (CurateT1W & DoT1W)
    I=1;
    FilesFolderNames=[]; Infos=[]; SeriesDescription=[];
    disp('Reading files')
    for Date=Dates       
        [FilesFolderNames{I},Infos{I},SeriesDescription{I}]=ReadPatientFiles([Folder PatientName filesep Date{1} filesep '**' filesep '*']); 
        I=I+1;
    end
    mkdir([WriteFolder PatientName])
    save([WriteFolder PatientName filesep 'FilesInfo.mat'],'FilesFolderNames','Infos','SeriesDescription'); 
else
    load([WriteFolder PatientName filesep 'FilesInfo.mat']);
end


if DoDWI
    % Curate DW-MRI    
    if CurateDWI
        try, rmdir([WriteFolder PatientName '\DWI\'],'s');  end         
        I=1;
        for Date=Dates
            CurateDataDWI(WriteFolder, Infos{I},SeriesDescription{I},PatientName,Date{1},FolderDWI);
            [Templates{I},TemplateDate{I},LocalRegTemplate{I}]=GetDW_RegistrationTemplates([WriteFolder PatientName '\DWI\' Date{1}]);
            I=I+1;
        end
        if exist([WriteFolder PatientName '\DWI\'])~=0
            save([WriteFolder PatientName '\DWI\Templates.mat'], 'LocalRegTemplate', 'TemplateDate', 'Templates');
        end
    end       
    % Global Register DW-MRI  
    if GlobalRegisterDWI  
        disp('Global Reg. ADC')
        I=1;
        if exist([WriteFolder PatientName '\DWI\\Templates.mat'])~=0
            load([WriteFolder PatientName '\DWI\\Templates.mat']);
            for Date=Dates  
                if ~isempty(TemplateDate{I})
                    obj = MCC_cPatient([WriteFolder PatientName '\DWI\' Date{1}]);
                    obj.Register_Date(Templates{I},TemplateDate{I}, {'none','rigid','similarity','affine'});
                    clear obj
                    OrganizeRegisteredDWI([WriteFolder PatientName '\DWI\' Date{1}],WriteFolder) 
                end
                I=I+1;
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
                            MakeAnnotationsTemplate([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\'],I,LocalRegTemplate{I});
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
        disp('Local Reg. ADC')
        I=1;
        if exist([WriteFolder PatientName '\DWI\Templates.mat'])~=0
            load([WriteFolder PatientName '\DWI\Templates.mat']);
            for Date=Dates   
                if exist([WriteFolder PatientName '\DWI\Registered\' Date{1} '\Global\'])~=0
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
        try, rmdir([WriteFolder PatientName filesep 'T1W' filesep],'s');  end
        I=1;
        for Date=Dates
            CurateT1WFun(WriteFolder,Infos{I},SeriesDescription{I},PatientName,Date{1},Folder_T1); 
            CurateDCEfun(Infos{I},SeriesDescription{I},PatientName,Date{1},Folder_DCE,WriteFolder);
            I=I+1;
        end
    end  
    % Annotate T1
    if AnnotateT1W 
        disp('Annotate T1 masks')
        Ind=cellfun(@(x) exist([WriteFolder PatientName filesep 'T1W' filesep x filesep])~=0, Dates);
        I=1;
        for Date=Dates(Ind)  
            if exist([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep])~=0 
                f=0;
                while f==0
                    try
                        MakeAnnotationsTemplate([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep],I);
                        f=1;
                    catch 
                        disp(' ');disp('Error, try again.')
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
        if exist([WriteFolder PatientName filesep 'T1W' filesep])~=0
            I=1;
            if ~exist([WriteFolder PatientName filesep 'T1W' filesep 'ImageTemplates.mat'])
                for Date=Dates
                    if exist([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep])~=0
                        ImageTemplate{I}=GetT1RegistrationTemplate([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep]);  
                    end
                    I=I+1;
                end
                ImageTemplate(cellfun(@isempty,ImageTemplate))=[];
                save([WriteFolder PatientName filesep 'T1W' filesep 'ImageTemplates.mat'], 'ImageTemplate');
            else
                load([WriteFolder PatientName filesep 'T1W' filesep 'ImageTemplates.mat']);
            end            
            try, rmdir([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep],'s');  end
            try, rmdir([WriteFolder PatientName filesep 'T1W' filesep 'DateRegistered' filesep],'s');  end
            try, rmdir([WriteFolder PatientName filesep 'T1W' filesep 'DateRegisteredLocal' filesep],'s');  end            
            obj = MCC_cPatient([WriteFolder PatientName filesep 'T1W' filesep]);
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
        fid=fopen([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep 'SpetialRegT1.txt'],'w');
        fprintf(fid, 'True');
        fclose(fid);
    end
    % Local registration 
    if LocalRegisterT1W        
        disp('Local registration T1w')
        load([WriteFolder PatientName filesep 'T1W' filesep 'ImageTemplates.mat']);
        I=1;
        for Date=Dates    
            if exist([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep])~=0
                VisualizeLocalRegistrationT1([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep 'Processed' filesep 'Registered'],...
                    ImageTemplate{I},roiFactorT1, IncludeArteryLR);
                 I=I+1;
            end           
        end
        fid=fopen([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep 'roiFactorT1.txt'],'w');
        fprintf(fid, ['Artery: ' num2str(IncludeArteryLR) char(10) num2str(roiFactorT1)]);
        fclose(fid);
    end
    % Obtain T1 maps
    if GetT1        
        disp('Computing T1 maps')
        for Date=Dates
            if exist([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep])~=0
                T1map([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep],[WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep],...
                    RefTiss, ScaleT1Ims, gammaT1, gammaM0, CorrectBiasBetweenSlices,AssumeT1s,AirTh,solver,Date{1},Dates,...
                    WriteFolder,PatientName,CorrectFA,WeightT1); 
            end
        end
        fid=fopen([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep 'T1mapppingParameters.txt'],'w');
        fprintf(fid, ['CorrectFA: ' num2str(CorrectFA) char(10) 'Scale: ' num2str(ScaleT1Ims) char(10)...
            'RefTiss: ' RefTiss char(10) 'Gammas: ' num2str(gammaT1) ', ' num2str(gammaM0) char(10)...
            'BiasCorrection: ' num2str(CorrectBiasBetweenSlices) char(10) 'AssumeT1: ' num2str(AssumeT1s) char(10)...
            'AirTh: ' num2str(AirTh) char(10) 'Solver: ' solver]);
        fclose(fid);
    end
    % Get Concentration maps
    if GetConc
        disp('Computing C(t)')
        load([WriteFolder PatientName filesep 'T1W' filesep 'ImageTemplates.mat']);
        I=1;
        for Date=Dates      
            if exist([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep 'DCE_t=1' filesep])~=0   
                [r,ContrastBrand]=Get_rContrast([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep],ContrastBrand);
                [Cons_NW_NL,Conds_NW_NL,Times]=GetConcentrations([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep],...
                                    [WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep],r,ImageTemplate{I}, ScaleMo, AssumeT1sC, AssumeT1Artery);                 
                save([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep 'Concentrations.mat'], 'Cons_NW_NL', 'Conds_NW_NL',...
                        'Times','r');
                I=I+1;
            end
        end
        save([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep 'AssumeT1sC.mat'],'AssumeT1sC')
        fid=fopen([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep 'ConcentrationParameteres.txt'],'w');
        fprintf(fid,['ContrastBrand: ' ContrastBrand char(10) 'r: ' num2str(r) char(10) 'AssumedT1: ' num2str(AssumeT1sC)...
            char(10) 'AssumedT1 Artery: ' num2str(AssumeT1Artery)]);
        fclose(fid);
    end
    % Get DCE parametric Maps
    if GetPar       
        disp('AUC maps')
        for Date=Dates 
            if exist([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep 'DCE_t=1' filesep])~=0
                AUC_NW=Get90sAUC([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep],...
                    [WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep]);   
                save([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep '90sAUCs.mat'],...
                    'AUC_NW');
                close all force
            end
        end
        disp('Tofts maps')
        for Date=Dates  
            if exist([WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep 'DCE_t=1' filesep])~=0
                Folder2=[WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep];
                [Vp_NW_NL,Ve_NW_NL,Kt_NW_NL]=GetKinetcisParameters(WriteFolder,Folder2,...
                    [WriteFolder PatientName filesep 'T1W' filesep Date{1} filesep]);   
                save([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Date{1} filesep 'Local' filesep 'Perfussion_Parameters_Maps.mat'],...
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




