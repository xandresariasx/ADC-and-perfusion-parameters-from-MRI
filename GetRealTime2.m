function TimesOut=GetRealTime2(Info)

try
    if contains(Info.Filename,'SMIL')
        TriggerTime=Info.TriggerTime/1000;
    else
        TriggerTime=Info.TriggerTime;
    end
catch
    TriggerTime=0;
end

if isempty(TriggerTime)
    TriggerTime=0;
end

AcquisitionTime=str2num(Info.AcquisitionTime);

TimesOut=floor(AcquisitionTime/10000)*60^2 +  floor(mod(AcquisitionTime,10000)/100)*60  +  mod(AcquisitionTime, 100);

TimesOut=TimesOut+TriggerTime;

