function Time=GetTriggerTime(Info)

try
    Time=Info.TriggerTime;
    if isempty(Time)
        Time=0;
    end
catch
    Time=0;
end