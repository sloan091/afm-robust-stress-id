function data_d = aggData(time,data,dt,method)

if strncmp(method,'mean',4)
    data_d = retime(timetable(time,data),'Daily',...
        @(x) nanmean(x));
elseif strncmp(method,'sum',3)
    data_d = retime(timetable(time,data),'Daily',...
        @(x) sum(x)*(dt/24));
else
    error('Failed to specify valid aggregation method')
end

end