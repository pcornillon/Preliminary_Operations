function get_nadir_info(start_date_time, end_date_time)
% get_nadir_info - read nadir information for all granules in the specified period - PCC
%
% This function will read MODIS Aqua SST 11um granules in sequence for the
%  given year and save the nadir track, the tracks of the two edges of each
%  swath, the detector number, the orbit number and calculate the time for
%  each. It strings them all together and saves them on a monthly basis. It
%  will reset the vectors for each at the beginning of the next month,
%  after saving.
%
% The function also records all granules that are missing, writing their
%  estimated name to a cell array and adding nans in place of the missing
%  scanline information.
%
% INPUT
%   start_date_time - year, month, day, hour minute vector for start time.
%   end_date_time - year, month, day, hour minute vector for end time.
%
% OUTPUT
%
%
% EXAMPLE
%   get_nadir_info([2009 12 31 22 0], [2011 1 1 2 0]) - to process all
%    granules from 10 PM on December 31 2009 through 2 AM January 1 2011.
%    Including the last 2 hours of 2009 and the first 2 hours of 2011 will
%    assure that the first and last orbit touching 2010 are captured.

nc_read = 0;

if computer == 'MACI64'
    sw3 = 0;
    base_dir_in = '/Volumes/Aqua-1/MODIS_R2019/';
    base_dir_out = '/Volumes/Aqua-1/MODIS_R2019/Test/';
else
    sw3 = 1;
    base_dir_in =  's3://podaac-ops-cumulus-protected/MODIS_A-JPL-L2P-v2019.0/';
    base_dir_out =  's3://podaac-ops-cumulus-protected/MODIS_A-JPL-L2P-v2019.0/';
%     s3://podaac-ops-cumulus-protected/MODIS_A-JPL-L2P-v2019.0/    20100619075507-JPL-L2P_GHRSST-SSTskin-MODIS_A-N-v02.0-fv01.0.nc
end

Diary_File = [base_dir_out, 'Orbits/Logs/' strrep(num2str(now), '.', '_') '.txt'];
diary(Diary_File)

% Initialize various arrays and scalars

plotem = 0;
iFig = 1;

latlim = -78;

if (length(start_date_time) ~= 6) | (length(end_date_time) ~= 6)
    disp(['Input start and end time vectors must be 6 elements long. start_date_time: ' ...
        num2str(start_date_time) ' and end_date_time: ' num2str(end_date_time)])
    keyboard
end

matlab_time_start = datenum(start_date_time);
matlab_time_end = datenum(end_date_time);

if matlab_time_end < matlab_time_start
    disp(['End time: ' datestr(matlab_time_end) ' comes before start time: ' datestr(matlab_time_start)])
    keyboard
end

iFilename = 0;
iMissing = 0;

clear filenames missing_granules

imatlab_time = matlab_time_start;
month_save = -100;
day_save = -100;

% Loop over granules

tic

while imatlab_time <= matlab_time_end
    
    imatlab_time = imatlab_time + 5 / (24 * 60);
    
    % Get the date and time for this granule and generate strings for use
    % in the names.
    
    [iyear, imonth, iday, ihour, iminute, isecond] = datevec(imatlab_time);
    
    iyears = convertStringsToChars(num2str(iyear));
    imonths = return_a_string(imonth);
    idays = return_a_string(iday);
    ihours = return_a_string(ihour);
    iminutes = return_a_string(iminute);
    
    if iday ~= day_save
        day_save = iday;
        disp(['Working on ' imonths '/' idays '/' iyears ' at time ' num2str(toc)])
    end
    
    % If this is a new month save the nadir track info and reinialize the vectors.
    
    if (imonth ~= month_save)
        
        if exist('nlon')
            % First find the start and end of each complete orbit. Note
            % that 100 scan lines will be added to each end to assure
            % overlap needed for gradient and fronts calculations.
            
            diff_nlat = diff(nlat);
            nn = find(abs(nlat(1:end-1)-latlim)<0.1 & detnum(1:end-1)==5 & diff_nlat<0);
            
            diff_nn = diff(nn);
            mm = find(diff_nn > 50);
            nn_start = [nn(mm)' nn(end)] - 100;
            nn_end = [nn_start(2:end)-1 length(nlat)-100] + 100;
            
            nn_start(1) = max(nn_start(1), 1);
            nn_end(end) = min([nn_end(end), length(nlon)]);
            
            orbit_info.scan_line_start = nn_start;
            orbit_info.scan_line_end = nn_end;
            
            orbit_info.nlons = nlon(nn_start);
            orbit_info.nlats = nlat(nn_start);
            
            orbit_info.nlone = nlon(nn_end);
            orbit_info.nlate = nlat(nn_end);
            
            orbit_info.matlab_times = matlab_time(nn_start);
            orbit_info.matlab_timee = matlab_time(nn_end);
            
            orbit_info.det_nums = detnum(nn_start);
            orbit_info.det_nume = detnum(nn_end);
            
            orbit_info.filenames = filenames(filename_index(nn_start));
            orbit_info.filenamee = filenames(filename_index(nn_end));
            
            orbit_info.filenames_index = filenames(filename_index(nn_start));
            orbit_info.filenamee_index = filenames(filename_index(nn_end));
            
            % Now save.
            
            % Subtract a year from the value used in the filename if the
            % previous month, month_save=12, since the year would have been
            % incremented when going from month 12 to month 1.
            
            year_for_name = iyears;
            if month_save == 12
                year_for_name = convertStringsToChars(num2str(iyear - 1));
            end
            
            save( [base_dir_out, 'Orbits/nadir_info_' year_for_name '_' return_a_string(month_save)], ...
                'filenames', 'orbit_number', 'filename_index', 'matlab_time', '*lon', '*lat', 'detnum', ...
                'nsol_z', 'scan_line_in_file', 'orbit_info', 'latlim');
            
            matlab_time = matlab_time(nn_start(end):end);
            nlat = nlat(nn_start(end):end);
            nlon = nlon(nn_start(end):end);
            slat = slat(nn_start(end):end);
            slon = slon(nn_start(end):end);
            elat = elat(nn_start(end):end);
            elon = elon(nn_start(end):end);
            nsol_z = nsol_z(nn_start(end):end);
            detnum = detnum(nn_start(end):end);
            scan_line_in_file = scan_line_in_file(nn_start(end):end);
            
            jgranule = 0;
            for igranule=filename_index(nn_start(end)):filename_index(end)
               jgranule = jgranule + 1;
               temp(jgranule) = filenames(igranule);
            end
            filenames = temp;
            filename_index = filename_index(nn_start(end):end) - filename_index(nn_start(end)) + 1;

            iFilename = jgranule;
                        
            month_save = imonth;
            
            if plotem; iFig = plot_orbits(iFig, nlon, nlat, orbit_info); end
        else
            matlab_time = [];
            nlat = [];
            nlon = [];
            slat = [];
            slon = [];
            elat = [];
            elon = [];
            nsol_z = [];
            detnum = [];
            
            filename_index = [];
            scan_line_in_file = [];
            
            iFilename = 0;
            
            clear filenames missing_granules
            
            month_save = imonth;
        end
    end
    
    iFilename = iFilename + 1;
    
    % See if this file exists. If not, nans for missing data.
    
    file_found = 1;
    
    if sw3
        % Need to put sw3 stuff here.
    else
        new_fi_list = dir([base_dir_in 'day/' iyears '/AQUA_MODIS.' iyears imonths idays 'T' ihours iminutes '*']);
        
        if isempty(new_fi_list)
            new_fi_list = dir([base_dir_in 'night/' iyears '/AQUA_MODIS.' iyears imonths idays 'T' ihours iminutes '*']);
            
            if isempty(new_fi_list)
                disp(['Could not file file: AQUA_MODIS.' iyears imonths idays 'T' ihours iminutes '*'])
                file_found = 0;
            end
        end
        
        fi = [new_fi_list(1).folder '/' new_fi_list(1).name];
    end
    
    if file_found
% % %         fi = [new_fi_list(1).folder '/' new_fi_list(1).name];
        
        % Get info about this file.
        
        if nc_read
            info = ncinfo(fi);
            
            if strcmp(info.Dimensions(1).Name, 'number_of_lines') ~= 1
                disp(['Wrong dimension: ' info.Dimensions(1).Name])
                keyboard
            end
            nscans = info.Dimensions(1).Length;
            
            filenames{iFilename} = fi;
            
            filename_index = [filename_index; int32(ones(nscans,1)*iFilename)];
            
            if strcmp(info.Attributes(5).Name, 'orbit_number') ~= 1
                disp(['Wrong attribute: ' info.Dimensions(2).Name])
                keyboard
            end
            orbit_number(iFilename) = int32(info.Attributes(5).Value);
            
            sl_year = ncread(fi, '/scan_line_attributes/year');
            sl_day = ncread(fi, '/scan_line_attributes/day');
            sl_msec = floor(ncread(fi, '/scan_line_attributes/msec'));
            matlab_time = [matlab_time; datenum(sl_year, 1, sl_day) + sl_msec/86400000];
            
            nlat = [nlat; single(ncread(fi, '/scan_line_attributes/clat'))];
            nlon = [nlon; single(ncread(fi, '/scan_line_attributes/clon'))];
            
            slat = [slat; single(ncread(fi, '/scan_line_attributes/slat'))];
            slon = [slon; single(ncread(fi, '/scan_line_attributes/slon'))];
            
            elat = [elat; single(ncread(fi, '/scan_line_attributes/elat'))];
            elon = [elon; single(ncread(fi, '/scan_line_attributes/elon'))];
            
            nsol_z = [nsol_z; int16(ncread(fi, '/scan_line_attributes/csol_z'))];
            
            detnum = [detnum; int8(ncread(fi, '/scan_line_attributes/detnum'))];
            
            scan_line_in_file = [scan_line_in_file int16([1:nscans])];
        else
            info = h5info(fi);
            
            if strcmp(info(1).Datasets(3).Name, 'number_of_lines') ~= 1
                disp(['Wrong dimension: ' info.Dimensions(1).Name])
                keyboard
            end
            nscans = info(1).Datasets(3).Dataspace.Size;
            
            filenames{iFilename} = fi;
            
            filename_index = [filename_index; int32(ones(nscans,1)*iFilename)];
            
            orbit_number(iFilename) = int32(h5readatt( fi, '/', 'orbit_number'));
            
            sl_year = double(h5read(fi, '/scan_line_attributes/year'));
            sl_day = double(h5read(fi, '/scan_line_attributes/day'));
            sl_msec = double(floor(h5read(fi, '/scan_line_attributes/msec')));
            matlab_time = [matlab_time; datenum(sl_year, 1, sl_day) + sl_msec/86400000];
            
            nlat = [nlat; single(h5read(fi, '/scan_line_attributes/clat'))];
            nlon = [nlon; single(h5read(fi, '/scan_line_attributes/clon'))];
            
            slat = [slat; single(h5read(fi, '/scan_line_attributes/slat'))];
            slon = [slon; single(h5read(fi, '/scan_line_attributes/slon'))];
            
            elat = [elat; single(h5read(fi, '/scan_line_attributes/elat'))];
            elon = [elon; single(h5read(fi, '/scan_line_attributes/elon'))];
            
            nsol_z = [nsol_z; int16(h5read(fi, '/scan_line_attributes/csol_z'))];
            
            detnum = [detnum; int8(h5read(fi, '/scan_line_attributes/detnum'))];
            
            scan_line_in_file = [scan_line_in_file int16([1:nscans])];
        end
    else
        iMissing = iMissing + 1;
        missing_granules{iMissing} = [base_dir_in 'night/' iyears '/AQUA_MODIS.' iyears imonths idays 'T' ihours iminutes];
        
        filenames{iFilename} = missing_granules{iMissing};
        
        orbit_number(iFilename) = int32(nan);
        
        filename_index = [filename_index; int32(nan(nscans,1))];
        
        matlab_time = [matlab_time; double(nan(nscans,1))];
        
        nlat = [nlat; single(nan(nscans,1))];
        nlon = [nlon; single(nan(nscans,1))];
        
        slat = [slat; single(nan(nscans,1))];
        slon = [slon; single(nan(nscans,1))];
        
        elat = [elat; single(nan(nscans,1))];
        elon = [elon; single(nan(nscans,1))];
        
        nsol_z = [nsol_z; int16(nan(nscans,1))];
        
        detnum = [detnum; int8(nan(nscans,1))];
        
        scan_line_in_file = [scan_line_in_file int16(1:nscans)];
    end
end

% Save the last month's worth of data. Start by getting the beginning and
% end info for each orbit.

diff_nlat = diff(nlat);
nn = find(abs(nlat(1:end-1)-latlim)<0.1 & detnum(1:end-1)==5 & diff_nlat>0);

diff_nn = diff(nn);
mm = find(diff_nn > 50);
if isempty(mm)
    disp(['The data for this month does not include the beginning of an orbit.'])
    orbit_info = [];
else
    
    nn_start = [nn(mm)' nn(end)];
    nn_end = [nn_start(2:end)-1 length(nlat)-100] + 100;
    
    nn_start(1) = max(nn_start(1), 1);
    nn_end(end) = min(nn_end(end), length(nlon));
    
    orbit_info.scan_line_start = nn_start;
    orbit_info.scan_line_end = nn_end;
    
    orbit_info.nlons = nlon(nn_start);
    orbit_info.nlats = nlat(nn_start);
    
    orbit_info.nlone = nlon(nn_end);
    orbit_info.nlate = nlat(nn_end);
    
    orbit_info.matlab_times = matlab_time(nn_start);
    orbit_info.matlab_timee = matlab_time(nn_end);
    
    orbit_info.det_nums = detnum(nn_start);
    orbit_info.det_nume = detnum(nn_end);
    
    orbit_info.filenames = filenames(filename_index(nn_start));
    orbit_info.filenamee = filenames(filename_index(nn_end));
    
    orbit_info.filenames_index = filenames(filename_index(nn_start));
    orbit_info.filenamee_index = filenames(filename_index(nn_end));
end

% Subtract a year from the value used in the filename if the
% previous month, month_save=12, since the year would have been
% incremented when going from month 12 to month 1.

year_for_name = iyears;
if month_save == 12
    year_for_name = convertStringsToChars(num2str(iyear - 1));
end

save( [base_dir_out 'Orbits/nadir_info_' year_for_name '_' return_a_string(month_save)], ...
    'filenames', 'orbit_number', 'filename_index', 'matlab_time', '*lon', '*lat', 'detnum', ...
    'nsol_z', 'scan_line_in_file', 'orbit_info', 'latlim');

toc

if plotem; iFig = plot_orbits(iFig, nlon, nlat, orbit_info); end

%% Functions.

function iFig = plot_orbits(iFig, nlon, nlat, orbit_info)
figure(iFig)
iFig = iFig + 1;
clf

% xx = nlon;
% xx(xx<0) = xx(xx<0) +360;
% plot(xx,nlat,'.k')
plot(nlon,nlat,'.k')
hold on

load coastlines.mat
% yy = coastlon;
% yy(yy>0) = yy(yy>0) - 360;
% diffyy = diff(yy);
% nn = find(abs(diffyy)>50);
% yy(nn+1) = nan;
% plot(yy,coastlat, 'r', linewidth=2)
plot(coastlon,coastlat,'c',linewidth=2)

plot(orbit_info.nlons, orbit_info.nlats, 'ok', markersize=10, markerfacecolor='g')
plot(orbit_info.nlone, orbit_info.nlate, 'ok', markersize=10, markerfacecolor='r')
plot(nlon(orbit_info.scan_line_start(1):orbit_info.scan_line_end(1)), nlat(orbit_info.scan_line_start(1):orbit_info.scan_line_end(1)), 'r.', linewidth=3)
for i=1:length(orbit_info.scan_line_start)
    plot(nlon(orbit_info.scan_line_start(i):orbit_info.scan_line_end(i)), nlat(orbit_info.scan_line_start(i):orbit_info.scan_line_end(i)), '.', linewidth=3)
end
