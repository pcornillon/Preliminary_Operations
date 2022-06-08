% submit_Run_Fix_MODIS_Masks_nccopy - this script will copy variables from one netCDF file to another - PCC
%
% Specifically, it will copy year, day, msec, slon, slat, clon, clat, elon,
%  elat, csol_z, sst, sstref, qual_sst, flags_sst and tilt from the
%  original file obtained from Goddard to a new file in the 'original'
%  directory on one of the large disks. It will also generate a diary file.
%
% It will ask for the years to copy and for the starting month and day if 
%  the year is entered with a minus sign.

Diary_File = ['/Volumes/Aqua-1/Fronts/MODIS_Aqua_L2/Logs/Rewrite_Original_' strrep(num2str(now), '.', '_') '.txt'];
diary(Diary_File)

% Get the name of the disk to which the data are to be written.

disk_out = input('Enter the name of the output disk: ', 's');
base_dir_out = ['/Volumes/' disk_out '/MODIS_Aqua_L2/'];
% base_dir_out = '/Volumes/Aqua-1/Fronts/MODIS_Aqua_L2/';

% Get the year, month and day at which to start processing.

year_list = input('Enter year(s) to process cell array (e.g., {''2003'' ''2004'' ''2005'' ''2006'' ''2007''}): ');

% If the first element in year list is entered with a - in front, then ask
%  for month and day to start. Otherwise, get the month and day of the most 
%  recently processed file for this year and determine month and day to 
%  start if requested

if strfind(year_list{1}, '-')
    year_list{1} = year_list{1}(2:end);
    
    month_start = input(['Enter the month to start with ', year_list{1} ' (1 for January or cr): ']);
    if isempty(month_start); month_start = 1; end
    
    day_start = input(['Enter the day to start with ', num2str(month_start) '/' year_list{1} ' (1 for first day or cr): ']);
    if isempty(day_start); day_start = 1; end
else    
    file_list = dir( ['/Volumes/Aqua-1/Fronts/MODIS_Aqua_L2/Original/' year_list{1} '/AQUA*']);
    if isempty(file_list)
        month_start = 1;
        day_start = 1;
    else
        date_list = cat(1,file_list.datenum);
        nn = find(max(date_list) == date_list);
        
        month_start = str2num(file_list(nn).name(16:17));
        day_start = str2num(file_list(nn).name(18:19));
    end
end

disp(['Successfully started job submitted for ' year_list{1} ' starting at month ' num2str(month_start) ' and day ' num2str(day_start)])

matlab_start = datenum(str2num(year_list{1}), month_start, day_start);

base_dir = '/Volumes/Aqua-1/Fronts/MODIS_Aqua_L2/';

for iYear=1:length(year_list) % Loop over years to process ......................................
    
    fprintf('Processing %5.0f \n', iYear)
    
    % Get the list of files to consider for processing.
    
    YearS = year_list{iYear};
    file_list = dir(['/Volumes/Aqua-1/MODIS_R2019/night/' YearS '/AQUA*.nc']);
    
    tic
     for iFile=1:length(file_list) % Loop over granules *******************************************
        
        file_in = [file_list(iFile).folder '/' file_list(iFile).name];
        
        if mod(iFile,1000) == 0
            fprintf('Processed file #%i - %c - Elapsed time %6.0f seconds. \n', iFile, fine_in, toc)
            tic
        end
        
        % Has this file already been processed?
        
        nn_name_start = strfind(file_in, '/');
        nn_name_end = strfind(file_in, '.nc');
        good_base_filename_out = strrep(file_in(nn_name_start(end)+1:nn_name_end-1), '.', '_');
        good_filename_out = [good_base_filename_out '_original.nc4'];
        
        % Get year and month to put this granule in the proper directory.
        
        nn_year = strfind(file_in, 'AQUA_MODIS.') + 11;
        YearS = file_in(nn_year:nn_year+3);
        MonthS = file_in(nn_year+4:nn_year+5);
        DayS = file_in(nn_year+6:nn_year+7);
        
        % Check to see if this file is for a pass on or after the  start
        % month and day specified for the 1st year. If so process, if not
        % go to the next pass.
        
        if datenum(str2num(YearS), str2num(MonthS), str2num(DayS)) >= matlab_start
            
            file_out = [base_dir_out 'Original/' YearS '/' MonthS '/' good_filename_out];
            
            if exist(file_out) == 2
                fprintf('%i: %c has already been processed; skipping to the next file. \n', iFile, file_list(iFile).name)
            else
                
                status = system(['/usr/local/bin/nccopy -w -V year,day,msec,slon,slat,clon,clat,elon,elat,csol_z,sst,sstref,qual_sst,flags_sst,tilt ' file_in ' ' file_out]);
%                 status = system(['/usr/local/bin/ncks -C -v year,day,msec,slon,slat,clon,clat,elon,elat,csol_z,sst,sstref,qual_sst,flags_sst,tilt ' file_in ' ' file_out]);
                if status ~= 0
                    fprintf('Problem writing fields for granule number %i named %c \n', iFile, file_out)
                end 
            end
        end
    end
end

toc
