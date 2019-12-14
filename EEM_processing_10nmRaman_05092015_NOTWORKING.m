%%% IMPORTANT PLEASE READ BEFORE USE!!!!!

%%% Updated to DIVIDE by correction factor
%%% Code is used to normalize Exitation Emission matrices to Raman unitsand correct for inner filter effects
%%% This script uses the flucut command (in PLS toolbox) to correct for
%%% inner filter effects. This toolbox must be installed on your machine
%%% Written by Karl Meingast(kmmeinga@mtu.edu) %%%%

%%% IMPORTANT PLEASE READ BEFORE USE!!!!!
%%% Naming conventions checklist
%%% Files MUST NOT start with "a","b" or "s". 
%%% Files can be named whatever the aqualog software will allow otherwise
%%% Sarna and blank are named sarna(date).csv and blank(data).csv
%%% very important that abs files get named abs_ in front of the name of the sample (all numbers)
%%% sarna and blank abs files are named abs_s(date) and abs_b(date)

%%% Last modified 4/17/15 Updates: Fixed file naming convention problem.
%%% Files may be named anything now (other than starting with a,b,s... 
%%% This was done by asking if name of the (kth) element of S (data
%%% structure) is a character. If it is not, the data in S(k) is empty
%%% becuae it corresponed with an absorption file which was included into
%%% the S(k) EEM element.
%%% Not maybe the tidiest way, but works well as a workaround.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc   % Clear command window
    clear % Clear workspace


%%% Take me to the spot %%
    foldername = uigetdir; %prompts you for folder to look in
    cd(foldername) % goes to that folder
%%% Make a list of all files
    list = dir('*.csv'); % lists all csv file in that folder. No extra csv files should be in there

%%%%%%%%%%%%%%%%%%READ IN THE DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find abs files added 8/4/2014


    for i=1:length(list); % start main loop through all csv files
        if (list(i,1).name(1,1)) == 'a';
        abslist(i).name = list(i,1).name; % Find absorption files

        elseif (list(i,1).name(1,1)) == 'b'; % find blank
%% Import the data from the blank
            [~, ~, raw] = xlsread(list(i,1).name); % read raw data
%% Replace non-numeric cells with 0.0
            R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
            raw(R) = {0.0}; % Replace non-numeric cells

%% Create output variable
            blank = cell2mat(raw); % change to matrix

%%Read in absorption file for blank
            filename = dir('abs_b*');
            [~, ~, alldata] = xlsread(filename.name); % read abs data file

            abs_spec = (rot90(cell2mat(alldata(2:end,10)))); % Takes abs column from entire output

   %%% flucut the blank

    %% Set up otions for flucut
   
        opts = [];
        opts.LowZero = 'on';
         opts.innerfilter.use = 'on';
        opts.innerfilter.spectra = abs_spec; % This is the only use of absorbance file
        opts.innerfilter.wavelengths = flipud(rot90(blank(1,2:end))); % Get wavelengths in correct orientation for flucut

    
        realsize = size(blank); %find size of data
        subblankdata = blank(2:realsize(1),2:realsize(2));% cut wavelenghts
        


   

        A = zeros(1,125,121); % Create dataframe to feed to flucut
        A(1,:,:) = subblankdata; %Fill with blank's data
        A = dataset(A);
        
        A.axisscale{2} = blank(2:end,1);
        A.axisscale{3} = flipud(rot90(blank(1,2:end)));
        
        blankflucutwithife_data = flucut(A, 15,15, opts);% Cut blank

        %% Clear temporary variables
       % clearvars raw R; % clear stuff to go fasterrrrr

%% Create Sarna water %%%%%%%%%
        elseif (list(i,1).name(1,1)) == 's'; % find starna named ramandate.csv
            [~, ~, raw] = xlsread(list(i,1).name); % READ DATA
            %% Replace non-numeric cells with 0.0
            R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
            raw(R) = {0.0}; % Replace non-numeric cells
            
            %% Create output variable
            sarna = cell2mat(raw); % cange to matrix
            
            %% flucut Sarna
            
            %%Read in abs filename
            filename = dir('abs_s*'); %Find abs file for starna
            [~, ~, alldata] = xlsread(filename.name); % read abs data file
            
            abs_spec = (rot90(cell2mat(alldata(2:end,10)))); % Takes abs column from entire output
        
            
            %% Set up otions for flucut
            
            opts = [];
            opts.LowZero = 'on';
            opts.innerfilter.use = 'on';
            opts.innerfilter.spectra = abs_spec;
            opts.innerfilter.wavelengths = flipud(rot90(sarna(1,2:end)));
            
            
            
            realsize = size(sarna); %find size of data
            subsarnadata = sarna(2:realsize(1),2:realsize(2));% cut wavelenghts
            
            
            
            
            
            A = zeros(1,125,121);
            A(1,:,:) = subsarnadata;
            A = dataset(A);
            
            A.axisscale{2} = sarna(2:end,1);% Assign axis scales to sarna dataset
            A.axisscale{3} = flipud(rot90(sarna(1,2:end)));
            
            
            sarna_flucutwithife_data = flucut((A), 15,15, opts); % No cutting of raman line from Sarna. We need this to normalize to.
            
            %% Clear temporary variables
           % clearvars raw R; %fasterrrrr
            
        else % this should find all csv files that are left ie not blank or starna
            
            splitname = strtok(list(i,1).name,'.'); % Extract name of data file
            ask = ['Enter dilution factor for: ',splitname , '    eg. THIS SHOULD BE BETWEEN 0 AND 1   ']; % This creates prompt asking for dilutin factor
            dilution_factor = input(ask);
         
     
            
            %% Import the data 
            [~, ~, raw] = xlsread(list(i,1).name); % read data file
            
            
            %% Replace non-numeric cells with 0.0
            R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
            raw(R) = {0.0}; % Replace non-numeric cells
            
            %% Create output variable
            RAW = cell2mat(raw); % change to matrix
            %eval([char(splitname),'=RAW;']);% rename raw but I don't think I need this
            S(i).data = RAW;% build S... a matalb structure that will have all data associated with each file. this puts the data in it These are really useful ways to handle data!!!! If you start using matlab, learn how to use these from the get go. they make life easier
            S(i).name = splitname; % this associates the file name to S
            
            absfilename = ['abs_',splitname,'.csv'];
            
            [~, ~, alldata] = xlsread(absfilename); % read abs data file
            
            S(i).absdata = alldata;
            %% Clear temporary variables
            clearvars raw R RAW splitname;
        end
    end
    
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Extract Raman Vlaues from sarna %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%% Integral value for starna %%%
    x=round((57/3)); %Not currently using this number but derived from LawaetZ and Stedmon 2008 as dx for raman integral
    y = sarna_flucutwithife_data.data(:,51:69,85);  %Following Lawaetz and Stedmon 2008; Ex=350, Em=371-428 by 3nm steps
    ydx = y*3; %Multiplies each y value on raman peak by 3 nm, since these values represent 3 nm incriments
    % calculating the area under the curve
    raman1 = nansum(ydx);
    clear y
    
    %%%%%%%%%%%%%%%%%%%%%% Find max value of blank looking near from Sarna index
    %%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
    normblanktosarna =  blankflucutwithife_data.data./sarna_flucutwithife_data.data; %normalixe blank to sarna
    
    
    %%% Integral value for blank %%%
    x=round((57/3)); 
    y = blankflucutwithife_data.data(:,51:69,85);  
    ydx = y*3;
    % calculating the area under the curve
    raman2 = nansum(ydx);
    raman = raman2/raman1; % correction factor decided 8/14/14 we may not need this
    
    
    %%%%%%%%%%%%%% Normalize samples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:length(S) % now loop through all structures that were built
        if ischar(S(k).name) == 0 % This is a workaround for empty data in S dataframe resulting when absorbtion file is found but not read. These abs files are merged into S where EEMs values exist
        print('nodata')
        else
        %% Build abs dataset for flucut
        abs_spec = (rot90(cell2mat(S(k).absdata(2:end,10)))); % Takes abs column from entire output
        
        %%%% Updated 5/9/2015 to convert abs data to absorption
        %%%% coefficients. Uses (2.303*Abs)/l where l is the pathlength of
        %%%% 1cm or 0.010 m
        
      abs_coeff_data = horzcat(S(k).absdata(:,1),S(k).absdata(:,10));
      absorption_coefficients(1,:) = {'absorption_coeffs'}; %STUCK HERE 5/9/2015
      absorption_coefficients = (((cell2mat(S(k).abs_coefficients(2:end,2))*2.303)/0.01)/dilution_factor);
      abs_coeff_data = [abs_coeff_data absorption_coefficients]
        %% Set up otions for flucut
        
        opts = [];
        opts.LowZero = 'on';
        opts.innerfilter.use = 'on';
        opts.innerfilter.spectra = abs_spec;
        opts.innerfilter.wavelengths = flipud(rot90(S(k).data(1,2:end)));
        opts.RamanCorrect = 'on';
        opts.RamanWidth = 10;
        
        
        realsize = size(S(k).data); %find size of data
        subdata = S(k).data(2:realsize(1),2:realsize(2));% cut wavelenghts
        A = zeros(1,125,121);
        A(1,:,:) = subdata;
        A = dataset(A);
        
        A.axisscale{2} = S(k).data(2:end,1);
        A.axisscale{3} = flipud(rot90(S(k).data(1,2:end)));
        
        
        S(k).flucutwithife_data = flucut(A, 15,15, opts);
        
        %% DO WE NEED TO SUBTRACT A BLANK????? CURRRENTLY NOT DUE TO LAEWAETZ AND STEDMON 2008
        
        normdata = (S(k).flucutwithife_data.data)./raman1; % normalizing to raman peak at 350 nm from sarna water % DO WE NEED TO SUBTRACT A BLANK -TJV & KMM 8/14/14
        
        S(k).normalized_flucut_dilution = normdata/dilution_factor; % 05/08/2015 KMM - DIVIDE BY DILUTION FACTOR INSTEAD OF MULTIPLY.
       S(k).normalized_flucut_dilution(S(k).normalized_flucut_dilution < 0) = NaN;
        % Updated 7/30/14 to incorporate flucut
        
        
        
        %figure(k+2) % open figure
        %plot(xdata,ydata); % plot figure same as blank. should be close to flat indicating max values were taken from close to raman line of sarna
    end
    end
    
    
    mkdir('matlab_norm_outputs')% make subfolder to put new files in
    cd ('matlab_norm_outputs') % go to that folderr
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  write data to files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk = 1:length(S)% loop throgh data frmaes
        if ischar(S(kk).name) == 0
            print('nodata')
        else
        output = zeros(125,121);
        output(:,:) =  S(kk).normalized_flucut_dilution(1,:,:);
        output = flipud(transpose(output));
        dlmwrite([S(kk).name,'_processed','.','csv'],output); % write a normalized data file for each file. output is comma seperated values with .MAT ending
        %% Write out abs file too
        output_abs = S(kk).abs_coefficients;
        dlmwrite([S(kk).name,'_processed_abs_coefficients','.','csv'],output_abs);
        end
        
    end
    
    cd .. %%% BACK TO WHERE WE STARTED %%%%%%%%%%%%%%%%%%
    %clear all %% Remove all figures
    
