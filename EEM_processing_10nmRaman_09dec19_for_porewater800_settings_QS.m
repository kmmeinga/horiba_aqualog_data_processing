%%% IMPORTANT PLEASE READ BEFORE USE!!!!!

%%% Updated to DIVIDE by correction factor
%%% Code is used to normalize Exitation Emission matrices to Raman units and correct for inner filter effects
%%% This script uses the flucut command (in PLS toolbox) to correct for
%%% inner filter effects. This toolbox must be installed on your machine
%%% Written by Karl Meingast(kmmeinga@mtu.edu) %%%%

%%% IMPORTANT PLEASE READ BEFORE USE!!!!!
%%% Naming conventions checklist
%%% Files MUST NOT start with "a","b" or "s". 
%%% Files can be named whatever the aqualog software will allow otherwise
%%% Sarna and blank are named sarna(date).csv and blank(data).csv
%%% very important that abs files get named abs_ in front of the name of the sample
%%% starna and blank files are named abs_s(date) and abs_b(date) #check
%%% this

%%% Last modified 4/17/15 Updates: Fixed file naming convention problem.
%%% Files may be named anything now (other than starting with a,b,s... 
%%% This was done by asking if name of the (kth) element of S (data
%%% structure) is a character. If it is not, the data in S(k) is empty
%%% becuae it corresponed with an absorption file which was included into
%%% the S(k) EEM element.
%%% Not maybe the tidiest way, but works well as a workaround.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc% Clear command window
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
        if (list(i,1).name(1,1:3)) == 'abs';
        abslist(i).name = list(i,1).name; % Find absorption files

        elseif (list(i,1).name(1,1:5)) == 'blank'; % find blank
%% Import the data from the blank
            [~, ~, raw] = xlsread(list(i,1).name); % read raw data
%% Replace non-numeric cells with 0.0
            R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
            raw(R) = {0.0}; % Replace non-numeric cells

%% Create output variable
            blank = cell2mat(raw); % change to matrix

%%Read in absorption file for blank
            filename = dir('abs_blank*');
            [~, ~, alldata] = xlsread(filename.name); % read abs data file

            abs_spec = (rot90(cell2mat(alldata(2:end,10)))); % Takes abs column from entire output

   %%% flucut the blank

    %% Set up options for flucut
   
        opts = [];
        opts.LowZero = 'on';
         opts.innerfilter.use = 'on';
        opts.innerfilter.spectra = abs_spec; % This is the only use of absorbance file
        opts.innerfilter.wavelengths = flipud(rot90(blank(1,2:end))); % Get wavelengths in correct orientation for flucut

    
        realsize = size(blank); %find size of data
        subblankdata = blank(2:realsize(1),2:realsize(2));% cut wavelenghts
        


   

        A = zeros(1,125,188); % Create dataframe to feed to flucut
        A(1,:,:) = subblankdata; %Fill with blank's data
        A = dataset(A);
        
        A.axisscale{2} = blank(2:end,1);
        A.axisscale{3} = flipud(rot90(blank(1,2:end)));
        
        blankflucutwithife_data = flucut(A, 15,15, opts);% Cut blank

        %% Clear temporary variables
       % clearvars raw R; % clear stuff to go fasterrrrr


            %% Create Sarna water %%%%%%%%%
        elseif (list(i,1).name(1,1:6)) == 'starna'; % find starna named ramandate.csv
            [~, ~, raw] = xlsread(list(i,1).name); % READ DATA
            %% Replace non-numeric cells with 0.0
            R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
            raw(R) = {0.0}; % Replace non-numeric cells
            
            %% Create output variable
            sarna = cell2mat(raw); % cange to matrix
            
            %% flucut Sarna
            
            %%Read in abs filename
            filename = dir('abs_starna*'); %Find abs file for starna
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
            
            
            
            
            
            A = zeros(1,125,188);
            A(1,:,:) = subsarnadata;
            A = dataset(A);
            
            A.axisscale{2} = sarna(2:end,1);% Assign axis scales to sarna dataset
            A.axisscale{3} = flipud(rot90(sarna(1,2:end)));
            
            
            sarna_flucutwithife_data = A % No cutting of raman line from Sarna. We need this to normalize to.
            
            
            %% Create Quinine sulfate 10/14/18 added KM %%%%%%%%%
        elseif (list(i,1).name(1,1:2)) == 'QS'; % find QS named ramandate.csv
            [~, ~, raw] = xlsread(list(i,1).name); % READ DATA
            %% Replace non-numeric cells with 0.0
            R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
            raw(R) = {0.0}; % Replace non-numeric cells
            
            %% Create output variable
            QS = cell2mat(raw); % cange to matrix
            
            %% flucut QS
            
            %%Read in abs filename
            filename = dir('abs_QS*'); %Find abs file for QS
            [~, ~, alldata] = xlsread(filename.name); % read abs data file
            
            abs_spec = (rot90(cell2mat(alldata(2:end,10)))); % Takes abs column from entire output
        
           
            
            %% Set up otions for flucut QS
            
            opts = [];
            opts.LowZero = 'on';
            opts.innerfilter.use = 'on';
            opts.innerfilter.spectra = abs_spec;
            opts.innerfilter.wavelengths = flipud(rot90(QS(1,2:end)));
            
            
            
            realsize = size(QS); %find size of data
            subQSdata = QS(2:realsize(1),2:realsize(2));% cut wavelenghts
            
            
            
            
            
            A = zeros(1,125,188);
            A(1,:,:) = subQSdata;
            A = dataset(A);
            
            A.axisscale{2} = QS(2:end,1);% Assign axis scales to QS dataset
            A.axisscale{3} = flipud(rot90(QS(1,2:end)));
            
            
          QS_flucutwithife_data = flucut((A), 15,15, opts); % No cutting of raman line from Sarna. We need this to normalize to.
            
            %% Clear temporary variables
           % clearvars raw R; %fasterrrrr
            
        else % this should find all csv files that are left ie not blank or starna
            
            splitname = strtok(list(i,1).name,'.'); % Extract name of data file
            ask = ['Enter dilution factor for: ',splitname , '    eg. THIS SHOULD BE BETWEEN 0 AND 1   ']; % This creates prompt asking for dilutin factor
            S(i).dilution_factor = input(ask);
         
     
            
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
            clearvars raw R RAW splitname alldata;
        end
    end
    
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Extract Raman Vlaues from sarna %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%% Integral value for QS %%%
   y = QS_flucutwithife_data.data(:,26:125,152) %Ex=350, Em=290-700 by 3nm Cartisano et al 2018

    % calculating the area under the curve
    QS = nansum(y)/10;
    clear y
    
    %%% iNT value starna for RU %%%
       y = sarna_flucutwithife_data.data(:,50:81,152);  %Following Lawaetz and Stedmon 2008; Ex=350, Em=371-428 by 3nm steps
    % calculating the area under the curve
    raman1 = nansum(y);
    clear y
    %%%%%%%%%%%%%%%%%%%%%% Find max value of blank looking near from Sarna index    
    
    %%%%%%%%%%%%%% Normalize samples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:length(S) % now loop through all structures that were built
        if ischar(S(k).name) == 0 % This is a workaround for empty data in S dataframe resulting when absorbtion file is found but not read. These abs files are merged into S where EEMs values exist
        ;
        else
        %% Build abs dataset for flucut
        abs_spec = (rot90(cell2mat(S(k).absdata(2:end,10)))); % Takes abs column from entire output
        
     
     
        %% Set up otions for flucut
        
        opts = [];
        opts.LowZero = 'on';
        opts.innerfilter.use = 'on';
        opts.innerfilter.spectra = abs_spec;
        opts.innerfilter.wavelengths = flipud(rot90(S(k).data(1,2:end)));
        opts.RamanCorrect = 'on';
        opts.RamanWidth = 10; %This adjusts raman width (nm) to cut from samples.
        
        
        realsize = size(S(k).data); %find size of data
        subdata = S(k).data(2:realsize(1),2:realsize(2));% cut wavelenghts
        A = zeros(1,125,188);
        A(1,:,:) = subdata;
        A = dataset(A);
        
        A.axisscale{2} = S(k).data(2:end,1);
        A.axisscale{3} = flipud(rot90(S(k).data(1,2:end)));
        
        
        S(k).flucutwithife_data = flucut(A, 15,15, opts);
        
        %% DO WE NEED TO SUBTRACT A BLANK????? CURRRENTLY NOT DUE TO LAEWAETZ AND STEDMON 2008
        
        normdata = (S(k).flucutwithife_data.data)./QS; % normalizing to raman peak at 350 nm from sarna water % DO WE NEED TO SUBTRACT A BLANK -TJV & KMM 8/14/14
        
        S(k).normalized_flucut_dilution = normdata/S(k).dilution_factor; % 05/08/2015 KMM - DIVIDE BY DILUTION FACTOR INSTEAD OF MULTIPLY.
       S(k).normalized_flucut_dilution(S(k).normalized_flucut_dilution < 0) = NaN;
        % Updated 7/30/14 to incorporate flucut
        
        normdataru = (S(k).flucutwithife_data.data)./raman1; % normalizing to raman peak at 350 nm from sarna water % DO WE NEED TO SUBTRACT A BLANK -TJV & KMM 8/14/14
        
        S(k).normalized_flucut_dilution_ru = normdataru/S(k).dilution_factor; % 05/08/2015 KMM - DIVIDE BY DILUTION FACTOR INSTEAD OF MULTIPLY.
       S(k).normalized_flucut_dilution_ru(S(k).normalized_flucut_dilution_ru < 0) = NaN;
        
        %figure(k+2) % open figure
        %plot(xdata,ydata); % plot figure same as blank. should be close to flat indicating max values were taken from close to raman line of sarna
    end
    end
    
    
    mkdir('matlab_norm_outputs')% make subfolder to put new files in
    cd ('matlab_norm_outputs') % go to that folderr
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  write data to files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk = 1:length(S)% loop through data frmaes
        if ischar(S(kk).name) == 0
            print('nodata')
        else
        output = zeros(125,188);
        output(:,:) =  S(kk).normalized_flucut_dilution(1,:,:);
        output = flipud(transpose(output));
        dlmwrite([S(kk).name,'_processed_EEM_QSE','.','csv'],output); % write a normalized data file for each file. output is comma seperated values with .MAT ending
           %%%% Updated 5/9/2015 to convert abs data to absorption
        %%%% coefficients. Uses (2).303*Abs)/l where l is the pathlength of
        %%%% 1cm or 0.010 m
       figure(100+kk)
       surf(rot90(fliplr(output)))
       title([S(kk).name,'_processed_EEM_QSE'])
       
      absorption_coefficients = (((cell2mat(S(kk).absdata(2:end,10))*2.303)/0.01)/S(kk).dilution_factor); % 
      abs_coeff_export = table(cell2mat(S(kk).absdata(2:end,1)),absorption_coefficients,'VariableNames',{'Wavelength_nm' 'Absoprtion_coeff_m_1'});
       %% Write out abs file too 
      writetable(abs_coeff_export,[S(kk).name,'_processed_abs_coefficients','.','csv'],'Delimiter',',');
        clear abs_coeff_export
        end
    end
        %%% Write out starna units %%%
          for kkk = 1:length(S)% loop through data frmaes
        if ischar(S(kkk).name) == 0
            print('nodata')
        else
        output = zeros(125,188);
        output(:,:) =  S(kkk).normalized_flucut_dilution_ru(1,:,:);
        output = flipud(transpose(output));
        dlmwrite([S(kkk).name,'_processed_EEM_Raman_Units','.','csv'],output); % write a normalized data file for each file. output is comma seperated values with .MAT ending
           %%%% Updated 5/9/2015 to convert abs data to absorption
        end
    end
    
    cd .. %%% BACK TO WHERE WE STARTED %%%%%%%%%%%%%%%%%%
    %clear all %% Remove all figures
    
