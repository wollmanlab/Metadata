classdef Metadata < handle
    % METADATA class stores all the metadata, i.e. acqusition and
    % experimental data, for a specific set of images.
    % It is used by the Scope class to store all the information (including
    % use supplied one) and therby is critical to acqusition process.
    %
    % There are three "types" of uses: 1. store metadata, 2. retrive metaedata
    % and 3. retrive images/stacks by metadata
    %
    % Storing Metadata:
    % 1. add a new image to metadata using method: addNewImage
    %
    %         MD.addNewImage(filename,'prop1',value1,'prop2',value2)
    %
    % method returns the index of the image which could be used to easily
    % add additional metadata to that image.
    % 2. add to existing image(s) using image index by calling addToImages
    %
    %         MD.addToImages(indexes,'prop1',value,'prop2',value2)
    %
    % Retriving Metadata:
    % 1. First opertaion is to get the indexes of the specific images using
    %
    %         indx = MD.getIndex('prop1',value1,'prop2',value2)
    %
    %    getIndex is a fundamental operation of Metadata that
    %    accepts a series of conditions that metadata should have and it
    %    returns the indxes of all images that obey these criteria.
    %
    % Note: getIndex has a "hack" where if the type has the ending _indx
    % that instead of Value the user can supply the index to the UnqValue
    % cell array - just a shorthand syntax.
    %
    % 2. given indexes, getting metadata is just a simple call to
    %
    %        Values = MD.getSpecificMetadataByIndex(indx,Types)
    %
    % Retriving images:
    % just use imread(MD,prop1,value1,...) / stkread(MD,prop1,value1) with conditions similar to MD.getIndex(...)
    %
    
    properties
        basepth='';
        acqname
        
        ImgFiles = cell(0,1); % a long list of filenames
        Values= {}; % Channels,Zslices,Timestamps_avg,Positins, ...
        Types = {}; % "Header"  - shows what is stored in each cell of Values
        OldValuesSize = 0;
        
        NonImageBasedData = {};
        
        Project
        Username
        Dataset
        
        Description
        
        dieOnReadError = false;
        
        defaultTypes = {  'Zindex' 'Z' 'Channel'    'Exposure'    'Fluorophore'    'Marker'    'PixelSize'    'PlateType'    'Position'    'Skip'    'TimestampFrame'    'TimestampImage'    'XY'    'acq'    'frame'    'group' 'AllInputs' 'XYbeforeTransform' 'driftTform'}'
    end
    
    properties (Transient = true)
        pth
        verbose = true;
    end
    
    properties (Dependent = true)
        sz
        NewTypes
    end
    
    methods
        function MD = Metadata(pth,acqname,matflag)%constructor
            if nargin==0
                return
            end
            
            %  pth=regexprep(pth,'data4','bigstore');
            %  pth=regexprep(pth,'data3','bigstore');
            
            % if a Metadata.mat file exist just load it
            if exist(fullfile(pth,'Metadata.mat'),'file')
                if exist(fullfile(pth,'Metadata.txt'),'file') && nargin<3
                    s.MD=MD.readDSV(pth);
                else
                    s=load(fullfile(pth,'Metadata.mat'));
                    MD = s.MD;
                end
                MD=s.MD;
                MD.pth = pth;
                
                %% we're done.
                return
            else
                % check to see if Metadata.mat files exist in multiple
                % subdirectories, if so just read them all and merge,
                files = rdir([pth filesep '**' filesep 'Metadata.mat']);
                if ~isempty(files)
                    pths = regexprep({files.name},[filesep 'Metadata.mat'],'');
                    for i=1:numel(pths)
                        MDs(i) = Metadata(pths{i}); %#ok<AGROW>
                    end
                    MD = merge(MDs);
                    return
                end
                if nargin==1
                    error('Couldn''t find Metadata.mat in path: %s\n please check',pth)
                end
            end
            % Creates an EMPTY Metadata object with only pth and acqname
            MD.pth = pth;
            MD.acqname = acqname;
            
            %% add the AllInputs field if not there already
            %             if ~ismember('AllInputs',MD.Types)
            %                 MD.mergeTypes(MD.NewTypes,'AllInputs');
            %             end
            
        end
        
        
        function new_md = readDSV(~, pth)
            delimiter = '\t';
            s = load(fullfile(pth, 'Metadata.mat'));
            M = readtable(fullfile(pth, 'Metadata.txt'),'delimiter',delimiter);
            
            if any(strcmp('XY',M.Properties.VariableNames))
                M.XY = cellfun(@str2num,M.XY,'UniformOutput', false);%Parse XY values
            end
            if any(strcmp('XYbeforeTransform',M.Properties.VariableNames))
                M.XYbeforeTransform = cellfun(@str2num,M.XYbeforeTransform,'UniformOutput', false);%Parse XY values
            end
            if any(strcmp('driftTform',M.Properties.VariableNames))
                M.driftTform = cellfun(@str2num,M.driftTform,'UniformOutput', false);%Parse XY values
            end
            
            types = M.Properties.VariableNames;
            values = table2cell(M);
            s.MD.Values = values(:, 1:end-1);
            s.MD.Types = types(1:end-1);
            s.MD.ImgFiles = values(:, end)';
            new_md = s.MD;
            %            NO NO! MD made by scope MUST match MD made by reload!
            %            new_md.convert_type_datatype('XYbeforeTransform', @num2str);
            %            new_md.convert_type_datatype('XY', @num2str);
        end
        
        function MD = convert_type_datatype(MD, type, type_func)
            idx = find(cellfun(@(x) strcmp(x, type), MD.Types));
            if ~any(idx)
                disp('Type not found so nothing will happen.')
                return
            end
            %             assert(any(idx), 'Type not found in Metadata')
            new_type = cellfun(type_func, MD.Values(:, idx), 'UniformOutput', false);
            %             new_type = cell2mat(new_type);
            %             new_type = mat2cell(new_type, size(new_type, 1));
            MD.Values(:, idx) = new_type;
        end
        
        
        function sz = get.sz(MD)
            sz = size(MD.Values); %number of types?
        end
        
        function NT = get.NewTypes(MD)
            NT = setdiff(MD.Types,MD.defaultTypes); %C = setdiff(A,B) returns the data in A that is not in B.
        end
        
        function MD = removeImagesByIndex(MD,indx)
            MD.ImgFiles(indx)=[];
            MD.Values(indx,:)=[];
        end
        
        
        function prj = project(MD,func,varargin)%what is func here?
            % same logic as stkread but performs functio func on all images
            T = varargin(1:2:end);
            V = varargin(2:2:end);
            indx = MD.getIndex(T,V);
            files = cellfun(@(f) fullfile(MD.pth,f),MD.ImgFiles(indx),'uniformoutput',0);
            files = regexprep(files,'\\','/');
            prj = stackProject(files,'func',func);
        end
        
        function out = stkfun(MD,func,varargin)
            
            % update verbose to the status of the publishing flag if it
            % exist
            vrb=getappdata(0,'publishing');
            if ~isempty(vrb)
                D.verbose = ~vrg;
            end
            
            % same logic as stkread but performs functio func on all images
            T = varargin(1:2:end);
            V = varargin(2:2:end);
            indx = MD.getIndex(T,V);
            n=0;
            out = cell(numel(indx),1);
            for i=1:numel(indx)
                filename = fullfile(MD.basepth,MD.pth,MD.ImgFiles{indx(i)});
                filename = regexprep(filename,'\\',filesep);
                MD.verbose && fprintf(repmat('\b',1,n)); %#ok<VUNUS>
                msg = sprintf('processing image %s, number %g out of %g\n',filename,i,numel(indx));
                n=numel(msg);
                msg = regexprep(msg,'\\','\\\\');
                MD.verbose && fprintf(msg); %#ok<VUNUS>
                try
                    tf = Tiff(filename,'r');
                    img = tf.read();
                    img = mat2gray(img,[0 2^16]);
                    out{i} = func(img);
                catch e
                    warning('Couldn''t read image %s, error message was: %s',filename,e.message);  %#ok<WNTAG>
                end
            end
        end
        
        function [stk,indx,filename] = stkread(MD,varargin)
            % reads a stack of images based on criteria.
            % criteria must be supplied in type,value pair
            % there are six speical cases for properties that are not
            % really types:
            %
            % stkread(MD,...,'sortby',prop,...)
            %       This will sort the stack using the property prop
            %       Default is to sort by TimestampFrame - to do otherwise,
            %       pass a different field or empty (...,'sortby','',...)
            %
            % stkread(MD,...'max',mx,...)
            %       Only read mx images
            %
            % stkread(MD,...,'specific',nm,...)
            %       Reads a specific plane, nm could be a number or
            %       'first','last','median' that will be converted to
            %       numbers.
            %
            % stkread(MD,...,'timefunc',@(t) t<datenum('June-26-2013 19:46'),...)
            %       Reads images up to specific timepoint, good for
            %       runaway
            %       experiments...
            %
            % stkread(MD,...,'resize',sz,...)
            %       Resizes images as they are read.
            %
            % stkread(MD,...,'groupby',grp,...)
            %       Loads multiple stacks and groups them by grp
            %
            % stkread(MD,...,'register',)
            %
            %
            % Function  assumes that all images have the same size (!)
            
            % update verbose to the status of the publishing flag if it
            % exist
            vrb=getappdata(0,'publishing');
            if ~isempty(vrb)
                MD.verbose = ~vrb;
            end
            
            if numel(varargin)==2 && iscell(varargin{1}) && iscell(varargin{2})
                T=varargin{1};
                V=varargin{2};
            else
                T = varargin(1:2:end);
                V = varargin(2:2:end);
            end
            
            %% find out if there is a groupby
            if ismember('groupby',T)
                groupby = V{ismember(T,'groupby')};
                V(ismember(T,'groupby'))=[];
                T(ismember(T,'groupby'))=[];
                Grp = unique(MD,groupby);
                if isnumeric(Grp)
                    Grp = num2cell(Grp);
                end
                stk = cell(size(Grp));
                for i=1:numel(Grp)
                    stk{i} = stkread(MD,[{groupby} T],[Grp{i} V]);
                end
                indx = unique(MD,groupby);
                indx(cellfun(@isempty,stk))=[];
                stk(cellfun(@isempty,stk))=[];
                return
            end
            
            %% figure out if I need to resize
            resize = V(ismember(T,'resize'));
            V(ismember(T,'resize'))=[];
            T(ismember(T,'resize'))=[];
            if isempty(resize)
                resize = 1;
            else
                resize = resize{1};
            end
            
            montage = V(ismember(T,'montage'));
            if isempty(montage)
                montage = false;
            else
                montage = montage{1};
            end
            
            V(ismember(T,'montage'))=[];
            T(ismember(T,'montage'))=[];
            
            func = V(ismember(T,'func'));
            if isempty(func)
                func = @(m) m;
            else
                func = func{1};
            end
            V(ismember(T,'func'))=[];
            T(ismember(T,'func'))=[];
            
            registerflag = V(ismember(T,'register'));
            if isempty(registerflag)
                registerflag=0;
            else
                registerflag=registerflag{1};
            end
            V(ismember(T,'register'))=[];
            T(ismember(T,'register'))=[];
            
            
            awtflag = V(ismember(T,'blindflatfield'));
            if isempty(awtflag)
                awtflag=0;
            else
                awtflag=awtflag{1};
            end
            V(ismember(T,'blindflatfield'))=[];
            T(ismember(T,'blindflatfield'))=[];
            
            flatfieldcorrection= V(ismember(T,'flatfieldcorrection'));
            if isempty(flatfieldcorrection)
                flatfieldcorrection=false; % Don't correct by default in stkread! (!) AOY
            else
                flatfieldcorrection=flatfieldcorrection{1};
            end
            V(ismember(T,'flatfieldcorrection'))=[];
            T(ismember(T,'flatfieldcorrection'))=[];
            
            %% get indexes
            indx = MD.getIndex(T,V);
            if isempty(indx)
                stk=[];
                return
            end
            
            %% get image size and init the stack
            try
                filename = MD.getImageFilename({'index'},{indx(1)});
                warning('off')
                info = imfinfo(filename);
                warning('on')
            catch  %#ok<CTCH>
                try
                    filename = MD.getImageFilename({'index'},{indx(2)});
                    warning('off')
                    info = imfinfo(filename);
                    warning('on')
                catch %#ok<CTCH>
                    if MD.dieOnReadError
                        error('Files not found to read stack')
                    else
                        warning('Files not found to read stack')
                        info.Height = 2048;
                        info.Width = 2064;
                    end
                end
            end
            
            stk = cell(numel(indx),1);
            
            % deal with blind ff
            if awtflag
                fltfieldnames = MD.getSpecificMetadataByIndex('Channel',indx);
                unqFltFieldNames = unique(fltfieldnames);
                FlatFields = zeros([[info.Height info.Width] numel(unqFltFieldNames)],'single');
                handledflatfields = {};
            end
            
            %% read the images needed for flat field correction
            if flatfieldcorrection
                fltfieldnames = MD.getSpecificMetadataByIndex('FlatField',indx);
                unqFltFieldNames = unique(fltfieldnames);
                %changed FlatFields init because siz has already been
                %resized, leading to problems downstream
                FlatFields = zeros([[info.Height info.Width] numel(unqFltFieldNames)],'uint16');
                for i=1:numel(unqFltFieldNames)
                    try
                        warning('off')
                        info1 = imfinfo(fullfile(MD.pth,['flt_' unqFltFieldNames{i} '.tif']));
                        warning('on')

                        FlatFields(:,:,i)=imread(fullfile(MD.pth,['flt_' unqFltFieldNames{i} '.tif']), 'Info', info1);
                    catch
                        [pth2,~]=fileparts(MD.pth);
                        FlatFields(:,:,i)=imread(fullfile(pth2,['flt_' unqFltFieldNames{i} '.tif']));
                    end
                end
            end
            %% read the stack
            n=0;
            filename=cell(numel(indx),1);
            for i=1:numel(indx)
                filename{i} = MD.getImageFilename({'index'},{indx(i)});
                MD.verbose && fprintf(repmat('\b',1,n));%#ok<VUNUS>
                msg = sprintf('reading image %s, number %g out of %g\n',filename{i},i,numel(indx));
                n=numel(msg);
                msg = regexprep(msg,'\\','\\\\');
                MD.verbose && fprintf(msg);%#ok<VUNUS>
                try
                    if exist([filename{i} '.bz2'],'file')
                        system(sprintf('bunzip2 %s',[filename{i} '.bz2']));
                    end
                    
                    %img = imread(filename);
                    %tf = Tiff(filename{i},'r');
                    %img = tf.read();
                    %img=single(img)/2^16;
                    warning('off')
                    info = imfinfo(filename{i});
                    warning('on')
                    num_images = numel(info);
                    
                    
                    blnk = zeros([info(1).Height info(1).Width]);
                    blnk = imresize(blnk,resize);
                    siz = [size(blnk),num_images];
                    img = zeros(siz,'single');
                    
                    for k = 1:num_images
                        img1 = imread(filename{i}, k, 'Info', info);
                        img1=single(img1)/2^16;
                        
                        
                        
                        if flatfieldcorrection
                            try
                                flt = FlatFields(:,:,ismember(unqFltFieldNames,fltfieldnames{i}));
                                img1 = doFlatFieldCorrection(MD,img1,flt);
                            catch
                                MD.dieOnReadError = 1;
                                error('Could not find flatfield files, to continue without it add flatfieldcorrection, false to stkread call');
                            end
                        end
                        
                        
                        if awtflag %blind flat field correction
                            if ~ismember(fltfieldnames{i},handledflatfields)
                                ind1=ismember(unqFltFieldNames,fltfieldnames{i});
                                handledflatfields{ind1} = fltfieldnames{i};
                                
                                disp(['Calculating flat field for channel ' fltfieldnames{i}])
                                
                                img2awt = mean(img1,3);
                                awtImage = awt2Dlite(img2awt,9);
                                FlatFields(:,:,ind1)=squeeze(awtImage(:,:,:,end));
                            end
                            flt = FlatFields(:,:,ismember(unqFltFieldNames,fltfieldnames{i}));
                            
                            %img1 = max(flt(:)).*(img1./repmat(flt,1,1,size(img1,3)));

                            img1 = (img1-repmat(flt,1,1,size(img1,3)));
                            img1 = img1.*(img1>0);
                            %img1 = max(flt(:))+img1;
                        end
                        
                        
                        
                        
                        if registerflag
                            fprintf('Registering... ')
                            TformT = MD.getSpecificMetadataByIndex('driftTform',indx(i));
                            TformT = TformT{1};
                            if size(TformT,2)==9
                                Tform = affine2d(reshape(TformT,3,3)');
                                img1 = imwarp(img1,Tform,'OutputView',imref2d(size(img1)));
                            else
                                fprintf('No drift correction transform found.')
                            end
                            fprintf('...done\n')
                        end
                        
                        
                        if resize~=1
                            img1 = imresize(img1,resize);
                        end
                        
                        
                        img1 = func(img1);
                        
                        img(:,:,k) = img1; %add to singlefilestack
                        
                    end
                    
                    stk{i}=img;
                    %stk = cat(3,stk,img);
                    %stk = permute(cell2mat(stk),[2 3 1]);
                catch e
                    if MD.dieOnReadError
                        error('Couldn''t read image %s, error message was: %s',filename{i},e.message);  %#ok<WNTAG>
                    else
                        warning('Couldn''t read image %s, error message was: %s',filename{i},e.message);  %#ok<WNTAG>
                    end
                end
            end
            stk = cat(3,stk{:});
            pos = MD.getSpecificMetadataByIndex('Position',indx);
            grp = MD.getSpecificMetadataByIndex('group',indx);
            if montage && ~isequal(pos,grp)
                %% make into a 2D montage
                splt = regexp(pos,'_','split');
                r = cellfun(@(s) str2double(s{4}),splt);
                c = cellfun(@(s) str2double(s{3}),splt);
                
                mntg = single(zeros(size(stk,1)*max(r),size(stk,2)*max(c)));
                for i=1:1:numel(r)
                    ixr=(r(i)-1)*size(stk,1)+(1:size(stk,1));
                    ixc=(c(i)-1)*size(stk,2)+(1:size(stk,2));
                    mntg(ixr,ixc)=flipud(stk(:,:,i));
                end
                stk=mntg;
            end
            
            
        end
        
        function  img = doFlatFieldCorrection(~,img,flt,varargin)
            % inputs
            arg.cameraoffset = 8974/2^16;
            arg = parseVarargin(varargin,arg);
            
            % the flt that is passed in uint16, convert...
            flt = double(flt)-3814/2^16;
            
            %flt = double(flt)-arg.cameraoffset/2^16;
            flt = flt./nanmean(flt(:));
            for i=1:size(img,3)
                img1 = img(:,:,i);
                img1 = double(img1-arg.cameraoffset)./flt+arg.cameraoffset;
                img1(flt<0.05) = prctile(img1(unidrnd(numel(img1),10000,1)),1); % to save time, look at random 10K pixels and not all of them...
                img(:,:,i) = img1;
            end
            % deal with artifacts
            img(img<0)=0;
            img(img>1)=1;
        end
        
        function Time = getTime(MD,varargin)
            if numel(varargin)==2 && iscell(varargin{1}) && iscell(varargin{2})
                T=varargin{1};
                V=varargin{2};
            else
                T = varargin(1:2:end);
                V = varargin(2:2:end);
            end
            units = V(ismember(T,'units'));
            V(ismember(T,'units'))=[];
            T(ismember(T,'units'))=[];
            if isempty(units)
                units = 'seconds';
            else
                units = units{1};
            end
            
            indx = MD.getIndex(T,V);
            Tstart = MD.getSpecificMetadata('TimestampFrame','specific','first');
            Time = MD.getSpecificMetadataByIndex('TimestampFrame',indx);
            Time=cat(1,Time{:});
            Time=Time-Tstart{1};
            switch units
                case {'d','days'}
                case {'h','hours'}
                    Time=Time*24;
                case {'m','min','minutes'}
                    Time=Time*24*60;
                case {'s','sec','seconds'}
                    Time=Time*24*60*60;
            end
            
        end
        
        % getIndex gets two cell arrays one for critieria and one for
        % values and it returns the indexes of images that are true for
        % these criteria
        function indx = getIndex(M,T,V)
            
            % 1. First opertaion is to get the indexes of the specific images using
            %
            %         indx = MD.getIndex('prop1',value1,'prop2',value2)
            %
            %    getIndex is a fundamental operation of Metadata that
            %    accepts a series of conditions that metadata should have and it
            %    returns the indxes of all images that obey these criteria.
            %
            % Note: getIndex has a "hack" where if the type has the ending _indx
            % that instead of Value the user can supply the index to the UnqValue
            % cell array - just a shorthand syntax.
            
            %% deal with three special cases, sortby, max and specific
            % to sort the stack and max to limit the number of images in
            % the stack and get a specific plane from the stack.
            
            sortby = V(ismember(T,'sortby'));
            if isempty(sortby)
                sortby='TimestampFrame';
            end
            V(ismember(T,'sortby'))=[];
            T(ismember(T,'sortby'))=[];
            
            mx = V(ismember(T,'max'));
            V(ismember(T,'max'))=[];
            T(ismember(T,'max'))=[];
            
            specific = V(ismember(T,'specific'));
            V(ismember(T,'specific'))=[];
            T(ismember(T,'specific'))=[];
            
            timefunc = V(ismember(T,'timefunc'));
            V(ismember(T,'timefunc'))=[];
            T(ismember(T,'timefunc'))=[];
            
            % get index via criteria
            tf = true(size(M.ImgFiles));
            tf=tf(:);
            % make Types and Vlaues into cells if needed
            if ~iscell(T)
                T = {T};
            end
            if ~iscell(V)
                V = {V};
            end
            for i=1:numel(T)
                indx  = strfind(T{i},'_indx');
                if ~isempty(indx)
                    T{i} = T{i}(1:(indx-1));
                    unq = unique(MD,T{i});
                    V{i} = unq{V{i}};
                end
                ixcol = ismember(M.Types,T(i));
                if ~any(ixcol)
                    error('Types requested are wrong - check for typos');
                end
                if isnumeric(V{i}(1))
                    vtmp = M.Values(:,ixcol);
                    vtmp(cellfun(@isempty,vtmp))={NaN};
                    assert(all(cellfun(@(v) numel(v)==1,vtmp)),'Error in metadata - there are vector values!')
                    va = cat(1,vtmp{:});
                else
                    va = M.Values(:,ixcol);
                end
                tf = tf & ismember(va,V{i});
            end
            indx = find(tf);
            if isempty(indx) % no point sorting an empty array...
                return
            end
            
            %% perform extra operations (sort, specific, max)
            if ~isempty(timefunc)
                T=M.getSpecificMetadataByIndex('TimestampFrame',indx);
                T=cat(1,T{:});
                timefunc=timefunc{1};
                indx = indx(timefunc(T));
            end
            
            if ~isempty(sortby)
                Vsort = M.Values(indx,ismember(M.Types,sortby));
                %adding this try-catch as timefunc filtering can create
                %an empty Vsort
                try
                    Vsort{1};
                catch
                    error('Error, well has no images due to timefunc filter')
                end
                if isnumeric(Vsort{1})
                    vtmp = Vsort;
                    vtmp(cellfun(@isempty,vtmp))={NaN};
                    Vsort = cat(1,vtmp{:});
                end
                [~,ordr]=sort(Vsort);
                indx = indx(ordr);
            end
            
            if ~isempty(mx)
                indx = indx(1:mx{1});
            end
            
            if ~isempty(specific)
                specific = specific{1};
                if ischar(specific)
                    switch specific
                        case 'last'
                            specific = numel(indx);
                        case 'first'
                            specific = 1;
                        case 'median'
                            specific = ceil(numel(indx/2));
                    end
                end
                indx = indx(specific);
            end
            
            
        end%end of getIndex function
        
        % addNewImages adds a specific image to the Metadata structure, it
        % also adds its properties / values pairs to the MD structure
        % creating new column if necessary.
        function ix = addNewImage(M,filename,varargin)
            T = varargin(1:2:end);
            V = varargin(2:2:end);
            if any(ismember(M.ImgFiles,filename))
                %msgbox('Warning - added image already exist in Metadata object!');
            end
            % add files to list
            M.ImgFiles{end+1} = filename;
            ix = numel(M.ImgFiles);
            for i=1:numel(T)
                ixcol = ismember(M.Types,T{i});
                if ~any(ixcol)
                    M.Types{end+1}=T{i};
                    ixcol = numel(M.Types);
                end
                M.Values(ix,ixcol)=V(i);
            end
        end
        
        % add metadata to existing image by index.
        function addToImages(M,ix,varargin)
            T = varargin(1:2:end);
            V = varargin(2:2:end);
            for i=1:numel(T)
                ixcol = ismember(M.Types,T{i});
                if any(ixcol)
                    M.Values(ix,ixcol)=repmat(V(i),numel(ix),1);
                else
                    M.Types{end+1}=T{i};
                    M.Values(ix,end+1)=repmat(V(i),numel(ix),1);
                end
            end
        end
        
        function mergeTypes(MD,Types,newname)
            V = MD.getSpecificMetadata(Types);
            for i=1:numel(Types)
                if isnumeric(V{1,i}) || islogical(V{1,i})
                    V(:,i) = cellfun(@(s) sprintf('%0.2g',s),V(:,i),'uniformoutput',0);
                end
            end
            for i=1:MD.sz(1)
                newV = V{i,1};
                for j=2:size(V,2)
                    newV = [newV '_' V{i,j}];  %#ok<AGROW>
                end
                MD.addToImages(i,newname,newV);
            end
        end
        
        % returns the image filename that could include any subfolders
        % that are down of the Metadata.mat file. In case that the
        % images where saved on Windows OS it will replace the \ with
        % appropriate filesep.
        function [filename,indx] = getImageFilename(M,Types,Values)
            if strcmp(Types{1},'index')
                indx = Values{1};
            else
                indx = M.getIndex(Types,Values);
            end
            
            if numel(indx) ~=1
                error('criteria should be such that only one image is returned - please recheck criteria');
            end
            filename = M.ImgFiles{indx};
            filename = fullfile(M.basepth,M.pth,filename);
            filename = regexprep(filename,'\\',filesep);
            
            filename = regexprep(filename,'data3','bigstore');
            filename = regexprep(filename,'data4','bigstore');
        end
        

        function [filename,indx] = getImageFilenameRelative(M,Types,Values)
            if strcmp(Types{1},'index')
                indx = Values{1};
            else
                indx = M.getIndex(Types,Values);
            end
            
            if numel(indx) ~=1
                error('criteria should be such that only one image is returned - please recheck criteria');
            end
            filename = M.ImgFiles{indx};
            filename = regexprep(filename,'\\',filesep);
        end
        
        % read an image using criteria, use the more allaborate stkread
        % for fancy options.
        function [img,indx] = imread(M,varargin)
            % gets an image based on criteria types and value pairs
            T = varargin(1:2:end);
            V = varargin(2:2:end);
            
            [imgfilename,indx] = getImageFilename(M,T,V);
            if ~exist(imgfilename,'file')
                keyboard;
            end
            tf = Tiff(imgfilename,'r');
            img = tf.read();
        end
        
        % allow the user to get data on an image from its index
        function Vmd = getSpecificMetadataByIndex(M,T,indx)
            ixcols = ismember(M.Types,T);
            Vmd = M.Values(indx,ixcols);
            %Values is A X B cell array, A is the total number of images, B is number of types
        end
        
        % allows the user to get few data types based on crietria pairs
        function [Vmd,indx] = getSpecificMetadata(MD,Ttoget,varargin)
            % This is an important methods, really need to annotated well
            %how varargin should be like, need an example
            if numel(varargin)==2 && iscell(varargin{1}) && iscell(varargin{2})
                T=varargin{1};
                V=varargin{2};
            else
                T = varargin(1:2:end);
                V = varargin(2:2:end);
            end
            
            %% find out if there is a groupby
            if ismember('groupby',T)
                groupby = V{ismember(T,'groupby')};
                V(ismember(T,'groupby'))=[];
                T(ismember(T,'groupby'))=[];
                Grp = unique(MD,groupby);
                if isnumeric(Grp)
                    Grp = num2cell(Grp);
                end
                Vmd = cell(size(Grp));
                for i=1:numel(Grp)
                    Vmd{i} = getSpecificMetadata(MD,Ttoget,[{groupby} T],[Grp{i} V]);
                end
                return
            end
            
            indx = MD.getIndex(T,V);
            Vmd = getSpecificMetadataByIndex(MD,Ttoget,indx);
        end
        
        % get the number of values of a specific type. Useful for for loops
        % i.e. for i=1:numOf(MD,'Position')
        function N = numOf(M,Type)
            if ~iscell(Type)
                Type = {Type};
            end
            N=cellfun(@(t) numel(unique(M,t)),Type);
        end
        
        % runs the matlab function grpstats on a metadata where user
        % spcifies a type to be a grouping, a type to calculate and a cell
        % array of stats to calculate.
        function varargout = grpstats(M,typetocalculate,groupingtype,whichstats)
            ixcol = ismember(M.Types,typetocalculate);
            X = M.Values(:,ixcol);
            if ~isnumeric(X{1})
                error('can only do grpstats to Types that are numeric!')
            end
            X=cat(1,X{:});
            ixcol = ismember(M.Types,groupingtype);
            grp = M.Values(:,ixcol);
            varargout = cell(size(whichstats));
            for i=1:numel(whichstats)
                varargout{i} = grpstats(X,grp,whichstats{i});
            end
        end
        
        % return the unique value of a specific type
        function [unq,Grp] = unique(M,Type,varargin)
            % method will create a unique list of values for a type (or
            % cell array of types)
            % default behavior is to remove any Nan, [],'',{}
            
            % create a cell array of Type (or default to all)
            if nargin==1
                Type = M.Types;
            end
            if ~iscell(Type)
                Type = {Type};
            end
            
            Grp={};
            if ~isempty(varargin)
                %% transform varargin into T/V pair
                if numel(varargin)==2 && iscell(varargin{1}) && iscell(varargin{2})
                    T=varargin{1};
                    V=varargin{2};
                else
                    T = varargin(1:2:end);
                    V = varargin(2:2:end);
                end
                
                %% find out if there is a groupby
                if ismember('groupby',T)
                    groupby = V{ismember(T,'groupby')};
                    V(ismember(T,'groupby'))=[];
                    T(ismember(T,'groupby'))=[];
                    Grp = unique(M,groupby);
                    if isnumeric(Grp)
                        Grp = num2cell(Grp);
                    end
                    unq = cell(size(Grp));
                    for i=1:numel(Grp)
                        unq{i} = unique(M,Type,[{groupby} T],[Grp{i} V]);
                    end
                    if all(cellfun(@iscell,unq)) && all(cellfun(@(m) numel(m)==1,unq))
                        unq = cellfun(@(x) x{1},unq,'uniformoutput',0);
                    end
                    return
                end
                indx = getIndex(M,T,V);
            else
                indx = 1:size(M.Values,1);
            end
            v=cell(numel(indx),numel(Type));
            for i=1:numel(Type)
                v(:,i) = M.Values(indx,ismember(M.Types,Type{i}));
                ix = find(cellfun(@numel,v(:,i))==1);
                for j=1:numel(ix)
                    %                     if ~iscell(vv{1})
                    %                         vv={vv(1)};
                    %                     end
                    v(ix(j),i) = v(ix(j),i);
                end
                
                if isnumeric(v{1,i}) || islogical(v{1,i})
                    ix=cellfun(@isempty,v(:,i));
                    v(ix,i)={nan};
                end
                
            end
            % unique set
            if size(v,2)==1 % only single type
                
                % "uncell" cell with a cell of size 1
                ix = find(cellfun(@(c) iscell(c) & numel(c)==1,v));
                for i=1:numel(ix)
                    v{ix(i)}=v{ix(i)}{1};
                end
                
                % assert that there are no more cells
                assert(~any(cellfun(@iscell,v)),'There are cells in the Metadata that are of size >1, please check')
                
                % remove [],'',Nan
                ix = cellfun(@(c) isempty(c) || (isnumeric(c) && isnan(c)),v);
                v(ix)=[];
                
                if isempty(v)
                    unq={};
                    return
                end
                
                % assrt that all values are numeric, char or logical;
                assert(all(cellfun(@(c) isnumeric(c) || islogical(c) || ischar(c),v)),'Error in propery, must be number, logical or string');
                
                % act on the nmeric and string portions of v seperately
                ix = cellfun(@(c) isnumeric(c) || islogical(c),v);
                vnumeric = cat(1,v{ix});
                unqnumeric = unique(vnumeric);
                ix = cellfun(@(c) ischar(c),v);
                unqchar = unique(v(ix));
                if isempty(unqnumeric)
                    unq = unqchar;
                elseif isempty(unqchar)
                    unq = unqnumeric;
                else
                    unqnumeric = num2cell(unqnumeric);
                    unq = [unqnumeric; unqchar];
                end
            else
                unq = uniqueRowsCA(v);
            end
            
        end
        
        % simple save to file.
        function saveMetadata(MD,pth)
            % saves the Metadata object to path pth
            if iscell(MD.pth) && numel(MD.pth)>1
                error('A composite Metadata can''t be saved!!');
            end
            if nargin==1
                pth = MD.pth;
            end
            % save the header as Metadata.mat and the rest as txt file
            V=MD.Values;
            MD.appendMetadataDSV(pth, 'Metadata.txt');
            MD.Values={};
            
            save(fullfile(pth,'Metadata.mat'),'MD')
            MD.Values=V;
            %MD.exportMetadata(pth);
        end
        
        
        function saveMetadataMat(MD,pth)
            % saves the Metadata object to path pth
            if iscell(MD.pth) && numel(MD.pth)>1
                error('A composite Metadata can''t be saved!!');
            end
            if nargin==1
                pth = MD.pth;
            end
            
            save(fullfile(pth,'Metadata.mat'),'MD')
        end
        
        function appendMetadataDSV(MD, pth, fname)
            delimiter = '\t';
            V=MD.Values;
            vcount = size(V, 1);
            oldvcount = MD.OldValuesSize;
            if vcount-oldvcount > 0
                
                if exist(fullfile(pth, fname))
                    newvalues = V(end-(vcount-oldvcount-1):end, :);
                    newfnames = MD.ImgFiles';
                    newfnames = strrep(newfnames, '\', '/');
                    newfnames = newfnames(end-(vcount-oldvcount-1):end);
                    md_export_csv = [newvalues newfnames];
                    
                    cell2csvAppend(fullfile(pth, fname), md_export_csv, delimiter);
                    MD.OldValuesSize = size(MD.Values, 1);
                else
                    
                    cell2csv(fullfile(pth, fname), cat(2, MD.Types, {'filename\n'}),  delimiter);
                    MD.appendMetadataDSV(pth, fname);
                end
            end
        end
        
        function exportMetadata(MD, pth, varargin)
            delimiter = '\t';
            fnames = MD.ImgFiles';
            fnames = fullfile(MD.pth, fnames);
            fnames = strrep(fnames, '\', '/');
            md_export_csv = [MD.Values fnames];
            md_types = cat(2, MD.Types, 'filename');
            md_export_csv = cat(1, md_types, md_export_csv);
            cell2csv(fullfile(pth, 'Metadata.txt'), md_export_csv, delimiter);
        end
        
        % simple disp override
        function disp(MD)
            if numel(MD)>1
                warning('Array of Metadata with %g elements - showing only the first one!\n',numel(MD));
            end
            % displays the Types and upto 10 rows
            MD(1).Types
            if size(MD(1).Values,1)>=10
                MD(1).Values(1:10,:)
            else % display all rows
                MD(1).Values
            end
        end
        
        % method will merge an array of Metadatas into a single Metadata
        % a "composite" metadata can't be saved!
        function MD = merge(MDs)
            
            %% start by finding the "base" path for all different MDs
            allpths = regexprep({MDs.pth},'//','/');
            prts = regexp(allpths,filesep,'split');
            [~,ordr]=cellfun(@sort,prts,'uniformoutput',false);
            N=cellfun(@(n) 1:numel(n),ordr,'uniformoutput',false);
            revordr = cell(size(ordr));
            for i=1:numel(ordr)
                revordr{i}(ordr{i})=N{i};
            end
            bs=prts{1};
            for i=2:numel(prts)
                bs = intersect(prts{i},bs);
            end
            ix = find(ismember(prts{1},bs));
            bspth = prts{1}{ix(1)};
            for i=2:numel(ix)
                bspth = [bspth filesep prts{1}{ix(i)}];  %#ok<AGROW>
            end
            bspth = regexprep(bspth,'//','/');
            rest = cellfun(@(p) p((numel(bspth)+1):end),allpths,'uniformoutput',0);
            
            %% create new Metadata
            MD = Metadata;
            MD.pth = bspth;
            for i=1:numel(MDs)
                %%
                V=MDs(i).Values;
                T=MDs(i).Types;
                
                % add new Types to MD
                [~,ixnew]=setdiff(T,MD.Types);
                [~,ixexisting,ixorderInMD]=intersect(T,MD.Types);
                
                % construct the cell to add
                Vnew = cell(size(V,1),numel(MD.Types) + numel(ixnew));
                Vnew(:,ixorderInMD)=V(:,ixexisting);
                Vnew(:,numel(MD.Types) +(1:numel(ixnew)))=V(:,ixnew);
                
                
                MD.Types = [MD.Types T(ixnew)];
                % add empty cols for all new Types i.e. no value exist in
                % current MD for them.
                MD.Values = [MD.Values cell(size(MD.Values,1),numel(ixnew))];
                MD.Values = [MD.Values; Vnew];
                
                % fix the filenames
                MD.ImgFiles = [MD.ImgFiles cellfun(@(f) fullfile(rest{i},f),MDs(i).ImgFiles,'uniformoutput',0)];
            end
            
            
        end
        
        % tabulate
        function tbl = tabulate(MD,T1,varargin)
            if numel(varargin)==2 && iscell(varargin{1}) && iscell(varargin{2})
                T=varargin{1};
                V=varargin{2};
            else
                T = varargin(1:2:end);
                V = varargin(2:2:end);
            end
            indx = MD.getIndex(T,V);
            Out = MD.getSpecificMetadataByIndex(T1,indx);
            tbl=tabulate(Out);
        end
        
        function plotMetadataHeatmap(MD,Type,varargin)
            
            arg.plate = Plate;
            arg.removeempty = false;
            arg.colormap = [0 0 0; jet(256)];
            arg.fig = 999;
            arg.default = 'all'; % what to do if not Type is sopecified. default is All, alternative is 'dialog' to get type.
            arg = parseVarargin(varargin,arg);
            
            if strcmp('Type','?')
                arg.default='dialog';
                Type='';
            end
            if nargin==1 || isempty(Type)
                switch arg.default
                    case 'dialog'
                        Type = listdlg('PromptString','Please choose:',...
                            'SelectionMode','single',...
                            'ListString',MD.NewTypes);
                    case 'all'
                        Type = 'AllInputs';
                        %% add the AllInputs field if not there already
                        if ~ismember('AllInputs',MD.Types)
                            MD.mergeTypes(MD.NewTypes,'AllInputs');
                        end
                end
            end
            
            Pos = unique(MD,'group');
            Val = unique(MD,Type,'groupby','group');
            
            
            arg.plate.x0y0=[0 0];
            
            if ~isempty(arg.fig)
                arg.plate.Fig=struct('fig',arg.fig,'Wells',{' '});
            end
            figure(arg.fig);
            clf
            
            %% remove empties from Val if needed
            if arg.removeempty
                if iscell(Val{1})
                    Val = cellfun(@(m) m(cellfun(@(x) ~isempty(x),m)),Val,'uniformoutput',0);
                else
                    Val = cellfun(@(m) m(~isnan(m)),Val,'uniformoutput',0);
                end
                Pos(cellfun(@isempty,Val))=[];
                Val(cellfun(@isempty,Val))=[];
            end
            
            %% if Val is only Char, make it into a cell array of cells
            if all(cellfun(@ischar,Val))
                Val = cellfun(@(m) {m},Val,'uniformoutput',0);
            end
            
            %% check to see that Val has single value per item:
            if ~all(cellfun(@(m) numel(m)==1,Val))
                errordlg('Cound not plot a heat map - need to have single value per well');
                error('Cound not plot a heat map - need to have single value per well');
            end
            
            %% if all values are numeric - draw a continous heatmap:
            Val = cat(1,Val{:});
            if isnumeric(Val)
                %%
                msk = nan(arg.plate.sz);
                for i=1:numel(Pos)
                    msk(ismember(arg.plate.Wells,Pos(i)))=Val(i);
                end
                msk=msk./max(msk(:));
                subplot('position',[0.1 0.1 0.7 0.8])
                arg.plate.plotHeatMap(msk,'colormap',arg.colormap);
                title(Type,'fontsize',13)
                subplot('position',[0.01 0.99 0.01 0.01])
                imagesc(unique(Val))
                set(gca,'xtick',[],'ytick',[])
                colorbar('position',[0.9 0.1 0.05 0.8])
            elseif all(cellfun(@ischar,Val))
                
                %%
                Val = regexprep(Val,'_',' ');
                msk = nan(arg.plate.sz);
                unq =unique(Val);
                for i=1:numel(unq)
                    ix = ismember(Val,unq{i});
                    msk(ismember(arg.plate.Wells,Pos(ix)))=i;
                end
                msk=msk./max(msk(:));
                msk_ix = gray2ind(msk,256);
                unq_ix = unique(msk_ix);
                unq_ix = setdiff(unq_ix,0);
                clr = arg.colormap;
                if ~isempty(arg.fig)
                    figure(arg.fig);
                end
                subplot('position',[0.1 0.1 0.7 0.8])
                arg.plate.plotHeatMap(msk,'colormap',arg.colormap);
                if strcmp(Type,'AllInputs')
                    alltypes = cellfun(@(m) [m ' '],MD.NewTypes,'Uniformoutput',0);
                    title(cat(2,alltypes{:}),'fontsize',13);
                else
                    title(Type,'fontsize',13)
                end
                subplot('position',[0.825 0.1 0.15 0.8])
                set(gca,'xtick',[],'ytick',[]);
                
                
                for i=1:numel(unq)
                    text(0.1,i/(numel(unq)+1),unq{i},'color',clr(unq_ix(i),:),'fontsize',25);
                end
            else
                errordlg('Cound not plot a heat map - must be either all numeric or all char');
                error('Cound not plot a heat map - must be either all numeric or all char');
            end
            
        end
        function possibleShiftingFrames = findAcqIndexes(MD, position, channel, varargin)
            arg.timefunc = @(t) true(size(t));
            arg = parseVarargin(varargin, arg);
            
            tbl=MD.tabulate('acq','Position',position,'Channel',channel,'timefunc',arg.timefunc);
            tbl=cat(1,tbl{:,2});
            possibleShiftingFrames=cumsum(tbl(1:end-1))+1;
        end
        
        % performs a crosstab operation to create a 2D table of counts for
        % type propertoies
        function [table,labels] = crosstab(MD,T1,T2,varargin)
            
            if ~isempty(varargin)
                %% transform varargin into T/V pair
                T = varargin(1:2:end);
                V = varargin(2:2:end);
                indx = getIndex(MD,T,V);
            else
                indx = 1:size(MD.Values,1);
            end
            
            
            V1 = MD.Values(indx,ismember(MD.Types,T1));
            if isnumeric(V1{1}) || islogical(V1{1})
                ix = cellfun(@isempty,V1);
                V1(ix)={NaN};
                V1=cat(1,V1{:});
            end
            V2 = MD.Values(indx,ismember(MD.Types,T2));
            if isnumeric(V2{1}) || islogical(V2{1})
                ix = cellfun(@isempty,V2);
                V2(ix)={NaN};
                V2=cat(1,V2{:});
            end
            [table,~,~,labels]=crosstab(V1,V2);
        end
        
        
        function makeTileConfig(MD)
            configFileName = fullfile(MD.pth,'TileConfiguration.txt');
            fid = fopen( configFileName, 'wt' );
            fprintf( fid, '%s\n', '# Define the number of dimensions we are working on');
            fprintf( fid, '%s\n', 'dim = 3');
            fprintf( fid, '%s\n', '# Define the image coordinates (in pixels)');
            
            Tiles = MD.unique('Tile');
            Channels = MD.unique('Channel');
            frames = MD.unique('frame');
            for indFrame = 1:numel(frames)
                for indCh=1:numel(Channels)
                    XYZ = [cell2mat(MD.getSpecificMetadata('XY','frame',frames(1),'Channel',Channels{indCh})), cell2mat(MD.getSpecificMetadata('Z','frame',frames(1),'Channel',Channels{indCh}))]./unique(cell2mat(MD.getSpecificMetadata('PixelSize')));
                    for indTile=1:numel(Tiles)
                        stkStr = [num2str((indCh-1)*numel(Tiles)+indTile-1) '; ' num2str(indFrame-1) '; (' num2str(XYZ(indTile,1)) ', ' num2str(XYZ(indTile,2)) ', ' num2str(XYZ(indTile,3)) ')' ];
                        fprintf( fid, '%s\n', stkStr);
                    end
                end
            end
            
            fclose(fid);
        end
        
    end
end
