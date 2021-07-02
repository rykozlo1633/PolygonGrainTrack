%% Particle tracking
%Track bidisperse grains (regular polygons or disks) using correlation template matching
%If shape or size distribution is complex/arbitrary, this will not be appropraite method.

%Directories of image, where to write new image with drawn grains, and where
% to write detected positions and orientations.
image_dir = 'ProcessedImages/LabFrame_';
imwrite_dir = 'ProcessedImages/TrackedGrains_'; %
poswrite_dir = 'Positions/';

%Set and image numbers to consider
startnum = 1;  endnum = startnum; %22000;
setnumber= 1; %can be range %%%*******1 = Pentagon image, 2 = Disk image in this example.
namesets = {'sample1','sample2'};
shapesets = {'Pentagon','Circle'};

%Size of frames
xsize = 1276; ysize = 700;
[X,Y] = meshgrid(1:xsize,1:ysize);

%Framerate of acquisition
fps = 120;
deltaT = 1/fps;

%This warning is unnecessary, will turn back on at end
warning('off','MATLAB:colon:nonIntegerIndex');

filt = (fspecial('gaussian', 7,1));

for index = setnumber
    
    name_dir = namesets{index};
    shape = shapesets{index};
    
    %Make directories for writing
    mkdir([imwrite_dir name_dir '/']);
    mkdir([poswrite_dir name_dir '/']);
    
    %Thresholds for large and small polygons/disks will depend on resolution
    % and size of grains; must test these numbers carefully.
    switch shape
        case 'Triangle'
            numsides = 3;
            largepolythresh = 520;
            smallpolythresh = 240;
        case 'Square'
            numsides = 4;
            largepolythresh = 600;
            smallpolythresh = 300;
        case 'Pentagon'
            numsides = 5;
            largepolythresh = 810;
            largepolythresh2 = 790;
            smallpolythresh = 420;
        case 'Hexagon'
            numsides = 6;
            largepolythresh = 915;
            smallpolythresh = 450;
        case 'Heptagon'
            numsides = 7;
            largepolythresh = 985;
            smallpolythresh = 490;
        case 'Circle'
            numsides = 8;
            largepolythresh = 1070;
            smallpolythresh = 650;
            showbars = true; %Find bars on disks to obtain rotations
    end
    
    %Load masks reprsenting the shapes of interest. A mask is binary.
    load(['ShapeMasks/' shape 'Masks.mat']); %SEE maskproduction.m for sample code to make masks.
    
    im_direc = dir([image_dir name_dir '/*.bmp']);
    
    for imindex = startnum:endnum
        
        
        %******************************************************
        %******************************************************
        %***Any pre-processing of images should happen here,***
        %***including binarization of image with grains set ***
        %***to 1 and non-grains to 0. Edge subtraction with ***
        %***e.g., Sobel, may help with edge-edge contacts******
        %******************************************************
        %******************************************************
        bin_n = imread([image_dir name_dir '/' im_direc(imindex).name]);
        % ^ obtained by image pre-processing, which could happen above.
        
        %If disks (or other shapes) marked with bar, either process image and obtain
        % binarized bar image, or import bar image from elsewhere (as done in this sample)
        if shape(1) == 'C' && showbars
            bars = imread(['ProcessedImages/ProcessedDiskImage_Bars.bmp']);
        end
        
        final_n = int8(bin_n);
        final_n(final_n==0) = -1; %final_n is binarized image; -1 assigned to 0, int8.
        
        %Convolve with large disks
        test_xc2 = conv2(double(final_n),double(maskdisk_l));
        test_xc2 = medfilt2(test_xc2,[3,3]);
        test_xc2 = conv2(test_xc2,filt,'same');
        padding = (size(maskdisk_l)-1)/2;
        test_xc2( 1:padding(1)+25,:) = 0; test_xc2( end-padding(1)-25+1:end,:) = 0;
        test_xc2( :,1:padding(2)+25) = 0; test_xc2( :,end-padding(2)-25+1:end) = 0;
        
        if numsides>0
            %Convolve with large polygon
            padding2 = (size(maskpoly_l(:,:,1))-1)/2;
            test_xc2_poly = single(zeros([size(final_n)+padding2*2 size(maskpoly_l,3)])); %Also a GB!!
            for layer1 = 1:size(maskpoly_l,3)
                test_xc2_poly(:,:,layer1) = conv2(final_n,maskpoly_l(:,:,layer1));
                test_xc2_poly(:,:,layer1) = medfilt2(test_xc2_poly(:,:,layer1),[3,3]);
                test_xc2_poly(:,:,layer1) = conv2((test_xc2_poly(:,:,layer1)),filt,'same');
            end
            test_xc2_poly( 1:padding2(1)+25,:,:) = 0; test_xc2_poly( end-padding2(1)-25+1:end,:,:) = 0;
            test_xc2_poly( :,1:padding2(2)+25,:) = 0; test_xc2_poly( :,end-padding2(2)-25+1:end,:) = 0;
            
        end
        
        %Arrays for large and small grains
        partcents_l = []; %must be grown because number of particles unknown
        partcents_s = []; %must be grown because number of particles unknown
        
        %Find peaks in large polygon convolution
        if numsides>0
            %Polygons
            while 1 %Could insert conditional: If number of grains is known, then
                        %stop search once this number has been detected
                        %(potentially within some error, using a threshold
                        %of correlation peak)
                        
                [maxs,ymaxs] = max(test_xc2_poly);
                [maxs2,xmaxs] = max(maxs);
                [testmax,zmax] = max(maxs2);
                xmax_orig = xmaxs(zmax);
                ymax_orig = ymaxs(1,xmax_orig,zmax);
                
                ex = test_xc2_poly(ymax_orig-center_l_poly+1:ymax_orig-center_l_poly+...
                    size_maskpoly_l,xmax_orig-center_l_poly+1:xmax_orig-center_l_poly+size_maskpoly_l,:);
                ex(ex<0) = 0;
                
                zmaxes = NaN(size(ex,3),1);
                for number = 1:size(ex,3)
                    try
                        zmaxes(number) = max(ex(:,:,number),[],'all');
                    catch %Error if running earlier than R2018b
                        szex = numel(ex(:,:,number));
                        zmaxes(number) = max( reshape(ex(:,:,number),[szex,1]) );
                    end
                        
                end
                zmaxes_sm = smoothdata(zmaxes,'Gaussian',5);
                %figure(1);plot(zmaxes_sm)
                [~,zmax2] = max(zmaxes_sm);
                zmax = zmax2;
                [maxs,ymaxs] = max(ex(:,:,zmax));
                [testmax,xmaxs] = max(maxs);
                xmax = xmaxs;
                ymax = ymaxs(1,xmax);
                
                %Option 2: find centroid around max point. (Other peak-finding method could be used here)
                shiftz = 0;
                if zmax == 2; shiftz = 1; elseif zmax == 1; shiftz= 2;end;
                if zmax == maxang+1; shiftz = -2; elseif zmax == maxang; shiftz= -1;end;
                shiftx = 0;
                if xmax == 2; shiftx = 1; elseif xmax == 1; shiftx= 2;end;
                if xmax == size(ex,2); shiftx = -2; elseif xmax == size(ex,2)-1; shiftx= -1;end;
                shifty = 0;
                if ymax == 2; shifty = 1; elseif ymax == 1; shifty= 2;end;
                if ymax == size(ex,1); shifty = -2; elseif ymax == size(ex,1)-1; shifty= -1;end;
                
                ex = circshift(ex,[shifty,shiftx,shiftz]);
                zmax = zmax+shiftz; xmax = xmax+shiftx; ymax = ymax+shifty;
                com_reg = ex(ymax-2:ymax+2,xmax-2:xmax+2,zmax-2:zmax+2);
                M = sum(com_reg(:));
                for dims = 1:3
                    shp = ones(1,3);
                    shp(dims) = size(com_reg,dims);
                    rep = size(com_reg);
                    rep(dims) = 1;
                    ind = repmat(reshape(1:size(com_reg,dims),shp),rep);
                    C(dims) = sum(ind(:).*com_reg(:))./M;
                end
                C = C+[ymax-2-1, xmax-2-1, zmax-2-1];
                ex = circshift(ex,-[shifty,shiftx,shiftz]);
                xmax = C(2)-shiftx; ymax = C(1)-shifty; zmax = C(3)-shiftz;
                
                ymax = ymax + ymax_orig - center_l_poly;
                xmax = xmax + xmax_orig - center_l_poly;
                
                %remove self-convolved mask from convolution image
                test_xc2_poly(ymax-center_l_poly+1:ymax-center_l_poly+size_maskpoly_l,xmax-center_l_poly+1:xmax-center_l_poly+size_maskpoly_l,:) = ...
                    test_xc2_poly(ymax-center_l_poly+1:ymax-center_l_poly+size_maskpoly_l,xmax-center_l_poly+1:xmax-center_l_poly+size_maskpoly_l,:)...
                    - selfmask.poly_l{ round(zmax) };
                
                %Once detected grain is below threshold, stop searching for large grains
                if testmax<= largepolythresh
                    break
                end
                partcents_l = [partcents_l; xmax-padding2(2) ,ymax-padding2(1), (zmax-1)*angular_res];
                
            end
            fprintf('     Finished with large polygon convolution! \n')
        else
            
            while 1 %Could insert conditional: If number of grains is known, then
                        %stop search once this number has been detected
                        %(potentially within some error, using a threshold
                        %of correlation peak)
                %Disks
                
                [maxs,ymaxs] = max(test_xc2);
                [testmax,xmax_orig] = max(maxs);
                ymax_orig = ymaxs(xmax_orig);
                ex = test_xc2(ymax_orig-center_l_disk+1:ymax_orig-center_l_disk+size_maskdisk_l,xmax_orig-center_l_disk+1:xmax_orig-center_l_disk+size_maskdisk_l,:);
                
                ex(ex<0) = 0;
                
                [maxs,ymaxs] = max(ex);
                [~,xmax] = max(maxs);
                ymax = ymaxs(xmax);
                %Option 2: find COM around max point. (Other peak-finding method could be used)
                shiftx = 0;
                if xmax == 2; shiftx = 1; elseif xmax == 1; shiftx= 2;end;
                if xmax == size(ex,2); shiftx = -2; elseif xmax == size(ex,2)-1; shiftx= -1;end;
                shifty = 0;
                if ymax == 2; shifty = 1; elseif ymax == 1; shifty= 2;end;
                if ymax == size(ex,1); shifty = -2; elseif ymax == size(ex,1)-1; shifty= -1;end;
                
                ex = circshift(ex,[shifty,shiftx]);
                xmax = xmax+shiftx; ymax = ymax+shifty;
                com_reg = ex(ymax-2:ymax+2,xmax-2:xmax+2);
                M = sum(com_reg(:));
                C = NaN(1,2);
                for dims = 1:2
                    shp = ones(1,2); shp(dims) = size(com_reg,dims);
                    rep = size(com_reg); rep(dims) = 1;
                    ind = repmat(reshape(1:size(com_reg,dims),shp),rep);
                    C(dims) = sum(ind(:).*com_reg(:))./M;
                end
                C = C+[ymax-2-1, xmax-2-1];
                %ex = circshift(ex,-[shifty,shiftx]);
                ymax = C(2) - shiftx + ymax_orig-center_l_disk;
                xmax = C(1) - shifty + xmax_orig-center_l_disk;
                
                test_xc2(ymax-center_l_disk+1:ymax-center_l_disk+size_maskdisk_l,xmax-center_l_disk+1:xmax-center_l_disk+size_maskdisk_l,:) = ...
                    test_xc2(ymax-center_l_disk+1:ymax-center_l_disk+size_maskdisk_l,xmax-center_l_disk+1:xmax-center_l_disk+size_maskdisk_l,:)...
                    - selfmask.disk_l;
                
                if testmax<=largepolythresh
                    break
                end
                partcents_l = [partcents_l; xmax-padding(2) ,ymax-padding(1)];
            end
            
            fprintf('     Finished with large disk convolution! \n')
        end
        
        %Draw image of detected large grains and remove from final_n before small grain detection
        poly_rad = sqrt(innerrad_l^2 + size_large^2/4);
        if numsides==3
            vertices = NaN(length(partcents_l),2*numsides);
            for ll = 1:length(partcents_l(:,1))
                v1 = [partcents_l(ll,1),partcents_l(ll,2)] + [0 , -poly_rad];
                v1 = [[cosd(partcents_l(ll,3)) -sind(partcents_l(ll,3)); sind(partcents_l(ll,3)) cosd(partcents_l(ll,3))] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v2 = [[cosd(120) -sind(120); sind(120) cosd(120)] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v3 = [[cosd(120) -sind(120); sind(120) cosd(120)] * [v2 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                vertices(ll,:) = [v1,v2,v3];
            end
        elseif numsides==4
            vertices = NaN(length(partcents_l),2*numsides);
            for ll = 1:length(partcents_l(:,1))
                v1 = [partcents_l(ll,1),partcents_l(ll,2)] + [poly_rad*cosd(45) , -poly_rad*sind(45)];
                v1 = [[cosd(partcents_l(ll,3)) -sind(partcents_l(ll,3)); sind(partcents_l(ll,3)) cosd(partcents_l(ll,3))] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v2 = [[cosd(90) -sind(90); sind(90) cosd(90)] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v3 = [[cosd(90) -sind(90); sind(90) cosd(90)] * [v2 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v4 = [[cosd(90) -sind(90); sind(90) cosd(90)] * [v3 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                vertices(ll,:) = [v1,v2,v3,v4];
            end
        elseif numsides==5
            vertices = NaN(length(partcents_l),2*numsides);
            for ll = 1:length(partcents_l(:,1))
                v1 = [partcents_l(ll,1),partcents_l(ll,2)] + [0 , -poly_rad];
                v1 = [[cosd(partcents_l(ll,3)) -sind(partcents_l(ll,3)); sind(partcents_l(ll,3)) cosd(partcents_l(ll,3))] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v2 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v3 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v2 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v4 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v3 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v5 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v4 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                vertices(ll,:) = [v1,v2,v3,v4,v5];
            end
        elseif numsides == 6
            vertices = NaN(length(partcents_l),2*numsides);
            for ll = 1:length(partcents_l(:,1))
                v1 = [partcents_l(ll,1),partcents_l(ll,2)] + [poly_rad*cosd(60) , -poly_rad*sind(60)];
                v1 = [[cosd(partcents_l(ll,3)) -sind(partcents_l(ll,3)); sind(partcents_l(ll,3)) cosd(partcents_l(ll,3))] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v2 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v3 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v2 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v4 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v3 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v5 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v4 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v6 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v5 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                vertices(ll,:) = [v1,v2,v3,v4,v5,v6];
            end
        elseif numsides == 7
            vertices = NaN(length(partcents_l),2*numsides);
            for ll = 1:length(partcents_l(:,1))
                v1 = [partcents_l(ll,1),partcents_l(ll,2)] + [0 , -poly_rad];
                v1 = [[cosd(partcents_l(ll,3)) -sind(partcents_l(ll,3)); sind(partcents_l(ll,3)) cosd(partcents_l(ll,3))] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v2 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v1 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v3 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v2 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v4 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v3 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v5 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v4 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v6 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v5 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                v7 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v6 - [partcents_l(ll,1),partcents_l(ll,2)]]']' + [partcents_l(ll,1),partcents_l(ll,2)];
                vertices(ll,:) = [v1,v2,v3,v4,v5,v6,v7];
            end
        end
        
        %Subtract large pentagons from original image before searching for small.
        if numsides>0
            %old method: use insertShape (requires Computer Vision toolbox and isn't great
            %  with subpixel centroids -- it's more for visual than quantitative analysis)
            %subtraction_image = insertShape(single(bin_n),'FilledPolygon',vertices,'Color',zeros(length(vertices),3),'Opacity',1);
            %new method: draw the polygons yourself!
            subtraction_image = single(bin_n);
            for dl = 1:size(partcents_l,1)
                dx = partcents_l(dl,1); dy = partcents_l(dl,2); dr = ceil(poly_rad);
                if floor(dy-dr)<=0; dy = dy-floor(dy-dr)+1; end
                if floor(dx-dr)<=0; dx = dx-floor(dx-dr)+1; end
                if ceil(dy+dr)>=size(bin_n,1); dy = dy-(ceil(dy+dr)-size(bin_n,1))-1; end
                if ceil(dx+dr)>=size(bin_n,2); dx = dx-(ceil(dx+dr)-size(bin_n,2))-1; end
                subsampx = X( floor(dy-dr):ceil(dy+dr), floor(dx-dr):ceil(dx+dr) );
                subsampy = Y( floor(dy-dr):ceil(dy+dr), floor(dx-dr):ceil(dx+dr) );
                subsampx2 = subsampx(:); subsampy2 = subsampy(:);
                in = inpolygon( subsampx2, subsampy2, vertices(dl,[1:2:end,1])', vertices(dl,[2:2:end,2])' );
                allin = [subsampx2(in), subsampy2(in)];
                for q = 1:length(allin)
                    subtraction_image( allin(q,2), allin(q,1)) = 0;
                end
            end
        else
            %subtraction_image = insertShape(single(bin_n),'FilledCircle',[partcents_l,ones(length(partcents_l),1)*(size_large)/2],'Color',zeros(length(partcents_l),3),'Opacity',1);
            subtraction_image = single(bin_n);
            for pl = 1:size(partcents_l,1)
                px = partcents_l(pl,1); py = partcents_l(pl,2); pr = size_large/2;
                for xx = floor(px-pr):ceil(px+pr)
                    for yy = floor(py-pr):ceil(py+pr)
                        if (xx-px)^2+(yy-py)^2<=pr^2
                            subtraction_image(yy,xx) = 0;
                        end
                    end
                end
            end
        end
        %subtraction_image = subtraction_image(:,:,1);
        subtraction_image(subtraction_image<1) = -1;
        
        %Convolve with small disk
        test_xc2 = conv2(double(subtraction_image),double(maskdisk_s));
        test_xc2 = medfilt2(test_xc2,[3,3]);
        test_xc2 = conv2(test_xc2,filt,'same');
        padding = (size(maskdisk_s)-1)/2;
        test_xc2( 1:padding(1)+25,:) = 0; test_xc2( end-padding(1)-25+1:end,:) = 0;
        test_xc2( :,1:padding(2)+25) = 0; test_xc2( :,end-padding(2)-25+1:end) = 0;
        
        if numsides>0
            %Convolve with small polygon
            padding2 = (size(maskpoly_s(:,:,1))-1)/2;
            test_xc2_poly = single(zeros([size(subtraction_image)+padding2*2 size(maskpoly_s,3)])); %Also a GB!!
            for layer1 = 1:size(maskpoly_s,3)
                test_xc2_poly(:,:,layer1) = conv2(subtraction_image,maskpoly_s(:,:,layer1));
                test_xc2_poly(:,:,layer1) = medfilt2(test_xc2_poly(:,:,layer1),[3,3]);
                test_xc2_poly(:,:,layer1) = conv2((test_xc2_poly(:,:,layer1)),filt,'same');
            end
            test_xc2_poly( 1:padding2(1)+25,:,:) = 0; test_xc2_poly( end-padding2(1)-25+1:end,:,:) = 0;
            test_xc2_poly( :,1:padding2(2)+25,:) = 0; test_xc2_poly( :,end-padding2(2)-25+1:end,:) = 0;
        end
        
        %Find peaks in small polygon convolution
        
        if numsides>0
            while 1 %Could insert conditional: If number of grains is known, then
                        %stop search once this number has been detected
                        %(potentially within some error, using a threshold
                        %of correlation peak)
                
                [maxs,ymaxs] = max(test_xc2_poly);
                [maxs2,xmaxs] = max(maxs);
                [testmax,zmax] = max(maxs2);
                xmax_orig = xmaxs(zmax);
                ymax_orig = ymaxs(1,xmax_orig,zmax);
                
                ex = test_xc2_poly(ymax_orig-center_s_poly+1:ymax_orig-center_s_poly+size_maskpoly_s,xmax_orig-center_s_poly+1:xmax_orig-center_s_poly+size_maskpoly_s,:);
                ex(ex<0) = 0;
                
                zmaxes = NaN(size(ex,3),1);
                for number = 1:size(ex,3)
                    try
                        zmaxes(number) = max(ex(:,:,number),[],'all');
                    catch %Error if running earlier than R2018b
                        szex = numel(ex(:,:,number));
                        zmaxes(number) = max( reshape(ex(:,:,number),[szex,1]) );
                    end
                end
                zmaxes_sm = smoothdata(zmaxes,'Gaussian',5);
                %figure(1);plot(zmaxes_sm)
                [~,zmax2] = max(zmaxes_sm);
                zmax = zmax2;
                [maxs,ymaxs] = max(ex(:,:,zmax));
                [testmax,xmaxs] = max(maxs);
                xmax = xmaxs;
                ymax = ymaxs(1,xmax);
                
                if testmax<=smallpolythresh
                    break
                end
                
                %Option 2: find COM around max point.
                shiftz = 0;
                if zmax == 2; shiftz = 1; elseif zmax == 1; shiftz= 2;end;
                if zmax == maxang+1; shiftz = -2; elseif zmax == maxang; shiftz= -1;end;
                shiftx = 0;
                if xmax == 2; shiftx = 1; elseif xmax == 1; shiftx= 2;end;
                if xmax == size(ex,2); shiftx = -2; elseif xmax == size(ex,2)-1; shiftx= -1;end;
                shifty = 0;
                if ymax == 2; shifty = 1; elseif ymax == 1; shifty= 2;end;
                if ymax == size(ex,1); shifty = -2; elseif ymax == size(ex,1)-1; shifty= -1;end;
                
                ex = circshift(ex,[shifty,shiftx,shiftz]);
                zmax = zmax+shiftz; xmax = xmax+shiftx; ymax = ymax+shifty;
                com_reg = ex(ymax-2:ymax+2,xmax-2:xmax+2,zmax-2:zmax+2);
                M = sum(com_reg(:));
                C = [0,0,0];
                for dims = 1:3
                    shp = ones(1,3);
                    shp(dims) = size(com_reg,dims);
                    rep = size(com_reg);
                    rep(dims) = 1;
                    ind = repmat(reshape(1:size(com_reg,dims),shp),rep);
                    C(dims) = sum(ind(:).*com_reg(:))./M;
                end
                C = C+[ymax-2-1, xmax-2-1, zmax-2-1];
                ex = circshift(ex,-[shifty,shiftx,shiftz]);
                xmax = C(2)-shiftx; ymax = C(1)-shifty; zmax = C(3)-shiftz;
                
                ymax = ymax + ymax_orig -center_s_poly;
                xmax = xmax + xmax_orig -center_s_poly;
                
                test_xc2_poly(ymax-center_s_poly+1:ymax-center_s_poly+size_maskpoly_s,xmax-center_s_poly+1:xmax-center_s_poly+size_maskpoly_s,:) = ...
                    test_xc2_poly(ymax-center_s_poly+1:ymax-center_s_poly+size_maskpoly_s,xmax-center_s_poly+1:xmax-center_s_poly+size_maskpoly_s,:)...
                    - selfmask.poly_s{ round(zmax) };
                
                
                partcents_s = [partcents_s; xmax-padding2(2) ,ymax-padding2(1), (zmax-1)*angular_res];
                
            end
            fprintf('     Finished with small polygon convolution! \n')
        else
            allmax = [];
            while 1 %Could insert conditional: If number of grains is known, then
                        %stop search once this number has been detected
                        %(potentially within some error, using a threshold
                        %of correlation peak)
                %Disks
                
                [maxs,ymaxs] = max(test_xc2);
                [testmax,xmax_orig] = max(maxs);
                ymax_orig = ymaxs(xmax_orig);
                allmax = [allmax;testmax];
                ex = test_xc2(ymax_orig-center_s_disk+1:ymax_orig-center_s_disk+size_maskdisk_s,xmax_orig-center_s_disk+1:xmax_orig-center_s_disk+size_maskdisk_s,:);
                ex(ex<0) = 0;
                
                [maxs,ymaxs] = max(ex);
                [~,xmax] = max(maxs);
                ymax = ymaxs(xmax);
                %Option 2: find COM around max point.
                shiftx = 0;
                if xmax == 2; shiftx = 1; elseif xmax == 1; shiftx= 2;end;
                if xmax == size(ex,2); shiftx = -2; elseif xmax == size(ex,2)-1; shiftx= -1;end;
                shifty = 0;
                if ymax == 2; shifty = 1; elseif ymax == 1; shifty= 2;end;
                if ymax == size(ex,1); shifty = -2; elseif ymax == size(ex,1)-1; shifty= -1;end;
                
                ex = circshift(ex,[shifty,shiftx]);
                xmax = xmax+shiftx; ymax = ymax+shifty;
                com_reg = ex(ymax-2:ymax+2,xmax-2:xmax+2);
                M = sum(com_reg(:));
                C = NaN(1,2);
                for dims = 1:2
                    shp = ones(1,2); shp(dims) = size(com_reg,dims);
                    rep = size(com_reg); rep(dims) = 1;
                    ind = repmat(reshape(1:size(com_reg,dims),shp),rep);
                    C(dims) = sum(ind(:).*com_reg(:))./M;
                end
                C = C+[ymax-2-1, xmax-2-1];
                %ex = circshift(ex,-[shifty,shiftx]);
                ymax = C(2) - shiftx + ymax_orig-center_s_disk;
                xmax = C(1) - shifty + xmax_orig-center_s_disk;
                
                test_xc2(ymax-center_s_disk+1:ymax-center_s_disk+size_maskdisk_s,xmax-center_s_disk+1:xmax-center_s_disk+size_maskdisk_s,:) = ...
                    test_xc2(ymax-center_s_disk+1:ymax-center_s_disk+size_maskdisk_s,xmax-center_s_disk+1:xmax-center_s_disk+size_maskdisk_s,:)...
                    - selfmask.disk_s;
                
                if testmax<=smallpolythresh
                    break
                end
                partcents_s = [partcents_s; xmax-padding(2) ,ymax-padding(1)];
            end
            
            fprintf('     Finished with small disk convolution! \n')
        end
        
        
        %Find vertices of small polygons for drawing on final image
        poly_rad = sqrt(innerrad_s^2 + size_small^2/4);
        if numsides==3
            vertices2 = NaN(length(partcents_s),2*numsides);
            for ll = 1:length(partcents_s(:,1))
                v1 = [partcents_s(ll,1),partcents_s(ll,2)] + [0 , -poly_rad];
                v1 = [[cosd(partcents_s(ll,3)) -sind(partcents_s(ll,3)); sind(partcents_s(ll,3)) cosd(partcents_s(ll,3))] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v2 = [[cosd(120) -sind(120); sind(120) cosd(120)] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v3 = [[cosd(120) -sind(120); sind(120) cosd(120)] * [v2 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                vertices2(ll,:) = [v1,v2,v3];
            end
        elseif numsides==4
            vertices2 = NaN(length(partcents_s),2*numsides);
            for ll = 1:length(partcents_s(:,1))
                v1 = [partcents_s(ll,1),partcents_s(ll,2)] + [poly_rad*cosd(45) , -poly_rad*sind(45)];
                v1 = [[cosd(partcents_s(ll,3)) -sind(partcents_s(ll,3)); sind(partcents_s(ll,3)) cosd(partcents_s(ll,3))] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v2 = [[cosd(90) -sind(90); sind(90) cosd(90)] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v3 = [[cosd(90) -sind(90); sind(90) cosd(90)] * [v2 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v4 = [[cosd(90) -sind(90); sind(90) cosd(90)] * [v3 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                vertices2(ll,:) = [v1,v2,v3,v4];
            end
        elseif numsides==5
            vertices2 = NaN(length(partcents_s),2*numsides);
            for ll = 1:length(partcents_s(:,1))
                v1 = [partcents_s(ll,1),partcents_s(ll,2)] + [0 , -poly_rad];
                v1 = [[cosd(partcents_s(ll,3)) -sind(partcents_s(ll,3)); sind(partcents_s(ll,3)) cosd(partcents_s(ll,3))] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v2 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v3 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v2 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v4 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v3 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v5 = [[cosd(72) -sind(72); sind(72) cosd(72)] * [v4 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                vertices2(ll,:) = [v1,v2,v3,v4,v5];
            end
        elseif numsides == 6
            vertices2 = NaN(length(partcents_s),2*numsides);
            for ll = 1:length(partcents_s(:,1))
                v1 = [partcents_s(ll,1),partcents_s(ll,2)] + [poly_rad*cosd(60) , -poly_rad*sind(60)];
                v1 = [[cosd(partcents_s(ll,3)) -sind(partcents_s(ll,3)); sind(partcents_s(ll,3)) cosd(partcents_s(ll,3))] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v2 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v3 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v2 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v4 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v3 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v5 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v4 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v6 = [[cosd(60) -sind(60); sind(60) cosd(60)] * [v5 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                vertices2(ll,:) = [v1,v2,v3,v4,v5,v6];
            end
        elseif numsides == 7
            vertices2 = NaN(length(partcents_s),2*numsides);
            for ll = 1:length(partcents_s(:,1))
                v1 = [partcents_s(ll,1),partcents_s(ll,2)] + [0 , -poly_rad];
                v1 = [[cosd(partcents_s(ll,3)) -sind(partcents_s(ll,3)); sind(partcents_s(ll,3)) cosd(partcents_s(ll,3))] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v2 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v1 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v3 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v2 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v4 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v3 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v5 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v4 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v6 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v5 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                v7 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v6 - [partcents_s(ll,1),partcents_s(ll,2)]]']' + [partcents_s(ll,1),partcents_s(ll,2)];
                vertices2(ll,:) = [v1,v2,v3,v4,v5,v6,v7];
            end
        end
        
        %***Next step uses insertShape to draw, requires Computer Vision toolbox.
        %***Alternative: use plot for the polygons, and viscircles or plot for the disks, then print figure
        if numsides>0
            final_n_poly = insertShape(uint8(cat(3,bin_n,bin_n,bin_n))*255, 'FilledPolygon', vertices,'opacity',0.5,'color',repmat([255,0,0],[size(vertices,1),1]));
            final_n_poly = insertShape(final_n_poly, 'FilledPolygon', vertices2,'opacity',0.5,'color',repmat([0,0,255],[size(vertices2,1),1]));
        elseif numsides==0
            %For disks, if bar drawn on surface, use bar to detect rotation
            
            [Xmesh,Ymesh] = meshgrid(1:size(bin_n,2),1:size(bin_n,1));
            r_s = size_small/2;
            r_l = size_large/2;
            if size(partcents_s,2)==2
                partcents_s = [partcents_s, NaN(length(partcents_s),1)];
                partcents_l = [partcents_l, NaN(length(partcents_l),1)];
            end
            
            for indsmall = 1:length(partcents_s)
                botx = floor(partcents_s(indsmall,1)-r_s);
                boty = floor(partcents_s(indsmall,2)-r_s);
                topx = ceil(partcents_s(indsmall,1)+r_s);
                topy = ceil(partcents_s(indsmall,2)+r_s);
                if botx<1; botx = 1; end
                if topx>size(bin_n,2); topx = size(bin_n,2); end
                if boty<1; boty = 1; end
                if topy>size(bin_n,1); topy = size(bin_n,1); end
                locX = Xmesh(boty:topy,botx:topx);
                locY = Ymesh(boty:topy,botx:topx);
                R = zeros(size(locX));
                R( (locX-partcents_s(indsmall,1)).^2 + (locY-partcents_s(indsmall,2)).^2 < (r_s-3)^2 ) = 1;
                mask_temp = NaN(size(R));
                mask_temp(logical(R)) = 1;
                locbin = double(bars(boty:topy,botx:topx)).*double(mask_temp);
                locbin(isnan(mask_temp)) = 0;
                locbin = bwareaopen(locbin,10);
                regions = regionprops(locbin,'Area','Orientation','Centroid');
                if ~isempty(regions)
                    if size(regions,1)>1
                        regionpix = regionprops(locbin,'PixelIdxList');
                        regionpix = extractfield( regionpix, 'PixelIdxList');
                        
                        [regionpix_row,regionpix_col] = ind2sub(size(locbin), regionpix');
                        centroid = [mean(regionpix_col), mean(regionpix_row)];
                        orientations = regorientation( regionpix_col, regionpix_row, centroid );
                        areas = numel(regionpix);
                        com = centroid;
                    else
                        areas = regions.Area;
                        orientations = regions.Orientation;
                        com = regions.Centroid;
                    end
                else
                    orientations = NaN;
                end
                partcents_s(indsmall,3) = -orientations; %ranges -90 to 90 with x-axis
                
            end
            for indlarge = 1:length(partcents_l)
                botx = floor(partcents_l(indlarge,1)-r_l);
                boty = floor(partcents_l(indlarge,2)-r_l);
                topx = ceil(partcents_l(indlarge,1)+r_l);
                topy = ceil(partcents_l(indlarge,2)+r_l);
                if botx<1; botx = 1; end
                if topx>size(bin_n,2); topx = size(bin_n,2); end
                if boty<1; boty = 1; end
                if topy>size(bin_n,1); topy = size(bin_n,1); end
                locX = Xmesh(boty:topy,botx:topx);
                locY = Ymesh(boty:topy,botx:topx);
                R = zeros(size(locX));
                R( (locX-partcents_l(indlarge,1)).^2 + (locY-partcents_l(indlarge,2)).^2 < (r_l-3)^2 ) = 1;
                mask_temp = NaN(size(R));
                mask_temp(logical(R)) = 1;
                locbin = double(bars(boty:topy,botx:topx)).*double(mask_temp);
                locbin(isnan(mask_temp)) = 0;
                locbin = bwareaopen(locbin,10);
                regions = regionprops(locbin,'Area','Orientation','Centroid');
                if ~isempty(regions)
                    if size(regions,1)>1
                        regionpix = regionprops(locbin,'PixelIdxList');
                        regionpix = extractfield( regionpix, 'PixelIdxList');
                        
                        [regionpix_row,regionpix_col] = ind2sub(size(locbin), regionpix');
                        centroid = [mean(regionpix_col), mean(regionpix_row)];
                        orientations = regorientation( regionpix_col, regionpix_row, centroid );
                        areas = numel(regionpix);
                    else
                        areas = regions.Area;
                        orientations = regions.Orientation;
                        com = regions.Centroid;
                    end
                else
                    orientations = NaN;
                end
                partcents_l(indlarge,3) = -orientations; %ranges -90 to 90 with x-axis
                
            end
            
            final_n_poly = insertShape(uint8(cat(3,bin_n,bin_n,bin_n))*255, 'FilledCircle', [partcents_l(:,1:2), ones(length(partcents_l),1)*size_large/2],'opacity',0.5,'color',repmat([255,0,0],[size(partcents_l,1),1]));
            final_n_poly = insertShape(final_n_poly, 'FilledCircle', [partcents_s(:,1:2), ones(length(partcents_s),1)*size_small/2],'opacity',0.5,'color',repmat([0,0,255],[size(partcents_s,1),1]));
            %Draw orientation lines
            lineends_l = [partcents_l(:,1)+size_large/2*cosd(partcents_l(:,3)), ...
                partcents_l(:,2)+size_large/2*sind(partcents_l(:,3)), ...
                partcents_l(:,1)-size_large/2*cosd(partcents_l(:,3)), ...
                partcents_l(:,2)-size_large/2*sind(partcents_l(:,3)) ];
            lineends_s = [partcents_s(:,1)+size_small/2*cosd(partcents_s(:,3)), ...
                partcents_s(:,2)+size_small/2*sind(partcents_s(:,3)), ...
                partcents_s(:,1)-size_small/2*cosd(partcents_s(:,3)), ...
                partcents_s(:,2)-size_small/2*sind(partcents_s(:,3)) ];
            lineends_l(isnan(lineends_l(:,1)),:) = [];
            lineends_s(isnan(lineends_s(:,1)),:) = [];
            final_n_poly = insertShape(final_n_poly,'line', lineends_l, 'opacity', 1,'color',repmat([1,1,0.3]*255,[size(lineends_l,1),1]),'linewidth',2);
            final_n_poly = insertShape(final_n_poly,'line', lineends_s, 'opacity', 1,'color',repmat([1,1,0.3]*255,[size(lineends_s,1),1]),'linewidth',2);
            
        end
        
        figure(101);imshow(final_n_poly);
        
        imwrite(final_n_poly,[imwrite_dir name_dir '/' im_direc(imindex).name(1:end-4) '.jpg'],'quality',50);
        
        dlmwrite([poswrite_dir name_dir '/' im_direc(imindex).name(1:end-4) '.txt'], ...
            [[partcents_s, (size_small)*ones(length(partcents_s),1)]; [partcents_l, ...
            (size_large)*ones(length(partcents_l),1)]],'delimiter', ' ','Precision',9);
        
    end
end

%clean up
warning('on','MATLAB:colon:nonIntegerIndex');
clearvars
%close all
clc
