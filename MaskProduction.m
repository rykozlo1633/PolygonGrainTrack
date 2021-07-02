%% Mask production
%Produce regular polygon masks and save to directory to be loaded in GrainTrack.m 
%Current code assume bidisperse packing, so 2 mask sizes.

clearvars selfmask
savemask_dir = 'ShapeMasks/';

%Pick shape
shape = 'Pentagon';

switch shape
    %All parameters below will depend on the mask you need given image resolution,
    %grain size, etc.
    case 'Triangle'
        numsides = 3;
        size_layer_l = 2; %outer -1 shell; must be in units of pixels
        size_layer_s = 2;
        angular_res = 1; %pixels 1
        maxang = 360/numsides-1;
        size_small = 30; %pixels (side length)
        size_large = 40; %pixels (side length)
    case 'Square'
        numsides = 4;
        size_layer_l = 2; %outer -1 shell; must be in units of pixels
        size_layer_s = 2;
        angular_res = 1; %pixels 1
        maxang = 360/numsides-1;
        size_small = 22; %pixels (side length)
        size_large = 30.5; %pixels (side length)
    case 'Pentagon'
        numsides = 5;
        size_layer_l = 2; %outer -1 shell; must be in units of pixels
        size_layer_s = 2;
        angular_res = 1; %pixels 1
        maxang = 360/numsides-1;
        size_small = 18; %pixels (side length)
        size_large = 24; %pixels (side length)
        size_small_poly = size_small;
        size_large_poly = size_large;
    case 'Hexagon'
        numsides = 6;
        size_layer_l = 2; %outer -1 shell; must be in units of pixels
        size_layer_s = 2;
        angular_res = 1; %pixels 1
        maxang = 360/numsides-1;
        size_small = 16; %pixels (side length)
        size_large = 20.5; %pixels (side length)
    case 'Heptagon'
        numsides = 7;
        size_layer_l = 2; %outer -1 shell; must be in units of pixels
        size_layer_s = 2;
        angular_res = 1; %pixels 1
        maxang = floor(360/numsides);
        size_small = 14.0; %pixels (side length)
        size_large = 18.0; %pixels (side length)
    case 'Circle'
        numsides = 0;
        size_layer_l = 2;% %outer -1 shell; must be in units of pixels
        size_layer_s = 2;%
        size_small = 33; %pixels (diameter) 
        size_large = 41; %pixels (diameter) 
end

%Obtain Inner radius and size of polygon mask
[maskpoly_s_temp, innerrad_s] = polygonmask(numsides,size_small,0);
[maskpoly_l_temp, innerrad_l] = polygonmask(numsides,size_large,0);
%figure;imshow(maskpoly_l_temp);

%Disk masks
maskdisk_s = polygonmask(0,ceil(2*innerrad_s),0); %if pixel units, the ceil not necessary
maskdisk_s = padarray(maskdisk_s,[size_layer_s,size_layer_s],0,'both');
maskdisk_l = polygonmask(0,ceil(2*innerrad_l),0);
maskdisk_l = padarray(maskdisk_l,[size_layer_l,size_layer_l],0,'both');

%Create mask for outer -1 boundary
maskdisk_s_out = polygonmask(0,ceil(2*innerrad_s)+2*size_layer_s,0); %if pixel units, the ceil not necessary
maskdisk_l_out = polygonmask(0,ceil(2*innerrad_l)+2*size_layer_l,0);

%Hold on to original mask
maskdisk_s_temp = maskdisk_s;
maskdisk_l_temp = maskdisk_l;

%Add -1 shell to disk masks
maskdisk_l = int8(maskdisk_l) + (int8(maskdisk_l)-int8(maskdisk_l_out));
maskdisk_s = int8(maskdisk_s) + (int8(maskdisk_s)-int8(maskdisk_s_out));

%Library of self-convolved disk masks:
selfmask.disk_s = conv2(maskdisk_s,maskdisk_s_temp);
size_maskdisk_s = size(selfmask.disk_s,1);
center_s_disk = (size_maskdisk_s+1)/2;
selfmask.disk_l = conv2(maskdisk_l,maskdisk_l_temp);
size_maskdisk_l = size(selfmask.disk_l,1);
center_l_disk = (size_maskdisk_l+1)/2;
%Note: The second term is the _temp mask because we want to subtract these
%self-convolved masks from the convolution image...we are essentially
%detecting peaks with the outer -1 shell but eliminated found particles
%without this shell. This prevents overlap in original image but means that
%calculations are approximate rather than exact.

if numsides~=0
    %Polygon masks
    
    counter = 1;
    maskpoly_s = NaN([size(maskpoly_s_temp)+[size_layer_s,size_layer_s]*2, length(0:angular_res:maxang)]);
    clear maskpoly_s_temp
    for theta = 0:angular_res:maxang
        [maskpoly_s_ang,~,outerrad_s] =  polygonmask(numsides,size_small,theta);
        maskpoly_s_ang = padarray(maskpoly_s_ang,[size_layer_s,size_layer_s],0,'both');
        tempmask = bwmorph(maskpoly_s_ang,'thicken',size_layer_s);
        tempmask = bwmorph(tempmask,'diag');
        tempmask2 = logical(tempmask-maskpoly_s_ang);
        maskpoly_s_ang = int8(maskpoly_s_ang);
        maskpoly_s_ang(tempmask2) = -1;
        maskpoly_s(:,:,counter) = maskpoly_s_ang;
        counter = counter+1;
    end
    clear maskpoly_s_ang
    
    counter = 1;
    maskpoly_l = NaN([size(maskpoly_l_temp)+[size_layer_l,size_layer_l]*2, length(0:angular_res:maxang)]);
    clear maskpoly_l_temp
    for theta = 0:angular_res:maxang
        [maskpoly_l_ang,~,outerrad_l] =  polygonmask(numsides,size_large,theta);
        maskpoly_l_ang = padarray(maskpoly_l_ang,[size_layer_l,size_layer_l],0,'both');
        tempmask = bwmorph(maskpoly_l_ang,'thicken',size_layer_l);
        tempmask = bwmorph(tempmask,'diag');
        tempmask2 = logical(tempmask-maskpoly_l_ang);
        maskpoly_l_ang = int8(maskpoly_l_ang);
        maskpoly_l_ang(tempmask2) = -1;
        maskpoly_l(:,:,counter) = maskpoly_l_ang;
        counter = counter+1;
    end
    clear maskpoly_l_ang tempmask counter
    
    maskpoly_s_temp = maskpoly_s;
    maskpoly_s_temp(maskpoly_s_temp<0) = 0;
    maskpoly_l_temp = maskpoly_l;
    maskpoly_l_temp(maskpoly_l_temp<0) = 0;
    selfmask.poly_s = cell(size(maskpoly_s_temp,3),1);
    for layer1 = 1:size(maskpoly_s_temp,3)
        for layer2 = 1:size(maskpoly_s_temp,3)
            selfmask.poly_s{layer1}(:,:,layer2) = single(conv2(rot90(maskpoly_s_temp(:,:,layer1),2),maskpoly_s(:,:,layer2)));
        end
    end
    selfmask.poly_l = cell(size(maskpoly_l_temp,3),1);
    for layer1 = 1:size(maskpoly_l_temp,3)
        for layer2 = 1:size(maskpoly_l_temp,3)
            selfmask.poly_l{layer1}(:,:,layer2) = single(conv2(rot90(maskpoly_l_temp(:,:,layer1),2),maskpoly_l(:,:,layer2)));
            % Why rot90? Because we're doing convolutions, which flip over
            % both axes! So when we want to subtract mask from original
            % image, what we actually want to subtract is the flipped image
            % (from convolution).
        end
    end
   
    size_maskpoly_s = size(selfmask.poly_s{1},1);
    size_maskpoly_l = size(selfmask.poly_l{1},1);
    center_s_poly = (size_maskpoly_s+1)/2;
    center_l_poly = (size_maskpoly_l+1)/2;
    
end

%Save mat files
if numsides ~= 0
    save([savemask_dir shape 'Masks.mat'],...
        'center_l_disk','center_s_disk','center_l_poly','center_s_poly','maskdisk_l',...
        'maskdisk_s','maskpoly_s','maskpoly_l','size_maskdisk_l','size_maskdisk_s',...
        'size_maskpoly_s','size_maskpoly_l','selfmask','size_small','size_large',...
        'size_layer_s','size_layer_l','numsides','angular_res','maxang','innerrad_s','innerrad_l')
else
    save([savemask_dir shape 'Masks.mat'],...
        'center_l_disk','center_s_disk','maskdisk_l','maskdisk_s','size_maskdisk_l',...
        'size_maskdisk_s','selfmask','size_small','size_large','size_layer_s','size_layer_l','numsides',...
        'innerrad_s','innerrad_l')
end
fprintf(['Masks for ' shape 's saved \n'])


%%
function [maskpoly, innerrad, distvorttocent] = polygonmask(numsides,sidelength,angrot)
% Creates mask of side length floor(distanceFromVortexToCenter*2)
% with polygon COM at the center and rotated at angle angrot (degrees)
% from 0 degrees, which is set as a flat face at y=0. Masks are periodic in
% angle from 0->360/numsides-1 degrees.
% innerrad = is optional output for inner circle radius
% NOTE: numsides = 0 is a circle mask with radius sidelength

% if numsides~=0
%     sidelength = sidelength-0.5; %adjusts for poly2mask
% end

if numsides==3
    tri_0degree = sidelength * [[0 1 cosd(60)]',[0 0 sind(60)]'];
    distvorttocent = tri_0degree(3,2) - sidelength/2*tand(30);
    %fit triangle into space of largest sweeping area of triangle by
    %considering triangle at 30 and 60 rotation, where vertices are in minimal
    %points for x and y range
    tri_0degree(:,1) = tri_0degree(:,1) + distvorttocent*(1 - cosd(30));
    tri_0degree(:,2) = tri_0degree(:,2) + distvorttocent*(1 - cosd(60));
    %find center
    tri_cent = [ tri_0degree(1,1)+sidelength/2 tri_0degree(1,2)+distvorttocent*sind(30) ];
    %unrotated, centered triangle
    tri_0degree = [ tri_0degree(:,1)-tri_cent(1) tri_0degree(:,2)-tri_cent(2)];
    %rotate triangle
    tri_rot = [cosd(angrot) -sind(angrot); sind(angrot) cosd(angrot)] * tri_0degree';
    tri_rot = [tri_rot(1,:)'+tri_cent(1) tri_rot(2,:)'+tri_cent(2)];
    %create mask
    maskpoly = poly2mask(tri_rot(:,1),tri_rot(:,2),floor(2*distvorttocent),floor(2*distvorttocent));
    
elseif numsides == 4
    sqr_0degree = sidelength* [[0 1 1 0]',[0 0 1 1]' ];
    distvorttocent = sidelength/sqrt(2);
    %fit square into space of largest sweeping area of square by
    %considering square at 45 rotation, where vertices are in minimal
    %points for x and y range
    sqr_0degree(:,1) = sqr_0degree(:,1) + distvorttocent*(1 - cosd(45));
    sqr_0degree(:,2) = sqr_0degree(:,2) + distvorttocent*(1 - cosd(45));
    %find center
    sqr_cent = [ sqr_0degree(1,1)+sidelength/2 sqr_0degree(1,2)+distvorttocent*sind(45) ];
    %unrotated, centered square
    sqr_0degree = [ sqr_0degree(:,1)-sqr_cent(1) sqr_0degree(:,2)-sqr_cent(2)];
    %rotate square
    sqr_rot = [cosd(angrot) -sind(angrot); sind(angrot) cosd(angrot)] * sqr_0degree';
    sqr_rot = [sqr_rot(1,:)'+sqr_cent(1) sqr_rot(2,:)'+sqr_cent(2)];
    %create mask
    maskpoly = poly2mask(sqr_rot(:,1),sqr_rot(:,2),floor(2*distvorttocent),floor(2*distvorttocent));
    
elseif numsides == 5
    %Create pentagonal mask of side length pentsidelength and angle angrot
    pent_0degree = sidelength* [[cosd(72) cosd(72)+1 cosd(72)+1+cosd(72) cosd(72)+1/2 0]' ,...
        [0 0 sind(72) sind(72)+cosd(54) sind(72)]'] ;
    distvorttocent = pent_0degree(4,2)-sidelength/2*tand(54);
    %fit pentagon into space of largest sweeping area of pentagon
    pent_0degree(:,1) = pent_0degree(:,1) + (distvorttocent - sidelength*(1/2+cosd(72)));
    pent_0degree(:,2) = pent_0degree(:,2) + (distvorttocent - sidelength/2*tand(54));
    %find center
    pent_cent = [sidelength/2+sidelength*cosd(72) + pent_0degree(5,1) , sidelength/2*sind(54)/cosd(54) + pent_0degree(1,2)];
    %unrotated, centered pentagon
    pent_0degree = [pent_0degree(:,1)-pent_cent(1), pent_0degree(:,2)-pent_cent(2)];
    %rotate pentagon
    pent_rot = [cosd(angrot) -sind(angrot); sind(angrot) cosd(angrot)] * pent_0degree';
    pent_rot = [pent_rot(1,:)' + pent_cent(1), pent_rot(2,:)' + pent_cent(2)] ;
    %createmask
    maskpoly = poly2mask(pent_rot(:,1),pent_rot(:,2),floor(2*distvorttocent),floor(2*distvorttocent));

elseif numsides==6
    hex_0degree = sidelength * [[cosd(60) 1+cosd(60) 1+2*cosd(60) 1+cosd(60) cosd(60) 0]',...
        [0 0 sind(60) 2*sind(60) 2*sind(60) sind(60)]'];
    distvorttocent = sidelength/2/cosd(60);
    %fit hexagon into space of largest sweeping area of hexagon by
    %considering hexagon at 30 rotation, where vertices are in minimal
    %points for y range
    hex_0degree(:,1) = hex_0degree(:,1) + 0;
    hex_0degree(:,2) = hex_0degree(:,2) + distvorttocent*(1 - cosd(30));
    %find center
    hex_cent = [ hex_0degree(6,1)+distvorttocent hex_0degree(6,2) ];
    %unrotated, centered hexagon
    hex_0degree = [ hex_0degree(:,1)-hex_cent(1) hex_0degree(:,2)-hex_cent(2)];
    %rotate hexagon
    hex_rot = [cosd(angrot) -sind(angrot); sind(angrot) cosd(angrot)] * hex_0degree';
    hex_rot = [hex_rot(1,:)'+hex_cent(1) hex_rot(2,:)'+hex_cent(2)];
    %create mask
    maskpoly = poly2mask(hex_rot(:,1),hex_rot(:,2),floor(2*distvorttocent),floor(2*distvorttocent));
    
elseif numsides==7
    hepang = (180-360/7)/2;
    distvorttocent = sidelength/2/cosd(hepang);
    v1 = sidelength * [cosd(360/7), 0] ;
    center = v1 + [sidelength/2, distvorttocent*sind(hepang)];
    v2 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v1 - center]']' + center;
    v3 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v2 - center]']' + center;
    v4 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v3 - center]']' + center;
    v5 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v4 - center]']' + center;
    v6 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v5 - center]']' + center;
    v7 = [[cosd(360/7) -sind(360/7); sind(360/7) cosd(360/7)] * [v6 - center]']' + center;
    
    hep_0degree = [v1;v2;v3;v4;v5;v6;v7];
    %fit heptagon into space of largest sweeping area of heptagon by
    %considering heptagon at 30 rotation, where vertices are in minimal
    %points for y range
    hep_0degree(:,1) = hep_0degree(:,1) + distvorttocent*(1-cosd(360/7/4));
    hep_0degree(:,2) = hep_0degree(:,2) + distvorttocent*(1 - cosd(360/7/2));
    %find center
    hep_cent = [ hep_0degree(7,1)+distvorttocent hep_0degree(1,2)+distvorttocent*sind(hepang) ];
    %unrotated, centered heptagon
    hep_0degree = [ hep_0degree(:,1)-hep_cent(1) hep_0degree(:,2)-hep_cent(2)];
    %rotate heptagon
    hep_rot = [cosd(angrot) -sind(angrot); sind(angrot) cosd(angrot)] * hep_0degree';
    hep_rot = [hep_rot(1,:)'+hep_cent(1) hep_rot(2,:)'+hep_cent(2)];
    %create mask
    maskpoly = poly2mask(hep_rot(:,1),hep_rot(:,2),floor(2*distvorttocent),floor(2*distvorttocent));
elseif numsides == 0
    maskpoly = zeros([sidelength,sidelength]);
    for xx = 1:size(maskpoly,2)
        for yy = 1:size(maskpoly,1)
            if (xx-0.5-size(maskpoly,2)/2)^2 + (yy-0.5-size(maskpoly,1)/2)^2 <= (sidelength/2)^2
                maskpoly(yy,xx) = 1;
            end
        end
    end
    %Using increments of 0.5 for xx and yy ensure that pixel space corresponds
    %to real space. This centers disk within each mask array, whether even or
    %odd size is used.
    innerrad = sidelength/2;
end

maskpoly = bwmorph(maskpoly,'diag'); %spurred corners are filled.
if numsides~=0
    innerrad = sqrt(distvorttocent^2 - sidelength^2/4);
end

end
