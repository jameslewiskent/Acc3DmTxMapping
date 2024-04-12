function [Maps,Rel_Images,cumtime,Rel_Maps]  = ReconB1TIAMOMaps(MIDs,Desired_ro,Desired_pe1,Desired_pe2,settings)
% Take two MIDs, the first one is the 8 relative GRE maps, the second is
% the 2 absolute shim setting modes
% This method uses TxLR to jointly reconstruct the 8 relative and 2
% absolute maps

if isfield(settings,'niters')
    niters = settings.niters;
else
    niters = 50;
end

if  isfield(settings,'RetroMask')
    RetroMask = settings.RetroMask;
    RetroUS = 1;
else
    RetroMask = 'NaN';
    RetroUS = 0;
end

if isfield(settings,'kshift')
    kshift = settings.kshift;
else
    kshift = [0 0 0];
end

if isfield(settings,'RedoRecon')
    RedoRecon = settings.RedoRecon;
else
    RedoRecon = 0;
end

if isfield(settings,'m')
    m = settings.m;
else
    m = 3;
end

if isfield(settings,'NCC')
    NCC = settings.NCC;
else
    NCC = 0;
end

if isfield(settings,'CCType') && (settings.CCType == 1 || strcmp(settings.CCType,'Geometric'))
    CCType = 'Geometric';
elseif isfield(settings,'CCType') && (settings.CCType == 2 ||strcmp(settings.CCType,'ESPIRiT'))
    CCType = 'ESPIRiT';
elseif isfield(settings,'CCType') && (settings.CCType == 3 ||strcmp(settings.CCType,'SingularValue'))
    CCType = 'SingularValue';
else
    CCType = 'Geometric';
end

if isfield(settings,'plotyn')
    plotyn = settings.plotyn;
else
    plotyn = 0;
end

if isfield(settings,'calcsens')
    calcsens = settings.calcsens;
else
    calcsens = 1;
end

if isfield(settings,'Shims')
    Shim_Setting1 = settings.Shims(1,:);
    Shim_Setting2 = settings.Shims(2,:);
else
    % Define encoding matrix used to acquire both complimentary shims
    [Enc_Matrix] = Gen_Enc_Mat(16,8);
    Shim_Setting1 = conj(Enc_Matrix(3,:)); % Assuming CP+
    Shim_Setting2 = conj(Enc_Matrix(5,:)); % Assuming CP2+
end

rxkernel = [5 5];
txkernel = [5 5];
txlrkernel = [5 5];
eig_thresh = 0.02;
nnWeights = [50 50];
calibsize = [24 24];
ccdim = 1;
tukey_val = 0.7;
se = strel('sphere',1);

for MID_n = 1:size(MIDs,2)
    MID = MIDs(MID_n);
    twix = mapVBVD(MID);
    twix.image.flagIgnoreSeg = true;
    twix.image.flagRemoveOS = true;
    opts    =   struct();
    if MID_n == 1 && isfield(twix,'noise')
        twix.noise.flagRemoveOS = true;
        noise = twix.noise('');
        opts.noise = cov(noise);
    end
    NRX = twix.image.NCha;
    NTX = 8;
    volts(MID_n) = twix.hdr.Spice.TransmitterReferenceAmplitude;
    scan_params{MID_n} = ['TransmitterReferenceAmplitude',twix.hdr.Meas.TransmitterReferenceAmplitude,'FlipAngle',twix.hdr.Meas.adFlipAngleDegree];
    
    NCol = twix.image.NCol/2; % /2 removes over-sampling in read direction
    kdata{MID_n} = squeeze(twix.image(''));
    Acceleration_Factor{MID_n} = round(numel(kdata{MID_n})/nnz(kdata{MID_n}));% Determine acceleration factor used
end

if twix.hdr.Config.NPar == 1 % 2D
    kdata{1} = permute(kdata{1},[5,2,1,3,4]);
    kdata{2} = permute(kdata{2},[6,2,1,3,4,5]);
end

if size(MIDs,2) == 3
    % Absolute shim maps in seperate twix files so combine.
    kdata{2} = cat(5,permute(kdata{2},[1 2 3 4 6 5]),permute(kdata{3},[1 2 3 4 6 5])); kdata{3} = [];
    if volts(2) ~= volts(3)
        error('different voltages!')
    end
end

if any(size(kdata{1},3:4) ~= size(kdata{2},3:4)) % Sizes don't match, crop down
    crop_size = min(size(kdata{1},3:4),size(kdata{2},3:4));
    
    % Crop kdata1
    centre = size(kdata{1},3:4)/2;
    crop1 = centre(1)-floor((crop_size(1)-1)/2):centre(1)+ceil((crop_size(1)-1)/2);
    crop2 = centre(2)-floor((crop_size(2)-1)/2):centre(2)+ceil((crop_size(2)-1)/2);
    kdata{1} = kdata{1}(:,:,crop1,crop2,:);
    
    % Crop kdata2
    centre = size(kdata{2},3:4)/2;
    crop1 = centre(1)-floor((crop_size(1)-1)/2):centre(1)+ceil((crop_size(1)-1)/2);
    crop2 = centre(2)-floor((crop_size(2)-1)/2):centre(2)+ceil((crop_size(2)-1)/2);
    kdata{2} = kdata{2}(:,:,crop1,crop2,:,:);
end

if NCC == 0
    NCC = NRX;
end
txlr_range = 1:(size(kdata{1},5)+size(kdata{2},5)*size(kdata{2},6));

recon_parameters = {'eig_thresh','Desired_pe1','Desired_pe2','kshift'};
upsample_parameters= {'MID','niters','nnWeights','rxkernel','txkernel','txlrkernel','NRX','txlr_range'};
filename = ['ReconData_MIDs',num2str(MIDs(1)),'and',num2str(MIDs(2)),'niters',num2str(niters),'vc',num2str(NCC),CCType,'.mat']; % fully-sampled k-space reconstructed from under-sampled k-space

if RetroUS == 1
    filename = strcat('RETRO_Acc',num2str(round(numel(RetroMask)/nnz(RetroMask))),'_',filename); % This is a retrospective undersampling
end
tic
cumtime(1) = 0;
if ~isfile(filename) || RetroUS == 1 % only if results not already saved or not fully sampled
    
    kdata_joint = cat(5,kdata{1},kdata{2}(:,:,:,:,:,1),kdata{2}(:,:,:,:,:,2)); % Currently [Kx Rx Ky Kz Tx]
    
    if RetroUS == 1
        kdata_joint = kdata_joint .*repmat(RetroMask,[size(kdata_joint,1) size(kdata_joint,2) 1 1 1]);
    end
    if plotyn == 1
        figure(); imagesc(imtile(abs(squeeze(kdata_joint(12,1,:,:,:))))); title('Undersampled kspace')
    end
    USMask = zeros(size(kdata_joint)); USMask(kdata_joint~=0)=1;
    if NCC ~= 0 && NRX > 8 && NCC ~= NRX
        % Perform Coil Compression
        clear CCDATA mtx
        
        DATA = permute(kdata_joint,[1 3 4 2 5]);
        USMask = zeros(size(DATA)); USMask(DATA~=0)=1;
        calib = sum(DATA(:,:,:,:,:),5); % Find a single calibration matrix for all transmit
        if strcmp(CCType,'Geometric')
            tempmtx =  calcGCCMtx(calib,ccdim,1); % calcGCCMtx or calcECCMtx or calcSCCMtx
            tempmtx = tempmtx(:,1:NCC,:,:);
            mtx = alignCCMtx(tempmtx,NCC);
            for Tx_n = txlr_range
                CCDATA(:,:,:,:,Tx_n) = CC(DATA(:,:,:,:,Tx_n),mtx,ccdim);
            end
        elseif strcmp(CCType,'ESPIRiT')
            tempmtx =  calcECCMtx(calib,ccdim,NCC); % calcGCCMtx or calcECCMtx or calcSCCMtx
            tempmtx = tempmtx(:,1:NCC,:,:);
            mtx = alignCCMtx(tempmtx,NCC);
            for Tx_n = txlr_range
                CCDATA(:,:,:,:,Tx_n) = CC(DATA(:,:,:,:,Tx_n),mtx,ccdim);
            end
        elseif strcmp(CCType,'SingularValue')
            tempmtx =  calcSCCMtx(calib); % calcGCCMtx or calcECCMtx or calcSCCMtx
            mtx = tempmtx(:,1:NCC,:,:);
            for Tx_n = txlr_range
                CCDATA(:,:,:,:,Tx_n) = CC(DATA(:,:,:,:,Tx_n),mtx);
            end
        end
        kdata_joint= permute(USMask(:,:,:,1:NCC,:).*CCDATA,[1 4 2 3 5]);
    end
    
    % Perform FFT in RO dimension prior to recon
    kdata_joint = ifftc(kdata_joint,1);
    
    % pre-allocate matrices
    kdata_full = kdata_joint;
    if Acceleration_Factor{1} ~= 1 || RetroUS == 1
        parfor layer_n = 1:size(kdata_joint,1)
            kdata_joint_layer = permute(squeeze(kdata_joint(layer_n,:,:,:,txlr_range)),[2 3 1 4]); % Concatenate transmit images and get into right dimensions ADMM wants under-sampled [kx,ky,Rx,Tx] k-space
            kdata_full(layer_n,:,:,:,txlr_range) = permute(admm_txlr(double(kdata_joint_layer), txlrkernel, niters, nnWeights,opts),[3 1 2 4]);
            disp(['Kx-Ky Layer ',num2str(layer_n),'/',num2str(NCol),' k-space fully up-sampled'])
        end
    else
        disp('Data is already fully sampled.')
    end
    
    % Undo FFT in RO
    kdata_full = fftc(kdata_full,1);
    
    cumtime(1) = toc;
    save(filename,'kdata_full','USMask',upsample_parameters{:})
    cumtime(2) = toc;
elseif isfile(filename)
    load(filename,'kdata_full',upsample_parameters{:})
    disp('Loaded in fully sampled dataset.')
    cumtime = 'NaN';
end

if ismember('Maps', who('-file', filename))
    load_check = load(filename,'eig_thresh','Desired_pe1','Desired_pe2','kshift');
    try
        recon_param_unchanged_bool = (load_check.eig_thresh == eig_thresh && load_check.Desired_pe2 == Desired_pe2 && load_check.Desired_pe1 == Desired_pe1 && all(load_check.kshift == kshift));
    catch
        disp('Could not load a parameter.');
        recon_param_unchanged_bool = 0;
    end
    if recon_param_unchanged_bool && RedoRecon == 0
        disp('Loading in previously reconstructed images.')
        load(filename,'Rel_Images','Ref_Shim1','Ref_Shim2','Prep_Shim1','Prep_Shim2','tx_sens')
    end
end

% recon results not saved, or recon values are different
if ~ismember('Maps', who('-file', filename)) || ~recon_param_unchanged_bool || RedoRecon
    disp('Reconstructing Images.')
    permute_set = cat(1,[0,2,4,6,3,1,7,5]'+1,(9:size(kdata_full,5))'); % Since excitation order is 0 5 1 4 2 7 3 6 for relative maps
    kdata_full = kdata_full(:,:,:,:,permute_set);
    
    % Phase roll data if required
    if any(kshift ~=0)
        kdata_full = kdata_full.*repmat(permute(exp(1i*linspace(0,kshift(1)*pi,size(kdata_full,1))),[2 1]),1,size(kdata_full,2),size(kdata_full,3),size(kdata_full,4),size(kdata_full,5));
        kdata_full = kdata_full.*repmat(permute(exp(1i*linspace(0,kshift(2)*pi,size(kdata_full,3))),[3 1 2]),size(kdata_full,1),size(kdata_full,2),1,size(kdata_full,4),size(kdata_full,5));
        kdata_full = kdata_full.*repmat(permute(exp(1i*linspace(0,kshift(3)*pi,size(kdata_full,4))),[4 3 1 2]),size(kdata_full,1),size(kdata_full,2),size(kdata_full,3),1,size(kdata_full,5));
    end
    if plotyn == 1
        figure(); imagesc(imtile(squeeze(sum(abs(kdata_full(12,:,:,:,:)),2)))); title('Reconstructed kspace (Summed Abs Rx)')
    end
    
    centre = size(kdata_full,3:4)/2;
    calib1 = centre(1)-floor((calibsize(1)-1)/2):centre(1)+ceil((calibsize(1)-1)/2);
    calib2 = centre(2)-floor((calibsize(2)-1)/2):centre(2)+ceil((calibsize(2)-1)/2);
    if any(calibsize < rxkernel)
        rxkernel = [calibsize(1) calibsize(2)];
        txkernel = [calibsize(1) calibsize(2)];
    end
    
    % Define windowing filters
    TukeyFilter = tukeywin(calibsize(1)+2,tukey_val); TukeyFilter = TukeyFilter(2:end-1,1);
    
    if Desired_ro ~= 1
    % Filter in RO
    kdata_full = kdata_full.*repmat(TukeyFilter,[1,size(kdata_full,2:5)]);
    end
    
    % Zero-pad in RO
    kdata_zp = padarray(padarray(kdata_full,ceil(([Desired_ro,0,0,0,0]-[size(kdata_full,1),0,0,0,0])/2),'pre'),floor(([Desired_ro,0,0,0,0]-[size(kdata_full,1),0,0,0,0])/2),'post');
    
    if plotyn == 1
        figure(); imagesc(imtile(abs(squeeze(kdata_zp(round(size(kdata_zp,1)/2),1,:,:,:))))); title('Zero padded kspace')
    end
    
    %   Perform fft in col dimension prior to rx_sens
    kdata_hybrid = ifftc(kdata_zp,1); clear kdata_zp
    % Calculate Tx and Rx sensitivities already been FFT'd in RO dimension
    rx_sens = zeros(size(kdata_hybrid,1),Desired_pe1,Desired_pe2,NCC);
    tx_sens = zeros(size(kdata_hybrid,1),Desired_pe1,Desired_pe2,NTX);
    parfor col_n = 1:size(kdata_hybrid,1) % Calculated sensitivites based solely of relative maps
        rx_sens(col_n,:,:,:) = tx_espirit(permute(squeeze(kdata_hybrid(col_n,:,calib1,calib2,1:8)),[2,3,4,1]), [Desired_pe1 Desired_pe2], rxkernel, eig_thresh);
        tx_sens(col_n,:,:,:) = tx_espirit(permute(squeeze(kdata_hybrid(col_n,:,calib1,calib2,1:8)),[2,3,1,4]), [Desired_pe1 Desired_pe2], txkernel, eig_thresh);
    end
    
    % Define windowing filters
    TukeyRot = RotArray(tukeywin(calibsize(1)+2,tukey_val)); TukeyRot = TukeyRot(2:end-1,2:end-1);
    
    % Filter in both phase dimensions
    kdata_hybrid = kdata_hybrid.*repmat(permute(TukeyRot,[3 4 1 2]),[size(kdata_hybrid,1),size(kdata_hybrid,2),1,1,size(kdata_hybrid,5)]);
    
    % Zero-pad in both phase dimensions
    kdata_zp = padarray(padarray(kdata_hybrid,ceil(([0,0,Desired_pe1,Desired_pe2,0]-[0,0,size(kdata_hybrid,3),size(kdata_hybrid,4),0])/2),'pre'),floor(([0,0,Desired_pe1,Desired_pe2,0]-[0,0,size(kdata_hybrid,3),size(kdata_hybrid,4),0])/2),'post'); clear kdata_hybrid
    
    % fft in lin and par dimension:
    fft_dims = [3 4];
    for f = fft_dims
        kdata_zp = ifftc(kdata_zp,f);
    end
    
    % Combine sensitivity maps with images, sum across receive channels
    images  = (sum(kdata_zp.*conj(permute(rx_sens,[1 4 2 3])),2));
    
    % Put singleton dimension at end of array
    images = permute(images,[1,3,4,5,2]);
    
    if plotyn == 1
        figure('color','w','Position',[595,195,575,417]); imagesc(imtile(squeeze(abs(images(size(images,1)/2,:,:,9:12))))); axis image off
        figure('color','w','Position',[485,252,555,431]); imagesc(imtile(squeeze(abs((180./pi).*acos(real(images(size(images,1)/2,:,:,11:12)./images(size(images,1)/2,:,:,9:10)))))),[0 180]); axis image off
        colormap(turbo); cb = colorbar; cb.Label.String = ['\alpha [',char(176),']'];
        
        figure('color','w'); imagesc(imtile(abs(squeeze(images(size(images,1)/2,:,:,:))))); title('Images'); axis image off
        figure('color','w'); imagesc(imtile(abs(squeeze(images(:,size(images,2)/2,:,:))))); title('Images'); axis image off
        figure('color','w'); imagesc(imtile(abs(squeeze(images(:,:,size(images,3)/2,:))))); title('Images'); axis image off
    end
    
    Rel_Images = images(:,:,:,1:8);
    cumtime(3) = toc;
    Ref_Shim1 = images(:,:,:,9);
    Ref_Shim2 = images(:,:,:,10);
    Prep_Shim1 = images(:,:,:,11);
    Prep_Shim2 = images(:,:,:,12);
    save(filename,'Rel_Images','Ref_Shim1','Ref_Shim2','Prep_Shim1','Prep_Shim2','tx_sens',recon_parameters{:},'-append')
    cumtime(4) = toc;
end

if plotyn == 1
    figure('color','w','Position',[556,13,644,631])
    tiledlayout('flow','TileSpacing','compact','Padding','compact')
    nexttile; imagesc(abs(squeeze(Ref_Shim1(size(Ref_Shim1,1)/2,:,:))),[0 max(abs(Ref_Shim1),[],'all')]); title('Ref Shim1'); axis image off;
    nexttile; imagesc(abs(squeeze(Ref_Shim2(size(Ref_Shim2,1)/2,:,:))),[0 max(abs(Ref_Shim1),[],'all')]); title('Ref Shim2'); axis image off;
    nexttile; imagesc(abs(squeeze(Prep_Shim1(size(Prep_Shim1,1)/2,:,:))),[0 max(abs(Ref_Shim1),[],'all')]); title('Prep Shim1'); axis image off;
    nexttile; imagesc(abs(squeeze(Prep_Shim2(size(Prep_Shim2,1)/2,:,:))),[0 max(abs(Ref_Shim1),[],'all')]); title('Prep Shim2'); axis image off;
end

disp('Calculating Maps.')

% Calculate relative maps from transmit sensitivities or relative images
if calcsens == 1
    Rel_Maps =  tx_sens./(sum(abs(tx_sens),4));
else
    Rel_Maps = Rel_Images./(sum(abs(Rel_Images),4));
end

if plotyn == 1
    figure();imagesc(imtile(abs(squeeze(Rel_Maps(size(Rel_Maps,1)/2,:,:,:)))),[0 1]); title('Relative Maps [Magnitude]')
    figure();imagesc(imtile(angle(squeeze(Rel_Maps(12,:,:,:)).*exp(-1i.*angle(squeeze(Rel_Maps(12,:,:,1)))))),[-pi pi]); title('Relative Maps [Phase]')
end

% Load lookup table
load('sandwich_lookup_table.mat','x_query','fx_interp');

% Remove ratios outside of dynamic range Shim 1
Image_Ratio_Shim1 = real(Prep_Shim1 ./ Ref_Shim1);
Image_Ratio_Shim1(Image_Ratio_Shim1 > max(fx_interp)) = max(fx_interp);
Image_Ratio_Shim1(Image_Ratio_Shim1 < min(fx_interp)) = min(fx_interp);

% Mask for Saturated values Shim 1
Saturated_Mask_Shim1 = ones(size(Image_Ratio_Shim1));
Saturated_Mask_Shim1(Image_Ratio_Shim1 == single(max(fx_interp))) = 0;
Saturated_Mask_Shim1(Image_Ratio_Shim1 == single(min(fx_interp))) = 0;
Saturated_Mask_Shim1 = imerode(Saturated_Mask_Shim1,se);

% Remove ratios outside of dynamic range Shim 2
Image_Ratio_Shim2 = real(Prep_Shim2 ./ Ref_Shim2);
Image_Ratio_Shim2(Image_Ratio_Shim2 > max(fx_interp)) = max(fx_interp);
Image_Ratio_Shim2(Image_Ratio_Shim2 < min(fx_interp)) = min(fx_interp);

% Mask for Saturated values Shim 2
Saturated_Mask_Shim2 = ones(size(Image_Ratio_Shim2));
Saturated_Mask_Shim2(Image_Ratio_Shim2 == single(max(fx_interp))) = 0;
Saturated_Mask_Shim2(Image_Ratio_Shim2 == single(min(fx_interp))) = 0;
Saturated_Mask_Shim2 = imerode(Saturated_Mask_Shim2,se);

% Create a NaN mask for lookup-table saturated values
Nan_Saturated_Mask_Shim1 = nan(size(Saturated_Mask_Shim1));
Nan_Saturated_Mask_Shim1(Saturated_Mask_Shim1 == 1) = 1;
Nan_Saturated_Mask_Shim2 = nan(size(Saturated_Mask_Shim2));
Nan_Saturated_Mask_Shim2(Saturated_Mask_Shim2 == 1) = 1;

% Calculate absolute maps from lookup table
FA_Shim1 = zeros(size(Image_Ratio_Shim1));
FA_Shim2 = zeros(size(Image_Ratio_Shim2));
for x = 1:size(Image_Ratio_Shim1,1)
    for y = 1:size(Image_Ratio_Shim1,2)
        for z = 1:size(Image_Ratio_Shim1,3)
            [~,min_ind] = min(abs(Image_Ratio_Shim1(x,y,z) - fx_interp));
            FA_Shim1(x,y,z) = (180/pi)*x_query(min_ind); % Measured FA in degrees
            [~,min_ind] = min(abs(Image_Ratio_Shim2(x,y,z) - fx_interp));
            FA_Shim2(x,y,z) = (180/pi)*x_query(min_ind); % Measured FA in degrees
        end
    end
end

% Now calculate relative shims for each mode
Rel_Shim1 = sum(bsxfun(@times,Rel_Maps,permute(Shim_Setting1,[1 3 4 2])),4);
Rel_Shim2 = sum(bsxfun(@times,Rel_Maps,permute(Shim_Setting2,[1 3 4 2])),4);

if plotyn == 1
    figure('color','w','units','normalized','Position',[0.3,0.1,0.383463541666667,0.717852684144819]); tiledlayout('flow','TileSpacing','compact','Padding','compact');
    nexttile; imagesc(imtile(abs(squeeze(FA_Shim1(size(FA_Shim1,1)/2,:,:)))),[0 180]); axis image off; title('Flip Angle Map Shim1');
    nexttile; imagesc(imtile(abs(squeeze(FA_Shim2(size(FA_Shim1,1)/2,:,:)))),[0 180]); axis image off; title('Flip Angle Map Shim2');
    nexttile; imagesc(imtile(abs(squeeze(Rel_Shim1(size(FA_Shim1,1)/2,:,:))))); axis image off; title('Relative Map Shim1');
    nexttile; imagesc(imtile(squeeze(abs(Rel_Shim2(size(FA_Shim1,1)/2,:,:))))); axis image off; title('Relative Map Shim2');
    
    figure('color','w','units','normalized','Position',[0.3,0.1,0.383463541666667,0.717852684144819]); tiledlayout('flow','TileSpacing','compact','Padding','compact');
    nexttile; imagesc(abs(squeeze(FA_Shim1(:,size(FA_Shim1,2)/2,:))),[0 180]); axis image off; title('Flip Angle Map Shim1');
    nexttile; imagesc(abs(squeeze(FA_Shim2(:,size(FA_Shim1,2)/2,:))),[0 180]); axis image off; title('Flip Angle Map Shim2');
    nexttile; imagesc(abs(squeeze(Rel_Shim1(:,size(FA_Shim1,2)/2,:)))); axis image off; title('Relative Map Shim1');
    nexttile; imagesc(squeeze(abs(Rel_Shim2(:,size(FA_Shim1,2)/2,:)))); axis image off; title('Relative Map Shim2')
    
    figure('color','w','units','normalized','Position',[0.3,0.1,0.383463541666667,0.717852684144819]); tiledlayout('flow','TileSpacing','compact','Padding','compact');
    nexttile; imagesc(imtile(abs(squeeze(FA_Shim1(:,:,size(FA_Shim1,3)/2)))),[0 180]); axis image off; title('Flip Angle Map Shim1');
    nexttile; imagesc(imtile(abs(squeeze(FA_Shim2(:,:,size(FA_Shim1,3)/2)))),[0 180]); axis image off; title('Flip Angle Map Shim2');
    nexttile; imagesc(imtile(abs(squeeze(Rel_Shim1(:,:,size(FA_Shim1,3)/2))))); axis image off; title('Relative Map Shim1');
    nexttile; imagesc(imtile(squeeze(abs(Rel_Shim2(:,:,size(FA_Shim1,3)/2))))); axis image off; title('Relative Map Shim2');
end

% Divide FA shims by relative shims to get absolute scaling factors
c1 = abs(FA_Shim1./Rel_Shim1);
c2 = abs(FA_Shim2./Rel_Shim2);

% Weighted combination of scaling factors
C = ((c1.*abs(Ref_Shim1).^m) + (c2.*abs(Ref_Shim2).^m))./ (abs(Ref_Shim1).^m + abs(Ref_Shim2).^m);

% Scale relative maps to absolute maps
Maps = C.*Rel_Maps;

% Interpolate any lookup-table saturated 'holes' in maps
Interp_Maps = zeros(size(Maps));
for Tx_n = 1:8
    for layer_n = 1:size(Maps,1)
        Interp_Maps(layer_n,:,:,Tx_n) = inpaint_nans(squeeze(abs(double(Maps(layer_n,:,:,Tx_n))).*Nan_Saturated_Mask_Shim1(layer_n,:,:).*Nan_Saturated_Mask_Shim2(layer_n,:,:)),1);
    end
end
Maps = Interp_Maps.*exp(-1i.*angle(Rel_Maps)); % Rephase interpolated maps

cumtime(5) = toc;

Maps = Maps.*exp(-1i.*angle(Rel_Maps(:,:,:,1))); % Phase relative to first channel

if plotyn == 1
    figure('color','w'); tiledlayout('flow','TileSpacing','compact','Padding','compact');
    nexttile; imagesc(imtile(abs(squeeze(sum(Maps(size(Maps,1)/2,:,:,:),4)))),[0 160]); title('Summed to CP'); axis image off
    nexttile; imagesc(imtile(abs(squeeze(sum(Maps(:,size(Maps,2)/2,:,:),4)))),[0 160]); title('Summed to CP'); axis image off
    nexttile; imagesc(imtile(abs(squeeze(sum(Maps(:,:,22,:),4)))),[0 160]); title('Summed to CP'); axis image off
end

save(filename,'Maps','Rel_Maps','-append')
end

