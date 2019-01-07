function psf = calcPsfMFM(psfSz, voxSz, lambdas, varargin)
% 2017.04.27_Xiang & Kuan modified
% psf_c: conventional psf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% after this run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MicroscopeMFM3DPSF.m on the
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% left by changing sampling_space_x = 0.05;
%% MFM DOE load
% a1 + a2 = circlular support
%% Setup Parameters
%% Grating Parameters
load('Diffractive-optical-element-data\GratingPhaseChangePart_25.mat'); %[0, -1], -1 means shift
load('Diffractive-optical-element-data\GratingPhaseUnchangePart_25.mat'); % [0, +1], +1 means no phase change,  
grating_pixel_size  = 1300; % nm

%% Setup Parameters
%% Parse positional Parameters
if ~exist('psfSz', 'var') || isempty(psfSz)    
    psfSz = [1024, 1024, 121];    
end
if ~exist('lambdas', 'var') || isempty(lambdas)    
    lambdas = 520;    
    %lambdas = linspace(515*(10^(-3)),525*(10^(-3)),11);
end
NLambda = length(lambdas);
%% Parse name-value parameters
p = inputParser;  % CaseSensitive=0, PartialMatching=1
p.addParameter('verbose',   true); % true/false, whether show result
p.addParameter('psfCSz',   [101, 101, 121]); % 1D numeric vector
p.addParameter('wLambdas',  ones(1,NLambda),    @(x) isempty(x) || (isnumeric(x) && isvector(x))); % 1D numeric vector
% Microcope Parameters
p.addParameter('M',         60); % maginification size; 20x/0.5 or 40x/0.95 
p.addParameter('NA',        1.27); % 
p.addParameter('f_tl',      200*1e6); % um, tube lens 200mm 
p.parse(varargin{:});
cellfun(@(x) evalin('caller', sprintf('%s=p.Results.%s;', x, x)), fields(p.Results)); % e.g. eval('fSz = p.Results.fSz'); Or use in.fSz after set in =p.Results;
% compute dependent parameters
if  ~exist('voxSz', 'var') || isempty(voxSz)
    voxSz(1:2) = 13*1e3/(400/200)/M; %camera_pixel_size/(f2/f1*M), camera_pixel_size=13um, f2=400mm, f1=200mm; Pixel size after Microscope ,why much less than  diff_limit_x=central_lambda/2/NA= 0.205um    
    voxSz(3) = 50; % 50nm, much less thant diff_limit_z = 2*central_lambda/(NA)^2 = 0.645um;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psf = zeros(psfSz);
psfLambdaZ = zeros([psfSz(1:2), NLambda]);
%% 3D PSF generation

[~, ~, csfLambda] = calcPsfMicroscope(psfCSz, voxSz, lambdas, 'verbose', verbose, 'wLambdas', wLambdas, 'M', M, 'NA', NA, 'f_tl', f_tl);

psfC = floor(psfSz(1:2)/2); %[512; 512]
for nZ = 1:psfSz(3)    
    nZ
    for nLambda = 1:NLambda
        nLambda
        lambda = lambdas(nLambda); % um
        nn = ceil(400*1e6*lambda/grating_pixel_size/(13*1e3)); % nn = ceil(f2*lambda/grating_pixel_size/camera_pixel_size);
        a2_1 = zeroPadding(a2, nn);
        a1_1 = zeroPadding(a1, nn);
        %x_out_resolution = lambda/central_lambda*sampling_space_x;
        grating = a2_1 + a1_1.*(-exp(1i*pi*lambda/520)); % grating = a2_1 + a1_1.*(-exp(1i*pi*lambda/central_lambda));
        center_grating = ceil(size(grating,1)/2);
       %% generate the field in native image plane u_i
        tic;        
        csfLambda_zeroPad = zeroPadding(csfLambda(:, :, nZ, nLambda), nn); 
        csfCam = fft2(fftshift(fft2(csfLambda_zeroPad)).*grating);  % psf on camera space
        csfCamCrop = csfCam(center_grating-(psfC(1)-1):center_grating+psfC(1), center_grating-(psfC(2)-1):center_grating+psfC(2));
        psfLambdaZ(:,:,nLambda) = csfCamCrop.*conj(csfCamCrop);
        toc;
    end
    psf(:,:,nZ) = sum(psfLambdaZ,3);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psf = psf./(max(psf(:)));
if verbose
    figure;show3dPointsWithProjections(psf, voxSz, 0.001);title('MFM 3D PSF')
end
OTF = (fftn(psf));
OTF = fftshift(abs(OTF));
OTF = OTF./max(OTF(:));
if verbose
    figure;show3dPointsWithProjections(OTF, 1./(voxSz.*psfSz), 0.01); title('MFM 3D OTF')
end


end

