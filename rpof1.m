%using 2D Riesz pahse optical flow
%based on:
% https://doi.org/10.1016/j.ymssp.2022.110044

%however, it works good especially if we apply temporal bandpass filtering
%for the extracted signal or for the phase difference
%however, this needs that the spatial center frequency be known since the
%theoretical base is that the intensity is a monocomponent which is not
%true in most of the cases and this put challenge on how to scale the
%displacement from the phase differences

close all
clear
% filename='D:/data/point1_new2_cropped.avi';
%   filename='D:/data/point2_2_cropped.avi';  
%   filename='D:/data/point3_cropped.avi';
filename='D:/data/forced2.mp4';
outDir = 'D:/data/Results/';
video_file = filename;

vr = VideoReader(video_file);
samplingRate = 500;% vr.FrameRate; %or set it manually
frameRange = [500 1600];
nF=frameRange(2)-frameRange(1)+1;
mmPerPixel = 0.124;
scale_factor  = 1;
%    pt = [293 630; %y ,x coordinates
%     292 624;
%     290 634];
% pt = [259 560; %y ,x coordinates
%     259 561;
%     259 562;
%     259 563;
%     259 564;
%     259 566;
%     259 567;
%     259 568;
%     259 569;
%     259 570;
%     259 602;
%     259 661];
%  pt=[179 504  %filename='D:/data/point2_2_cropped.avi';
%      179 505
%      179 506
%      179 507
%      179 508
%      179 509
%      179 510
%      179 520];
% pt=[218 555  %filename='D:/data/point3_cropped.avi';
%     218 556
%     218 557
%     218 558
%     218 559
%     218 560
%     218 561
%     218 562];
pt=[492 914  %'D:/data/forced2.mp4';
    492 915
    492 916
    492 917
    492 918
    492 919
    492 920];
readFrame = @(k) imresize(rgb2y(im2single(vr.read(frameRange(1)+k-1))), scale_factor);

loCutoff = 20;
hiCutoff = 26;
sigma = 3;

% Points to plot the motion in (y, x) format
[h, w, ~] = size(readFrame(1));



% Compute the maximal frequency response and modulation for each
% filter response.
spaceBlurExtent = 2*ceil(sigma);
x = -spaceBlurExtent:spaceBlurExtent;
blurKernel = exp(-(x.*x)/(2*sigma.^2+eps));
blurKernel = blurKernel./sum(blurKernel);
extra = (numel(blurKernel)-1)/2;
if extra < 2
    extra = 2;
end

%% point extraction

Np = size(pt,1);
phix = zeros(Np,1);
phiy = zeros(Np,1);
pIDx = zeros(Np,1);
pIDy = zeros(Np,1);
motion = zeros(Np,2, nF,'single');
tic


%inline is slower but much less memory is required

frameIDX = 24;
Icf = readFrame(frameIDX);
figure
imagesc(Icf);
Icf = Icf - mean(mean(Icf));
Ift = fft2(Icf);
[~,idx]=max(abs(Ift(:)));
[iy,ix]=ind2sub(size(Ift),idx);
wy=2*pi*iy/h;
wx=2*pi*ix/w;
fs=sqrt(wx^2+wy^2);
Ift = fftshift(Ift);
% figure
% mesh(abs(Ift))

d1=4;
d2=6;
Hf=zeros(h,w);
for u=1:w
    for v=1:h
        d=sqrt((u-w/2)^2+(v-h/2)^2);
        Hf(v,u)=1/(1+(d/d2)^10)-1/(1+(d/d1)^10);
    end
end
Ift = Ift .* Hf;
%Ift=ifftshift(Ift);
if mod(w,2) == 0
    w1 = -w/2:w/2-1;
else
    w1=-(w-1)/2:(w-1)/2;
end
if mod(h,2) == 0
    h1 = -h/2:h/2-1;
else
    h1=-(h-1)/2:(h-1)/2;
end
[u,v] = meshgrid(w1,h1);
%[u,v] = meshgrid(0:w-1,0:h-1);
ur = 1i*u ./ (0.00001+sqrt(u.^2+v.^2));
vr = 1i*v ./ (0.00001+sqrt(u.^2+v.^2));



%Ift=ifftshift(Ift);
%Ift = 0.000000001*Ift;
Irec = abs(ifft2(ifftshift(Ift)));
Rs1= abs(ifft2(ifftshift(Ift .* -ur)));%x-motion
Rs2= abs(ifft2(ifftshift(Ift .* -vr)));%y-motion
% figure
% imagesc(Rs1);
%         return
q0=zeros(Np,1);
q1=q0;
q2=q0;
r0=q0;
r1=q0;
r2=q0;

for frameIDX = 1:nF
    fprintf('Detecting motion for frame %d\n', frameIDX);
    Icf = readFrame(frameIDX); %this is luminance data
    Ift = fft2(Icf-mean(Icf));%-mean(Icf)
    Ift = fftshift(Ift);
    Ift = Ift .* Hf;
    Irec = abs(ifft2(ifftshift(Ift)));
    Rs1= abs(ifft2(ifftshift(Ift .* -ur)));
    Rs2= abs(ifft2(ifftshift(Ift .* -vr)));


    %             figure(3)
    %             imagesc(abs(Irec));
    %             Icf = Icf - mean(mean(Icf));
    for ptIdx=1:Np
        m = pt(ptIdx,1); %y coordinate
        n = pt(ptIdx,2); %x coordinate
        r0(ptIdx) = Irec(m,n);
        r1(ptIdx) = Rs1(m,n);
        r2(ptIdx) = Rs2(m,n);

        if frameIDX > 1
            p0=r0(ptIdx)*q0(ptIdx)+r1(ptIdx)*q1(ptIdx)+r2(ptIdx)*q2(ptIdx);
            p1=-r0(ptIdx)*q1(ptIdx)+r1(ptIdx)*q0(ptIdx);
            p2=-r0(ptIdx)*q2(ptIdx)+r2(ptIdx)*q0(ptIdx);
            pphix=(p1/sqrt(p1^2+p2^2))*acos(p0/sqrt(p0^2+p1^2+p2^2));
            pphiy=(p2/sqrt(p1^2+p2^2))*acos(p0/sqrt(p0^2+p1^2+p2^2));

            motion(ptIdx,1,frameIDX) = motion(ptIdx,1,frameIDX)+pphix;
            motion(ptIdx,2,frameIDX) = motion(ptIdx,2,frameIDX)+pphiy;
        end
        q0(ptIdx) = r0(ptIdx);
        q1(ptIdx) = r1(ptIdx);
        q2(ptIdx) = r2(ptIdx);

    end






end
%clear currentPhase previousPhase referencePyramid;


tr=toc
tt=(0:nF-1)*1/samplingRate;
for ptIDX = 1:Np
    % Extract vertical motion at points of interest
    %         yp = squeeze(motion(pt(ptIDX,1), pt(ptIDX, 2), 2, :));
    %         xp = squeeze(motion(pt(ptIDX,1), pt(ptIDX, 2), 1, :));
    xp = squeeze(motion(ptIDX, 1, :));
    yp = squeeze(motion(ptIDX, 2, :));
    xp = mmPerPixel * xp/fs;
    yp = -mmPerPixel * yp/fs;



    if ptIDX == 1
        figure()
        plot(tt,xp);
        ycap=yp;
        figure()
        plot(tt,yp);
        %             xxp=yp.*hanning(nF);
        %             yf=fft(xxp);
        %             fa=(0:nF/2-1)*samplingRate/nF;
        %             figure()
        %             plot(fa,abs(yf(1:nF/2)))
    end
end
% return
yys0=(mmPerPixel/wy)*squeeze(sum(motion(:, 2, :),1))/Np;%averaged value
[B_band, A_band] = butter(2, [loCutoff hiCutoff]/(samplingRate/2));
figure()
plot(tt,yys0);
yys0= filter(B_band, A_band, yys0, []);
figure()
plot(tt,yys0);
stp=511;
yys=yys0(stp:stp+499);

n=length(yys);
tt2=(0:n-1)*1/samplingRate;
figure()
plot(tt2,yys);
xxp=yys.*hanning(n);
yf=fft(xxp);
fa=(0:n/2-1)*samplingRate/n;
figure()
plot(fa,abs(yf(1:n/2)))
% [mxv,i1] = max(abs(yf(1:nF/2)))
pathFolderResults = 'd:/data/simple/forced2yrpof.txt';
%   writematrix(yys0,pathFolderResults,'Delimiter','tab')



