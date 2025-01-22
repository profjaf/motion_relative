%using differential pahse optical flow
%based on:
% https://doi.org/10.1016/j.ymssp.2023.111089
%however, it is not working porperly since we do not know the center
%frequency or wave length, also it is very noisy
%the method need much work and fine tuning in order to be applicable

close all
clear
%   filename='D:/data/point1_new2_cropped.avi';  
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
%    pt = [259 560; %y ,x coordinates
%        259 561;
%        259 562;
%        259 563;
%        259 564;
%        259 566;
%        259 567;
%        259 568;
%        259 569;
%        259 570;
%        259 602;
%        259 661];
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
       [max_value,idx]=max(abs(Ift(:)));
       [iy,ix]=ind2sub(size(Ift),idx);
       wy=2*pi*iy/h;
       wx=2*pi*ix/w;
       
       for ptIdx=1:Np
           m = pt(ptIdx,1); %y coordinate
           n = pt(ptIdx,2); %x coordinate
           pIDx(ptIdx)=Icf(m-1,n+1)+2*Icf(m,n+1)+Icf(m+1,n+1)-(Icf(m-1,n-1)+2*Icf(m,n-1)+Icf(m+1,n-1));
           pIDy(ptIdx)=Icf(m+1,n-1)+2*Icf(m+1,n)+Icf(m+1,n+1)-(Icf(m-1,n-1)+2*Icf(m-1,n)+Icf(m-1,n+1));
       end

        for frameIDX = 1:nF
            fprintf('Detecting motion for frame %d\n', frameIDX);
            Icf = readFrame(frameIDX); %this is luminance data
%             Icf = Icf - mean(mean(Icf));
            for ptIdx=1:Np
                m = pt(ptIdx,1); %y coordinate
                n = pt(ptIdx,2); %x coordinate 
                Icdx =Icf(m,:);
                Icf2 = Icf - mean(Icdx);
                IDx=Icf2(m-1,n+1)+2*Icf2(m,n+1)+Icf2(m+1,n+1)-(Icf2(m-1,n-1)+2*Icf2(m,n-1)+Icf2(m+1,n-1));
                IDAx=Icf2(m,n)-1i*0.37*IDx;
                pphix = phix(ptIdx);
                phix(ptIdx)=angle(IDAx);
               
                Icdy = Icf(:,n);
%                 Icft = fft((Icdy- mean(Icdy)).*hanning(h));
%                 [mxv,i1] = max(abs(Icft(1:h/2)));
%                 w1=2*pi*i1/h
                Icf2 = Icf - mean(Icdy);
                IDy=Icf2(m+1,n-1)+2*Icf2(m+1,n)+Icf2(m+1,n+1)-(Icf2(m-1,n-1)+2*Icf2(m-1,n)+Icf2(m-1,n+1));
                IDAy=Icf2(m,n)-1i*(0.37)*IDy;
                pphiy = phiy(ptIdx);
                phiy(ptIdx) = angle(IDAy);
                

                if frameIDX > 1
                    motion(ptIdx,1,frameIDX) = motion(ptIdx,1,frameIDX)+mod(pi+phix(ptIdx)-pphix, 2*pi)-pi;
                    motion(ptIdx,2,frameIDX) = motion(ptIdx,2,frameIDX)+mod(pi+phiy(ptIdx)-pphiy, 2*pi)-pi;
                end
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
        xp = mmPerPixel * xp;
        yp = mmPerPixel * yp;
%         figure()
%         plot(tt,xp);
      
       
        if ptIDX == 1
            figure()
            plot(tt,yp);
%             xxp=yp.*hanning(nF);
%             yf=fft(xxp);
%             fa=(0:nF/2-1)*samplingRate/nF;
%             figure()
%             plot(fa,abs(yf(1:nF/2)))
        end
    end
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
pathFolderResults = 'd:/data/simple/forced2ydpof.txt';
%   writematrix(yys0,pathFolderResults,'Delimiter','tab')




