%using hilbert pahse optical flow
%based on:
% https://doi.org/10.1016/j.ymssp.2021.108418

%however, it works good especially if we apply temporal bandpass filtering
%for the extracted signal or for the phase difference
%however, this needs that the spatial center frequency be known since the
%theoretical base is that the intensity is a monocomponent which is not
%true in most of the cases and this put challenge on how to scale the
%displacement from the phase differences

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

%    pt = [260 20; %y ,x coordinates 259, 560 filename='D:/data/point2_2_cropped.avi';
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
       Icdy=Icf(:,20);
       figure
       plot(Icdy)
%        return
       Icf = Icf - mean(mean(Icf));
       Ift = fft2(Icf);
       [max_value,idx]=max(abs(Ift(:)));
       [iy,ix]=ind2sub(size(Ift),idx);
       wy=2*pi*iy/h;
       wx=2*pi*ix/w;
       Ift = fftshift(Ift);
       figure
       mesh(abs(Ift))
       
       d1=2;
       d2=8;
       Hf=zeros(h,w);
       for u=1:w
           for v=1:h
               d=sqrt((u-w/2)^2+(v-h/2)^2);
               Hf(v,u)=1/(1+(d/d2)^10)-1/(1+(d/d1)^10);
           end           
       end
       Ift = Ift .* Hf;
       Irec = ifft2(ifftshift(Ift));
       figure
       imagesc(abs(Irec));
        % return
       pIcf = Icf;
       for ptIdx=1:Np
           m = pt(ptIdx,1); %y coordinate
           n = pt(ptIdx,2); %x coordinate
           pIDx(ptIdx)=Icf(m-1,n+1)+2*Icf(m,n+1)+Icf(m+1,n+1)-(Icf(m-1,n-1)+2*Icf(m,n-1)+Icf(m+1,n-1));
           pIDy(ptIdx)=Icf(m+1,n-1)+2*Icf(m+1,n)+Icf(m+1,n+1)-(Icf(m-1,n-1)+2*Icf(m-1,n)+Icf(m-1,n+1));
       end

        for frameIDX = 1:nF
            fprintf('Detecting motion for frame %d\n', frameIDX);
            Icf = readFrame(frameIDX); %this is luminance data
            Ift = fft2(Icf-mean(Icf));%-mean(Icf)
            Ift = fftshift(Ift);
            Ift = Ift .* Hf;
            Irec = abs(ifft2(ifftshift(Ift)));
%             figure(3)
%             imagesc(abs(Irec));
%             Icf = Icf - mean(mean(Icf));
            for ptIdx=1:Np
                m = pt(ptIdx,1); %y coordinate
                n = pt(ptIdx,2); %x coordinate 
                Icdx = Irec(m,:);
                Icdx = Icdx - mean(Icdx);
                %applying Hilbert transform
                xh = fft(Icdx);                
                xh=2.0*xh;
                xh(floor(h/2+1):h)=0;                
                xa=ifft(xh);

                pphix = phix(ptIdx);
                phix(ptIdx)=angle(xa(n));
               
                %Icdy = Irec(:,n);
                Icdy = Icf(:,n);
                Icdy = Icdy - mean(Icdy);
                %applying Hilbert transform
                yh = fft(Icdy);                
                yh=2.0*yh;
                yh(floor(h/2+1):h)=0;                
                ya=ifft(yh);

                pphiy=phiy(ptIdx) ;            
                phiy(ptIdx) = angle(ya(m));
                

                if frameIDX > 1
                    motion(ptIdx,1,frameIDX) = motion(ptIdx,1,frameIDX)+mod(pi+phix(ptIdx)-pphix, 2*pi)-pi;
                    motion(ptIdx,2,frameIDX) = motion(ptIdx,2,frameIDX)+mod(pi+phiy(ptIdx)-pphiy, 2*pi)-pi;
                end

            end
            pIcf = Icf;

            

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
        xp = mmPerPixel * xp/wx;
        yp = mmPerPixel * yp/wy;
%          figure()
%          plot(tt,xp);
      
       
        if ptIDX == 1
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
yys0=(mmPerPixel/wy)*squeeze(sum(motion(:, 2, :),1))/Np;%averaged value
% [B_band, A_band] = butter(2, [10 15]/250);
[B_band, A_band] = butter(2, [loCutoff hiCutoff]/(samplingRate/2));
figure()
plot(tt,yys0);
% return
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
pathFolderResults = 'd:/data/simple/forced2yhpof.txt';
%  writematrix(yys0,pathFolderResults,'Delimiter','tab')



