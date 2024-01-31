%------Long Nguyen---------------------------------------------------------
%------Grant-Funded Researcher---------------------------------------------
%------The University of Adelaide------------------------------------------
%------North Terrace SA 5005-----------------------------------------------

clear all

clc

mFile_name = mfilename('fullpath');

nameCell = regexp(mFile_name,'[\d]+','match');

iTer = str2double(nameCell{1,length(nameCell)});

nTimeStamp = str2double(nameCell{1,length(nameCell)-1});

tHzon = str2double(nameCell{1,length(nameCell)-2});

nVessel = str2double(nameCell{1,length(nameCell)-3});

nCam = str2double(nameCell{1,length(nameCell)-4});

%--------STEP 1: INITIALIZATION--------------------------------------------

% ground truth of vessel position

file_name = sprintf('scenario_%dcams_%dtargets_%dhorizons_%dTimeStamps_%d.mat',nCam,nVessel,tHzon,nTimeStamp,iTer);

folder_name = sprintf('%dcams_%dtargets_%dhorizons_%dTimeStamps',nCam,nVessel,tHzon,nTimeStamp);

load(['/home/lnguyen/rb61/lnguyen/T_ITS_2023/jobs_monarch/' folder_name '/scenarios/' file_name],"posV","iPosV")

%---------------------------------

t = 1; % seconds

hCam = 20; % meters height of the camera from the sea level

Dsafe = 100; % meters

R = 16e3; % meters

pG = 99.97/100; % gate probability

pD = 0.8; % detection probability

lbda = 10^-7; % /m^2

gma = 16;

alpha_1=1;

alpha_2=1;

alpha_3=1;

alpha_4=1;

alpha_5=1;

Dmin=80; %tan((90-20-10.5/2)*pi/180)*hCam; % meters

Dmax=R;

% eThrshld = 7*10^5;

eThrshld = 1*10^4;

rng(iTer,"twister");

%--------------------------------------------------------------------------
% dynamic model of vessels
F_v=cell(1,nVessel);
Q_v=cell(1,nVessel);
for i=1:nVessel
    F_v{1,i}=[1 t 0 0;...
              0 1 0 0;...
              0 0 1 t;...
              0 0 0 1];

    Q_v{1,i}=[36   0     0     0;... 
              0    0     0     0;...
              0    0     36    0;...
              0    0     0     0];

end
%--------------------------------------------------------------------------
% measurement model of static radars (linear and Gaussian distribution)
posR=[95000;0;-95000;0];

H_r=[1 0 0 0;...
     0 0 1 0];

p_r=13;

X_r_pr=cell(1,nTimeStamp);
for i=1:nTimeStamp
    X_r_pr{1,i}=zeros(4,1,nVessel);
end

X_r_up=cell(1,nTimeStamp);
for i=1:nTimeStamp
    X_r_up{1,i}=zeros(4,1,nVessel);
end

P_r_pr=cell(1,nTimeStamp);
for i=1:nTimeStamp
    P_r_pr{1,i}=zeros(4,4,nVessel);
end

P_r_up=cell(1,nTimeStamp);
for i=1:nTimeStamp
    P_r_up{1,i}=zeros(4,4,nVessel);
end

K_r=cell(1,nTimeStamp);
for i=1:nTimeStamp
    K_r{1,i}=zeros(4,2,nVessel);
end

M_r=cell(1,nTimeStamp);
for i=1:nTimeStamp
    M_r{1,i}=zeros(2,1,nVessel);
end
%--------------------------------------------------------------------------
% measurement model of static AIS (linear and Gaussian distribution)
% posA=[-95000;0;95000;0];

H_A=[1 0 0 0;...
     0 0 1 0];

p_A=10^3;

X_A_pr=cell(1,nTimeStamp);
for i=1:nTimeStamp
    X_A_pr{1,i}=zeros(4,1,nVessel);
end

X_A_up=cell(1,nTimeStamp);
for i=1:nTimeStamp
    X_A_up{1,i}=zeros(4,1,nVessel);
end

P_A_pr=cell(1,nTimeStamp);
for i=1:nTimeStamp
    P_A_pr{1,i}=zeros(4,4,nVessel);
end

P_A_up=cell(1,nTimeStamp);
for i=1:nTimeStamp
    P_A_up{1,i}=zeros(4,4,nVessel);
end

K_A=cell(1,nTimeStamp);
for i=1:nTimeStamp
    K_A{1,i}=zeros(4,2,nVessel);
end

M_A=cell(1,nTimeStamp);
for i=1:nTimeStamp
    M_A{1,i}=zeros(2,1,nVessel);
end
%--------------------------------------------------------------------------
% measurement model of moveable IRST (linear and Gaussian distribution)
H_C=cell(1,nCam);

X_C_pr=cell(1,nCam);

X_C_up=cell(1,nCam);

P_C_pr=cell(1,nCam);

P_C_up=cell(1,nCam);

optP_C_up=cell(1,nCam);

K_C=cell(1,nCam);

M_C=cell(1,nCam);

p_C=13;

x_Cv=cell(1,nCam);

distCV=cell(1,nCam);

for i=1:nCam
    H_C{1,i}=[1 0 0 0;...
              0 0 1 0]; % observation matrix

    X_C_pr{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        X_C_pr{1,i}{1,ii}=zeros(4,1,nVessel);
    end
    
    X_C_up{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        X_C_up{1,i}{1,ii}=zeros(4,1,nVessel);
    end

    P_C_pr{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        P_C_pr{1,i}{1,ii}=zeros(4,4,nVessel);
    end
    
    P_C_up{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        P_C_up{1,i}{1,ii}=zeros(4,4,nVessel);
    end
    
    optP_C_up{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        optP_C_up{1,i}{1,ii}=zeros(4,4,nVessel);
    end
    
    K_C{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        K_C{1,i}{1,ii}=zeros(4,2,nVessel);
    end
    
    M_C{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        M_C{1,i}{1,ii}=zeros(2,1,nVessel);
    end
    % decision variables for assigning cameras to target vessels
    x_Cv{1,i}=zeros(nTimeStamp,nVessel);

    distCV{1,i}=zeros(nVessel,nTimeStamp);

end    

beCV=cell(1,nTimeStamp);
for i=1:nTimeStamp
   beCV{1,i}=zeros(nCam,nVessel); 
end
%--------------------------------------------------------------------------
% sVessel = 0; % m/s
% hVessel = 90; % degrees

nAngleLv = 5;
nSpeedLv = 2;

sCamMin = 1*0; % (m/s) = 0 knot in cruise ship
sCamMax = 1*9; % (m/s) = 18 knots in cruise ship

aCamMax = pi; % (rad)

%--------------------------------------------------------------------------

obTmV=zeros(nVessel,nTimeStamp);

optVal=cell(1,nVessel);

optEstVal=cell(1,nVessel);

fnEstVal=cell(1,nVessel);

%--------------------------------------------------------------------------
% Initialize the positions, speeds and heading angles of radars

iPosCam=cell(1,nCam);

posCam=cell(1,nCam);

OrtCam=cell(1,nCam);

cPosCam=cell(1,nCam);

optPosCam=cell(1,nCam);

for i=1:nCam   
    if nCam>=4    
        if i>0&&i<=round(nCam/4)
            iPosCam{1,i}=[0*randi([0 100],1,1);...
                      0*randi([0 100],1,1);...
                      0;...
                      0];
        elseif i>round(nCam/4)&&i<=2*round(nCam/4)
            iPosCam{1,i}=[0*randi([-100 0],1,1);...
                      0*randi([0 100],1,1);...
                      0;...
                      0];
        elseif i>2*round(nCam/4)&&i<=3*round(nCam/4)
            iPosCam{1,i}=[0*randi([-100 0],1,1);...
                      0*randi([-100 0],1,1);...
                      0;...
                      0];
        elseif i>3*round(nCam/4)&&i<=nCam
            iPosCam{1,i}=[0*randi([0 100],1,1);...
                      0*randi([-100 0],1,1);...
                      0;...
                      0];
        end    
    elseif nCam==1 
            iPosCam{1,i}=[0*randi([0 100],1,1);...
                                  0*randi([0 100],1,1);...
                                  0;...
                                  0];
    elseif nCam==2
        if i==1
            iPosCam{1,i}=[0*randi([0 100],1,1);...
                      0*randi([0 100],1,1);...
                      0;...
                      0];
        elseif i==2
            iPosCam{1,i}=[0*randi([-100 0],1,1);...
                      0*randi([0 100],1,1);...
                      0;...
                      0];
        end
    elseif nCam==3
        if i==1
            iPosCam{1,i}=[0*randi([0 100],1,1);...
                      0*randi([0 100],1,1);...
                      0;...
                      0];
        elseif i==2
            iPosCam{1,i}=[0*randi([-100 0],1,1);...
                      0*randi([0 100],1,1);...
                      0;...
                      0];
        elseif i==3
            iPosCam{1,i}=[0*randi([-100 0],1,1);...
                      0*randi([-100 0],1,1);...
                      0;...
                      0];            
        end         
    end

    posCam{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        posCam{1,i}{1,ii}=zeros(4,1); %[x;y;speed;angle]
        if ii==1
            posCam{1,i}{1,ii}=iPosCam{1,i};
        end
    end

    OrtCam{1,i}=zeros(1,nTimeStamp);

    cPosCam{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        cPosCam{1,i}{1,ii}=zeros(4,1); %[x;y;speed;angle]
        if ii==1
            cPosCam{1,i}{1,ii}=iPosCam{1,i};
        end
    end

    optPosCam{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        optPosCam{1,i}{1,ii}=zeros(4,1); %[x;y;speed;angle]
        if ii==1
            optPosCam{1,i}{1,ii}=iPosCam{1,i};
        end
    end

end
%-----STEP 3: THE PROPOSED ALGORITHM---------------------------------------
% tHzon = 10;

exeTime=zeros(1,nTimeStamp-tHzon);

for i=(tHzon+1):nTimeStamp

    % start counting time
    tStart = tic; 

    for ii=(i-tHzon):i
        if ii==i-tHzon
            %--------------------------------------------------------------
            % present
            if ii==1
                for iii=1:nVessel                    
                    % radar
                    % estimate of radars based on Kalman filter
                    % predicted estimate of radars
                    X_r_pr{1,ii}(:,:,iii)=F_v{1,iii}*iPosV{1,iii};
                    P_r_pr{1,ii}(:,:,iii)=F_v{1,iii}*[5 5 0 0;...
                                                      5 5 0 0;...
                                                      0 0 5 5;...
                                                      0 0 5 5]*F_v{1,iii}'...
                                                      +Q_v{1,iii};
                    % AIS
                    % predicted estimate of AIS
                    X_A_pr{1,ii}(:,:,iii)=F_v{1,iii}*iPosV{1,iii};
                    P_A_pr{1,ii}(:,:,iii)=F_v{1,iii}*[5 5 0 0;...
                                                      5 5 0 0;...
                                                      0 0 5 5;...
                                                      0 0 5 5]*F_v{1,iii}'...
                                                      +Q_v{1,iii};
                    % camera
                    % predicted estimate of camera
                    for iiii=1:nCam
                        X_C_pr{1,iiii}{1,ii}(:,:,iii)=F_v{1,iii}*iPosV{1,iii};
                        P_C_pr{1,iiii}{1,ii}(:,:,iii)=F_v{1,iii}*[5 5 0 0;...
                                                                  5 5 0 0;...
                                                                  0 0 5 5;...
                                                                  0 0 5 5]*F_v{1,iii}'...
                                                                  +Q_v{1,iii};
    
                        OrtCam{1,iiii}(1,ii)=0;
                    end
                end
            else
                for iii=1:nVessel
                    % radar
                    % estimate of radars based on Kalman filter
                    % predicted estimate of radars
                    X_r_pr{1,ii}(:,:,iii)=F_v{1,iii}*X_r_up{1,ii-1}(:,:,iii);
                    P_r_pr{1,ii}(:,:,iii)=F_v{1,iii}*P_r_up{1,ii-1}(:,:,iii)...
                                          *F_v{1,iii}'+Q_v{1,iii};
                    % AIS
                    % predicted estimate of AIS
                    X_A_pr{1,ii}(:,:,iii)=F_v{1,iii}*X_A_up{1,ii-1}(:,:,iii);
                    P_A_pr{1,ii}(:,:,iii)=F_v{1,iii}*P_A_up{1,ii-1}(:,:,iii)...
                                          *F_v{1,iii}'+Q_v{1,iii};
                    % camera
                    % predicted estimate of camera
                    for iiii=1:nCam
                        X_C_pr{1,iiii}{1,ii}(:,:,iii)=F_v{1,iii}*X_C_up{1,iiii}{1,ii-1}(:,:,iii);               
                        P_C_pr{1,iiii}{1,ii}(:,:,iii)=F_v{1,iii}*P_C_up{1,iiii}{1,ii-1}(:,:,iii)*F_v{1,iii}'...
                                                      +Q_v{1,iii};
                    end
                end
            end
            %--------------------------------------------------------------
            % radar
            R_r=zeros(2,2,nVessel);
            for iii=1:nVessel
                d_rv=sqrt((posV{1,iii}{1,ii}(1,1)-posR(1,1))^2+(posV{1,iii}{1,ii}(3,1)-posR(3,1))^2);
                sigma_r=p_r/100*d_rv;
                R_r(:,:,iii)=sigma_r^2.*[1 0;0 1];
                % measurement
                M_r{1,ii}(:,:,iii)=mvnrnd([posV{1,iii}{1,ii}(1,1) posV{1,iii}{1,ii}(3,1)],(R_r(:,:,iii)))';
                % kalman gain
                K_r{1,ii}(:,:,iii)=P_r_pr{1,ii}(:,:,iii)*H_r'*(H_r*P_r_pr{1,ii}(:,:,iii)*H_r'+R_r(:,:,iii))^-1;
            end
            %--------------------------------------------------------------
            % AIS
            R_A=zeros(2,2,nVessel);
            for iii=1:nVessel
                %d_Av=sqrt((posV{1,iii}{1,ii}(1,1)-posA(1,1))^2+(posV{1,iii}{1,ii}(3,1)-posA(3,1))^2);
                
                sigma_A=p_A; %/100*d_Av;
                R_A(:,:,iii)=sigma_A^2.*[1 0;0 1];
                % measurement
                M_A{1,ii}(:,:,iii)=mvnrnd([posV{1,iii}{1,ii}(1,1) posV{1,iii}{1,ii}(3,1)],(R_A(:,:,iii)))';
                % kalman gain
                K_A{1,ii}(:,:,iii)=P_A_pr{1,ii}(:,:,iii)*H_A'*(H_A*P_A_pr{1,ii}(:,:,iii)*H_A'+R_A(:,:,iii))^-1;
            end
            %--------------------------------------------------------------
            % camera
            R_C=cell(1,nCam);
            tempCnt=zeros(1,nVessel);

            for iiii=1:nCam
                R_C{1,iiii}=zeros(2,2,nVessel);
                for iii=1:nVessel
                    d_Cv=sqrt((posV{1,iii}{1,ii}(1,1)-posCam{1,iiii}{1,ii}(1,1))^2+(posV{1,iii}{1,ii}(3,1)-posCam{1,iiii}{1,ii}(2,1))^2);
                    sigma_C=p_C/100*d_Cv;
                    R_C{1,iiii}(:,:,iii)=sigma_C^2.*[1 0;0 1]; 
        
                    % measurement
                    if posCam{1,iiii}{1,ii}(1,1)<posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)==posV{1,iii}{1,ii}(3,1)
                        angleCV=0;    

                    elseif posCam{1,iiii}{1,ii}(1,1)<posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)<posV{1,iii}{1,ii}(3,1)
                        angleCV=atan(abs(posCam{1,iiii}{1,ii}(2,1)-posV{1,iii}{1,ii}(3,1))/abs(posCam{1,iiii}{1,ii}(1,1)-posV{1,iii}{1,ii}(1,1)))*180/pi;            
                    
                    elseif posCam{1,iiii}{1,ii}(1,1)==posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)<posV{1,iii}{1,ii}(3,1)
                        angleCV=90;            
                    
                    elseif posCam{1,iiii}{1,ii}(1,1)>posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)<posV{1,iii}{1,ii}(3,1)
                        angleCV=90+atan(abs(posCam{1,iiii}{1,ii}(1,1)-posV{1,iii}{1,ii}(1,1))/abs(posCam{1,iiii}{1,ii}(2,1)-posV{1,iii}{1,ii}(3,1)))*180/pi;            
                    
                    elseif posCam{1,iiii}{1,ii}(1,1)>posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)==posV{1,iii}{1,ii}(3,1)
                        angleCV=180;            
                    
                    elseif posCam{1,iiii}{1,ii}(1,1)>posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)>posV{1,iii}{1,ii}(3,1)
                        angleCV=180+atan(abs(posCam{1,iiii}{1,ii}(2,1)-posV{1,iii}{1,ii}(3,1))/abs(posCam{1,iiii}{1,ii}(1,1)-posV{1,iii}{1,ii}(1,1)))*180/pi;            
                    
                    elseif posCam{1,iiii}{1,ii}(1,1)==posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)>posV{1,iii}{1,ii}(3,1)
                        angleCV=270;            
                    
                    elseif posCam{1,iiii}{1,ii}(1,1)<posV{1,iii}{1,ii}(1,1)&&posCam{1,iiii}{1,ii}(2,1)>posV{1,iii}{1,ii}(3,1)
                        angleCV=270+atan(abs(posCam{1,iiii}{1,ii}(1,1)-posV{1,iii}{1,ii}(1,1))/abs(posCam{1,iiii}{1,ii}(2,1)-posV{1,iii}{1,ii}(3,1)))*180/pi;
                    
                    end
         
                    if Dmin<=d_Cv&&d_Cv<=Dmax&&(OrtCam{1,iiii}(1,ii)-8.4/2<=angleCV&&angleCV<=OrtCam{1,iiii}(1,ii)+8.4/2)
                        M_C{1,iiii}{1,ii}(:,:,iii)=mvnrnd([posV{1,iii}{1,ii}(1,1) posV{1,iii}{1,ii}(3,1)],(R_C{1,iiii}(:,:,iii)))';
                    else                        
                        M_C{1,iiii}{1,ii}(:,:,iii)=inf;

                        tempCnt(1,iii)=tempCnt(1,iii)+1;
                    end

                    % kalman gain
                    K_C{1,iiii}{1,ii}(:,:,iii)=P_C_pr{1,iiii}{1,ii}(:,:,iii)*H_C{1,iiii}'*(H_C{1,iiii}*P_C_pr{1,iiii}{1,ii}(:,:,iii)*H_C{1,iiii}'+R_C{1,iiii}(:,:,iii))^-1;
      
                end
            end
            
            for iii=1:nVessel
                temp=0;
                if tempCnt(1,iii)==nCam
                    if ii==1  
                        temp=0;
                    else
                        for iiii=1:nCam
                            if max(x_Cv{1,iiii}(1:ii-1,iii))==1
                                tempMtrx=find(obTmV(iii,1:ii-1)==inf);
                                if ~isempty(tempMtrx)
                                    temp=1;
                                else
                                    temp=0;
                                end
                            end
                        end
                    end

                    if ii==1
                        if temp==0
                            obTmV(iii,ii) = 1; %randi([100 1000],1,1);

%                             obTmV(iii,ii)=200;
                        else
                            obTmV(iii,ii)=inf;
                        end
                    else  
                        if temp==0
                            obTmV(iii,ii)=obTmV(iii,ii-1)+1;   
                        else
                            obTmV(iii,ii)=inf;
                        end
                    end

                elseif tempCnt(1,iii)<nCam
                    if fnEstVal{1,iii}(1,ii)<=eThrshld                  
                        obTmV(iii,ii)=inf;
                    else
                        obTmV(iii,ii)=obTmV(iii,ii-1)+1;
                    end
                end
            end
            %--------------------------------------------------------------
            allMea=cell(1,nVessel);
            for iii=1:nVessel       
                allMea{1,iii}{1,1}=M_r{1,ii}(:,:,iii);
                allMea{1,iii}{1,2}=M_A{1,ii}(:,:,iii);
                for iiii=1:nCam
                    allMea{1,iii}{1,2+iiii}=M_C{1,iiii}{1,ii}(:,:,iii);
                end
            end
            %--------------------------------------------------------------
            % radars
            for iiiii=1:nVessel
                % predicted measurement
                M_pr_r=H_r*X_r_pr{1,ii}(:,:,iiiii);
                % innovation covariance
                S_r=H_r*P_r_pr{1,ii}(:,:,iiiii)*H_r'+R_r(:,:,iiiii);

                S_r=(S_r+S_r.')/2;

                % validated region
                vr_r=cell(1,length(allMea{1,iiiii})); % radar, AIS, and cameras
                for iii=1:length(allMea{1,iiiii})
                    if isequal(allMea{1,iiiii}{1,iii},[inf;inf])
                         vr_r{1,iii}=[inf;inf];
                    else
                        if (allMea{1,iiiii}{1,iii}-M_pr_r)'*S_r^-1*(allMea{1,iiiii}{1,iii}-M_pr_r)<=gma
                            vr_r{1,iii}=allMea{1,iiiii}{1,iii};
                        else
                            vr_r{1,iii}=[inf;inf];
                        end
                    end
                end
                
                lRaMea=zeros(1,length(vr_r));
                for iii=1:length(lRaMea)
                    if isequal(vr_r{1,iii},[inf;inf])
                        lRaMea(1,iii)=inf;
                    else
                        lRaMea(1,iii)=mvnpdf(vr_r{1,iii}',M_pr_r',S_r)/lbda;
                    end
                end
                
                sum=0;
                for iii=1:length(lRaMea)
                    if lRaMea(1,iii)==inf
                    else
                        sum=sum+lRaMea(1,iii);
                    end
                end
                
                bt_0=(1-pD*pG)/(1-pD*pG+sum);
                bt=zeros(1,length(lRaMea));
                for iii=1:length(bt)
                    if lRaMea(1,iii)==inf
                        bt(1,iii)=inf;
                    else
                        bt(1,iii)=lRaMea(1,iii)/(1-pD*pG+sum);
                    end
                end
                % updated estimate mean
                vSum=zeros(2,1);
                for iii=1:length(bt)
                    if bt(1,iii)==inf
                    else
                        vSum=vSum+bt(1,iii)*(vr_r{1,iii}-H_r*X_r_pr{1,ii}(:,:,iiiii));
                    end
                end
                X_r_up{1,ii}(:,:,iiiii)=X_r_pr{1,ii}(:,:,iiiii)+K_r{1,ii}(:,:,iiiii)*vSum;

                % updated estimate covariance
                pC=P_r_pr{1,ii}(:,:,iiiii)-K_r{1,ii}(:,:,iiiii)*S_r*K_r{1,ii}(:,:,iiiii)';
                
                vSumS=zeros(2,1);
                for iii=1:length(vr_r)
                    if bt(1,iii)==inf
                    else
%                         vSumS=vSumS+bt(1,iii)*vr_r{1,iii}*vr_r{1,iii}';

                        vSumS=vSumS+bt(1,iii)*(vr_r{1,iii}-H_r*X_r_pr{1,ii}(:,:,iiiii))*(vr_r{1,iii}-H_r*X_r_pr{1,ii}(:,:,iiiii))';

                    end
                end
                pS=K_r{1,ii}(:,:,iiiii)*(vSumS-vSum*vSum')*K_r{1,ii}(:,:,iiiii)';  
    
                P_r_up{1,ii}(:,:,iiiii)=bt_0*P_r_pr{1,ii}(:,:,iiiii)+(1-bt_0)*pC+pS;
            end
            %--------------------------------------------------------------
            % AIS
            for iiiii=1:nVessel
                % predicted measurement
                M_pr_A=H_A*X_A_pr{1,ii}(:,:,iiiii);
                % innovation covariance
                S_A=H_A*P_A_pr{1,ii}(:,:,iiiii)*H_A'+R_A(:,:,iiiii);

                S_A=(S_A+S_A.')/2;

                % validated region
                vr_A=cell(1,length(allMea{1,iiiii})); % radar, AIS, and cameras
                for iii=1:length(allMea{1,iiiii})
                    if isequal(allMea{1,iiiii}{1,iii},[inf;inf])
                        vr_A{1,iii}=[inf;inf];
                    else
                        if (allMea{1,iiiii}{1,iii}-M_pr_A)'*S_A^-1*(allMea{1,iiiii}{1,iii}-M_pr_A)<=gma
                            vr_A{1,iii}=allMea{1,iiiii}{1,iii};
                        else
                            vr_A{1,iii}=[inf;inf];
                        end
                    end
                end
                
                lRaMea=zeros(1,length(vr_A));
                for iii=1:length(lRaMea)
                    if isequal(vr_A{1,iii},[inf;inf])
                        lRaMea(1,iii)=inf;
                    else
                        lRaMea(1,iii)=mvnpdf(vr_A{1,iii}',M_pr_A',S_A)/lbda;
                    end
                end
                
                sum=0;
                for iii=1:length(lRaMea)
                    if lRaMea(1,iii)==inf
                    else
                        sum=sum+lRaMea(1,iii);
                    end
                end
                
                bt_0=(1-pD*pG)/(1-pD*pG+sum);
                bt=zeros(1,length(lRaMea));
                for iii=1:length(bt)
                    if lRaMea(1,iii)==inf
                        bt(1,iii)=inf;
                    else
                        bt(1,iii)=lRaMea(1,iii)/(1-pD*pG+sum);
                    end
                end

                % updated estimate mean
                vSum=zeros(2,1);
                for iii=1:length(bt)
                    if bt(1,iii)==inf
                    else
                        vSum=vSum+bt(1,iii)*(vr_A{1,iii}-H_A*X_A_pr{1,ii}(:,:,iiiii));
                    end
                end
                X_A_up{1,ii}(:,:,iiiii)=X_A_pr{1,ii}(:,:,iiiii)+K_A{1,ii}(:,:,iiiii)*vSum;

                % updated estimate covariance
                pC=P_A_pr{1,ii}(:,:,iiiii)-K_A{1,ii}(:,:,iiiii)*S_A*K_A{1,ii}(:,:,iiiii)';
                
                vSumS=zeros(2,1);
                for iii=1:length(vr_A)
                    if bt(1,iii)==inf
                    else
%                         vSumS=vSumS+bt(1,iii)*vr_A{1,iii}*vr_A{1,iii}';

                        vSumS=vSumS+bt(1,iii)*(vr_A{1,iii}-H_A*X_A_pr{1,ii}(:,:,iiiii))*(vr_A{1,iii}-H_A*X_A_pr{1,ii}(:,:,iiiii))';
                        
                    end
                end
                pS=K_A{1,ii}(:,:,iiiii)*(vSumS-vSum*vSum')*K_A{1,ii}(:,:,iiiii)';  
                
                P_A_up{1,ii}(:,:,iiiii)=bt_0*P_A_pr{1,ii}(:,:,iiiii)+(1-bt_0)*pC+pS;
            end
            %--------------------------------------------------------------
            for iiiii=1:nVessel
                for iiii=1:nCam
                    % camera
                    % predicted measurement
                    M_pr_C=H_C{1,iiii}*X_C_pr{1,iiii}{1,ii}(:,:,iiiii);
                    % innovation covariance
                    S_C=H_C{1,iiii}*P_C_pr{1,iiii}{1,ii}(:,:,iiiii)*H_C{1,iiii}'+R_C{1,iiii}(:,:,iiiii);

                    S_C=(S_C+S_C.')/2;
 
                    % validated region
                    vr_C=cell(1,length(allMea{1,iiiii})); % radar, AIS, and camera
                    for iii=1:length(allMea{1,iiiii})
                        if isequal(allMea{1,iiiii}{1,iii},[inf;inf])
                            vr_C{1,iii}=[inf;inf];
                        else
                            if (allMea{1,iiiii}{1,iii}-M_pr_C)'*S_C^-1*(allMea{1,iiiii}{1,iii}-M_pr_C)<=gma
                                vr_C{1,iii}=allMea{1,iiiii}{1,iii};
                            else
                                vr_C{1,iii}=[inf;inf];
                            end
                        end
                    end
                    
                    lRaMea=zeros(1,length(vr_C));
                    for iii=1:length(lRaMea)
                        if isequal(vr_C{1,iii},[inf;inf])
                            lRaMea(1,iii)=inf;
                        else
                            lRaMea(1,iii)=mvnpdf(vr_C{1,iii}',M_pr_C',S_C)/lbda;
                        end
                    end
                    
                    sum=0;
                    for iii=1:length(lRaMea)
                        if lRaMea(1,iii)==inf
                        else
                            sum=sum+lRaMea(1,iii);
                        end
                    end
                    
                    bt_0=(1-pD*pG)/(1-pD*pG+sum);
                    bt=zeros(1,length(lRaMea));
                    for iii=1:length(bt)
                        if lRaMea(1,iii)==inf
                            bt(1,iii)=inf;
                        else
                            bt(1,iii)=lRaMea(1,iii)/(1-pD*pG+sum);
                        end
                    end
                    % updated estimate mean
                    vSum=zeros(2,1);
                    for iii=1:length(bt)
                        if bt(1,iii)==inf
                        else
                            vSum=vSum+bt(1,iii)*(vr_C{1,iii}-H_C{1,iiii}*X_C_pr{1,iiii}{1,ii}(:,:,iiiii));
                        end
                    end
                    X_C_up{1,iiii}{1,ii}(:,:,iiiii)=X_C_pr{1,iiii}{1,ii}(:,:,iiiii)+K_C{1,iiii}{1,ii}(:,:,iiiii)*vSum;
                    % updated estimate covariance
                    pC=P_C_pr{1,iiii}{1,ii}(:,:,iiiii)-K_C{1,iiii}{1,ii}(:,:,iiiii)*S_C*K_C{1,iiii}{1,ii}(:,:,iiiii)';
                    
                    vSumS=zeros(2,1);
                    for iii=1:length(vr_C)
                        if bt(1,iii)==inf
                        else
%                             vSumS=vSumS+bt(1,iii)*vr_C{1,iii}*vr_C{1,iii}';

                            vSumS=vSumS+bt(1,iii)*(vr_C{1,iii}-H_C{1,iiii}*X_C_pr{1,iiii}{1,ii}(:,:,iiiii))*(vr_C{1,iii}-H_C{1,iiii}*X_C_pr{1,iiii}{1,ii}(:,:,iiiii))';
                       
                        end
                    end
                    pS=K_C{1,iiii}{1,ii}(:,:,iiiii)*(vSumS-vSum*vSum')*K_C{1,iiii}{1,ii}(:,:,iiiii)';  
                    
                    P_C_up{1,iiii}{1,ii}(:,:,iiiii)=bt_0*P_C_pr{1,iiii}{1,ii}(:,:,iiiii)+(1-bt_0)*pC+pS;

                end
            end
            % assign a camera to a target vessel
            for iiiii=1:nVessel
                for iiii=1:nCam
                    if ii == 1
                        distCV{1,iiii}(iiiii,ii)=sqrt((X_C_up{1,iiii}{1,ii}(1,1,iiiii)-posCam{1,iiii}{1,ii}(1,1))^2 ...
                            +(X_C_up{1,iiii}{1,ii}(3,1,iiiii)-posCam{1,iiii}{1,ii}(2,1))^2);
                       
                        beCV{1,ii}(iiii,iiiii) = 0/distCV{1,iiii}(iiiii,ii); %obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii); 

                    else
                        if nCam==nVessel
                            if x_Cv{1,iiii}(ii-1,iiiii)==1
                                if obTmV(iiiii,ii)~=inf
                                    beCV{1,ii}(iiii,iiiii)=10^20;
                                else
                                    beCV{1,ii}(iiii,iiiii)=-10;
                                end            
                            else
                                distCV{1,iiii}(iiiii,ii)=sqrt((X_C_up{1,iiii}{1,ii}(1,1,iiiii)-posCam{1,iiii}{1,ii}(1,1))^2 ...
                                    +(X_C_up{1,iiii}{1,ii}(3,1,iiiii)-posCam{1,iiii}{1,ii}(2,1))^2);

                                if obTmV(iiiii,ii)~=inf
                                    if ii < 1000*0
                                        beCV{1,ii}(iiii,iiiii) = randi([8 20],1,1)*10^5/distCV{1,iiii}(iiiii,ii);
                                    else
                                        beCV{1,ii}(iiii,iiiii)=fnEstVal{1,iiiii}(1,ii)/distCV{1,iiii}(iiiii,ii); %obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii); 
                                    end
                                else
                                    beCV{1,ii}(iiii,iiiii)=-10;
                                end
                                
                            end
                        else
                            if x_Cv{1,iiii}(ii-1,iiiii)==1
                                if obTmV(iiiii,ii)~=inf
                                    beCV{1,ii}(iiii,iiiii)=10^20;
                                else
                                    beCV{1,ii}(iiii,iiiii)=-10;
                                end    
                            else
                                distCV{1,iiii}(iiiii,ii)=sqrt((X_C_up{1,iiii}{1,ii}(1,1,iiiii)-posCam{1,iiii}{1,ii}(1,1))^2 ...
                                    +(X_C_up{1,iiii}{1,ii}(3,1,iiiii)-posCam{1,iiii}{1,ii}(2,1))^2);
            
                                if iiiii==1
                                    if obTmV(iiiii,ii)~=inf
                                        if ii < 1000*0
                                            beCV{1,ii}(iiii,iiiii) = randi([8 20],1,1)*10^5/distCV{1,iiii}(iiiii,ii);
                                        else
                                            beCV{1,ii}(iiii,iiiii)=fnEstVal{1,iiiii}(1,ii)/distCV{1,iiii}(iiiii,ii); %obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii); 
                                        end                                        
                                    else
                                        beCV{1,ii}(iiii,iiiii)=-10;
                                    end
                                else
                                    for iiiiii=1:iiiii-1
                                        if x_Cv{1,iiii}(ii,iiiiii)==1
                                            beCV{1,ii}(iiii,iiiii)=-10;                                            
                                            break
                                        else
                                            if obTmV(iiiii,ii)~=inf
                                                if ii < 1000*0
                                                    beCV{1,ii}(iiii,iiiii) = randi([8 20],1,1)*10^5/distCV{1,iiii}(iiiii,ii);
                                                else
                                                    beCV{1,ii}(iiii,iiiii)=fnEstVal{1,iiiii}(1,ii)/distCV{1,iiii}(iiiii,ii); %obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii); 
                                                end
                                            else
                                                beCV{1,ii}(iiii,iiiii)=-10;
                                            end                                            
                                        end
                                    end
                                end 
                            end
                        end
                    end
                end               
            end
            %----------------------------------------------------------
            tempBeCV=beCV{1,ii};
            for iiiiii=1:nCam
                tempMax=-1;
                tempCam=0;
                tempVsl=0;
                for iiii=1:nCam
                    for iiiii=1:nVessel
                        if tempBeCV(iiii,iiiii)>tempMax
                            tempMax=tempBeCV(iiii,iiiii);
                            tempCam=iiii;
                            tempVsl=iiiii;
                        end
                    end
                end
                %----------------------------------------------------------
                if tempMax~=-1
                    x_Cv{1,tempCam}(ii,tempVsl)=1;                    
                    for iiii=1:nCam
                        for iiiii=1:nVessel
                            if iiii==tempCam||iiiii==tempVsl
                                tempBeCV(iiii,iiiii)=-10;
                            end
                        end
                    end
                else
           
                end
            end

        else % fromm ii

            %--------------------------------------------------------------
            % future
            % radar
            for iiiii=1:nVessel

                if ii==i-tHzon+1
                else
                    X_r_up{1,ii-1}(:,:,iiiii) = X_r_pr{1,ii-1}(:,:,iiiii);
                end
                
                X_r_pr{1,ii}(:,:,iiiii)=F_v{1,iiiii}*X_r_up{1,ii-1}(:,:,iiiii);

                P_r_pr{1,ii}(:,:,iiiii)=F_v{1,iiiii}*P_r_up{1,ii-1}(:,:,iiiii)*F_v{1,iiiii}'+Q_v{1,iiiii};

                d_rv=sqrt((X_r_pr{1,ii}(1,1,iiiii)-posR(1,1))^2+(X_r_pr{1,ii}(3,1,iiiii)-posR(3,1))^2);

                sigma_r=p_r/100*d_rv;

                R_r=sigma_r^2.*[1 0;0 1];

                % updated estimate of radars
                S_r=H_r'*(R_r)^-1*H_r;

                P_r_up{1,ii}(:,:,iiiii)=((P_r_pr{1,ii}(:,:,iiiii))^-1+S_r)^-1;
            end
            %--------------------------------------------------------------
            % AIS
            for iiiii=1:nVessel

                if ii==i-tHzon+1
                else
                    X_A_up{1,ii-1}(:,:,iiiii) = X_A_pr{1,ii-1}(:,:,iiiii);
                end

                X_A_pr{1,ii}(:,:,iiiii)=F_v{1,iiiii}*X_A_up{1,ii-1}(:,:,iiiii);

                P_A_pr{1,ii}(:,:,iiiii)=F_v{1,iiiii}*P_A_up{1,ii-1}(:,:,iiiii)*F_v{1,iiiii}'+Q_v{1,iiiii};

                %d_Av=sqrt((X_A_pr{1,ii}(1,1,iiiii)-posA(1,1))^2+(X_A_pr{1,ii}(3,1,iiiii)-posA(3,1))^2);

                sigma_A=p_A; %/100*d_Av;

                R_A=sigma_A^2.*[1 0;0 1];

                % updated estimate of AIS
                S_A=H_A'*(R_A)^-1*H_A;

                P_A_up{1,ii}(:,:,iiiii)=((P_A_pr{1,ii}(:,:,iiiii))^-1+S_A)^-1;
            end
            %--------------------------------------------------------------
            % camera           
            for iiiiii=1:nCam
                for iiiii=1:nVessel
                    if x_Cv{1,iiiiii}(i-tHzon,iiiii)==1 %x_Cv{1,iiiiii}(ii-1,iiiii)==1

                        if ii==i-tHzon+1
                        else
                            X_C_up{1,iiiiii}{1,ii-1}(:,:,iiiii) = X_C_pr{1,iiiiii}{1,ii-1}(:,:,iiiii);                        
                        end

                        X_C_pr{1,iiiiii}{1,ii}(:,:,iiiii)=F_v{1,iiiii}*X_C_up{1,iiiiii}{1,ii-1}(:,:,iiiii);

                        P_C_pr{1,iiiiii}{1,ii}(:,:,iiiii)=F_v{1,iiiii}*P_C_up{1,iiiiii}{1,ii-1}(:,:,iiiii)*F_v{1,iiiii}'+Q_v{1,iiiii};
            
                        maxObjVal=-10e20;
            
                        for iii=1:2*nAngleLv %iii=1:2*(nAngleLv+1)                
                            if iii<=nAngleLv %iii<=nAngleLv+1
                                if ii==2
%                                         cPosCam{1,iiiiii}{1,ii-1}(4,1)=0-(iii-1)*aCamMax/nAngleLv;
                                      cPosCam{1,iiiiii}{1,ii-1}(4,1)=0-(iii-1)*aCamMax/(nAngleLv-1);

                                else
%                                         cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)-(iii-1)*aCamMax/nAngleLv;
                                      cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)-(iii-1)*aCamMax/(nAngleLv-1);

                                end
                            else
                                if ii==2
%                                         cPosCam{1,iiiiii}{1,ii-1}(4,1)=0+(iii-nAngleLv-2)*aCamMax/nAngleLv;
                                      cPosCam{1,iiiiii}{1,ii-1}(4,1)=0+(iii-nAngleLv-1)*aCamMax/(nAngleLv-1);

                                else
%                                         cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)+(iii-nAngleLv-2)*aCamMax/nAngleLv;
                                      cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)+(iii-nAngleLv-1)*aCamMax/(nAngleLv-1);

                                end
                            end
            
                            %------------------------------------------------------------------
                            for iiii=1:nSpeedLv %iiii=1:(nSpeedLv+1)                
                                %cPosCam{1,iiiiii}{1,ii-1}(3,1)=sCamMin+(iiii-1)*(sCamMax-sCamMin)/nSpeedLv; 

                                cPosCam{1,iiiiii}{1,ii-1}(3,1)=sCamMin+(iiii-1)*(sCamMax-sCamMin)/(nSpeedLv-1);
            
                                cPosCam{1,iiiiii}{1,ii}(1,1)=cPosCam{1,iiiiii}{1,ii-1}(1,1)...
                                    +t*cPosCam{1,iiiiii}{1,ii-1}(3,1)*cos(cPosCam{1,iiiiii}{1,ii-1}(4,1));
            
                                cPosCam{1,iiiiii}{1,ii}(2,1)=cPosCam{1,iiiiii}{1,ii-1}(2,1)...
                                    +t*cPosCam{1,iiiiii}{1,ii-1}(3,1)*sin(cPosCam{1,iiiiii}{1,ii-1}(4,1));                
                    
                                d_Cv=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii))^2 ...
                                    +(cPosCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii))^2);  
                                %------------------------------------------
                                % term 1 for detection
                                if d_Cv>=Dmin
                                    dTerm1=1;

                                else
                                    dTerm1=Dmin-d_Cv;

                                end

                                % term 2 for detection
                                if d_Cv<=Dmax
                                    dTerm2=1;

                                else
                                    dTerm2=d_Cv-Dmax;

                                end
                                %------------------------------------------
                                % term 3 for collision avoidance between
                                % cameras and vessels
                                dTerm3=0;

                                for iiiiiii=1:nVessel
                                    d_Cv_2=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiiiii))^2 ...
                                        +(cPosCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiiiii))^2);  

                                    if d_Cv_2<Dsafe
                                        dTerm3=dTerm3+(Dsafe-d_Cv_2);

                                    else
                                        dTerm3=dTerm3+0;

                                    end
                                end
                                %------------------------------------------
                                % term 4 for collision avoidance between
                                % cameras and another cameras
                                dTerm4=0;

                                if iiiiii==1
                                    dTerm4=0;

                                else
                                    for iiiiiii=1:iiiiii-1
                                        dCC=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-cPosCam{1,iiiiiii}{1,ii}(1,1))^2 ...
                                            +(cPosCam{1,iiiiii}{1,ii}(2,1)-cPosCam{1,iiiiiii}{1,ii}(2,1))^2);

                                        if dCC<Dsafe
                                            dTerm4=dTerm4+(Dsafe-dCC);

                                        else
                                            dTerm4=dTerm4+0;

                                        end
                                    end        
                                end    
                                %------------------------------------------
                                if d_Cv>=Dmin&&d_Cv<=Dmax
                                    sigma_C=p_C/100*d_Cv;
                                    R_C=sigma_C^2.*[1 0;0 1]; 
                                    S_C_up=H_C{1,iiiiii}'*(R_C)^-1*H_C{1,iiiiii};
                                    P_C_up{1,iiiiii}{1,ii}(:,:,iiiii)=((P_C_pr{1,iiiiii}{1,ii}(:,:,iiiii))^-1+S_C_up)^-1;   

                                elseif d_Cv<Dmin
                                    P_C_up{1,iiiiii}{1,ii}(:,:,iiiii)=P_C_pr{1,iiiiii}{1,ii}(:,:,iiiii);

                                elseif d_Cv>Dmax
                                    P_C_up{1,iiiiii}{1,ii}(:,:,iiiii)=P_C_pr{1,iiiiii}{1,ii}(:,:,iiiii);    
                                end  

                                dTerm = trace(((P_r_up{1,ii}(:,:,iiiii))^-1+(P_A_up{1,ii}(:,:,iiiii))^-1 ...
                                    +(P_C_up{1,iiiiii}{1,ii}(:,:,iiiii))^-1)^-1);

                                if dTerm  <= eThrshld
                                    dTerm = 0;
                                else
                                    dTerm = dTerm-eThrshld;  
                                end
                                                               
                                objVal=-(alpha_1*dTerm...
                                    +alpha_3*dTerm1*dTerm2...
                                    +alpha_4*dTerm3...
                                    +alpha_5*dTerm4);

                                estVal=trace(((P_r_up{1,ii}(:,:,iiiii))^-1+(P_A_up{1,ii}(:,:,iiiii))^-1 ...
                                    +(P_C_up{1,iiiiii}{1,ii}(:,:,iiiii))^-1)^-1);
                                
                                if maxObjVal<objVal
            
                                    maxObjVal=objVal;

%                                     maxObs = dTerm;

                                    maxEstVal=estVal;
            
                                    optPosCam{1,iiiiii}{1,ii}(1,1)=cPosCam{1,iiiiii}{1,ii}(1,1);
                                    optPosCam{1,iiiiii}{1,ii}(2,1)=cPosCam{1,iiiiii}{1,ii}(2,1);
                                    optPosCam{1,iiiiii}{1,ii-1}(3,1)=cPosCam{1,iiiiii}{1,ii-1}(3,1);
                                    optPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-1}(4,1);
            
                                    optP_C_up{1,iiiiii}{1,ii}(:,:,iiiii)=P_C_up{1,iiiiii}{1,ii}(:,:,iiiii);
                                end
                            end
                        end 
            
                        posCam{1,iiiiii}{1,ii}(1,1)=optPosCam{1,iiiiii}{1,ii}(1,1);
                        posCam{1,iiiiii}{1,ii}(2,1)=optPosCam{1,iiiiii}{1,ii}(2,1);
                        posCam{1,iiiiii}{1,ii-1}(3,1)=optPosCam{1,iiiiii}{1,ii-1}(3,1);
                        posCam{1,iiiiii}{1,ii-1}(4,1)=optPosCam{1,iiiiii}{1,ii-1}(4,1);
            
                        P_C_up{1,iiiiii}{1,ii}=optP_C_up{1,iiiiii}{1,ii};
            
                        if posCam{1,iiiiii}{1,ii}(1,1)<X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)==X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=0;

                        elseif posCam{1,iiiiii}{1,ii}(1,1)<X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)<X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=atan(abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii))/abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)))*180/pi;
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)==X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)<X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=90;
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)>X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)<X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=90+atan(abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii))/abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)))*180/pi;
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)>X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)==X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=180;
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)>X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)>X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=180+atan(abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii))/abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)))*180/pi;
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)==X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)>X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=270;
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)<X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii)&&posCam{1,iiiiii}{1,ii}(2,1)>X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)
                            OrtCam{1,iiiiii}(1,ii)=270+atan(abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiii))/abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiii)))*180/pi;
                        
                        end
            
                        cPosCam{1,iiiiii}{1,ii}(1,1)=optPosCam{1,iiiiii}{1,ii}(1,1);
                        cPosCam{1,iiiiii}{1,ii}(2,1)=optPosCam{1,iiiiii}{1,ii}(2,1);
                        cPosCam{1,iiiiii}{1,ii-1}(3,1)=optPosCam{1,iiiiii}{1,ii-1}(3,1);
                        cPosCam{1,iiiiii}{1,ii-1}(4,1)=optPosCam{1,iiiiii}{1,ii-1}(4,1);
            
                        optVal{1,iiiii}(1,ii)=maxObjVal;

                        optEstVal{1,iiiii}(1,ii)=maxEstVal;
                    else

                    end
                end
            end

            for iiiii=1:nVessel
                for iiiiii=1:nCam
                    if x_Cv{1,iiiiii}(ii-1,iiiii)==1                       
                        fnEstVal{1,iiiii}(1,ii)=optEstVal{1,iiiii}(1,ii);
                        break
                    else
                        fnEstVal{1,iiiii}(1,ii)=trace(((P_r_up{1,ii}(:,:,iiiii))^-1+(P_A_up{1,ii}(:,:,iiiii))^-1)^-1);
                    end
                end
            end

        end
    end

    % stop counting time
    tempTime = toc(tStart);
    
    exeTime(1,i-tHzon)=exeTime(1,i-tHzon)+tempTime;
end

%plot(fnEstVal{1,1})

finExeTime = mean(exeTime);

% tempObTmV=zeros(1,nVessel);
% for i=1:nVessel
%     temp=find(obTmV(i,:)==inf);
%     if ~isempty(temp)
%         tempObTmV(1,i)=temp(1,1);
%     else
%         tempObTmV(1,i)=nTimeStamp-1;
%     end
% end

tempObTmV=zeros(1,nVessel);
for i=1:nVessel
    temp=find(obTmV(i,:)==inf);
    if ~isempty(temp)
        tempObTmV(1,i)=obTmV(i,temp(1,1)-1);
    else
        tempObTmV(1,i)=obTmV(i,nTimeStamp-1);
    end
end

finObTmV=max(tempObTmV(1,:));

disp('Done!')

%--------------------------------------------------------------------------
for i=1:nCam
    for ii=2:nTimeStamp
        posCam{1,i}{1,ii}(1,1)=posCam{1,i}{1,ii-1}(1,1)+t*posCam{1,i}{1,ii-1}(3,1)*cos(posCam{1,i}{1,ii-1}(4,1));
        posCam{1,i}{1,ii}(2,1)=posCam{1,i}{1,ii-1}(2,1)+t*posCam{1,i}{1,ii-1}(3,1)*sin(posCam{1,i}{1,ii-1}(4,1));   
    end
end
%--------------------------------------------------------------------------
file_name = sprintf('output_GD_%dcams_%dtargets_%dhorizons_%dTimeStamps_%d.mat',nCam,nVessel,tHzon,nTimeStamp,iTer);

save(['/home/lnguyen/rb61/lnguyen/T_ITS_2023/jobs_monarch/' folder_name '/output/GD/' file_name],"posCam","obTmV","OrtCam","fnEstVal","exeTime","finExeTime","finObTmV");
