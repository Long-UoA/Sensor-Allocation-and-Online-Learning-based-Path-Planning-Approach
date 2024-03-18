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

%--------STEP 1: INITIALIZATION---------------------------------------------
% ground truth of vessel position

file_name = sprintf('scenario_%dcams_%dtargets_%dhorizons_%dTimeStamps_%d.mat',nCam,nVessel,tHzon,nTimeStamp,iTer);

folder_name = sprintf('%dcams_%dtargets_%dhorizons_%dTimeStamps',nCam,nVessel,tHzon,nTimeStamp);

load(file_name,"posV","iPosV")

%----------------------------------------------------

t = 1; % second

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
% posA=[-1000000000;0;-1000000000;0];

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
% sVessel = 1; % m/s
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

tempPosCam=cell(1,nCam);

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

    tempPosCam{1,i}=cell(1,nTimeStamp);
    for ii=1:nTimeStamp
        tempPosCam{1,i}{1,ii}=zeros(4,1); %[x;y;speed;angle]
        if ii==1
            tempPosCam{1,i}{1,ii}=iPosCam{1,i};
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
%--------------------------------------------------------------------------
nameCamAct=cell(nSpeedLv*(2*(nAngleLv)),1);
for i=1:length(nameCamAct)
    nameCamAct{i,1}=[0,0];
end

for ii=1:nSpeedLv
    for i=1:length(nameCamAct)
        if i>(ii-1)*(2*(nAngleLv))&&i<=ii*(2*(nAngleLv))
            nameCamAct{i,1}(1,1)=ii;
            nameCamAct{i,1}(1,2)=i-(ii-1)*(2*(nAngleLv));            
        end       
    end
end

probCamAct=cell(nCam,nTimeStamp);
for i=1:nCam
    for ii=1:nTimeStamp
        probCamAct{i,ii}=zeros(1,length(nameCamAct));
        if ii==1
            for iii=1:length(probCamAct{i,ii})
                probCamAct{i,ii}(1,iii)=1/length(probCamAct{i,ii});
            end
        end
    end
end

camAct=cell(nCam,nTimeStamp);
camActSpdAng=cell(nCam,nTimeStamp);
for i=1:nCam
    for ii=1:nTimeStamp
        camAct{i,ii}=[0,0]; % speed,heading angle
        camActSpdAng{i,ii}=[0,0];
    end
end
%========================
nIter = 20;
%========================
probCamActIter=cell(nCam,nIter+1);
for i=1:nCam
    for ii=1:nIter
        probCamActIter{i,ii}=zeros(1,length(nameCamAct));
        if ii==1
            for iii=1:length(probCamActIter{i,ii})
                probCamActIter{i,ii}(1,iii)=1/length(probCamActIter{i,ii});
            end
        end
    end
end

camActIter=cell(nCam,nIter);
for i=1:nCam
    for ii=1:nIter
        camActIter{i,ii}=[0,0]; % speed,heading angle
    end
end

camUtl=cell(nCam,nIter);
for i=1:nCam
    for ii=1:nIter
        camUtl{i,ii}=zeros(1,length(nameCamAct));
    end
end

sumRegret=cell(1,nIter);
for i=1:nIter
    sumRegret{1,i}=zeros(nCam,length(nameCamAct),length(nameCamAct));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
end

epsilon = 0.01;

lambda=0.5;

%-----STEP 3: THE PROPOSED ALGORITHM---------------------------------------
% tHzon = 1;

exeTime=zeros(1,nTimeStamp- tHzon);

fnMSE=cell(1,nVessel);
for ii=1:nVessel
    fnMSE{1,ii}(1,nTimeStamp)=inf;
end

fnEst=cell(1,nVessel);
for ii=1:nVessel
    fnEst{1,ii}(1,nTimeStamp)=inf;
end

for i=(tHzon+1):nTimeStamp

    i

    % starting time
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
                    P_r_pr{1,ii}(:,:,iii)=F_v{1,iii}*(1/25)*[5 5 0 0;...
                                                      5 5 0 0;...
                                                      0 0 5 5;...
                                                      0 0 5 5]*F_v{1,iii}'...
                                                      +Q_v{1,iii};
                    % AIS
                    % predicted estimate of AIS
                    X_A_pr{1,ii}(:,:,iii)=F_v{1,iii}*iPosV{1,iii};
                    P_A_pr{1,ii}(:,:,iii)=F_v{1,iii}*(1/25)*[5 5 0 0;...
                                                      5 5 0 0;...
                                                      0 0 5 5;...
                                                      0 0 5 5]*F_v{1,iii}'...
                                                      +Q_v{1,iii};
                    % camera
                    % predicted estimate of camera
                    for iiii=1:nCam
                        X_C_pr{1,iiii}{1,ii}(:,:,iii)=F_v{1,iii}*iPosV{1,iii};
                        P_C_pr{1,iiii}{1,ii}(:,:,iii)=F_v{1,iii}*(1/25)*[5 5 0 0;...
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
%                 d_Av=sqrt((posV{1,iii}{1,ii}(1,1)-posA(1,1))^2+(posV{1,iii}{1,ii}(3,1)-posA(3,1))^2);
                
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
                            
%                             obTmV(iii,ii)=obTmV(iii,ii-1)+1; 
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
                    if fnEst{1,iii}(1,ii) <= eThrshld %1.5*10^5                   
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

                S_r=((S_r+S_r.')/2);
%                 S_r(1,1) = abs(S_r(1,1));
%                 S_r(2,2) = abs(S_r(2,2));

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

                S_A=((S_A+S_A.')/2);                
                
                %S_A(1,1) = abs(S_A(1,1));
                
                %S_A(2,2) = abs(S_A(2,2));

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

%                     S_C(1,1) = abs(S_C(1,1));
%                     S_C(2,2) = abs(S_C(2,2));
 
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
                    if ii==1
                        distCV{1,iiii}(iiiii,ii)=sqrt((X_C_up{1,iiii}{1,ii}(1,1,iiiii)-posCam{1,iiii}{1,ii}(1,1))^2 ...
                            +(X_C_up{1,iiii}{1,ii}(3,1,iiiii)-posCam{1,iiii}{1,ii}(2,1))^2);
                       
                        beCV{1,ii}(iiii,iiiii)=obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii); 

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
                                    beCV{1,ii}(iiii,iiiii)=obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii); 
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
                                        beCV{1,ii}(iiii,iiiii)=obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii); 
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
                                                beCV{1,ii}(iiii,iiiii)=obTmV(iiiii,ii)/distCV{1,iiii}(iiiii,ii);
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

        else % from ii to the future
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
            % REGRET-MATCHING-BASED ONLINE LEARNING ALGORITHM
            % the position of each camera in the future after playing an
            % action
            %--------------------------------------------------------------                        
            % REGRET-MATCHING-BASED ONLINE LEARNING ALGORITHM
            % a cam select its action at time t
            for iii=1:nCam % the number of cameras (players)
                for iiiiiii=1:nVessel
                    if x_Cv{1,iii}(i-tHzon,iiiiiii)==1
                        currTrg=iiiiiii;
                        break
                    else
                        currTrg=0;
                    end
                end

                if currTrg~=0
                    temp=rand;%tempRef;
                    iiii=1;    
                    iiiii=0; 
                    while (1)
                       if iiii>length(probCamAct{iii,ii-1})
                           break;
                       end
                       iiiii=iiiii+probCamAct{iii,ii-1}(1,iiii);
                       if temp<iiiii
                          break; 
                       end
                       iiii=iiii+1;
                    end
                    camAct{iii,ii-1}=nameCamAct{iiii,1};  % player 'i' selects action 'ii' at iteration 'it'
                                            
                    currSpLv=camAct{iii,ii-1}(1,1);
    
                    camActSpdAng{iii,ii-1}(1,1)=sCamMin+(currSpLv-1)*(sCamMax-sCamMin)/(nSpeedLv-1); 
    
                    currHaLv=camAct{iii,ii-1}(1,2);
    
                    if currHaLv<=nAngleLv
                        if ii==2
                            camActSpdAng{iii,ii-1}(1,2)=(0-(currHaLv-1)*aCamMax/(nAngleLv-1))*180/pi;
                        else
                            camActSpdAng{iii,ii-1}(1,2)=(camActSpdAng{iii,ii-2}(1,2)*pi/180-(currHaLv-1)*aCamMax/(nAngleLv-1))*180/pi;
                        end
                    else
                        if ii==2
                            camActSpdAng{iii,ii-1}(1,2)=(0+(currHaLv-nAngleLv-1)*aCamMax/(nAngleLv-1))*180/pi;
                        else
                            camActSpdAng{iii,ii-1}(1,2)=(camActSpdAng{iii,ii-2}(1,2)*pi/180+(currHaLv-nAngleLv-1)*aCamMax/(nAngleLv-1))*180/pi;
                        end
                    end 
                                
                end
            end
            %--------------------------------------------------------------
            for iii=1:nCam
                camActIter{iii,1} = camAct{iii,ii-1};
            end 

            for iiiiiii=1:nVessel

                fnEst{1,iiiiiii}(1,ii)=trace(((P_r_up{1,ii}(:,:,iiiiiii))^-1+(P_A_up{1,ii}(:,:,iiiiiii))^-1 ...
                                )^-1);  

                fnMSE{1,iiiiiii}(1,ii)=trace(((P_r_up{1,ii}(:,:,iiiiiii))^-1+(P_A_up{1,ii}(:,:,iiiiiii))^-1 ...
                                )^-1) - eThrshld;  

            end                      

            for iter=1:nIter   
                
                for iiiiii=1:nCam    
                    for iiiiiii=1:nVessel
                        if x_Cv{1,iiiiii}(i-tHzon,iiiiiii)==1
                            currTrg=iiiiiii;
                            break
                        else
                            currTrg=0;
                        end
                    end
                    
                    if currTrg~=0
                        actIdx=0;
                        for iiiiiii=1:length(nameCamAct)
                            if isequal(nameCamAct{iiiiiii,1},camActIter{iiiiii,iter})
                                actIdx=iiiiiii;
                                break
                            end
                        end
        
                        currSpLv=camActIter{iiiiii,iter}(1,1);
        
                        cPosCam{1,iiiiii}{1,ii-1}(3,1)=sCamMin+(currSpLv-1)*(sCamMax-sCamMin)/(nSpeedLv-1); 
        
                        currHaLv=camActIter{iiiiii,iter}(1,2);
        
                        if currHaLv<=nAngleLv
                            if ii==2
                                cPosCam{1,iiiiii}{1,ii-1}(4,1)=0-(currHaLv-1)*aCamMax/(nAngleLv-1);
                            else
                                cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)-(currHaLv-1)*aCamMax/(nAngleLv-1);
                            end
                        else
                            if ii==2
                                cPosCam{1,iiiiii}{1,ii-1}(4,1)=0+(currHaLv-nAngleLv-1)*aCamMax/(nAngleLv-1);
                            else
                                cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)+(currHaLv-nAngleLv-1)*aCamMax/(nAngleLv-1);
                            end
                        end   
        
                        cPosCam{1,iiiiii}{1,ii}(1,1)=cPosCam{1,iiiiii}{1,ii-1}(1,1)...
                            +t*cPosCam{1,iiiiii}{1,ii-1}(3,1)*cos(cPosCam{1,iiiiii}{1,ii-1}(4,1));
        
                        cPosCam{1,iiiiii}{1,ii}(2,1)=cPosCam{1,iiiiii}{1,ii-1}(2,1)...
                            +t*cPosCam{1,iiiiii}{1,ii-1}(3,1)*sin(cPosCam{1,iiiiii}{1,ii-1}(4,1)); 
        
                        tempPosCam{1,iiiiii}{1,ii}(1,1) = cPosCam{1,iiiiii}{1,ii}(1,1);
        
                        tempPosCam{1,iiiiii}{1,ii}(2,1) = cPosCam{1,iiiiii}{1,ii}(2,1);
        
                        % predict the position of target in t+1    
                        if ii==i-tHzon+1
                        else
                            X_C_up{1,iiiiii}{1,ii-1}(:,:,currTrg) = X_C_pr{1,iiiiii}{1,ii-1}(:,:,currTrg);
                        end
    
                        X_C_pr{1,iiiiii}{1,ii}(:,:,currTrg)=F_v{1,currTrg}*X_C_up{1,iiiiii}{1,ii-1}(:,:,currTrg);
            
                        P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg)=F_v{1,currTrg}*P_C_up{1,iiiiii}{1,ii-1}(:,:,currTrg)*F_v{1,currTrg}'+Q_v{1,currTrg};
        
                        d_Cv=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg))^2 ...
                            +(cPosCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg))^2);                        
        
                        % term 1  
                        if d_Cv>=Dmin&&d_Cv<=Dmax
                            sigma_C=p_C/100*d_Cv;
                            R_C=sigma_C^2.*[1 0;0 1]; 
                            S_C_up=H_C{1,iiiiii}'*(R_C)^-1*H_C{1,iiiiii};
                            P_C_up{1,iiiiii}{1,ii}(:,:,currTrg)=((P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg))^-1+S_C_up)^-1;   
        
                        elseif d_Cv<Dmin
                            P_C_up{1,iiiiii}{1,ii}(:,:,currTrg)=P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg);
        
                        elseif d_Cv>Dmax
                            P_C_up{1,iiiiii}{1,ii}(:,:,currTrg)=P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg);    
                        end  
                       
                        dTerm1=trace(((P_r_up{1,ii}(:,:,currTrg))^-1+(P_A_up{1,ii}(:,:,currTrg))^-1 ...
                        +(P_C_up{1,iiiiii}{1,ii}(:,:,currTrg))^-1)^-1);      

                        fnEst{1,currTrg}(1,ii)=dTerm1;

                        if dTerm1 <= eThrshld
                            dTerm1 = 0;                
                        else
                            dTerm1 = dTerm1-eThrshld;
                        end

                        fnMSE{1,currTrg}(1,ii)=dTerm1;                        
        
                        % term 2
                        dTerm2=0; %obTmV(currTrg,ii-1)/d_Cv;
        
                        % term 3
                        if d_Cv>=Dmin
                            dTerm31=1;                
                        else
                            dTerm31=Dmin-d_Cv;                
                        end                
                        
                        if d_Cv<=Dmax
                            dTerm32=1;                
                        else
                            dTerm32=d_Cv-Dmax;                
                        end
        
                        dTerm3=dTerm31*dTerm32;
        
                        % term 4   
                        dTerm4=0;
                        for iiiiiii=1:nVessel                          

                            d_Cv_2=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiiiii))^2 ...
                                        +(cPosCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiiiii))^2);
            
                            if d_Cv_2<Dsafe
                                dTerm4=dTerm4+(Dsafe-d_Cv_2);    
                            else
                                dTerm4=dTerm4+0;    
                            end
                        end
        
                        % term 5
                        dTerm5=0;
                        for iiiiiii=1:nCam
                            if iiiiiii~=iiiiii
                                dCC=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-cPosCam{1,iiiiiii}{1,ii}(1,1))^2 ...
                                    +(cPosCam{1,iiiiii}{1,ii}(2,1)-cPosCam{1,iiiiiii}{1,ii}(2,1))^2);
                
                                if dCC<Dsafe
                                    dTerm5=dTerm5+(Dsafe-dCC);        
                                else
                                    dTerm5=dTerm5+0;        
                                end 
                            end
                        end    
        
                        % real utility 
                        camUtl{iiiiii,iter}(1,actIdx)=-(alpha_1*dTerm1+alpha_2*dTerm2+alpha_3*dTerm3+alpha_4*dTerm4+alpha_5*dTerm5);        
                       
                        %----------------------------------------------------------
                        posCam{1,iiiiii}{1,ii}(1,1)=cPosCam{1,iiiiii}{1,ii}(1,1);
                        posCam{1,iiiiii}{1,ii}(2,1)=cPosCam{1,iiiiii}{1,ii}(2,1);
                        posCam{1,iiiiii}{1,ii-1}(3,1)=cPosCam{1,iiiiii}{1,ii-1}(3,1);
                        posCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-1}(4,1);
           
                        if posCam{1,iiiiii}{1,ii}(1,1)<X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)==X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=0;    
        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)<X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)<X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=atan(abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg))/abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)))*180/pi;                
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)==X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)<X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=90;  
        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)>X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)<X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=90+atan(abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg))/abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)))*180/pi;      
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)>X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)==X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=180;                
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)>X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)>X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=180+atan(abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg))/abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)))*180/pi;             
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)==X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)>X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=270;                
                        
                        elseif posCam{1,iiiiii}{1,ii}(1,1)<X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg)&&posCam{1,iiiiii}{1,ii}(2,1)>X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)
                            OrtCam{1,iiiiii}(1,ii)=270+atan(abs(posCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg))/abs(posCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg)))*180/pi;             
                        
                        end    
        
                        %----------------------------------------------------------
                        % expected utility
                        for iiiiiii=1:length(nameCamAct)
        
                            for iiiiiiii=1:nVessel
                                if x_Cv{1,iiiiii}(i-tHzon,iiiiiiii)==1
                                    currTrg=iiiiiiii;
                                    break
                                end
                            end
        
                            if iiiiiii~=actIdx
        
                                currSpLv=nameCamAct{iiiiiii,1}(1,1);
        
                                cPosCam{1,iiiiii}{1,ii-1}(3,1)=sCamMin+(currSpLv-1)*(sCamMax-sCamMin)/(nSpeedLv-1); 
                
                                currHaLv=nameCamAct{iiiiiii,1}(1,2);
                
                                if currHaLv<=nAngleLv
                                    if ii==2
                                        cPosCam{1,iiiiii}{1,ii-1}(4,1)=0-(currHaLv-1)*aCamMax/(nAngleLv-1);
                                    else
                                        cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)-(currHaLv-1)*aCamMax/(nAngleLv-1);
                                    end
                                else
                                    if ii==2
                                        cPosCam{1,iiiiii}{1,ii-1}(4,1)=0+(currHaLv-nAngleLv-1)*aCamMax/(nAngleLv-1);
                                    else
                                        cPosCam{1,iiiiii}{1,ii-1}(4,1)=cPosCam{1,iiiiii}{1,ii-2}(4,1)+(currHaLv-nAngleLv-1)*aCamMax/(nAngleLv-1);
                                    end
                                end   
                
                                cPosCam{1,iiiiii}{1,ii}(1,1)=cPosCam{1,iiiiii}{1,ii-1}(1,1)...
                                    +t*cPosCam{1,iiiiii}{1,ii-1}(3,1)*cos(cPosCam{1,iiiiii}{1,ii-1}(4,1));
                
                                cPosCam{1,iiiiii}{1,ii}(2,1)=cPosCam{1,iiiiii}{1,ii-1}(2,1)...
                                    +t*cPosCam{1,iiiiii}{1,ii-1}(3,1)*sin(cPosCam{1,iiiiii}{1,ii-1}(4,1));  
        
                                % predict the position of target in t+1 
                                if ii==i-tHzon+1
                                else
                                    X_C_up{1,iiiiii}{1,ii-1}(:,:,currTrg) = X_C_pr{1,iiiiii}{1,ii-1}(:,:,currTrg);
                                end

                                X_C_pr{1,iiiiii}{1,ii}(:,:,currTrg)=F_v{1,currTrg}*X_C_up{1,iiiiii}{1,ii-1}(:,:,currTrg);
                    
                                P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg)=F_v{1,currTrg}*P_C_up{1,iiiiii}{1,ii-1}(:,:,currTrg)*F_v{1,currTrg}'+Q_v{1,currTrg};
                
                                d_Cv=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,currTrg))^2 ...
                                    +(cPosCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,currTrg))^2); 
                
                                % term 1        
                                if d_Cv>=Dmin&&d_Cv<=Dmax
                                    sigma_C=p_C/100*d_Cv;
                                    R_C=sigma_C^2.*[1 0;0 1]; 
                                    S_C_up=H_C{1,iiiiii}'*(R_C)^-1*H_C{1,iiiiii};
                                    P_C_up{1,iiiiii}{1,ii}(:,:,currTrg)=((P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg))^-1+S_C_up)^-1;   
            
                                elseif d_Cv<Dmin
                                    P_C_up{1,iiiiii}{1,ii}(:,:,currTrg)=P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg);
            
                                elseif d_Cv>Dmax
                                    P_C_up{1,iiiiii}{1,ii}(:,:,currTrg)=P_C_pr{1,iiiiii}{1,ii}(:,:,currTrg);    
                                end  
                                
                                dTerm1 = trace(((P_r_up{1,ii}(:,:,currTrg))^-1+(P_A_up{1,ii}(:,:,currTrg))^-1 ...
                                +(P_C_up{1,iiiiii}{1,ii}(:,:,currTrg))^-1)^-1);    

                                if dTerm1 <= eThrshld
                                    dTerm1 = 0;
                                else
                                    dTerm1 = dTerm1 - eThrshld;
                                end
                
                                % term 2
                                dTerm2=0; %obTmV(currTrg,ii-1)/d_Cv; %=0;
                
                                % term 3
                                if d_Cv>=Dmin
                                    dTerm31=1;                
                                else
                                    dTerm31=Dmin-d_Cv;                
                                end                
                                
                                if d_Cv<=Dmax
                                    dTerm32=1;                
                                else
                                    dTerm32=d_Cv-Dmax;                
                                end
                
                                dTerm3=dTerm31*dTerm32;
                
                                % term 4 
                                dTerm4 = 0;
                                for iiiiiiii=1:nVessel
                                    d_Cv_2=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-X_C_pr{1,iiiiii}{1,ii}(1,1,iiiiiiii))^2 ...
                                        +(cPosCam{1,iiiiii}{1,ii}(2,1)-X_C_pr{1,iiiiii}{1,ii}(3,1,iiiiiiii))^2);  
                    
                                    if d_Cv_2<Dsafe
                                        dTerm4=dTerm4+(Dsafe-d_Cv_2);    
                                    else
                                        dTerm4=dTerm4+0;    
                                    end
                                end
                
                                % term 5
                                dTerm5 = 0;                     
                                for iiiiiiii=1:nCam
                                    if iiiiiiii~=iiiiii
                                        dCC=sqrt((cPosCam{1,iiiiii}{1,ii}(1,1)-cPosCam{1,iiiiiiii}{1,ii}(1,1))^2 ...
                                            +(cPosCam{1,iiiiii}{1,ii}(2,1)-cPosCam{1,iiiiiiii}{1,ii}(2,1))^2);
                        
                                        if dCC<Dsafe
                                            dTerm5=dTerm5+(Dsafe-dCC);        
                                        else
                                            dTerm5=dTerm5+0;        
                                        end 
                                    end
                                end

                                camUtl{iiiiii,iter}(1,iiiiiii) = -(alpha_1*dTerm1+alpha_2*dTerm2+alpha_3*dTerm3+alpha_4*dTerm4+alpha_5*dTerm5);                           
                            
                            end
                        end
                        %----------------------------------------------------------
                        cPosCam{1,iiiiii}{1,ii}(1,1) = tempPosCam{1,iiiiii}{1,ii}(1,1);
        
                        cPosCam{1,iiiiii}{1,ii}(2,1) = tempPosCam{1,iiiiii}{1,ii}(2,1);

                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%% Regret Matching %%%%%%%%%%%%%%%%%%%%%%%%%% 
                for iiiiii = 1:nCam
                    for iiiiiii=1:nVessel
                        if x_Cv{1,iiiiii}(i-tHzon,iiiiiii)==1
                            currTrg=iiiiiii;
                            break
                        else
                            currTrg=0;
                        end
                    end
                    
                    if currTrg~=0        
                        actIdx=0;
                        for iiiiiii=1:length(nameCamAct)
                            if (isequal(camActIter{iiiiii,iter},nameCamAct{iiiiiii,1})) %camAct{iiiiii,ii-1}==nameCamAct{iiiiiii,1}
                                actIdx=iiiiiii;
                                break
                            end
                        end
        
                        mu=0;
                        for iiiiiii= 1:length(nameCamAct)
                            if (~isequal(camActIter{iiiiii,iter},nameCamAct{iiiiiii,1})) %camAct{iiiiii,ii-1}~=nameCamAct{iiiiiii,1}
                                 if iter==1
                                      sumRegret{1,iter}(iiiiii,actIdx,iiiiiii)=...
                                      lambda*0 ...
                                      +(1-lambda)*(camUtl{iiiiii,iter}(1,iiiiiii)-camUtl{iiiiii,iter}(1,actIdx));
                                 else
                                      sumRegret{1,iter}(iiiiii,actIdx,iiiiiii)=...
                                      lambda*sumRegret{1,iter-1}(iiiiii,actIdx,iiiiiii)...
                                      +(1-lambda)*(camUtl{iiiiii,iter}(1,iiiiiii)-camUtl{iiiiii,iter}(1,actIdx));
                                 end
                                
                                 mu=mu+max(sumRegret{1,iter}(iiiiii,actIdx,iiiiiii),0)+epsilon;
        
                            end
                        end
        
                        for iiiiiii= 1:length(nameCamAct)
                            if (~isequal(camActIter{iiiiii,iter},nameCamAct{iiiiiii,1})) %camAct{iiiiii,ii-1}~=nameCamAct{iiiiiii,1}
                                probCamActIter{iiiiii,iter+1}(1,iiiiiii)=(1/mu)*max(sumRegret{1,iter}(iiiiii,actIdx,iiiiiii),0);
                            end
                        end
                       
                        sumProb=0;
                        for iiiiiii=1:length(probCamActIter{iiiiii,iter+1})
                            if iiiiiii~=actIdx
                                sumProb=sumProb+probCamActIter{iiiiii,iter+1}(1,iiiiiii);
                            end
                        end
        
                        probCamActIter{iiiiii,iter+1}(1,actIdx)=1-sumProb;    
                    end    
                end
                %----------------------------------------------------------
                for iii=1:nCam % the number of cameras (players) 
                    for iiiiiii=1:nVessel
                        if x_Cv{1,iii}(i-tHzon,iiiiiii)==1
                            currTrg=iiiiiii;
                            break
                        else
                            currTrg=0;
                        end
                    end

                    if currTrg~=0
                        temp=rand;%tempRef;
                        iiii=1;    
                        iiiii=0; 
                        while (1)
                           if iiii>length(probCamActIter{iii,iter+1})
                               break;
                           end
                           iiiii=iiiii+probCamActIter{iii,iter+1}(1,iiii);
                           if temp<iiiii
                              break; 
                           end
                           iiii=iiii+1;
                        end
                        camActIter{iii,iter+1}=nameCamAct{iiii,1};  % player 'i' selects action 'ii' at iteration 'it'
                    end
                end
            end

            for iii=1:nCam
                probCamAct{iii,ii} = probCamActIter{iii,iter+1};
            end

        end
    end

    tempTime = toc(tStart);

    exeTime(1,i-tHzon) = exeTime(1,i-tHzon)+tempTime;
end

finExeTime = mean(exeTime); %/nCam;

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

figure

for i=1:nVessel
    plot(fnEst{1,i})
    hold on 
end

%--------------------------------------------------------------------------

file_name = sprintf('output_OL_%dcams_%dtargets_%dhorizons_%dTimeStamps_%d.mat',nCam,nVessel,tHzon,nTimeStamp,iTer);

save(file_name,"posCam","obTmV","OrtCam","fnEst","exeTime","finExeTime","finObTmV");








