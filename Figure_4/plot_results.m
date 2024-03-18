clear all

clc

nVessel = 3;

nCam = 2;

tHzon = 1;

nTimeStamp = 10000;

iTer = 0;

R = 16e3;

file_name = sprintf('output_OL_%dcams_%dtargets_%dhorizons_%dTimeStamps_%d.mat',nCam,nVessel,tHzon,nTimeStamp,iTer);

load(file_name,"posCam","obTmV","OrtCam","fnEst","exeTime","finExeTime","finObTmV");

file_name = sprintf('scenario_%dcams_%dtargets_%dhorizons_%dTimeStamps_%d.mat',nCam,nVessel,tHzon,nTimeStamp,iTer);

load(file_name,"posV","iPosV");

%-----------------------------

bep = birdsEyePlot('XLim',[-50000 50000],'YLim',[-50000 50000]);

caPlotter(1,1)=coverageAreaPlotter(bep,'FaceColor','red','EdgeColor','red');

caPlotter(1,2)=coverageAreaPlotter(bep,'FaceColor','red','EdgeColor','red');

fieldOfView = 3;

camroll(-90);

for i=1:nTimeStamp

    if i>1
        if (posCam{1,1}{1,i}(1,1)==posCam{1,1}{1,i-1}(1,1)&&posCam{1,1}{1,i}(2,1)==posCam{1,1}{1,i-1}(2,1))...
          &&(posCam{1,2}{1,i}(1,1)==posCam{1,2}{1,i-1}(1,1)&&posCam{1,2}{1,i}(2,1)==posCam{1,2}{1,i-1}(2,1))
           break
        end

        if posCam{1,1}{1,i}(1,1)==posCam{1,1}{1,i-1}(1,1)&&posCam{1,1}{1,i}(2,1)==posCam{1,1}{1,i-1}(2,1)
        else
            plotCoverageArea(caPlotter(1,1),[posCam{1,1}{1,i}(1,1) posCam{1,1}{1,i}(2,1)],R,OrtCam{1,1}(1,i),fieldOfView);
        
            hold on
        
            plot(posCam{1,1}{1,i}(1,1),posCam{1,1}{1,i}(2,1),'>r','MarkerSize',2)    
        end

        if posCam{1,2}{1,i}(1,1)==posCam{1,2}{1,i-1}(1,1)&&posCam{1,2}{1,i}(2,1)==posCam{1,2}{1,i-1}(2,1)
        else
            plotCoverageArea(caPlotter(1,2),[posCam{1,2}{1,i}(1,1) posCam{1,2}{1,i}(2,1)],R,OrtCam{1,2}(1,i),fieldOfView);

            hold on

            plot(posCam{1,2}{1,i}(1,1),posCam{1,2}{1,i}(2,1),'>r','MarkerSize',2)
        end

    else
        plotCoverageArea(caPlotter(1,1),[posCam{1,1}{1,i}(1,1) posCam{1,1}{1,i}(2,1)],R,OrtCam{1,1}(1,i),fieldOfView);
    
        hold on
    
        plot(posCam{1,1}{1,i}(1,1),posCam{1,1}{1,i}(2,1),'>r','MarkerSize',2)

        plotCoverageArea(caPlotter(1,2),[posCam{1,2}{1,i}(1,1) posCam{1,2}{1,i}(2,1)],R,OrtCam{1,2}(1,i),fieldOfView);

        plot(posCam{1,2}{1,i}(1,1),posCam{1,2}{1,i}(2,1),'>r','MarkerSize',2)
    end

    plot(posV{1,1}{1,i}(1,1),posV{1,1}{1,i}(3,1),'bo','MarkerSize',2)

    plot(posV{1,2}{1,i}(1,1),posV{1,2}{1,i}(3,1),'bo','MarkerSize',2)

    plot(posV{1,3}{1,i}(1,1),posV{1,3}{1,i}(3,1),'bo','MarkerSize',2)

    % pause(0.5)

end


    
