%Destination adress to output the analyzed data
adress = 'C:\Users\mrash\Desktop\EBSD\Scan Rotation\';

%EBSD data main location
pname = 'C:\Users\mrash\OneDrive - Illinois Institute of Technology\LPBF\IN718\Project 2 - Scan rotation\EBSD data\';

%Samples file location and samples names
samples = {'/16/16.crc', '/17/17.crc', '/18/18.crc', '/19/19.crc','/20/20.crc', '/21/21.crc', '/22/22.crc', '/23/23.crc', '/24/24.crc', '/25/25.crc', '/26/26.crc', '/27/27.crc', '/28/28.crc', '/29/29.crc'};
sampleName = {'16','17','18','19','20','21','22','23','24','25','26','27','28','29'};

%Analyzing and extracting different materials properties and store them in
%lists
TP = [];
TwinPercent = [];
SchmidF = [];
grainSize = [];
TI = [];
ShapeF = [];
mergedGrainSize = [];
diameters = table;

%Change the material axis for IPF coloring reference
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

%iterating between samples to perform analysis
for num = 1:length(samples)
    sample = sampleName{num};
    fname = [pname samples{num}];
    
    %loading the EBSD data
    ebsd = EBSD.load(fname,CS,'interface','crc',...
        'convertEuler2SpatialReferenceFrame');
    
    %Extracting the grains data from the EBSD data
    [grains,ebsd.grainId] = calcGrains(ebsd('indexed'));
    ebsd(grains(grains.grainSize <= 20)) = [];
    %ebsd(grains(grains.equivalentRadius <= 2*min(grains.equivalentRadius))) = [];
    [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',10*degree,'removeQuadruplePoints');
    CS = grains.CS;
    gB = grains.boundary;
    gB_NiNi = gB('Ni-superalloy','Ni-superalloy');
    
    %Extracting and plotting twin grain boundaries' misorientation angle
    close all
    histogram(gB_NiNi.misorientation.angle./degree)
    ylim([0,20000]);
    xlabel('misorientation angle (degree)')
    fullFileName = fullfile(adress,[sample 'MisHist.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %Change the coloring reference of the IPF maps
    oM = ipdfHSVOrientationMapping(ebsd);
    %define the direction of the ipf
    oM.inversePoleFigureDirection = yvector;
    %convert the ebsd map orientations to a color based on the IPF
    color = oM.orientation2color(ebsd.orientations);
    
    %Extracting and plotting different types of twin grain bouandaries and 
    %triple points over band contrast image
    gB3 = gB_NiNi(angle(gB_NiNi.misorientation,CSL(3,ebsd.CS)) < 5*degree);
    isCSL3 = grains.boundary.isTwinning(CSL(3,ebsd.CS),5*degree);
    isCSL9 = grains.boundary.isTwinning(CSL(9,ebsd.CS),5*degree);
    isCSL27 = grains.boundary.isTwinning(CSL(27,ebsd.CS),5*degree);
    tPid = sum(isCSL3(grains.triplePoints.boundaryId),2)>=2;

    plot(ebsd,ebsd.bc,'micronbar', 'off')
    colormap gray
    hold on
    plot(gB, 'lineWidth', 2, 'lineColor', 'black')
    hold on
    plot(gB3,'lineColor','red','linewidth',3)
    hold off
    
    fullFileName = fullfile(adress,[sample 'Merged.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    

    % denoise the orientation map
    ebsdS = fill(ebsd('indexed'));
    
    %calculate and plot the ODF maps with different axis directions
    psi = calcKernel(grains.meanOrientation);
    odf = calcDensity(ebsdS.orientations,'kernel',psi);
    h = Miller(1,1,0,odf.CS);
    plotPDF(odf,h,'antipodal','silent');
    caxis([0,12]);
    fullFileName = fullfile(adress,[sample 'PF-110.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    close all
    h = Miller(1,1,1,odf.CS);
    plotPDF(odf,h,'antipodal','silent');
    caxis([0,12]);
    fullFileName = fullfile(adress,[sample 'PF-111.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %calculate and plot the KAM map
    ebsd = ebsd.gridify;
    kam = ebsd.KAM / degree;
    F = halfQuadraticFilter;
    F.alpha = 0.5;
    plot(ebsdS,ebsdS.KAM('threshold',2.5*degree) ./ degree,'micronbar','off')
    caxis([0,3])
    mtexColorbar
    mtexColorMap jet
    hold on
    plot(gB,'lineWidth',1.5)
    hold off
    fullFileName = fullfile(adress,[sample 'KAM.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %plotting the triple point junctions over the band contrast image
    plot(ebsd,ebsd.bc,'micronbar', 'off')
    colormap gray
    hold on
    plot(gB, 'lineWidth', 2, 'lineColor', 'black')
    hold on
    plot(gB3,'lineColor','white','linewidth',3)
    hold on
    plot(grains.triplePoints(tPid),'color','red','linewidth',2,'MarkerSize',8)
    hold off
    fullFileName = fullfile(adress,[sample 'TP.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %plotting the IPF maps
    oM.inversePoleFigureDirection = yvector;
    color = oM.orientation2color(ebsd.orientations);
    plot(ebsd, color,'micronbar', 'off');
    fullFileName = fullfile(adress,[sample 'Grains.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %calculate and plot misorientation angle distribution
    odf = calcDensity(ebsd.orientations);
    mdf = calcMDF(odf);
    close all
    plotAngleDistribution(mdf)
    hold all
    plotAngleDistribution(ebsd.CS,ebsd.CS)
    hold off
    legend('model ODF','uniform ODF')
    fullFileName = fullfile(adress,[sample 'MisAngle.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    odf = calcDensity(ebsd.orientations);
    mdf = calcMDF(odf);
    plotAxisDistribution(gB.misorientation, 'smooth');
    fullFileName = fullfile(adress,[sample 'MisAxis1.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    ind = gB.misorientation.angle>59*degree & gB.misorientation.angle<61*degree;
    mori = gB.misorientation(ind);
    scatter(mori);
    fullFileName = fullfile(adress,[sample 'MisAxis3D.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %Calculate Schmidt factor for each grain and plot it
    sS = slipSystem.fcc(ebsd.CS);
    sS = sS.symmetrise;
    sSLocal = grains.meanOrientation * sS;
    sigma = stressTensor.uniaxial(vector3d.Y);
    SF = sSLocal.SchmidFactor(sigma);

    % take the maxium allong the rows
    [SFMax,active] = max(SF,[],2);
    SchmidF(end+1) = mean(abs(SFMax));
    
    % plot the maximum Schmid factor
    plot(grains,SFMax,'micronbar','off','linewidth',2);
    mtexColorbar
    mtexColorMap jet
    caxis([0.35,0.5])
    fullFileName = fullfile(adress,[sample 'SF.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %Calculate and plot the Bingham distribution
    ebsd = ebsd('indexed').gridify;
    kappa = ebsd.curvature;
    alpha = kappa.dislocationDensity;
    dS = dislocationSystem.fcc(ebsd.CS);
    a = norm(ebsd.CS.aAxis);
    nu = 0.3;
    dS(dS.isEdge).u = 1;
    dS(dS.isScrew).u = 1 - 0.3;
    dSRot = ebsd.orientations * dS;
    [rho,factor] = fitDislocationSystems(kappa,dSRot);
    alpha = sum(dSRot.tensor .* rho,2);
    alpha.opt.unit = '1/um';
    kappa = alpha.curvature;
    newMtexFigure('nrows',3,'ncols',3);
    for i = 1:3
        for j = 1:3

        nextAxis(i,j)
        plot(ebsd,kappa{i,j},'micronBar','off')
        hold on; plot(grains.boundary,'linewidth',2); hold off

        end
    end
    setColorRange([-0.005,0.005])
    drawNow(gcm,'figSize','large');
    close all
    plot(ebsd,factor*sum(abs(rho .* dSRot.u),2),'micronbar','off')
    mtexColorMap('hot')
    mtexColorbar
    set(gca,'ColorScale','log');
    set(gca,'CLim',[1e12 1e16]);
    hold on
    plot(grains.boundary,'linewidth',2)
    hold off
    fullFileName = fullfile(adress,[sample 'DislocationDensity.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %calculate the number of triple points, the twin grain boundaries'
    %percentage, average grain size, texture index, and grain shape factor
    grainSize(end+1) = mean((2*grains.equivalentRadius));
    ShapeF(end+1) = mean((grains.shapeFactor));
    TP(end+1) = length(grains.triplePoints(tPid));
    TwinPercent(end+1) = (sum(isCSL3)+sum(isCSL9)+sum(isCSL27))/length(isCSL3);
    TI(end+1) = textureindex(odf);
    diameter = 2*grains(grains.equivalentRadius>2.5).equivalentRadius;
    close all
    h1 = histogram(diameter);
    h1.BinWidth = 5;
    ylim([0,200]);
    xlim([0,200]);
    xlabel('Equivalent Diameter (microns^2)')
    fullFileName = fullfile(adress,[sample 'GrainSizeDist.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    odf = calcDensity(ebsd.orientations);
    
    %merging twin grain boundaries with parent grains and calculating the
    %parent grain sizes, grain orientation spread, 
    CS = grains.CS;
    gB = grains.boundary;
    gB_NiNi = gB('Ni-superalloy','Ni-superalloy');
    gB3 = gB_NiNi(angle(gB_NiNi.misorientation,CSL(3,ebsd.CS)) < 3*degree);
    [mergedGrains,parentId] = merge(grains,gB3);
    mergedGrainSize(end+1) = mean((2*mergedGrains.equivalentRadius));
    mergedGrains.prop.GOS = accumarray(parentId,grains.GOS,size(mergedGrains),@nanmean);
    plot(mergedGrains,mergedGrains.GOS  ./ degree,'micronbar', 'off');
    mtexColorMap jet
    caxis([0,5])
    fullFileName = fullfile(adress,[sample 'GOS-merged.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    plot(grains,grains.GOS  ./ degree,'micronbar', 'off');
    mtexColorMap jet
    caxis([0,5])
    fullFileName = fullfile(adress,[sample 'GOS.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    
    %outputing the grain size data
    fileName = fullfile(adress,[sample 'GrainSize.xlsx']);
    writematrix(diameter,fileName);
    
    %calculating and plotting different ODF sections
    odf = calcDensity(ebsd.orientations);
    plotSection(odf, 'phi2', 89*degree, 'micronbar', 'off');
    mtexColorbar
    caxis([0,15])
    plot3d(odf);
    xlim([0,360]);
    ylim([0,90]);
    zlim([0,90]);
    xticks([0,90,180,270,360]);
    yticks([0,45,90]);
    zticks([0,45,90]);
    mtexColorMap('jet')
    
    close all
    odf = calcDensity(ebsdS.orientations);
    odf.SS = specimenSymmetry('222');
    plot(odf,'sections',18,'layout',[5 4],...
        'coordinates','off','xlabel','','ylabel','')
    plot3d(odf);
    mtexColorMap('jet')
    fullFileName = fullfile(adress,[sample 'ODF-3D.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
    plotSection(odf, 'phi2', 0*degree, 'micronbar', 'off');
    caxis([0,15])
    fullFileName = fullfile(adress,[sample 'ODF-0.tiff']);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300);
end
