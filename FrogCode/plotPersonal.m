function plots = plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
%points or lines must be organized in columns to properly plot using this
%function
%x and yError must be the same length as the number of points if each point
%has a differrnt error, otherwise it may stay 1 size if it will apply to
%all
hold on
%plot each point without transparency
f1 = scatter(x,y,MarkerSize,RGBa(:,1:3),'Filled');
if RGBa(:,4) == 0.8
    f1.MarkerFaceAlpha = 0.8;
end
%plot the line connecting data
plot(x,y,'Color',RGBa,'LineWidth',LineWidth);
%plot the x error bars, individual points must be in seperate columns
plot([(x - xError);(x + xError)],[y;y],'Color',RGBa,'LineWidth',LineWidth);
%plot left "caps"
scaledYCapSize = range(ylimRange)*CapSize*aspectRatio;
plot([(x - xError);(x - xError)],[y-scaledYCapSize;y+scaledYCapSize],'Color',RGBa,'LineWidth',LineWidth);
%plot right "caps"
plot([(x + xError);(x + xError)],[y-scaledYCapSize;y+scaledYCapSize],'Color',RGBa,'LineWidth',LineWidth);
%plot the y error bars
plot([x;x],[(y - yError);(y + yError)],'Color',RGBa,'LineWidth',LineWidth);
%plot lower "caps"
scaledXCapSize = range(xlimRange)*CapSize;
plot([x-scaledXCapSize;x+scaledXCapSize],[(y - yError);(y - yError)],'Color',RGBa,'LineWidth',LineWidth);
%plot upper "caps"
plot([x-scaledXCapSize;x+scaledXCapSize],[(y + yError);(y + yError)],'Color',RGBa,'LineWidth',LineWidth);
end
