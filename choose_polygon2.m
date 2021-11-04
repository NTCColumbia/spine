function [roi,Xsp,Ysp] = choose_polygon2(pixelsx, pixelsy)
% choose polygon from the current figure
% The function returns a vector of all the pixels indices in a chosen polygon

gcf;

count=1;
[Xsp,Ysp,cc] = ginput(1);
hold on
plot_sp = plot(Xsp(1),Ysp(1),'r','linewidth',1);

while cc~=3
    count=count+1;
    [Xsp(count),Ysp(count),cc]=ginput(1);
    set(plot_sp,'xdata',Xsp,'ydata',Ysp);
end

Xind = ones(pixelsy,1)*(1:pixelsx);
Yind = (1:pixelsy)'*ones(1,pixelsx);
in_pol = inpolygon(Xind,Yind,Xsp,Ysp);
% in_pol = in_pol';
vec_2D = double(in_pol);
roi = find(vec_2D==1);
% roi_2D = findn(vec_2D==1);

% vec = reshape(vec_2D,pixels^2,1);
% roi = find(vec==1);

%% Check
% mat = zeros(100);
% mat(roi) = 1;
% mat = reshape(mat,10000,1);
% figure;mimg(mat,100,100);

