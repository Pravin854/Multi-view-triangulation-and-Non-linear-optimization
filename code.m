close all;
clc;

load('./MatFilesQues1/cube_imgs.mat'); % load image coordinates
pts2D_view3 = squeeze(image_pts);
pts2D_view3(:,3,:) = ones(56,8)'; % homogenize matrix
pts2D_view3 = permute(pts2D_view3, [2 3 1]); 

load('./MatFilesQues1/projMatrices.mat'); % load projection matrices

% read all projection matrices into one
proj = cat(3, projMatrices{1}, projMatrices{2}, projMatrices{3}, projMatrices{4}, projMatrices{5}, projMatrices{6}, projMatrices{7}, projMatrices{8});
% disp(proj);
out = [];

for i = 1: 56

t = [];
p = [];

for j = 1: 8
    
    X1 = (pts2D_view3(:,:,j));
    tp = X1(3,:);
    X1_hom(1,:) = X1(1,:) ./ tp; % homogenize
    X1_hom(2,:) = X1(2,:) ./ tp; % homogenize
    X1 = (X1_hom(1, i));
    Y1 = (X1_hom(2, i));
    Pro1 = proj(: , : , j) ; % projection matrix of single image
    
    for k = 1:3
        tp1(k) = Pro1(3,k);
    end
    
    s1 = [X1*Pro1(3,4) - Pro1(1 , 4); Y1*Pro1(3,4) - Pro1(2 , 4)];
    
    t1 = [Pro1(1,1) - X1*tp1(1), Pro1(1,2) - X1*tp1(2), Pro1(1,3) - X1*tp1(3);
          Pro1(2,1) - Y1*tp1(1), Pro1(2,2) - Y1*tp1(2), Pro1(2,3) - Y1*tp1(3)];
        
    p = [p ; s1];  
    t = [t ; t1];

end

    out = [out ; (t\p)'];
        
end

Points = [out' ; ones(1, 56)]; % points after triangulation

figure;
plot3(Points(1, :),Points(2, :),Points(3, :),'*');
title('Cube after triangulation');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
legend('Cube point');