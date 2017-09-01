clear;
%close all;
clc;

[pts, tris, norms] = objread('p51_mustang.obj');
%[pts, tris, norms] = objread('box.obj');

numPts=size(pts,1);                 %points to projects

% v1=rand(3,1);p51_mustang
% v2=rand(3,1);
% v1=v1/norm(v1);
% v2=v2/norm(v2);
%
% v2=v2-(v2'*v1)*v1; %makes v1 and v2 perpendicular
% v2= v2/norm(v2);
% v3 = cross(v1,v2);
meanVal = mean(pts);
pts_centered = pts-ones(numPts,1)*meanVal; %centered the points about origin
theta = pi/4;
Ry = [cos(theta) 0 sin(theta); 0 1 0;-sin(theta) 0 cos(theta)];
Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R = Ry*Rx;
%R = [v1 v2 v3];
pts = (R*pts_centered')';

pts = pts+ones(numPts,1)*meanVal;
aabb = [min(pts);max(pts)];

resolution = [100,100];      %pixels of sensor %make0.5-0.15
pixel_size = [9.8,9.8]*10^(-6);     %meters/pixel size of each pixel
center_pt = 0.5*resolution;         %pixels center of image
f=10.5*10^(-3);                        %meter focal length

% resolution = round(0.01*[2048, 1536]);      %pixels of sensor %make0.5-0.15
% pixel_size = [3.2,3.2]*10^(-6);     %meters/pixel size of each pixel
% center_pt = 0.5*resolution;         %pixels center of image
% f=.05*6*10^(-3);

pts(:,3) = pts (:,3)+55;
I=zeros(resolution(2), resolution(1));
% for ptIdx=1:numPts
%     P=pts(ptIdx,:);
%     x_im = round((-f*P(1)/(pixel_size(1)*P(3)))+center_pt(1));
%     y_im = round((-f*P(2)/(pixel_size(2)*P(3)))+center_pt(2));
%     if(x_im > 0 && x_im < resolution(1) && ...
%             y_im > 0 && y_im < resolution(2))
%         I(y_im, x_im)=1;
%     end
% end
%
% imshow(I,[])
ka = 15;
kd = 200;
ks =10;
alpha = 20;
light_pos=[0 0 1]';

numtris = size(tris,1);

for x_im=1:resolution(1)
    %fprintf(1,'processing %f done...\n',(x_im));
    for y_im=1:resolution(2)
        p1=[(x_im-center_pt(1))*pixel_size(1)...
            (y_im-center_pt(2))*pixel_size(2)...
            -f]';
        p2=[0 0 0]';
        closestPt = [Inf, Inf, Inf];
        if(p1(2)>aabb(3) && p1(2)<aabb(4) && p1(1)>aabb(1) && p1(1)<aabb(2));
            for triIdx = 1: numtris,
                if(norms(triIdx,1)<0 || norms(triIdx,3)>0 || norms(triIdx,2)>0)
                    tripts = pts(tris(triIdx,:),:);
                    [hasIntersection, iPt,n] = IntersectRayPolygon(p1,p2,tripts(1,:)',tripts(2,:)',tripts(3,:)',0);
                    if(hasIntersection)
                        if(iPt(3)<closestPt(3) && iPt(3)>0)
                            closestPt=iPt;
                            triNormal = n/norm(n);
                        end
                    end
                else continue;
                end
            end
        end
        if(isfinite(closestPt(3)))
            
            lightVec = (light_pos-closestPt)/norm(light_pos - closestPt);
            %triNormal
            I(y_im , x_im) = ka + kd*lightVec'*triNormal;
            %I(y_im,x_im)=1;
        end
    end
end

figure, imshow(I,[])
