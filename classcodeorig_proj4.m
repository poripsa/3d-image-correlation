clear;
% close all;
clc;

%[pts, tris, norms] = objread('box.obj');
%resolution = round(0.05*[2048,1536]);      %pixels of sensor
resolution = [100,100];
pixel_size = [0.05,0.05]*10^(-6);     %meters/pixel size of each pixel
center_pt = 0.5*resolution;         %pixels center of image
f = 5.25*10^(-6);
%f=0.03*6*10^(-3);                        %meter focal length
%numPts=size(pts,1);                 %points to projects
camera_positions = [0.25 0 0; -0.25 0 0];
camera_y_angle = [3 ,-3];
numCameras = size(camera_positions,1);

%pts(:,3) = pts (:,3)+15;
I=zeros(resolution(2), resolution(1));
%define Sphere
%sphere_depth = 1;
%sphere_radii = [1 1];
% spherechar = double(uint8('s'));
% planechar = double(uint8('p'));
% sphere_centers =[[spherechar 000 000 sphere_depth]' [spherechar 400 300 sphere_depth]'];
ka =1;
kd=1;
light_pos=[0 0 1]';
p2=[0 0 0]';
sphere_centers=[1 0 3; -1 0 3; 0 1 3; 0 -1 3];
%sphere_centers
planes = [-1 0 0 1.25; 1 0 0 1.25; 0 -1 0 1.25; 0 1 0 1.25; 0 0 1 1.25; 0 0 -1 4.25];
sphere_radii=[0.5 0.5 0.5 0.5];
%sphere_centers = [0 0 2];
light_positions = [0 1 1; 1 -1 1];
light_power = [50,20];
ambient_light_power = 10;

numPixels = prod(resolution);
camera1_intr_points = zeros(numPixels,3);
camera2_intr_points = zeros(numPixels,3);
camera_intr_points = zeros(numPixels,3);
for cameraIdx = 1:numCameras,
    cam_translation = camera_positions(cameraIdx,:);
    cam_angle_rad = camera_y_angle(cameraIdx)*pi/180;
    cam_rotation = [cos(cam_angle_rad) 0 -sin(cam_angle_rad);0 1 0; sin(cam_angle_rad) 0 cos(cam_angle_rad)];
    %     if (cameraIdx ==1)
    %         %camera_intr_points = camera1_intr_points;
    %         camera1_intr_points = camera_intr_points;
    %     else
    %         %camera_intr_points = camera2_intr_points;
    %         camera2_intr_points = camera_intr_points;
    %     end
    %camera_int_points =
    for x_im=1:resolution(1)
        for y_im=1:resolution(2)
            %for a 3D ray emanating from pixel (x_im, y_im) into 3D space
            p1=[-(x_im-center_pt(1))*pixel_size(1)...
                (y_im-center_pt(2))*pixel_size(2)...
                -f]';
            p2=[0 0 0]';
            p1 = cam_rotation*p1+cam_translation';
            p2 = cam_rotation*p2+cam_translation';
            closest_iPt = [Inf,Inf,Inf];
            closest_norm = [Inf,Inf,Inf];
            for sphereIdx=1:size(sphere_centers,1),
                sphere_ctr =sphere_centers(sphereIdx,:);
                sphere_radius=sphere_radii(sphereIdx);
                [iPt1,iPt2]=LineIntersectSphere(p1,p2,sphere_ctr',sphere_radius);
                if(isfinite(iPt1)==1 & isfinite(iPt2)==1)
                    %if (iPt1(3)>0 && iPt2(3)>0)
                    if (iPt1(3) < iPt2(3))
                        iPt = iPt1;
                    else
                        iPt = iPt2;
                    end
                    %end
                    if(iPt(3)<closest_iPt(3))
                        closest_iPt = iPt';
                        normalVec =(iPt'-sphere_ctr)/norm(iPt'-sphere_ctr);
                        closest_norm = normalVec;
                    end
                end
            end
            for planeIdx= 1:size(planes,1),
                plane_eq=planes(planeIdx,:);
                V=p2-p1;
                p0=p2;
                normalVec=plane_eq(1:3);
                t=-(p0'*normalVec'+plane_eq(4))/(V'*normalVec');
                iPt=p0+t*V;
                if(iPt(3)<closest_iPt(3) & iPt(3)>0)
                    closest_iPt=iPt';
                    closest_norm = normalVec;
                end
            end
            
            if (isfinite(closest_iPt(3)))
                %strore intersection points
                xy_2_idx = (y_im-1)*resolution(1)+x_im;%xy loc to unique index where 1st val =1 and last is res
                camera_intr_points(xy_2_idx,:)=closest_iPt;
                I(y_im,x_im)=ka*ambient_light_power;
                
                for lightIdx = 1:size(light_positions,1),
                    light_pos = light_positions(lightIdx,:);
                    this_light_power = light_power(lightIdx);
                    %normalVec=(iPt-sphere_ctr)/norm(iPt-sphere_ctr);
                    lightVec=(light_pos-closest_iPt)/norm(light_pos-closest_iPt);
                    
                    I(y_im,x_im)=kd*this_light_power*lightVec*closest_norm';
                end
            end
        end
    end
    if (cameraIdx ==1)
        %camera_intr_points = camera1_intr_points;
        camera1_intr_points = camera_intr_points;
        label_string = sprintf('Right Camera');
    else
        %camera_intr_points = camera2_intr_points;
        camera2_intr_points = camera_intr_points;
        label_string = sprintf('Left Camera');
    end
    subplot(1,numCameras, cameraIdx),imshow(I,[]),xlabel(label_string);
    
end
%,cameraIdximshow(I,[])
%figure,imshow(I,[])
camera1_to_camera2_correspondence= zeros(resolution(2), resolution(1));
%compute correspondence
for cam1_y_im=1:resolution(2)
    for cam1_x_im=1:resolution(1)
        
        %retreive the 3d surface loaction for camera1 pixel (y_im,x_im)
        cam1_xy_2_idx = (cam1_y_im-1)*resolution(1)+cam1_x_im;%xy loc to unique index where 1st val =1 and last is res
        %cam2_y_im = floor((cam2_xy_2_idx -1)/resolution(1))+1;
        cam1_3d_pt = camera1_intr_points(cam1_xy_2_idx,:);
        %cam1_3d_pt = camera_intr_points(cam1_xy_2_idx,:);
        closest_pt_idx = -1;
        closest_pt = [Inf Inf Inf];
        closest_dist3d = Inf;
        
        for cam2_y_im=cam1_y_im:cam1_y_im %1:
            %for cam2_y_im=1:resolution(2)
            %cam2_y_im = cam1_y_im;
            for cam2_x_im=1:resolution(1)
                
                %retreive the 3d surface loaction for camera1 pixel (y_im,x_im)
                cam2_xy_2_idx = (cam2_y_im-1)*resolution(1)+cam2_x_im;%xy loc to unique index where 1st val =1 and last is res
                %cam2_xy_2_idx
                %cam2_y_im2 = floor((cam2_xy_2_idx -1)/resolution(1))+1;
                
                cam2_3d_pt = camera2_intr_points(cam2_xy_2_idx,:);%camera_intr_points
                %cam2_3d_pt = camera_intr_points(cam2_xy_2_idx,:);
                dist3d = norm(cam1_3d_pt -cam2_3d_pt);
                %is this correspondence the best
                if (dist3d < closest_dist3d)
                    closest_pt_Idx = cam2_xy_2_idx;
                    %closest_pt_Idx
                    closest_pt = cam2_3d_pt;
                    closest_dist3d = dist3d;
                end
                
            end
        end
        if (isfinite(closest_dist3d))
            camera1_to_camera2_correspondence(cam1_y_im,cam1_x_im)=closest_pt_Idx;
            %closest_pt_Idx
        end
    end
end
reconstructed_points = zeros(numPixels,3);
for cam1_y_im=1:resolution(2)
    for cam1_x_im=1:resolution(1)
        
        cam2_xy_2_idx= camera1_to_camera2_correspondence(cam1_y_im,cam1_x_im);
        cam2_y_im = floor((cam2_xy_2_idx -1)/resolution(1))+1;
        %cam2_x_im = floor((cam2_xy_2_idx)/resolution(2))+1;
        
        cam2_x_im = cam2_xy_2_idx-(cam2_y_im-1)*resolution(1);
        %cam2_y_im = cam2_xy_2_idx - (cam2_x_im-1)*resolution(2);
        %cam2_x_im = cam2_xy_2_idx
        %cam2_x_im = floor((cam2_xy_2_idx -cam2_y_im)/resolution(2))+1;
%         fprintf(1,'camera1 pixel (%d,%d)',...
%             cam1_x_im,cam1_y_im);
%         fprintf(1, 'corresponds to camera 2 pixel (%d,%d).\n',...
%             cam2_x_im,cam2_y_im);
        
        %compute v_r,v_l,
        
        cameraIdx=1;
        cam_translation = camera_positions(cameraIdx,:);
        cam_angle_rad = camera_y_angle(cameraIdx)*pi/180;
        cam_rotation = [cos(cam_angle_rad) 0 -sin(cam_angle_rad);0 1 0; sin(cam_angle_rad) 0 cos(cam_angle_rad)];
        %
        p1=[-(cam1_x_im-center_pt(1))*pixel_size(1)...
            (cam1_y_im-center_pt(2))*pixel_size(2)...
            -f]';
        p2=[0 0 0]';
        cam1_p1 = cam_rotation*p1+cam_translation';
        cam1_p2 = cam_rotation*p2+cam_translation';
        
        cameraIdx =2;
        cam_translation = camera_positions(cameraIdx,:);
        cam_angle_rad = camera_y_angle(cameraIdx)*pi/180;
        cam_rotation = [cos(cam_angle_rad) 0 -sin(cam_angle_rad);0 1 0; sin(cam_angle_rad) 0 cos(cam_angle_rad)];
        %
        p1=[-(cam2_x_im-center_pt(1))*pixel_size(1)...
            (cam2_y_im-center_pt(2))*pixel_size(2)...
            -f]';
        p2=[0 0 0]'; 
        cam2_p1 = cam_rotation*p1+cam_translation';
        cam2_p2 = cam_rotation*p2+cam_translation';
        
        vr = cam1_p2-cam1_p1;
        pr = cam1_p2;
        vl = cam2_p2-cam2_p1;
        pl = cam2_p2;
        
        lamda_vec = inv([vr -vl]'*[vr -vl])*[vr -vl]'*[pr-pl];
        p_reconstructed = 0.5*(lamda_vec(1)*vr+pr+lamda_vec(2)*vl+pl);
        aaa=1;
        cam1_xy_2_idx = (cam1_y_im-1)*resolution(1)+cam1_x_im;
        reconstructed_points(cam1_xy_2_idx,:) = p_reconstructed';
    end
end

figure, subplot(3,1,1), plot3(reconstructed_points(:,1),reconstructed_points(:,2),reconstructed_points(:,3),'b.','MarkerSize',5),xlabel('(1)Looking down from positive Z');
subplot(3,1,2), plot3(reconstructed_points(:,1),reconstructed_points(:,2),reconstructed_points(:,3),'b.','MarkerSize',5),xlabel('(2)Looking down from positive Y');
subplot(3,1,3), plot3(reconstructed_points(:,1),reconstructed_points(:,2),reconstructed_points(:,3),'b.','MarkerSize',5),ylabel('(3)Looking down on positive x-axis');

%figure, plot3(reconstructed_points(:,1),reconstructed_points(:,2),reconstructed_points(:,3),'b.','MarkerSize',20),xlabel('(1)Looking down from positive z');%,ylabel('(3)Looking down on positive x-axis');
axis equal;
ax = gca;
ax.Clipping = 'off';
%bundle adjustment wiki








