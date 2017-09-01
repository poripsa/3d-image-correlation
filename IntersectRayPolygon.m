function [intersect, ip,n] = IntersectRayPolygon( rp0, rp1, tp1, tp2, tp3, backface_culling)  
  %
  % This function intersects a line passing through the 3D points (rp0,rp1) with
  % the triangle formed by the three sequential 3D points (tp1,tp2,tp3). Note that
  % the implied triangular surface is defined by the 3 line segments ending with
  % endpoints (tp1,tp2), (tp2,tp3), and (tp3,tp1). The order of the points defines
  % orientation of the polynomial in terms of the "inward" pointing surface normal.
  %
  % (rp0, rp1)       - the 3D point pair through which the ray passes.
  %
  % (tp1,tp2,tp3)    - the 3D point ordered triplet specifying the vertices of the triangle
  %                    of the triangle to be intersected by the ray.
  %
  % backface_culling - this variable effects what is considered to be a valid polygon-ray
  %                    intersection. Specifically, when backface culling is set to 1 only 
  %                    triangles with normals pointing opposite the ray rp1-rp0 are 
  %                    considered for intersection. When backface culling is set to 0 
  %                    all triangles are considered for intersection.
  % 
  % 
  if (exist('backface_culling')==0)
    backface_culling = 1;
  end
  intersect = 0;
  ip = [inf; inf; inf];
  t = [tp1; tp2; tp3];
  e1 = (tp2-tp1);
  e2 = (tp3-tp1);
  n = cross(e1,e2);
  denominator = ((rp1-rp0)'*n);
  % backface culling
  if (backface_culling == 0)
    if (abs(denominator) < 1E-10)
      return
    end
  elseif (denominator < 1E-10)
    return
  end
  % no backface culling
  tparam = ((tp1-rp0)'*n)/denominator;
  planePoint = rp0 + tparam*(rp1-rp0);

  wVec = planePoint - tp1;
  vVec = e1;
  uVec = e2;
  ulengthsq = uVec'*uVec;
  vlengthsq = vVec'*vVec;
  uProjv = uVec'*vVec;
  uProjw = uVec'*wVec;
  vProjw = vVec'*wVec;
	denominator = (uProjv^2-ulengthsq*vlengthsq);
  sparam = (uProjv*vProjw - vlengthsq*uProjw)/denominator;
  tparam = (uProjv*uProjw - ulengthsq*vProjw)/denominator;
  if (sparam >= 0 && sparam <= 1 && tparam >= 0 && tparam <= 1 && sparam + tparam < 1)
    %ip = tp1 + sparam*uVec + tparam*vVec;
    %planePoint;
    %tp1
    %tp2
    %tp3
    ip=planePoint;
    intersect = 1;
  else 
    ip = [inf; inf; inf];
    intersect = 0;
  end
