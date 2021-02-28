# ZPinchArray

 This class is used to visualize the magnetic field and
%     j x B force on unit j ez for Z pinch cylindical arrays in the 2D x-y plane
%     The current flows in the z-direction
%
%     Functions:
%
%       ZpinchArray(N,R0,d0,I,isExploder) % Creates a Z pinch array 
%
%       [Bx_t,By_t,B] = getMagField(obj,pts_x, pts_y)
%             % Returns the total magnetic field from teh superposition of
%             % magentic field around each wire
%
%        out = getth(~,x,y) % Returns azmuthal angle for a point(x,y)
%
%        showBField(obj,M,L) % shows the total B field for the grid
%
%        [fx,fy] = getJxBDirection(obj,Bx,By,B) % returns the of the j x B force on a unit j 
