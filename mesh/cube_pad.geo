// measurement
dz  = 0.02;
dxy = 2;
z   = 1.5;
nxy = 20;

Point(6) = {0.5-dxy/2, 0.5-dxy/2, z};
Point(7) = {0.5-dxy/2, 0.5+dxy/2, z};
Point(8) = {0.5+dxy/2, 0.5-dxy/2, z};
Point(9) = {0.5+dxy/2, 0.5+dxy/2, z};

Line(6) = {6, 7};
Line(7) = {7, 9};
Line(8) = {9, 8};
Line(9) = {8, 6};

Transfinite Line {6, 7, 8, 9} = nxy;
Line Loop(10) = {6, 7, 8, 9};

Plane Surface(2) = {10};
Transfinite Surface {2};

h[] = Extrude{0, 0, dz} {Surface{2}; Layers{1};};
Physical Volume(2) = {h[1]};

// magnet
n = 10;

Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {0, 1, 0};
Point(4) = {1, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Transfinite Line {1, 2, 3, 4} = n;
Line Loop(5) = {1, 2, 3, 4};

Plane Surface(1) = {5};
Transfinite Surface {1};

m[] = Extrude{0, 0, 1} {Surface{1}; Layers{n};};
Physical Volume(1) = {m[1]};

