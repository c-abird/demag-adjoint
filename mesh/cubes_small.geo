n = 1;

// cube 1
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {0, 1, 0};
Point(4) = {1, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Transfinite Line {1, 2, 3, 4} = n;
Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};
Transfinite Surface {1};


// cube 2
Point(11) = {0, 0, 3};
Point(12) = {1, 0, 3};
Point(13) = {0, 1, 3};
Point(14) = {1, 1, 3};

Line(11) = {11, 12};
Line(12) = {12, 14};
Line(13) = {14, 13};
Line(14) = {13, 11};

Transfinite Line {11, 12, 13, 14} = n;
Line Loop(11) = {11, 12, 13, 14};

Plane Surface(11) = {11};
Transfinite Surface {11};

m1[] = Extrude{0, 0, 1} {Surface{1}; Layers{n};};
Physical Volume(1) = {m1[1]};

m2[] = Extrude{0, 0, 1} {Surface{11}; Layers{n};};
Physical Volume(2) = {m2[1]};
