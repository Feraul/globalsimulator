cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {-0.0167, 0.0289, 0, cl1};
Point(4) = {0.9833, 0.0289, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Physical Point(7) = {1, 2, 3, 4};
Physical Line(8) = {1, 2, 3, 4};
Physical Surface(9) = {6};

Transfinite Line {1,3} = 12 Using Progression 1.000000;
Transfinite Line {2,4} = 12 Using Progression 1.000000;
Transfinite Surface {6} = {1,2,3,4};
