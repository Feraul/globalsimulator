cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {0.5, 0, 0, cl1};
Point(3) = {0.5, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(8) = {2, 3, 4, 1};
Plane Surface(8) = {8};
Physical Point(9) = {1, 2, 3, 4};
Physical Line(10) = {1};
Physical Line(11) = {3};
Physical Line(12) = {2, 4};
Physical Surface(13) = {8};
