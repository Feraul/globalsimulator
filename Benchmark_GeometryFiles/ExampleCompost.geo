cl1 = 0.1;
cl2 = 0.005;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Point(5) = {0, 0.7, 0, cl1};
Point(6) = {1, 0.3, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 1};
Line(7) = {5, 6};
Line Loop(8) = {1, 2, -7, 6};
Plane Surface(8) = {8};
Line Loop(9) = {7, 3, 4, 5};
Plane Surface(9) = {9};

Physical Point(10) = {1, 2, 3, 4, 5, 6};
Physical Line(11) = {1, 2, 3, 4, 5, 6};
Physical Line(12) = {7};
Physical Surface(13) = {8};
Physical Surface(14) = {9};

Point(7) = {0.3, 0.1, 0, cl1};
Point(8) = {0.8, 0.5, 0, cl2};
Line(15) = {7, 8};

Line{15} In Surface{8};
Line{15} In Surface{9};

Physical Point(1000) = {7, 8};
Physical Line(1000) = {15};