cl__1 = 0.1;
Point(1) = {0, 0, 0, 0.1};
Point(2) = {1, 0, 0, 0.1};
Point(3) = {1, 1, 0, 0.1};
Point(4) = {0, 1, 0, 0.1};
Point(5) = {0.5, 0, 0, 0.1};
Point(6) = {0.5, 1, 0, 0.1};
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 6};
Line(5) = {6, 4};
Line(6) = {4, 1};
Line(7) = {5, 6};
Line Loop(9) = {1, 7, 5, 6};
Plane Surface(9) = {9};
Line Loop(11) = {-7, 2, 3, 4};
Plane Surface(11) = {11};
Physical Point(201) = {1, 2, 3, 4, 5, 6};
Physical Line(201) = {1, 2, 3, 4, 5, 6};
Physical Line(14) = {7};
Physical Surface(1) = {9};
Physical Surface(2) = {11};

Transfinite Line {1,5} = 6 Using Progression 1.000000;
Transfinite Line {2,4} = 6 Using Progression 1.000000;
Transfinite Surface {9} = {1,5,6,4} Right;

Transfinite Line {6,7} = 11 Using Progression 1.000000;
Transfinite Line {7,3} = 11 Using Progression 1.000000;
Transfinite Surface {11} = {5,2,3,6} Right;

