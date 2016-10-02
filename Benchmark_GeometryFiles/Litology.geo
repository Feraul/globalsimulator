cl1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {0.5, 0, 0, 1};
Point(3) = {1, 0, 0, 1};
Point(4) = {1, 0.5, 0, 1};
Point(5) = {1, 1, 0, 1};
Point(6) = {0.5, 1, 0, 1};
Point(7) = {0, 1, 0, 1};
Point(8) = {0, 0.5, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 8};
Line(10) = {3, 7};
Line(11) = {4, 6};
Line Loop(13) = {1, 9, 8};
Plane Surface(13) = {13};
Line Loop(15) = {-9, 2, 10, 7};
Plane Surface(15) = {15};
Line Loop(17) = {-10, 3, 11, 6};
Plane Surface(17) = {17};
Line Loop(19) = {-11, 4, 5};
Plane Surface(19) = {19};
Physical Point(20) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Line(21) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Line(22) = {9, 10, 11};
Physical Surface(23) = {13};
Physical Surface(24) = {15};
Physical Surface(25) = {17};
Physical Surface(26) = {19};

Transfinite Line {1,9} = 10 Using Progression 1.000000;
Transfinite Line {8,9} = 10 Using Progression 1.000000;
Transfinite Surface {13} = {1,2,8};


Transfinite Line {2,7} = 10 Using Progression 1.000000;
Transfinite Line {9,10} = 14 Using Progression 1.000000;
Transfinite Surface {15} = {2,3,7,8};


Transfinite Line {6,3} = 10 Using Progression 1.000000;
Transfinite Line {11,10} = 14 Using Progression 1.000000;
Transfinite Surface {17} = {7,6,4,3};


Transfinite Line {11,5} = 10 Using Progression 1.000000;
Transfinite Line {4,5} = 10 Using Progression 1.000000;
Transfinite Surface {19} = {4,5,6};

