Point(1) = {0, 0, 0, 1};
Point(2) = {16, 0, 0, 1};
Point(3) = {0, 1.1, 0, 1};
Point(4) = {16, 1.1, 0, 1};
Point(5) = {5, 1.1, 0, 0.2};
Point(6) = {6.5, 1.1, 0, 0.2};
Point(7) = {9.5, 1.1, 0, 0.2};
Point(8) = {11, 1.1, 0, 0.2};
Point(9) = {5, 0, 0, 0.2};
Point(10) = {6.5, 0, 0, 0.2};
Point(11) = {9.5, 0, 0, 0.2};
Point(12) = {11, 0, 0, 0.2};

Line(21) = {1, 3};
Line(22) = {2, 4};
Line(23) = {1, 9};
Line(24) = {9, 10};
Line(25) = {10, 11};
Line(26) = {11, 12};
Line(27) = {12, 2};
Line(28) = {3, 5};
Line(29) = {5, 6};
Line(30) = {6, 7};
Line(31) = {7, 8};
Line(32) = {8, 4};

Physical Line(1) = {21};
Physical Line(2) = {22};
Physical Line(3) = {23, 24, 25, 26, 27};
Physical Line(4) = {28, 29, 30, 31, 32};

Line Loop(45) = {21, 28, 29, 30, 31, 32, -22, -27, -26, -25, -24, -23};
Plane Surface(46) = {45};
Physical Surface(47) = {46};
