// Gmsh project created on Thu Mar 30 16:48:37 2023
SetFactory("OpenCASCADE");
//+

raff = 0.8;

raff2 = 0.1;

Point(1) = {0, 0, 0, raff};
//+
Point(2) = {-1, -1, 0,raff};
//+
Point(3) = {1, 1, 0, raff};
//+
Point(4) = {-1, 1, 0, raff};
//+
Point(5) = {0, -1, 0, raff};
//+
Point(6) = {-0, -0.5, 0, raff2};
//+
Point(7) = {1, 0, 0, raff};
//+
Point(8) = {0.5, -0, 0, raff2};
//+
Point(9) = {0, 0.5, 0, raff2};
//+
Point(10) = {-0.5, -0, 0, raff2};
//+
Line(1) = {4, 3};
//+
Line(2) = {7, 3};
//+
Line(3) = {7, 8};
//+
Line(4) = {8, 1};
//+
Line(5) = {1, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 2};
//+
Line(8) = {2, 4};
//+
Circle(9) = {6, 1, 10};
//+
Circle(10) = {10, 1, 9};
//+
Circle(11) = {8, 1, 9};


Physical Curve("circle_boundary", 19) = {5, 9, 10, 11, 4};
//+
Physical Curve("full_boundary", 20) = {8, 7, 6, 5, 4, 3, 2, 1};
//+
Physical Curve("circle_only", 21) = {9, 10, 11};
//+
Curve Loop(1) = {8, 1, -2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};

Curve{9} In Surface{1};
Curve{10} In Surface{1};
Curve{11} In Surface{1};
//+
Curve Loop(2) = {10, -11, 4, 5, 9};
//+
//Plane Surface(2) = {2};

//+
Physical Surface("domain", 22) = {1};

//+
Physical Surface("circle_domain", 23) = {3};
//+
Curve Loop(3) = {9, 10, -11, 4, 5};
