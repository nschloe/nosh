Point(1) = {-6, -2, -1, 1.0};
Point(2) = {-6, 2, -1, 1.0};
Point(3) = {6, 2, -1, 1.0};
Point(4) = {6, -2, -1, 1.0};
Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Point(5) = {-4, 0, -1, 1.0};
Point(6) = {-4, 1, -1, 1.0};
Point(7) = {-4, -1, -1, 1.0};
Point(8) = {-3, 0, -1, 1.0};
Point(9) = {-5, 0, -1, 1.0};
Circle(5) = {6, 5, 9};
Circle(6) = {9, 5, 7};
Circle(7) = {7, 5, 8};
Circle(8) = {8, 5, 6};
Point(10) = {11, 1, -1, 1.0};
Delete {
  Point{10};
}
Point(10) = {5, 1, -1, 1.0};
Point(11) = {3, 1, -1, 1.0};
Point(12) = {3, -1, -1, 1.0};
Point(13) = {5, -1, -1, 1.0};
Line(9) = {11, 10};
Line(10) = {13, 12};
Line(11) = {12, 11};
Line(12) = {10, 13};
Point(14) = {0, 1, -1, 1.0};
Point(15) = {-0.86602540378443904, -0.5, -1, 1.0};
Point(16) = {0.86602540378443904, -0.5, -1, 1.0};
Line(13) = {14, 15};
Line(14) = {15, 16};
Line(15) = {16, 14};
Line Loop(16) = {4, 1, 2, 3};
Plane Surface(17) = {16};
Line Loop(18) = {11, 9, 12, 10};
Line Loop(19) = {15, 13, 14};
Line Loop(20) = {8, 5, 6, 7};
Plane Surface(21) = {16, 18, 19, 20};
Plane Surface(22) = {16, 18, 19, 20};
Plane Surface(23) = {16, 18, 19, 20};
Delete {
  Surface{17};
}
Delete {
  Surface{21};
}
Delete {
  Surface{22};
}
Delete {
  Surface{23};
}
Plane Surface(23) = {16, 18, 19, 20};
Translate {0, 0, 1} {
  Surface{23};
}
Translate {0, 0, 1} {
  Surface{23};
}
Extrude {0, 0, 2} {
  Surface{23};
}
