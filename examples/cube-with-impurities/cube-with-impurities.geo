// Gmsh file for a cube of edge length 10 with a ball-shaped cavity in the
// the middle of radius 1. The object is centered around (0,0,0) and aligned
// with the axes.

// Define characteristic lengths of the edges to be used for point definition.
// Different characteristic lengths could be used for the corner points of
// the cube and the around the cavity.
lcar1 = 0.65;

// If we wanted to change these mesh sizes globally (without changing
// the above definitions), we could give a global scaling factor for
// all characteristic lengths on the command line with the `-clscale'
// option (or with `Mesh.CharacteristicLengthFactor' in an option
// file). For example, with:
//
// > gmsh t5.geo -clscale 1
//
// this input file produces a mesh of approximately 1,300 nodes and
// 11,000 tetrahedra. With
//
// > gmsh t5.geo -clscale 0.2
//
// the mesh counts approximately 350,000 nodes and 2.1 million
// tetrahedra. You can check mesh statistics in the graphical user
// interface with the `Tools->Statistics' menu.

// Define the corner points of the cube.
Point(1) = {5.0,5.0,5.0,lcar1};
Point(2) = {5.0,5.0,-5.0,lcar1};
Point(3) = {5.0,-5.0,5.0,lcar1};
Point(4) = {5.0,-5.0,-5.0,lcar1};
Point(5) = {-5.0,5.0,5.0,lcar1};
Point(6) = {-5.0,5.0,-5.0,lcar1};
Point(7) = {-5.0,-5.0,5.0,lcar1};
Point(8) = {-5.0,-5.0,-5.0,lcar1};

// Define cube lines.
Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,5};
Line(4) = {2,4};
Line(5) = {2,6};
Line(6) = {3,4};
Line(7) = {3,7};
Line(8) = {4,8};
Line(9) = {5,6};
Line(10) = {5,7};
Line(11) = {6,8};
Line(12) = {7,8};

// Define cube surfaces 
Line Loop(22) = {1,4,-6,-2};    Plane Surface(31) = {22};
Line Loop(23) = {1,5,-9,-3};    Plane Surface(32) = {23};
Line Loop(24) = {2,7,-10,-3};   Plane Surface(33) = {24};
Line Loop(25) = {4,8,-11,-5};   Plane Surface(34) = {25};
Line Loop(26) = {6,8,-12,-7};   Plane Surface(35) = {26};
Line Loop(27) = {9,11,-12,-10}; Plane Surface(36) = {27};

// Instead of using included files, we now use a user-defined function
// in order to carve some holes in the cube.
Function CheeseHole 
  // In the following commands we use the reserved variable name
  // `newp', which automatically selects a new point number. This
  // number is chosen as the highest current point number, plus
  // one. (Note that, analogously to `newp', the variables `newc',
  // `news', `newv' and `newreg' select the highest number amongst
  // currently defined curves, surfaces, volumes and `any entities
  // other than points', respectively.)
  p1 = newp; Point(p1) = {x,  y,  z,  lcar1} ;
  p2 = newp; Point(p2) = {x+r,y,  z,  lcar1} ;
  p3 = newp; Point(p3) = {x,  y+r,z,  lcar1} ;
  p4 = newp; Point(p4) = {x,  y,  z+r,lcar1} ;
  p5 = newp; Point(p5) = {x-r,y,  z,  lcar1} ;
  p6 = newp; Point(p6) = {x,  y-r,z,  lcar1} ;
  p7 = newp; Point(p7) = {x,  y,  z-r,lcar1} ;

  c1 = newreg; Circle(c1) = {p2,p1,p7};
  c2 = newreg; Circle(c2) = {p7,p1,p5};
  c3 = newreg; Circle(c3) = {p5,p1,p4};
  c4 = newreg; Circle(c4) = {p4,p1,p2};
  c5 = newreg; Circle(c5) = {p2,p1,p3};
  c6 = newreg; Circle(c6) = {p3,p1,p5};
  c7 = newreg; Circle(c7) = {p5,p1,p6};
  c8 = newreg; Circle(c8) = {p6,p1,p2};
  c9 = newreg; Circle(c9) = {p7,p1,p3};
  c10 = newreg; Circle(c10) = {p3,p1,p4};
  c11 = newreg; Circle(c11) = {p4,p1,p6};
  c12 = newreg; Circle(c12) = {p6,p1,p7};

  // We need non-plane surfaces to define the spherical holes. Here we
  // use ruled surfaces, which can have 3 or 4 sides:
  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1};
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2};
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3};
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4};
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5};
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6};
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7};
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8};

  // We then store the surface loops identification numbers in a list
  // for later reference (we will need these to define the final
  // volume):
  theloops[t] = newreg;
  Surface Loop(theloops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  thehole = newreg;
  // Don't define a volume for the ball -- it's supposed to be a cavity.
  Volume(thehole) = theloops[t] ;

Return

// We can use a `For' loop to generate five holes in the cube:

// Make cavity.
x = 4.15522 ; y = 1.72953 ; z = 1.85558; r = 0.44343;
t = 1;
Call CheeseHole;

x = 3.89734 ; y = -4.2627 ; z = -0.374816; r = 0.521064;
t = 2;
Call CheeseHole;

x = 1.82816 ; y = -2.36904 ; z = 3.06606; r = 0.352634;
t = 3;
Call CheeseHole;

x = 0.0199112 ; y = 1.36431 ; z = -2.02116; r = 0.342777;
t = 4;
Call CheeseHole;

x = 1.82075 ; y = 3.23302 ; z = 3.86981; r = 0.363792;
t = 5;
Call CheeseHole;

x = -3.83152 ; y = 1.71341 ; z = 1.09494; r = 0.322216;
t = 6;
Call CheeseHole;

x = -2.92572 ; y = -2.39699 ; z = 3.50972; r = 0.518325;
t = 7;
Call CheeseHole;

x = 3.37505 ; y = 4.13882 ; z = 4.24564; r = 0.361813;
t = 8;
Call CheeseHole;

x = 0.490033 ; y = 1.80707 ; z = 3.51843; r = 0.47904;
t = 9;
Call CheeseHole;

x = -0.636882 ; y = 2.14228 ; z = -2.88407; r = 0.323772;
t = 10;
Call CheeseHole;

x = -2.49579 ; y = 2.26235 ; z = -0.737592; r = 0.372971;
t = 11;
Call CheeseHole;

x = -0.963893 ; y = -2.91916 ; z = 3.71599; r = 0.344159;
t = 12;
Call CheeseHole;

x = 0.194181 ; y = -2.97932 ; z = -4.67614; r = 0.294657;
t = 13;
Call CheeseHole;

x = -0.605209 ; y = -4.17493 ; z = -2.21023; r = 0.483708;
t = 14;
Call CheeseHole;

x = -2.6323 ; y = 3.75534 ; z = -0.00831381; r = 0.209197;
t = 15;
Call CheeseHole;

x = 2.87222 ; y = -3.61846 ; z = -2.44607; r = 0.450271;
t = 16;
Call CheeseHole;

x = -2.55712 ; y = -3.97016 ; z = -0.777584; r = 0.607413;
t = 17;
Call CheeseHole;

x = -2.99567 ; y = -2.56708 ; z = -1.2204; r = 0.332564;
t = 18;
Call CheeseHole;

x = 3.98235 ; y = -2.31374 ; z = -1.90802; r = 0.26381;
t = 19;
Call CheeseHole;

x = -2.63074 ; y = 0.843452 ; z = 2.23671; r = 0.516049;
t = 20;
Call CheeseHole;

// We define a physical volume for each hole:
//  Physical Volume (t) = thehole ;
//  Printf("Hole %g (center = {%g,%g,%g}, radius = %g) has number %g!",
//	 t, x, y, z, r, thehole) ;

// Define outer surface of the cube.
theloops[0] = newreg;
Surface Loop(theloops[0]) = {31,32,33,34,35,36};

// The volume of the cube, without the 5 holes, is now defined by 6
// surface loops: the first surface loop defines the exterior surface;
// the surface loops other than the first one define holes.  (Again,
// to reference an array of variables, its identifier is followed by
// square brackets):
v = newreg;
Volume(v) = {theloops[]};

// We finally define a physical volume for the elements discretizing
// the cube, without the holes (whose elements were already tagged
// with numbers 1 to 5 in the `For' loop):
Physical Volume (newv) = v;
