#include "epot_bicgstabsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "geomplotter.hpp"
#include "config.h"
#include "gtkplotter.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
using namespace std;

double getRandom(double min, double max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

double vec2ang(Vec3D v1, Vec3D v2){
  return acos((v1 * v2) /(v1.norm2() * v2.norm2()));
}

int main( int argc, char **argv )
{
  double N = 100000000;
  double th_lim = M_PI/64;
  double accept = 0.004;
  double u,v,th,ph,ang,ph2, tan1;
  double b = 0.035;
  double g = 1/sqrt(1 - b * b);
  double alpha = 1.11074;
  double amax = alpha + 0.05;
  double amin = alpha - 0.05;
  double sim_reg = 2 * th_lim *(cos(amin) - cos(amax));
  double in_sphere = N * (4 * M_PI)/(sim_reg);
  Vec3D ref = Vec3D(sin(alpha), 0, cos(alpha));
  Vec3D temp;
  double in_range = 0;
  double in_shift = 0;
  for(int i = 0; i < N; i++){
    u = getRandom(0,1);
    v = getRandom(cos(amax),cos(amin));
    th = (u - 0.5) * 2 * th_lim;
    ph = acos(v);
    temp = Vec3D(sin(ph)*cos(th), sin(ph)*sin(th), cos(ph));
    ang = vec2ang(temp, ref);
    if(abs(ang) <= accept){
      in_range++;
    }
    tan1 = sin(ph)/(g*(cos(ph) + b));
    ph2 = atan(tan1);
    temp = Vec3D(sin(ph2)*cos(th), sin(ph2)*sin(th), cos(ph2));
    ang = vec2ang(temp, ref);
    if(abs(ang) <= accept){
      in_shift++;
    }
    
  }
  cout << "slice: " << N << "\t sphere: " << in_sphere << "\n";
  cout << "Counted: " << in_range << "\t fraction: " << in_range/in_sphere << "\n";
  cout << "Shifted: " << in_shift << "\t fraction: " << in_shift/in_sphere << "\n";
  cout << "effect" << (in_shift - in_range)/(in_range) << "\n";
  cout << "Thin ang ref: " << accept * accept/4 <<"\n";
     
}
