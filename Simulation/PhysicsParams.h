//---------------------------------------------------------------------------

#ifndef PhysicsParamsH
#define PhysicsParamsH
//---------------------------------------------------------------------------

#include <stdio.h>

#define SIMULATE_DYNAMIC 0
#define SIMULATE_STATIC  1
#define SIMULATE_RIGID   2
#define SIMULATE_NONE    3


struct PhysicsParams {
  bool warpedStiffness;

  float timeStep;                // s
  int   solverIterations;
  float gravityConst;            // m/s^2
  float grabForceRadius;         // m
  float dynamicDamping;
  float grabSpringConst;         // N/m

  float youngModulus;            // N/m^2
  float poissonRatio;
  float density;                 // kg / m^3
};

extern PhysicsParams physicsParams;

void physicsParamsInit();

#endif
