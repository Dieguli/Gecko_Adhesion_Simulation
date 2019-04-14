#include "PhysicsParams.h"
#include "Utils.h"

//---------------------------------------------------------------------------

PhysicsParams physicsParams;

//---------------------------------------------------------------------------
void physicsParamsInit()
//---------------------------------------------------------------------------
{
  physicsParams.warpedStiffness = true;

  physicsParams.timeStep = 0.01;          // s
  physicsParams.solverIterations = 20;
  physicsParams.gravityConst = 9.81;      // m/s^2
  physicsParams.grabForceRadius = 0.05;   // m
  physicsParams.dynamicDamping = 2.0;
  physicsParams.grabSpringConst = 100.0;  // N/m

  physicsParams.youngModulus = 0.02e6; // N/m^2
  physicsParams.poissonRatio = 0.40;
  physicsParams.density = 1000.0; // kg / m^3
}



