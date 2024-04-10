#include "system.h"
#include "FreeRTOS.h"
#include "task.h"
#include "commander.h"
#include "relative_localization.h"
#include "num.h"
#include "param.h"
#include "debug.h"
#include <stdlib.h> // random
#include "lpsTwrTag.h" // UWBNum
#include "configblock.h"
#include "uart2.h"
#include "log.h"
#include <math.h>
#include "arm_math.h"
#include "estimator_kalman.h"
#include "estimator.h"

#define USE_MONOCAM 0

static bool isInit;
static bool keepFlying = true;
static setpoint_t setpoint;
static float relaVarInCtrl[NumUWB][STATE_DIM_rl];
static float inputVarInCtrl[NumUWB][STATE_DIM_rl];
static uint8_t selfID;

static uint8_t in_sight = true;
static float datum_x = 0.0f;
static float datum_y = 0.0f;
static float datum_z = 0.0f;
static float datum_psi = 0.0f;
static uint8_t datum_id = 0;

static float xAB;
static float yAB;
static float zAB;
static float psiAB;

static float xAbs;
static float yAbs;
static float zAbs;
static float psiAbs;

static positionMeasurement_t absPos;

// static float relaCtrl_p = 2.0f;
// static float relaCtrl_i = 0.0001f;
// static float relaCtrl_d = 0.01f;


static void setHoverSetpoint(setpoint_t *setpoint, float vx, float vy, float z, float yawrate) {
  setpoint->mode.z = modeAbs;
  setpoint->position.z = z;
  setpoint->mode.yaw = modeVelocity;
  setpoint->attitudeRate.yaw = yawrate;
  setpoint->mode.x = modeVelocity;
  setpoint->mode.y = modeVelocity;
  setpoint->velocity.x = vx;
  setpoint->velocity.y = vy;
  setpoint->velocity_body = true;
  commanderSetSetpoint(setpoint, 3);
}

// static void set3DVel(setpoint_t *setpoint, float vx, float vy, float vz, float yawrate) {
//   setpoint->mode.z = modeVelocity;
//   setpoint->mode.yaw = modeVelocity;
//   setpoint->attitudeRate.yaw = yawrate;
//   setpoint->mode.x = modeVelocity;
//   setpoint->mode.y = modeVelocity;
//   setpoint->velocity.x = vx;
//   setpoint->velocity.y = vy;
//   setpoint->velocity.z = vz;
//   setpoint->velocity_body = true;
//   commanderSetSetpoint(setpoint, 3);
// }


// #define SIGN(a) ((a>=0)?1:-1)

// static float PreErr_x = 0;
// static float PreErr_y = 0;
// static float PreErr_z = 0;
// static float IntErr_x = 0;
// static float IntErr_y = 0;
// static float IntErr_z = 0;
// static uint32_t PreTime;


// static void formation0asCenter(float tarX, float tarY, float tarZ) {
//   float dt = (float)(xTaskGetTickCount()-PreTime)/configTICK_RATE_HZ;
//   PreTime = xTaskGetTickCount();
//   if(dt > 1) // skip the first run of the EKF
//     return;
//   // pid control for formation flight
//   float err_x = -(tarX - relaVarInCtrl[0][STATE_rlX]);
//   float err_y = -(tarY - relaVarInCtrl[0][STATE_rlY]);
//   float err_z = -(tarZ - relaVarInCtrl[0][STATE_rlZ]);
//   float pid_vx = relaCtrl_p * err_x;
//   float pid_vy = relaCtrl_p * err_y;
//   float pid_vz = relaCtrl_p * err_z;

//   // Derivative gain
//   float dx = (err_x - PreErr_x) / dt;
//   float dy = (err_y - PreErr_y) / dt;
//   float dz = (err_z - PreErr_z) / dt;
//   PreErr_x = err_x;
//   PreErr_y = err_y;
//   PreErr_z = err_z;
//   pid_vx += relaCtrl_d * dx;
//   pid_vy += relaCtrl_d * dy;
//   pid_vz += relaCtrl_d * dz;

//   // Integral gain
//   IntErr_x += err_x * dt;
//   IntErr_y += err_y * dt;
//   IntErr_z += err_z * dt;
//   pid_vx += relaCtrl_i * constrain(IntErr_x, -0.5, 0.5);
//   pid_vy += relaCtrl_i * constrain(IntErr_y, -0.5, 0.5);
//   pid_vz += relaCtrl_i * constrain(IntErr_z, -0.5, 0.5);
//   pid_vx = constrain(pid_vx, -1.5f, 1.5f);
//   pid_vy = constrain(pid_vy, -1.5f, 1.5f);
//   pid_vz = constrain(pid_vz, -0.3f, 0.3f);

//   set3DVel(&setpoint, pid_vx, pid_vy, pid_vz, 0);
// }


void relativeControlTask(void* arg) {
  systemWaitStart();

  // Main Loop
  while(1) {
    vTaskDelay(10);
    DEBUG_PRINT("Control loop\n");
    if(in_sight){
      ///TODO: pass position estimate to relative positioning
      continue; // Do not do anything if in sight of AR headset
    }



    // If we should fly and information can be read from the relative UWB localisation
    if(relativeInfoRead((float *)relaVarInCtrl, (float *)inputVarInCtrl) && keepFlying) {

      // Hold relative Position to leader drone. CONVERGENCE ASSUMED
      xAB = relaVarInCtrl[datum_id][STATE_rlX];
      yAB = relaVarInCtrl[datum_id][STATE_rlY];
      zAB = relaVarInCtrl[datum_id][STATE_rlZ];
      psiAB = relaVarInCtrl[datum_id][STATE_rlYaw];

      /// TODO: Investigate if minus should be used after datum_xyz
      psiAbs = datum_psi - psiAB;
      xAbs = datum_x + xAB * arm_cos_f32(psiAbs) + yAB * arm_sin_f32(psiAbs);
      yAbs = datum_y + xAB * arm_sin_f32(psiAbs) + yAB * arm_cos_f32(psiAbs);
      zAbs = datum_z + zAB;

      absPos.x = xAbs;
      absPos.y = yAbs;
      absPos.z = zAbs;

      estimatorEnqueuePosition(&absPos);



      //formation0asCenter(xAB, yAB, zAB);

      
        
          //-cosf(relaVarInCtrl[0][STATE_rlYaw])*relaXof2in1 + sinf(relaVarInCtrl[0][STATE_rlYaw])*relaYof2in1;
          //-sinf(relaVarInCtrl[0][STATE_rlYaw])*relaXof2in1 - cosf(relaVarInCtrl[0][STATE_rlYaw])*relaYof2in1;

    // If we should not fly: LANDING procedure
    } else {
      for (int i=1; i<5; i++) {
        if(selfID!=0){
        setHoverSetpoint(&setpoint, 0, 0, 0.3f-(float)i*0.05f, 0);
        vTaskDelay(M2T(10));
        }
      }
    }

  }
}


void relativeControlInit(void) {
  if (isInit)
    return;

  selfID = (uint8_t)(((configblockGetRadioAddress()) & 0x000000000f) - 1);

  xTaskCreate(relativeControlTask,"relative_Control",configMINIMAL_STACK_SIZE, NULL, 3, NULL);
  isInit = true;
}




// Logging and Parameters

PARAM_GROUP_START(arswarm)
PARAM_ADD(PARAM_UINT8, keepFlying, &keepFlying)
// PARAM_ADD(PARAM_FLOAT, relaCtrl_p, &relaCtrl_p)
// PARAM_ADD(PARAM_FLOAT, relaCtrl_i, &relaCtrl_i)
// PARAM_ADD(PARAM_FLOAT, relaCtrl_d, &relaCtrl_d)

PARAM_ADD(PARAM_UINT8, in_sight, &in_sight)
PARAM_ADD(PARAM_FLOAT, datum_x, &datum_x)
PARAM_ADD(PARAM_FLOAT, datum_y, &datum_y)
PARAM_ADD(PARAM_FLOAT, datum_z, &datum_z)
PARAM_ADD(PARAM_FLOAT, datum_psi, &datum_psi)
PARAM_ADD(PARAM_UINT8, datum_id, &datum_id)

PARAM_GROUP_STOP(arswarm)


LOG_GROUP_START(arswarm)
LOG_ADD(LOG_FLOAT, xAbs, &xAbs)
LOG_ADD(LOG_FLOAT, yAbs, &yAbs)
LOG_ADD(LOG_FLOAT, zAbs, &zAbs)
LOG_ADD(LOG_FLOAT, psiAbs, &psiAbs)

LOG_GROUP_STOP(arswarm)




