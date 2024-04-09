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
#include "estimator_kalman.h"
#define USE_MONOCAM 0

static bool isInit;
static bool onGround = true;
static bool keepFlying = false;
static setpoint_t setpoint;
static float relaVarInCtrl[NumUWB][STATE_DIM_rl];
static float inputVarInCtrl[NumUWB][STATE_DIM_rl];
static uint8_t selfID;
static float height;
static float initial_hover_height;

static float kf_init_layer_distance = 0.30f; // [m]
static uint8_t layer_number;
static uint8_t layer_map[5] = {2, 4, 3, 1, 0};

static uint8_t in_sight = true;
static float datum_x = 0.0f;
static float datum_y = 0.0f;
static float datum_z = 0.0f;
static uint8_t datum_id = 0;


static float relaCtrl_p = 2.0f;
static float relaCtrl_i = 0.0001f;
static float relaCtrl_d = 0.01f;
// static float NDI_k = 2.0f;
static char c = 0; // monoCam


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

static void set3DVel(setpoint_t *setpoint, float vx, float vy, float vz, float yawrate)
{
  setpoint->mode.z = modeVelocity;
  setpoint->mode.yaw = modeVelocity;
  setpoint->attitudeRate.yaw = yawrate;
  setpoint->mode.x = modeVelocity;
  setpoint->mode.y = modeVelocity;
  setpoint->velocity.x = vx;
  setpoint->velocity.y = vy;
  setpoint->velocity.z = vz;
  setpoint->velocity_body = true;
  commanderSetSetpoint(setpoint, 3);
}


static void flyRandomIn1meter(float vel){
  float randomYaw = (rand() / (float)RAND_MAX) * 6.28f; // 0-2pi rad
  float randomVel = vel*(rand() / (float)RAND_MAX); // 0-1 m/s
  float vxBody = randomVel * cosf(randomYaw);
  float vyBody = randomVel * sinf(randomYaw);
  for (int i=1; i<100; i++) {
    setHoverSetpoint(&setpoint, vxBody, vyBody, height, 0);
    vTaskDelay(M2T(10));
  }
  for (int i=1; i<100; i++) {
    setHoverSetpoint(&setpoint, -vxBody, -vyBody, height, 0);
    vTaskDelay(M2T(10));
  }
}



#define SIGN(a) ((a>=0)?1:-1)
static float targetX;
static float targetY;
static float targetZ;
static float PreErr_x = 0;
static float PreErr_y = 0;
static float PreErr_z = 0;
static float IntErr_x = 0;
static float IntErr_y = 0;
static float IntErr_z = 0;
static uint32_t PreTime;


static void formation0asCenter(float tarX, float tarY, float tarZ) {
  float dt = (float)(xTaskGetTickCount()-PreTime)/configTICK_RATE_HZ;
  PreTime = xTaskGetTickCount();
  if(dt > 1) // skip the first run of the EKF
    return;
  // pid control for formation flight
  float err_x = -(tarX - relaVarInCtrl[0][STATE_rlX]);
  float err_y = -(tarY - relaVarInCtrl[0][STATE_rlY]);
  float err_z = -(tarZ - relaVarInCtrl[0][STATE_rlZ]);
  float pid_vx = relaCtrl_p * err_x;
  float pid_vy = relaCtrl_p * err_y;
  float pid_vz = relaCtrl_p * err_z;

  // Derivative gain
  float dx = (err_x - PreErr_x) / dt;
  float dy = (err_y - PreErr_y) / dt;
  float dz = (err_z - PreErr_z) / dt;
  PreErr_x = err_x;
  PreErr_y = err_y;
  PreErr_z = err_z;
  pid_vx += relaCtrl_d * dx;
  pid_vy += relaCtrl_d * dy;
  pid_vz += relaCtrl_d * dz;

  // Integral gain
  IntErr_x += err_x * dt;
  IntErr_y += err_y * dt;
  IntErr_z += err_z * dt;
  pid_vx += relaCtrl_i * constrain(IntErr_x, -0.5, 0.5);
  pid_vy += relaCtrl_i * constrain(IntErr_y, -0.5, 0.5);
  pid_vz += relaCtrl_i * constrain(IntErr_z, -0.5, 0.5);
  pid_vx = constrain(pid_vx, -1.5f, 1.5f);
  pid_vy = constrain(pid_vy, -1.5f, 1.5f);
  pid_vz = constrain(pid_vz, -0.3f, 0.3f);

  // float rep_x = 0.0f;
  // float rep_y = 0.0f;
  // for(uint8_t i=0; i<NumUWB; i++){
  //   if(i!=selfID){
  //     float dist = relaVarInCtrl[i][STATE_rlX]*relaVarInCtrl[i][STATE_rlX] + relaVarInCtrl[i][STATE_rlY]*relaVarInCtrl[i][STATE_rlY];
  //     dist = sqrtf(dist);
  //     rep_x += -0.5f * (SIGN(0.5f - dist) + 1) / (abs(relaVarInCtrl[i][STATE_rlX]) + 0.001f) * SIGN(relaVarInCtrl[i][STATE_rlX]);
  //     rep_y += -0.5f * (SIGN(0.5f - dist) + 1) / (abs(relaVarInCtrl[i][STATE_rlY]) + 0.001f) * SIGN(relaVarInCtrl[i][STATE_rlY]);
  //   }
  // }
  // rep_x = constrain(rep_x, -1.5f, 1.5f);
  // rep_y = constrain(rep_y, -1.5f, 1.5f);

  // pid_vx = constrain(pid_vx + rep_x, -1.5f, 1.5f);
  // pid_vy = constrain(pid_vy + rep_y, -1.5f, 1.5f);

  set3DVel(&setpoint, pid_vx, pid_vy, pid_vz, 0);
}

// static void NDI_formation0asCenter(float tarX, float tarY){
//   float err_x = -(tarX - relaVarInCtrl[0][STATE_rlX]);
//   float err_y = -(tarY - relaVarInCtrl[0][STATE_rlY]);
//   float rela_yaw = relaVarInCtrl[0][STATE_rlYaw];
//   float Ru_x = cosf(rela_yaw)*inputVarInCtrl[0][STATE_rlX] - sinf(rela_yaw)*inputVarInCtrl[0][STATE_rlY];
//   float Ru_y = sinf(rela_yaw)*inputVarInCtrl[0][STATE_rlX] + cosf(rela_yaw)*inputVarInCtrl[0][STATE_rlY];
//   float ndi_vx = NDI_k*err_x + 0*Ru_x;
//   float ndi_vy = NDI_k*err_y + 0*Ru_y;
//   ndi_vx = constrain(ndi_vx, -1.5f, 1.5f);
//   ndi_vy = constrain(ndi_vy, -1.5f, 1.5f);  
//   setHoverSetpoint(&setpoint, ndi_vx, ndi_vy, height, 0);
// }

void relativeControlTask(void* arg) {
  static uint32_t takeoff_tick, height_change_tick;
  systemWaitStart();
  static logVarId_t logIdStateIsFlying;
  logIdStateIsFlying = logGetVarId("kalman", "inFlight");


  // Main Loop
  while(1) {

    //Check if we should fly
    vTaskDelay(10);
    if(selfID==0){
      keepFlying = logGetUint(logIdStateIsFlying);
      keepFlying = command_share(selfID, keepFlying);
      continue; // Do not send commands to leader drone, as it is manually controlled
    }

    keepFlying = command_share(selfID, keepFlying);


    // If we should fly and information can be read from the localisation
    if(relativeInfoRead((float *)relaVarInCtrl, (float *)inputVarInCtrl) && keepFlying) {


      // TAKE OFF if on ground
      if(onGround) {
        estimatorKalmanInit(); // reseting kalman filter
        vTaskDelay(M2T(2000));
        for (int i=0; i<50; i++) {
          setHoverSetpoint(&setpoint, 0, 0, initial_hover_height, 0);
          vTaskDelay(M2T(100));
        }
        onGround = false;

        // Define time datum
        takeoff_tick = xTaskGetTickCount();
        height_change_tick = takeoff_tick;
        // Initialise RNG differently for each drone
        srand(selfID * takeoff_tick);
        layer_number = selfID;
        height = initial_hover_height + layer_number * kf_init_layer_distance;
      }


      // Start two timers after takeoff: one for absolute time, and one for height changes
      uint32_t ticks_since_takeoff = xTaskGetTickCount() - takeoff_tick;
      uint32_t ticks_since_height_change = xTaskGetTickCount() - height_change_tick;

      // CONVERGENCE FLIGHT for 21s
      if(ticks_since_takeoff < 21000) {
        
        // ID1 starts on layer 1 = 1.10m.
        // ID2 starts on layer 2 = 1.40m;
        // ID3 starts on layer 3 = 1.70m;
        // ID4 starts on layer 4 = 2.00m;

        // Change hovering layer every 3s
        if (ticks_since_height_change > 3000) {
          layer_number++;
          layer_number = layer_number % 5;
          height = initial_hover_height + layer_map[layer_number] * kf_init_layer_distance;
          // Reset timer after height change
          height_change_tick = xTaskGetTickCount(); 
        }

        flyRandomIn1meter(0.7f);

        targetX = relaVarInCtrl[0][STATE_rlX];
        targetY = relaVarInCtrl[0][STATE_rlY];
        targetZ = relaVarInCtrl[0][STATE_rlZ];


      } else {
        // CONVERGENCE PROCEDURE FINISHED
        // Hold relative positions for 5s.
        if ( (ticks_since_takeoff > 21000)){ 
          srand((unsigned int) relaVarInCtrl[0][STATE_rlX]*100);
          formation0asCenter(targetX, targetY, targetZ);
          // NDI_formation0asCenter(targetX, targetY);
        }
       
          //-cosf(relaVarInCtrl[0][STATE_rlYaw])*relaXof2in1 + sinf(relaVarInCtrl[0][STATE_rlYaw])*relaYof2in1;
          //-sinf(relaVarInCtrl[0][STATE_rlYaw])*relaXof2in1 - cosf(relaVarInCtrl[0][STATE_rlYaw])*relaYof2in1;
        
      }



    // If we should not fly: LANDING procedure
    } else {
      if(!onGround){
        for (int i=1; i<5; i++) {
          if(selfID!=0){
          setHoverSetpoint(&setpoint, 0, 0, 0.3f-(float)i*0.05f, 0);
          vTaskDelay(M2T(10));
          }
        }
        onGround = true;
      } 
    }

  }
}


void relativeControlInit(void) {
  if (isInit)
    return;

  selfID = (uint8_t)(((configblockGetRadioAddress()) & 0x000000000f) - 5);

  xTaskCreate(relativeControlTask,"relative_Control",configMINIMAL_STACK_SIZE, NULL, 3, NULL);
  height = 1.0f;
  initial_hover_height = 0.8f;

 

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
PARAM_ADD(PARAM_UINT8, datum_id, &datum_id)


PARAM_GROUP_STOP(arswarm)

LOG_GROUP_START(mono_cam)
LOG_ADD(LOG_UINT8, charCam, &c)
LOG_GROUP_STOP(mono_cam)