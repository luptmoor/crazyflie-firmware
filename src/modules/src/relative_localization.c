#include "relative_localization.h"
#include "debug.h"

#include <string.h>
#include <stdint.h>
#include <math.h>
#include "arm_math.h"
#include "FreeRTOS.h"
#include "task.h"

#include "param.h"

#include "log.h"
#include "system.h"
#include "led.h"

#include <radiolink.h>
#include "estimator.h"
#include "estimator_kalman.h"
#include "lpsTwrTag.h"

static bool isInit;

// static float Qvx = 1.0f; // velocity deviation
// static float Qvy = 1.0f; // velocity deviation
// static float Qvz = 1.0f; // velocity deviation
// static float Qr = 0.7f; // yaw rate deviation
// static float Ruwb = 2.0f; // ranging deviation

static float procNoise_velX = 0.2f; //0.1f;  // velocity deviation
static float procNoise_velY = 0.2f; //0.15f; // velocity deviation
static float procNoise_velZ = 0.15f; //0.2f;  // velocity deviation
static float procNoise_ryaw = 0.12f; //0.12f; // yaw rate deviation
static float measNoise_uwb = 0.15f;  //0.06;  // ranging deviation

static float InitCovPos = 0.001f;
static float InitCovYaw = 0.0001f;

static float led_thresh = 0.2f;
static float pdop = 1.0f;  // Position dilution of precision / standard decviation of distance estimate

static uint8_t in_sight = 1;
static float datum_x = 0.0f;
static float datum_y = 0.0f;
static float datum_z = 0.0f;
static float datum_psi = 0.0f;
static uint8_t datum_id = 0;

static float xAB = 0.0f;
static float yAB = 0.0f;
static float zAB = 0.0f;
static float psiAB = 0.0f;

static float xAbs = 0.0f;
static float yAbs = 0.0f;
static float zAbs = 0.0f;
static float psiAbs = 0.0f;

point_t kalman_pos;

static positionMeasurement_t absPos;



static relaVariable_t relaVar[NumUWB]; // Array of "structs" that contain P, ...
static float inputVar[NumUWB][STATE_DIM_rl];

static float A[STATE_DIM_rl][STATE_DIM_rl];
static float h[STATE_DIM_rl] = {0};
static arm_matrix_instance_f32 H = {1, STATE_DIM_rl, h};
static arm_matrix_instance_f32 Am = { STATE_DIM_rl, STATE_DIM_rl, (float *)A};

// Temporary matrices for the covariance updates
static float tmpNN1d[STATE_DIM_rl * STATE_DIM_rl];
static arm_matrix_instance_f32 tmpNN1m = { STATE_DIM_rl, STATE_DIM_rl, tmpNN1d};
static float tmpNN2d[STATE_DIM_rl * STATE_DIM_rl];
static arm_matrix_instance_f32 tmpNN2m = { STATE_DIM_rl, STATE_DIM_rl, tmpNN2d};
static float K[STATE_DIM_rl];
static arm_matrix_instance_f32 Km = {STATE_DIM_rl, 1, (float *)K};
static float tmpNN3d[STATE_DIM_rl * STATE_DIM_rl];
static arm_matrix_instance_f32 tmpNN3m = {STATE_DIM_rl, STATE_DIM_rl, tmpNN3d};
static float HTd[STATE_DIM_rl * 1];
static arm_matrix_instance_f32 HTm = {STATE_DIM_rl, 1, HTd};
static float PHTd[STATE_DIM_rl * 1];
static arm_matrix_instance_f32 PHTm = {STATE_DIM_rl, 1, PHTd};

static bool fullConnect = false; // a flag for control (fly or not)

static float vxj, vyj, vzj, rj; // receive vx, vy, vz, gz
static float vxi, vyi, vzi, ri; // self vx, vy, vz, gz
static uint16_t dij; // distance between self i and other j
static float hi, hj; // height of robot i and j
static float hij;
static float dist[NumUWB];

static inline void mat_trans(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_trans_f32(pSrc, pDst)); }

static inline void mat_inv(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_inverse_f32(pSrc, pDst)); }

static inline void mat_mult(const arm_matrix_instance_f32 * pSrcA, const arm_matrix_instance_f32 * pSrcB, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_mult_f32(pSrcA, pSrcB, pDst)); }

static inline float arm_sqrt(float32_t in)
{ float pOut = 0; arm_status result = arm_sqrt_f32(in, &pOut); configASSERT(ARM_MATH_SUCCESS == result); return pOut; }

void relativeLocoInit(void) {
  if (isInit)
    return;
  xTaskCreate(relativeLocoTask,"relative_Localization",ZRANGER_TASK_STACKSIZE, NULL,ZRANGER_TASK_PRI,NULL );
  isInit = true;
}

void relativeLocoTask(void* arg) {
  systemWaitStart();
  for (int n=0; n<NumUWB; n++) {
    for (int i=0; i<STATE_DIM_rl; i++) {
      for (int j=0; j<STATE_DIM_rl; j++) {
        relaVar[n].P[i][j] = 0;
      }
    }

    // Initialise relative state as zero in case the AR update fails
    relaVar[n].S[STATE_rlX] = 0.0f;
    relaVar[n].S[STATE_rlY] = 1.0f;
    relaVar[n].S[STATE_rlZ] = 0;
    relaVar[n].S[STATE_rlYaw] = 0;

    // Initialise error covariance matrix with variances on diagonal
    relaVar[n].P[STATE_rlX][STATE_rlX] = InitCovPos;
    relaVar[n].P[STATE_rlY][STATE_rlY] = InitCovPos;
    relaVar[n].P[STATE_rlZ][STATE_rlZ] = InitCovPos;
    relaVar[n].P[STATE_rlYaw][STATE_rlYaw] = InitCovYaw;  

    relaVar[n].receiveFlag = false;
  }


  // Main loop
  while(1) {
    
    vTaskDelay(10);
    for (int n=0; n<NumUWB; n++) {

      // AR Swarm part
      if (in_sight==1) {
        continue;

        // 1. Get current absolute position estimate from KF
        estimatorKalmanGetEstimatedPos(&kalman_pos);
        xAbs = kalman_pos.x;
        yAbs = kalman_pos.y;
        zAbs = kalman_pos.z;

        // 2. Find relative position of datum drone in this drone's body frame
        psiAbs = 0;
        xAB = (xAbs - datum_x) * arm_sin_f32(psiAbs) + (yAbs - datum_y) * arm_cos_f32(psiAbs);
        yAB = -(xAbs - datum_x) * arm_cos_f32(psiAbs) + (yAbs - datum_y) * arm_sin_f32(psiAbs);
        zAB = zAbs - datum_z;

        // 3. Initialise relatie EKF with this info
        relaVar[n].S[STATE_rlX] = xAB;
        relaVar[n].S[STATE_rlY] = yAB;
        relaVar[n].S[STATE_rlZ] = zAB;


        /// TODO: Initialise relative kalman filter with estimates from AR headset
        //  relaVar[n].S[STATE_rlX] = 0;
        // relaVar[n].S[STATE_rlY] = 0;
        // relaVar[n].S[STATE_rlZ] = 0;
        // relaVar[n].S[STATE_rlYaw] = 0;

      // Not in sight
      } else {

        // If data received from peer drone n: fetch own data
        if (twrGetSwarmInfo(n, &dij, &vxj, &vyj, &vzj, &rj, &hj)){
          estimatorKalmanGetSwarmInfo(&vxi, &vyi, &vzi, &ri, &hi);

          if(relaVar[n].receiveFlag){
            uint32_t osTick = xTaskGetTickCount();
            float dtEKF = (float)(osTick - relaVar[n].oldTimetick)/configTICK_RATE_HZ;
            relaVar[n].oldTimetick = osTick;

            // 1. Estimate relative positions
            relativeEKF(n, vxi, vyi, vzi, ri, hi, vxj, vyj, vzj, rj, hj, dij, dtEKF);
            if(n==1){hij = hj-hi;}

            // Update Kalman model inputs
            inputVar[n][STATE_rlX] = vxj;
            inputVar[n][STATE_rlY] = vyj;
            inputVar[n][STATE_rlZ] = vzj;
            inputVar[n][STATE_rlYaw] = 0.0f;


            // 2. Estimate absolute position of this drone and send to estimator
            xAB = relaVar[datum_id].S[STATE_rlX];
            yAB = relaVar[datum_id].S[STATE_rlY];
            zAB = relaVar[datum_id].S[STATE_rlZ];
            psiAB = relaVar[datum_id].S[STATE_rlYaw];

            /// TODO: Investigate if minus should be used after datum_xyz
            //psiAbs = datum_psi - psiAB;
            

            // Verified geometry
            // Py:     xAbs = x_datum + xAB * arm_sin_f32(psiAbs) - yAB * arm_cos_f32(psiAbs);
            // Py:     yAbs = y_datum + xAB * arm_cos_f32(psiAbs) + yAB * arm_sin_f32(psiAbs);

            
            psiAbs = 0.0f;
            xAbs = datum_x - yAB;
            yAbs = datum_y + xAB;
            zAbs = datum_z + zAB; // works

            absPos.x = xAbs;
            absPos.y = yAbs;
            absPos.z = zAbs;

            estimatorEnqueuePosition(&absPos); 
          } 

        } else {
          relaVar[n].oldTimetick = xTaskGetTickCount();
          relaVar[n].receiveFlag = true;
          fullConnect = true;
        }

      }
    }


    // Light up all LEDs if position of leader drone [0] is very certain
    pdop = arm_sqrt(relaVar[0].P[STATE_rlX][STATE_rlX] + relaVar[0].P[STATE_rlY][STATE_rlY] + relaVar[0].P[STATE_rlZ][STATE_rlZ]);
    if (pdop < led_thresh) {
      ledSetAll();
    } else {
      ledClearAll();
    
    }

  }
}

void relativeEKF(int n, float vxi, float vyi, float vzi, float ri, float hi, float vxj, float vyj, float vzj, float rj, float hj, uint16_t dij, float dt) {
  // some preprocessing
  arm_matrix_instance_f32 Pm = {STATE_DIM_rl, STATE_DIM_rl, (float *)relaVar[n].P};
  float cyaw = arm_cos_f32(relaVar[n].S[STATE_rlYaw]);
  float syaw = arm_sin_f32(relaVar[n].S[STATE_rlYaw]);
  float xij = relaVar[n].S[STATE_rlX];
  float yij = relaVar[n].S[STATE_rlY];
  float zij = relaVar[n].S[STATE_rlZ];

  // 1. Prediction step: forward Euler
  relaVar[n].S[STATE_rlX] = xij + (cyaw*vxj-syaw*vyj-vxi+ri*yij)*dt;
  relaVar[n].S[STATE_rlY] = yij + (syaw*vxj+cyaw*vyj-vyi-ri*xij)*dt;
  relaVar[n].S[STATE_rlZ] = zij + (vzj-vzi)*dt;
  relaVar[n].S[STATE_rlYaw] = relaVar[n].S[STATE_rlYaw] + (rj-ri)*dt;


  // 2. Calculate Error Covariance Matrix P_k+1_k

  // Jacobian Fx (4x4)
  A[0][0] = 1;
  A[0][1] = ri*dt;
  A[0][2] = 0;
  A[0][3] = (-syaw*vxj-cyaw*vyj)*dt;
  A[1][0] = -ri*dt;
  A[1][1] = 1;
  A[1][2] = 0;
  A[1][3] = (cyaw*vxj-syaw*vyj)*dt;
  A[2][0] = 0;
  A[2][1] = 0;
  A[2][2] = 1;
  A[2][3] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = 0;
  A[3][3] = 1;

  mat_mult(&Am, &Pm, &tmpNN1m); // A P
  mat_trans(&Am, &tmpNN2m); // A'
  mat_mult(&tmpNN1m, &tmpNN2m, &Pm); // P = A P A'      misses BQB' (process noise)

  // BQB' = [ Qv*c^2 + Qv*s^2 + Qr*y^2 + Qv,                       -Qr*x*y, -Qr*y]
  //        [                       -Qr*x*y, Qv*c^2 + Qv*s^2 + Qr*x^2 + Qv,  Qr*x]
  //        [                         -Qr*y,                          Qr*x,  2*Qr]*dt^2
// [ Qv + Qr*yij*yij + Qv, -Qr*yij*xij,          0,    -Qr*yij]
// [ -Qr*xij*yij,          Qv + Qr*xij*xij + Qv, 0,    Qr*xij]
// [ 0,                    0,                    2*Qv, 0]
// [ -Qr*yij,              Qr*xij,               0,    2*Qr]

 // Addition of Process noise: APA' += BQB'
  cyaw = arm_cos_f32(relaVar[n].S[STATE_rlYaw]);
  syaw = arm_sin_f32(relaVar[n].S[STATE_rlYaw]);
  float dt2 = dt*dt;

  relaVar[n].P[0][0] += dt2*((1+cyaw*cyaw)*procNoise_velX + syaw*syaw*procNoise_velY + procNoise_ryaw*yij*yij);
  relaVar[n].P[0][1] += dt2*(-procNoise_ryaw*yij*xij + cyaw*syaw*(procNoise_velX-procNoise_velY));
  relaVar[n].P[0][3] += dt2*(-procNoise_ryaw*yij);
  relaVar[n].P[1][0] += dt2*(-procNoise_ryaw*xij*yij + cyaw*syaw*(procNoise_velX-procNoise_velY));
  relaVar[n].P[1][1] += dt2*((1+cyaw*cyaw)*procNoise_velY + syaw*syaw*procNoise_velX + procNoise_ryaw*xij*xij);
  relaVar[n].P[1][3] += dt2*(procNoise_ryaw*xij);
  relaVar[n].P[2][2] += dt2*(2*procNoise_velZ);
  relaVar[n].P[3][0] += dt2*(-procNoise_ryaw*yij);
  relaVar[n].P[3][1] += dt2*(procNoise_ryaw*xij);
  relaVar[n].P[3][3] += dt2*(2*procNoise_ryaw);


  // Measurements matrix H (1x4)
  xij = relaVar[n].S[STATE_rlX];
  yij = relaVar[n].S[STATE_rlY];
  zij = relaVar[n].S[STATE_rlZ];
  float distPred = arm_sqrt(xij*xij+yij*yij+zij*zij)+0.0001f;
  float distMeas = (float)(dij/1000.0f);
  distMeas = distMeas - (0.048f*distMeas + 0.65f); // UWB bias model

  dist[n] = distMeas;
  h[0] = xij/distPred;
  h[1] = yij/distPred;
  h[2] = zij/distPred;
  h[3] = 0;



  // 3. Kalman Gain and 4. State Update
  mat_trans(&H, &HTm); // H'
  mat_mult(&Pm, &HTm, &PHTm); // PH'
  float HPHR = powf(measNoise_uwb, 2);// HPH' + R
  for (int i=0; i<STATE_DIM_rl; i++) { // Add the element of HPH' to the above
    HPHR += H.pData[i]*PHTd[i]; // this obviously only works if the update is scalar (as in this function)
  }

  for (int i=0; i<STATE_DIM_rl; i++) {
    K[i] = PHTd[i]/HPHR; // kalman gain = (PH' (HPH' + R )^-1)
    relaVar[n].S[i] = relaVar[n].S[i] + K[i] * (distMeas - distPred); // state update
  }     

  // 5. Error Covariance Matrix P_k+1_k+1
  mat_mult(&Km, &H, &tmpNN1m); // KH
  for (int i=0; i<STATE_DIM_rl; i++) { tmpNN1d[STATE_DIM_rl*i+i] -= 1; } // KH - I
  mat_trans(&tmpNN1m, &tmpNN2m); // (KH - I)'
  mat_mult(&tmpNN1m, &Pm, &tmpNN3m); // (KH - I)*P
  mat_mult(&tmpNN3m, &tmpNN2m, &Pm); // (KH - I)*P*(KH - I)'
}

bool relativeInfoRead(float* relaVarParam, float* inputVarParam){
  if(fullConnect){
    for(int i=0; i<NumUWB; i++){

      // States and inputs can be accessed from controller. TODO: send error covariance too to judge convergence
      *(relaVarParam + i*STATE_DIM_rl + STATE_rlX) = relaVar[i].S[STATE_rlX];
      *(relaVarParam + i*STATE_DIM_rl + STATE_rlY) = relaVar[i].S[STATE_rlY];
      *(relaVarParam + i*STATE_DIM_rl + STATE_rlZ) = relaVar[i].S[STATE_rlZ];
      *(relaVarParam + i*STATE_DIM_rl + STATE_rlYaw) = relaVar[i].S[STATE_rlYaw];
      *(inputVarParam + i*STATE_DIM_rl + 0) = inputVar[i][STATE_rlX];
      *(inputVarParam + i*STATE_DIM_rl + 1) = inputVar[i][STATE_rlY];
      *(inputVarParam + i*STATE_DIM_rl + 2) = inputVar[i][STATE_rlYaw];
    }   
    return true;
  }
  else
    return false;    
}







// Logging and Parameters

// LOG_GROUP_START(relativePosition)
// LOG_ADD(LOG_FLOAT, rlX0, &relaVar[0].S[STATE_rlX])
// LOG_ADD(LOG_FLOAT, rlY0, &relaVar[0].S[STATE_rlY])
// LOG_ADD(LOG_FLOAT, rlZ0, &relaVar[0].S[STATE_rlZ])
// LOG_ADD(LOG_FLOAT, rlYaw0, &relaVar[0].S[STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, PX0, &relaVar[0].P[STATE_rlX][STATE_rlX])
// LOG_ADD(LOG_FLOAT, PY0, &relaVar[0].P[STATE_rlY][STATE_rlY])
// LOG_ADD(LOG_FLOAT, PZ0, &relaVar[0].P[STATE_rlZ][STATE_rlZ])
// LOG_ADD(LOG_FLOAT, PYaw0, &relaVar[0].P[STATE_rlYaw][STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, rlX1, &relaVar[1].S[STATE_rlX])
// LOG_ADD(LOG_FLOAT, rlY1, &relaVar[1].S[STATE_rlY])
// LOG_ADD(LOG_FLOAT, rlZ1, &relaVar[1].S[STATE_rlZ])
// LOG_ADD(LOG_FLOAT, rlYaw1, &relaVar[1].S[STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, PX1, &relaVar[1].P[STATE_rlX][STATE_rlX])
// LOG_ADD(LOG_FLOAT, PY1, &relaVar[1].P[STATE_rlY][STATE_rlY])
// LOG_ADD(LOG_FLOAT, PZ1, &relaVar[1].P[STATE_rlZ][STATE_rlZ])
// LOG_ADD(LOG_FLOAT, PYaw1, &relaVar[1].P[STATE_rlYaw][STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, dist1, &dist)
// LOG_GROUP_STOP(relativePosition)

// LOG_GROUP_START(swarmInput)
// LOG_ADD(LOG_FLOAT, vx0, &vxi)
// LOG_ADD(LOG_FLOAT, vy0, &vyi)
// LOG_ADD(LOG_FLOAT, vz0, &vzi)
// LOG_ADD(LOG_FLOAT, r0, &ri)
// LOG_ADD(LOG_FLOAT, vx1, &inputVar[1][STATE_rlX])
// LOG_ADD(LOG_FLOAT, vy1, &inputVar[1][STATE_rlY])
// LOG_ADD(LOG_FLOAT, vz1, &inputVar[1][STATE_rlZ])
// LOG_ADD(LOG_FLOAT, r1, &inputVar[1][STATE_rlYaw])
// LOG_GROUP_STOP(swarmInput)

// LOG_GROUP_START(relative_pos)
// LOG_ADD(LOG_FLOAT, rlX0, &relaVar[0].S[STATE_rlX])
// LOG_ADD(LOG_FLOAT, rlY0, &relaVar[0].S[STATE_rlY])
// LOG_ADD(LOG_FLOAT, rlZ0, &relaVar[0].S[STATE_rlZ])
// LOG_ADD(LOG_FLOAT, rlYaw0, &relaVar[0].S[STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, dist0, &dist[0])
// LOG_ADD(LOG_FLOAT, rlX1, &relaVar[1].S[STATE_rlX])
// LOG_ADD(LOG_FLOAT, rlY1, &relaVar[1].S[STATE_rlY])
// LOG_ADD(LOG_FLOAT, rlYaw1, &relaVar[1].S[STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, rlZ1, &relaVar[1].S[STATE_rlZ])
// LOG_ADD(LOG_FLOAT, dist1, &dist[1])
// LOG_ADD(LOG_FLOAT, rlX2, &relaVar[2].S[STATE_rlX])
// LOG_ADD(LOG_FLOAT, rlY2, &relaVar[2].S[STATE_rlY])
// LOG_ADD(LOG_FLOAT, rlZ2, &relaVar[2].S[STATE_rlZ])
// LOG_ADD(LOG_FLOAT, rlYaw2, &relaVar[2].S[STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, dist2, &dist[2])
// LOG_ADD(LOG_FLOAT, rlX3, &relaVar[3].S[STATE_rlX])  // What happens if NumUWB = 2? 
// LOG_ADD(LOG_FLOAT, rlY3, &relaVar[3].S[STATE_rlY])
// LOG_ADD(LOG_FLOAT, rlZ3, &relaVar[3].S[STATE_rlZ])
// LOG_ADD(LOG_FLOAT, rlYaw3, &relaVar[3].S[STATE_rlYaw])
// LOG_ADD(LOG_FLOAT, dist3, &dist[3])
// LOG_ADD(LOG_FLOAT, uncert0x, &relaVar[0].P[STATE_rlX][STATE_rlX])
// LOG_ADD(LOG_FLOAT, uncert0y, &relaVar[0].P[STATE_rlY][STATE_rlY])
// LOG_ADD(LOG_FLOAT, uncert0z, &relaVar[0].P[STATE_rlZ][STATE_rlZ])
// LOG_GROUP_STOP(relative_pos)

LOG_GROUP_START(arswarm)
LOG_ADD(LOG_FLOAT, xAbs, &xAbs)
LOG_ADD(LOG_FLOAT, yAbs, &yAbs)
LOG_ADD(LOG_FLOAT, zAbs, &zAbs)
LOG_ADD(LOG_FLOAT, psiAbs, &psiAbs)
LOG_GROUP_STOP(arswarm)



// PARAM_GROUP_START(relative_pos)
// PARAM_ADD(PARAM_FLOAT, noiFlowX, &procNoise_velX) // make sure the name is not too long
// PARAM_ADD(PARAM_FLOAT, noiFlowY, &procNoise_velY) // make sure the name is not too long
// PARAM_ADD(PARAM_FLOAT, noiFlowZ, &procNoise_velZ) // make sure the name is not too long
// PARAM_ADD(PARAM_FLOAT, noiGyroZ, &procNoise_ryaw)
// PARAM_ADD(PARAM_FLOAT, noiUWB, &measNoise_uwb)
// PARAM_ADD(PARAM_FLOAT, Ppos, &InitCovPos)
// PARAM_ADD(PARAM_FLOAT, Pyaw, &InitCovYaw)
// PARAM_ADD(PARAM_FLOAT, led_thresh, &led_thresh)
// PARAM_GROUP_STOP(relative_pos)

PARAM_GROUP_START(arswarm)
PARAM_ADD(PARAM_UINT8, in_sight, &in_sight)
PARAM_ADD(PARAM_FLOAT, datum_x, &datum_x)
PARAM_ADD(PARAM_FLOAT, datum_y, &datum_y)
PARAM_ADD(PARAM_FLOAT, datum_z, &datum_z)
PARAM_ADD(PARAM_FLOAT, datum_psi, &datum_psi)
PARAM_ADD(PARAM_UINT8, datum_id, &datum_id)
PARAM_GROUP_STOP(arswarm)

