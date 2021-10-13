#include "deck.h"
#include "FreeRTOS.h"
#include "task.h"
#include "system.h"
#include "uart2.h"
#include "log.h"
#include "debug.h"

#define DEBUG_MODULE "CURVACE"

#define SERIAL_HEADER 'C'

static bool isInit;

typedef struct __attribute__((packed)) curvaceData_s {
  uint8_t header;
  uint8_t v1;
  uint8_t v2;
  uint8_t v3;
} curvaceData_t;


void curvaceTask(void *param)
{
  //int velXid;
  //int velYid;
  //int velZid;

  char c = 'A';

  curvaceData_t packet;

  systemWaitStart();

  // Setup data to transfer
  //velXid = logGetVarId("posCtl", "targetVX");
  //velYid = logGetVarId("posCtl", "targetVY");
  //velZid = logGetVarId("posCtl", "targetVZ");

  TickType_t lastWakeTime = xTaskGetTickCount();

  while (1)
  {
    // Set the loop unlock time in ms
    vTaskDelayUntil(&lastWakeTime, M2T(10));

    // Assemble the data
    packet.header = SERIAL_HEADER;
    packet.v1 = c++;
    if (c > 'Z') {
      c = 'A';
    }
    packet.v2 = 10;
    packet.v3 = 13;
    
    //packet.targetVX = logGetFloat(velXid);
    //packet.targetVY = logGetFloat(velYid);
    //packet.targetVZ = logGetFloat(velZid);

    // Send the three floats, byte by byte, to UART2
    uart2SendDataDmaBlocking(sizeof(curvaceData_t), (uint8_t *)(&packet));
  }
}


static void curvaceInit()
{
  DEBUG_PRINT("Starting CurvACE task: reading flow from  UART2\n");

  xTaskCreate(curvaceTask, "CURVACE", configMINIMAL_STACK_SIZE, NULL, 1, NULL);

  // Configure uart and set the baud rate
  /**
   * TODO: here or in cf-board.mk as flag?
   */
  uart2Init(115200);

  isInit = true;
}


static bool curvaceTest()
{
  return isInit;
}


static const DeckDriver curvaceDriver = {
  .name = "CurvACE",
  .init = curvaceInit,
  .test = curvaceTest,
};


DECK_DRIVER(curvaceDriver);
