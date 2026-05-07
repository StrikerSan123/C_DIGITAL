#include <ArduinoRS485.h>
#include <ArduinoModbus.h>
#include "Controllino.h"

// Variables para el LED 1 (D0)
bool estadoHabilitado1 = true; //empieza o no encendido.
bool ultimoEstadoBtn1 = LOW; 

// Variables para el LED 2 (CAMBIADO A D12)
bool estadoHabilitado2 = true;
bool ultimoEstadoBtn2 = LOW;

// Definición de pines
const int pinLED1 = CONTROLLINO_D0;
const int pinLED2 = CONTROLLINO_D6; 

void setup() {
  ModbusRTUServer.begin(1, 19200);
  ModbusRTUServer.configureHoldingRegisters(0x00, 2); //dice cuantos registros se reservan
 
  // Valores iniciales HMI
  ModbusRTUServer.holdingRegisterWrite(0, 50); // LED 1 al 50%
  ModbusRTUServer.holdingRegisterWrite(1, 50); // LED 2 al 50%
  //salida para los leds
  pinMode(pinLED1, OUTPUT);
  pinMode(pinLED2, OUTPUT);
  
  // Entradas para los pulsantes
  pinMode(CONTROLLINO_I16, INPUT);
  pinMode(CONTROLLINO_I17, INPUT);
}

void loop() {
  ModbusRTUServer.poll();

  // --- LÓGICA BOTÓN 1 (I16) LED 1
  bool lectura1 = digitalRead(CONTROLLINO_I16);
  if (lectura1 == HIGH && ultimoEstadoBtn1 == LOW) {
    estadoHabilitado1 = !estadoHabilitado1;
    delay(50); // antirrebote
  }
  ultimoEstadoBtn1 = lectura1;

  if (estadoHabilitado1) {
    int val1 = ModbusRTUServer.holdingRegisterRead(0);
    analogWrite(pinLED1, map(val1, 0, 100, 0, 255));
  } else {
    analogWrite(pinLED1, 0);
  }


  // --- LÓGICA BOTÓN 2 (I17) LED 2
  bool lectura2 = digitalRead(CONTROLLINO_I17);
  if (lectura2 == HIGH && ultimoEstadoBtn2 == LOW) {
    estadoHabilitado2 = !estadoHabilitado2; 
    delay(50);  //antirebote
  }
  ultimoEstadoBtn2 = lectura2;

  if (estadoHabilitado2) {
    // Sigue leyendo el registro 40002 de la HMI
    int val2 = ModbusRTUServer.holdingRegisterRead(1); 
    analogWrite(pinLED2, map(val2, 0, 100, 0, 255));
  } else {
    analogWrite(pinLED2, 0);
  }
}