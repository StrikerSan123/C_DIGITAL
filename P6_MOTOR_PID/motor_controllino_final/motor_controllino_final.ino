#include <ArduinoRS485.h>
#include <ArduinoModbus.h>
#include "Controllino.h"

// --- Pines de hardware ---
const int pin_motor = CONTROLLINO_D0;
const int entrada   = CONTROLLINO_IN1;
const uint16_t PULSOS_POR_REV = 36;

// --- Variables de Medición ---
volatile unsigned long conteo_pulsos = 0;
double rpm_actual = 0;

// --- Variables PID Discreto ---
float e[3] = {0, 0, 0}; // Error actual, anterior, tras-anterior
float u = 0;            // Acción de control acumulada (PWM)
float ref = 0;          // Setpoint

// Constantes (Tiempos en segundos)
float Kp = 0.1; 
float Ti = 0.150; 
float Td = 0.1; 
const float ts = 0.05;   // Tiempo de muestreo: 50ms

void setup() {
  Serial.begin(115200);

  if (!ModbusRTUServer.begin(1, 19200)) {
    while (1); 
  }

  // HR0:RPM | HR1:Ref | HR2:PWM | HR3:Kp | HR4:Ti | HR5:Td
  ModbusRTUServer.configureHoldingRegisters(0x00, 6);

  // Inicializar valores en HMI (Ti y Td escalados x1000 para precisión)
  ModbusRTUServer.holdingRegisterWrite(3, (uint16_t)(Kp * 100));
  ModbusRTUServer.holdingRegisterWrite(4, (uint16_t)(Ti * 1000));
  ModbusRTUServer.holdingRegisterWrite(5, (uint16_t)(Td * 1000));

  pinMode(pin_motor, OUTPUT);
  pinMode(entrada, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(entrada), contarPulso, FALLING);

  // Timer1 para 50ms
  noInterrupts();
  TCCR1A = 0; TCCR1B = 0;
  TCCR1B |= (1 << WGM12) | (1 << CS12); 
  OCR1A = 3125; 
  TIMSK1 |= (1 << OCIE1A);
  interrupts();
}

void loop() {
  ModbusRTUServer.poll();

  // 1. Leer parámetros del HMI
  ref = (float)ModbusRTUServer.holdingRegisterRead(1);
  
  // Actualizar constantes dinámicamente
  Kp = (float)ModbusRTUServer.holdingRegisterRead(3) / 100.0;
  Ti = (float)ModbusRTUServer.holdingRegisterRead(4) / 1000.0;
  Td = (float)ModbusRTUServer.holdingRegisterRead(5) / 1000.0;
 Serial.print("PWM: ");
Serial.print(u);

Serial.print(" | REF: ");
Serial.print(ref);

Serial.print(" | RPM: ");
Serial.println(rpm_actual);
  // 2. Aplicar salida física (calculada en la ISR)
  analogWrite(pin_motor, (int)u);
//analogWrite(pin_motor, 254);
  // 3. Feedback al HMI
  ModbusRTUServer.holdingRegisterWrite(0, (uint16_t)rpm_actual);
  ModbusRTUServer.holdingRegisterWrite(2, (uint16_t)u);
}

void contarPulso() {
  conteo_pulsos++;
}

// ============================================================
// CÁLCULO DEL CONTROLADOR (Sincronizado a 100ms)
// ============================================================
ISR(TIMER1_COMPA_vect) {
  // A. Calcular velocidad actual
  unsigned long p = conteo_pulsos;
  conteo_pulsos = 0;
 float instant_rpm = (60.0 * p) / (PULSOS_POR_REV * ts);
  rpm_actual = (rpm_actual * 0.7) + (instant_rpm * 0.3); // Filtro

  // B. Algoritmo PID Discreto
  e[2] = e[1];
  e[1] = e[0];
  e[0] = ref - (float)rpm_actual;

  if (ref > 0) {
    // Ecuación de recurrencia que proporcionaste
    float du = Kp * ((e[0] - e[1]) + (ts / Ti) * e[0] + (Td / ts) * (e[0] - 2 * e[1] + e[2]));
    u += du;

    // Saturación de la salida (0-255 para PWM)
    if (u > 255) u = 255;
    if (u < 0)   u = 0;
  } else {
    u = 0;
    e[0] = 0; e[1] = 0; e[2] = 0; // Reset de errores
  }
}