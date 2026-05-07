/*
  ============================================================
  CONTROL DE MOTOR DC + MEDICIÓN DE RPM + HMI COOLMAY (MODBUS)
  ============================================================

  Objetivo:
  - El HMI (Coolmay) envía el Duty Cycle (0 a 100 %) por Modbus.
  - El Controllino convierte ese Duty a PWM (0 a 255) y lo aplica en D0.
  - Se mide la velocidad del motor contando pulsos (encoder).
  - Cada 1 segundo se calcula la velocidad en RPM.
  - Se envían las RPM al HMI por Modbus para verlas en un Trend Graph.

  Comunicación:
  - Modbus RTU por RS485
  - Slave ID = 1
  - Baudrate = 19200

  Holding Registers (16 bits):
  - HR0: RPM (0..6000 aprox.)  -> para mostrar en el HMI
  - HR1: DutyCycle (0..100)    -> leído desde el HMI (slider/scrollbar)
*/

#include <ArduinoRS485.h>
#include <ArduinoModbus.h>
#include "Controllino.h"

// ------------------------
// Pines de hardware
// ------------------------
const int pin_motor = CONTROLLINO_D0;     // Salida PWM al driver del motor
const int entrada   = CONTROLLINO_IN1;    // Entrada de pulsos (encoder)

// ------------------------
// Medición por pulsos
// ------------------------
volatile unsigned long conteo_pulsos = 0; // se incrementa en interrupción externa

const uint16_t PULSOS_POR_REV = 36;       // Pulsos por cada vuelta del eje
const float fs = 1.0f;                    // Frecuencia de muestreo (1 Hz => cada 1 s)

// RPM como entero (se calcula en el ISR del timer y se lee en loop)
volatile uint16_t rpm = 0;

// Prototipo de ISR externa
void contarPulso();

void setup() {
  Serial.begin(115200);

  // ------------------------
  // Iniciar Modbus RTU Server
  // ------------------------
  if (!ModbusRTUServer.begin(1, 19200)) {
    while (1); // si falla, se queda aquí
  }

  // -----------------------------------------
  // Configurar Holding Registers del servidor
  // -----------------------------------------
  // 2 registros desde la dirección 0:
  // HR0 -> RPM
  // HR1 -> DutyCycle
  ModbusRTUServer.configureHoldingRegisters(0x00, 2);

  // Valores iniciales
  ModbusRTUServer.holdingRegisterWrite(0, 0); // Holding register 0 con RPM
  ModbusRTUServer.holdingRegisterWrite(1, 0); // Holding register 1 con DutyCycle

  // ------------------------
  // Configurar pines
  // ------------------------
  pinMode(pin_motor, OUTPUT);

  // Pull-up interno 
  pinMode(entrada, INPUT_PULLUP);

  // Interrupción externa: contar pulsos por flanco de bajada
  attachInterrupt(digitalPinToInterrupt(entrada), contarPulso, FALLING);

  // --------------------------------------------
  // Timer1 en modo CTC para interrumpir cada 1 s
  // --------------------------------------------
  // Nota: estos registros son del ATmega2560 (Controllino Mega).
  noInterrupts();

  TCCR1A = 0;               // Modo normal en A
  TCCR1B = 0;               // Apagado para configurar
  TCCR1B |= B00000100;      // Prescaler = 256 (CS12 = 1)
  TIMSK1 |= B00000010;      // Habilitar interrupción por Compare Match A
  OCR1A = 62500 / fs;       // 1 segundo si F_CPU=16MHz y prescaler=256
  TCNT1 = 0;                // Reset contador

  interrupts();
}

void loop() {
  // Atender peticiones Modbus del HMI
  ModbusRTUServer.poll();

  // --------------------------------------------
  // 1) Leer DutyCycle desde el HMI (HR1)
  // --------------------------------------------
  uint16_t duty = ModbusRTUServer.holdingRegisterRead(1);

  // Seguridad: limitar a 0..100
  duty = constrain(duty, 0, 100);

  // Convertir duty(0..100) a PWM(0..255)
  uint8_t pwmValue = (uint8_t)((uint32_t)duty * 255 / 100);

  // Aplicar PWM al motor
  analogWrite(pin_motor, pwmValue);

  // --------------------------------------------
  // 2) Enviar RPM al HMI (HR0)
  // --------------------------------------------
  // Copia “atómica” porque rpm se actualiza en una ISR:
  uint16_t rpm_local;
  noInterrupts();
  rpm_local = rpm;
  interrupts();

  ModbusRTUServer.holdingRegisterWrite(0, rpm_local);

  // --------------------------------------------
  // 3) Mostrar por Serial (opcional para depurar)
  // --------------------------------------------
  Serial.print("Duty(%): ");
  Serial.print(duty);
  Serial.print("  PWM: ");
  Serial.print(pwmValue);
  Serial.print("  RPM: ");
  Serial.println(rpm_local);
}

// ============================================================
// ISR Timer1: se ejecuta cada 1 segundo y calcula las RPM
// ============================================================
ISR(TIMER1_COMPA_vect) {
  TCNT1 = 0; // reset del timer

  // Copiar y resetear pulsos
  unsigned long p = conteo_pulsos;
  conteo_pulsos = 0;

  // RPM = (pulsos/segundo) / (pulsos_por_rev) * 60
  float rpm_calc = (float)p * 60.0f * fs / (float)PULSOS_POR_REV;

  // Limitar para que siempre quepa en uint16_t y sea razonable
  if (rpm_calc < 0) rpm_calc = 0;
  if (rpm_calc > 6000) rpm_calc = 6000;

  rpm = (uint16_t)rpm_calc;
}

// ============================================================
// ISR externa: se ejecuta en cada pulso del sensor/encoder
// ============================================================
void contarPulso() {
  conteo_pulsos++;
}