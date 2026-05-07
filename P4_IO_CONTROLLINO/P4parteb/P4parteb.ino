#include <Controllino.h>
//PARTE B SEMAFORO
// 1. Definición de los estados del sistema (Restricción: enum)
enum EstadoSemaforo {
  A_VERDE_B_ROJO,
  A_AMARILLO_B_ROJO,
  A_ROJO_B_VERDE,
  A_ROJO_B_AMARILLO
};

// 2. Estructura para definir el comportamiento de cada estado (Restricción: struct)
struct ConfiguracionEstado {
  EstadoSemaforo estadoActual;
  unsigned long duracion;        // Tiempo en milisegundos
  EstadoSemaforo siguienteEstado;
};

// 3. Asignación de Pines (Semáforo A: D0-D2 | Semáforo B: D12-D14)
const int verdeA = CONTROLLINO_D0;
const int amarA  = CONTROLLINO_D1;
const int rojoA  = CONTROLLINO_D2;

const int verdeB = CONTROLLINO_D12;
const int amarB  = CONTROLLINO_D13;
const int rojoB  = CONTROLLINO_D14;

// 4. Tabla de la Máquina de Estados (FSM)
ConfiguracionEstado fsm[] = {
  {A_VERDE_B_ROJO,    2000, A_AMARILLO_B_ROJO}, // A fluye, B espera (2s)
  {A_AMARILLO_B_ROJO, 1000, A_ROJO_B_VERDE},    // A frena, B espera (1s)
  {A_ROJO_B_VERDE,    2000, A_ROJO_B_AMARILLO}, // A espera, B fluye (2s)
  {A_ROJO_B_AMARILLO, 1000, A_VERDE_B_ROJO}     // A espera, B frena (1s)
};

int indiceActual = 0;
unsigned long tiempoReferencia = 0;

void setup() {
  pinMode(verdeA, OUTPUT); 
  pinMode(amarA, OUTPUT); 
  pinMode(rojoA, OUTPUT);
  pinMode(verdeB, OUTPUT); 
  pinMode(amarB, OUTPUT); 
  pinMode(rojoB, OUTPUT);
  
  tiempoReferencia = millis(); // Inicializar cronómetro
}

void loop() {
  unsigned long tiempoActual = millis();

  // Lógica de transición de la FSM (Restricción: no bloqueante)
  if (tiempoActual - tiempoReferencia >= fsm[indiceActual].duracion) {
    
    // Buscamos el índice del siguiente estado
    EstadoSemaforo proximo = fsm[indiceActual].siguienteEstado;
    for(int i = 0; i < 4; i++) {
      if(fsm[i].estadoActual == proximo) {
        indiceActual = i;
        break;
      }
    }
    
    tiempoReferencia = tiempoActual; // Reiniciar cronómetro para el nuevo estado
  }

  // Actualización de las salidas físicas
  actualizarSalidas(fsm[indiceActual].estadoActual);
}

// Función para manejar el encendido/apagado de los LEDs
void actualizarSalidas(EstadoSemaforo estado) {
  // Primero apagamos todo por seguridad
  digitalWrite(verdeA, LOW); digitalWrite(amarA, LOW); digitalWrite(rojoA, LOW);
  digitalWrite(verdeB, LOW); digitalWrite(amarB, LOW); digitalWrite(rojoB, LOW);

  // Encendemos según el estado actual de la FSM
  switch (estado) {
    case A_VERDE_B_ROJO:
      digitalWrite(verdeA, HIGH);
      digitalWrite(rojoB, HIGH);
      break;
    case A_AMARILLO_B_ROJO:
      digitalWrite(amarA, HIGH);
      digitalWrite(rojoB, HIGH);
      break;
    case A_ROJO_B_VERDE:
      digitalWrite(rojoA, HIGH);
      digitalWrite(verdeB, HIGH);
      break;
    case A_ROJO_B_AMARILLO:
      digitalWrite(rojoA, HIGH);
      digitalWrite(amarB, HIGH);
      break;
  }
}
