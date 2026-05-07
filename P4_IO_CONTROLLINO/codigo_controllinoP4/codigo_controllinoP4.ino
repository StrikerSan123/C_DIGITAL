#include <Controllino.h>
//PARTE A SECUENCIA LEDS

// Asignacion Leds
const int led1 = CONTROLLINO_D0;
const int led2 = CONTROLLINO_D6;
const int led3 = CONTROLLINO_D12;
const int led4 = CONTROLLINO_D13;
const int led5 = CONTROLLINO_D14;
const int led6 = CONTROLLINO_D8;
const int led7 = CONTROLLINO_D2;
const int led8 = CONTROLLINO_D1;
const int led9 = CONTROLLINO_D7;

// Asignacion Botones
const int boton1 = CONTROLLINO_I16;
const int boton2 = CONTROLLINO_I17;
const int boton3 = CONTROLLINO_I18;

// Secuencias
const int espiralNormal[] = {led1, led2, led3, led4, led5, led6, led7, led8, led9};
const int espiralInversa[] = {led9, led8, led7, led6, led5, led4, led3, led2, led1};

// Puntero para la secuencia
const int* secuenciaActiva = nullptr; 
int pasoActual = 0;
const int totalPasos = 9;

// Variables para millis();
unsigned long tiempoAnterior = 0;
const long intervalo = 250; 
bool ejecutando = false;

void setup() {
  pinMode(led1, OUTPUT);
  pinMode(led2, OUTPUT);
  pinMode(led3, OUTPUT);
  pinMode(led4, OUTPUT);
  pinMode(led5, OUTPUT);
  pinMode(led6, OUTPUT);
  pinMode(led7, OUTPUT);
  pinMode(led8, OUTPUT);
  pinMode(led9, OUTPUT);
  
  pinMode(boton1, INPUT);
  pinMode(boton2, INPUT);
  pinMode(boton3, INPUT);
}


void loop() {
  // Lectura de botones
  if (digitalRead(boton1) == HIGH) {
    //Logica para conservar la "localizacion" al cambiar el sentido de giro
    if (secuenciaActiva == espiralInversa) {
    pasoActual = (totalPasos - 1) - pasoActual;
    }
    secuenciaActiva = espiralNormal;
    ejecutando = true;
    tiempoAnterior = millis(); // Reinicia el tiempo para que empiece de inmediato
  }

  if (digitalRead(boton2) == HIGH) {
    // Similar a la del boton 1, invierte el sentido manteniendo posicion
    if (secuenciaActiva == espiralNormal) {
    pasoActual = (totalPasos - 1) - pasoActual;
    }
    secuenciaActiva = espiralInversa;
    ejecutando = true;
    tiempoAnterior = millis();
  }

  if (digitalRead(boton3) == HIGH) {
    apagarTodo();
    ejecutando = false;
    secuenciaActiva = nullptr;
  }

  // Lógica de encendido no bloqueante
  if (ejecutando && secuenciaActiva != nullptr) {
    if (millis() - tiempoAnterior >= intervalo) {
      tiempoAnterior = millis();

      apagarTodo(); // Apaga el anterior
      
      // Acceso al pin mediante el puntero
      digitalWrite(*(secuenciaActiva + pasoActual), HIGH);
      
      pasoActual++;
      if (pasoActual >= totalPasos) {
        pasoActual = 0; // Reinicia la espiral
      }
    }
  }
}

void apagarTodo() {
  for (int i = 0; i < totalPasos; i++) {
    digitalWrite(espiralNormal[i], LOW);
  }
}
