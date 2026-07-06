int las2 = 4;
int signal = 0;
int cam = 2;
int start = 3;

int pot = 0;
int tid = 0;
int i =0;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(las2, INPUT);                               // Ingång för Laser2 Pin4 Plint 3
  pinMode(cam, OUTPUT);                               // Utgång till Reläingång 1 till kamera Plint 5 och 6
  pinMode(start, OUTPUT);                             // Utgång till Reläingång 2 till footswitch plint 4
  digitalWrite(start, HIGH);                          // Relä 2 ej aktiv
  digitalWrite(cam, HIGH);                            // Relå 1 ej aktiv

}

void loop() {
  i = 1;
  pot = analogRead(A0);
  tid = pot * (18000 / 1023.0);  
 //Serial.println(tid);

  //delay(1500);
  digitalWrite(start, LOW);
  delay(500);
  digitalWrite(start, HIGH);  

  for (i = 1; i <= tid; i=i+1) {                //fördröjning i antal sekunder som läst in Potentiometer
    Serial.println(i);
    delay(200);
    signal = digitalRead(las2);                          //Läs in aktivitet av Laser2 Plint 3
    if (signal == HIGH) {
      digitalWrite(cam, LOW);
      delay(100);
      digitalWrite(cam, HIGH);

do {
    signal = digitalRead(las2); 
} while (signal == HIGH);

      
      
    }

  }      
}
