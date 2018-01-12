#include <Wire.h>
#include <LCD.h>
#include <LiquidCrystal_I2C.h>
#include <math.h>

#define I2C_ADDR    0x3F
#define BACKLIGHT_PIN     3
#define En_pin  2
#define Rw_pin  1
#define Rs_pin  0
#define D4_pin  4
#define D5_pin  5
#define D6_pin  6
#define D7_pin  7
#define SCAN_SPEED 8 // scans per millisecond
#define SCAN_TIME 8

LiquidCrystal_I2C  lcd(I2C_ADDR, En_pin, Rw_pin, Rs_pin, D4_pin, D5_pin, D6_pin, D7_pin);

byte smiley[8] = {
  B00000,
  B10001,
  B00000,
  B00100,
  B10001,
  B01110,
  B00000,
};
byte waveR1[8] {
  B00000,
  B00110,
  B11001,
  B00000,
  B00000,
  B00000,
  B10000,
  B01000,
};
byte waveR2[8] {
  B00000,
  B00000,
  B00000,
  B00000,
  B11001,
  B00110,
  B10000,
  B10000,
};
byte waveL1[8] {
  B00000,
  B00000,
  B00000,
  B00000,
  B10011,
  B01100,
  B00001,
  B00001,
};
byte waveL2[8] {
  B00000,
  B01100,
  B10011,
  B00000,
  B00000,
  B00000,
  B00001,
  B00010,
};
byte note[8] {
  B00001,
  B00011,
  B00101,
  B01001,
  B01001,
  B01011,
  B11011,
  B11000,

};

double** DFT(double frequencies[][2], int* wave, int N) { //scanTime is in millis
  double real = 0;
  double imag = 0;

  for (int k = 0; k < N; k++) {
    Serial.print(wave[k]);
    for (int i = 0; i < N; i++) {
      real += wave[i] * cos(-2 * PI * k * i / N);
      imag += wave[i] * sin(-2 * PI * k * i / N);
    }

    frequencies[k][0] = real;
    frequencies[k][1] = imag;
    real = 0;
    imag = 0;
  }
  //return frequencies;
}


double* iDFT(double frequencies[][2], int timeDomain, int N) {
  double *samples = malloc(sizeof(double) * (timeDomain));
  double real = 0;
  double imag = 0;

  for (int k = 0; k < timeDomain; k++) {
    for (int i = 0; i < N; i++) {
      //calculate real and imaginary parts separately
      real = real + frequencies[i][0] * cos(2 * PI * k * i / N) - frequencies[i][1] * sin(2 * PI * k * i / N);
      //imag = imag + frequencies[i][1]*cos(2*PI*k*i/N) + frequencies[i][0]*sin(2*PI*k*i/N);   //save time since imaginary is not used
    }
    real /= N;
    //imag /= N;   //imaginary is not used
    samples[k] = real;
    real = 0;
    //imag = 0;   //imaginary is not used
  }
  return samples;
}

int findFrequency(int *wave, int N) {
  double temp[N * 2][2];

  //Apply FFT, set first element to 0
  DFT(temp, wave, N * 2);
  temp[0][0] = 0;
  temp[0][1] = 0;

  //Take magnitudes of all elements
  for (int i = 0; i < N * 2; i++) {
    temp[i][0] = temp[i][0] * temp[i][0] + temp[i][1] * temp[i][1];
    temp[i][1] = 0;
  //Serial.print(temp[0][0]);
  //Serial.print(" ");
  //Serial.print(temp[0][1]);
  }

  free(wave);       //must free here, or else insufficient memory

  double* autoCorrSamples = iDFT(temp, N, N * 2);
  for (int i = 1; i < N; i++) {
    autoCorrSamples[i] = autoCorrSamples[i] / autoCorrSamples[0];
    printf("%f\n", autoCorrSamples[i]);
  }

  //Find suitable peak - distance from highest peak to 0 is used to calculate frequency
  double peak = 0;
  int peakIndex = 0;

  for (int i = 1; i < N - 1; i++) {
    if (peak <= autoCorrSamples[i] && (autoCorrSamples[i] > autoCorrSamples[i - 1]) && (autoCorrSamples[i] > autoCorrSamples[i + 1])) {
      peakIndex = i;
      peak = autoCorrSamples[i];
    }
  }

  free(autoCorrSamples);
  return round(SCAN_SPEED * 1000 / peakIndex);
}

void setup() {
  // put your setup code here, to run once:
  lcd.begin (20, 4);
  lcd.clear();
  Serial.begin(9600);

  // Switch on the backlight
  lcd.setBacklightPin(BACKLIGHT_PIN, POSITIVE);
  lcd.setBacklight(HIGH);
  lcd.setCursor(0, 0);
}

void loop() {
  // put your main code here, to run repeatedly:
  int N = SCAN_SPEED * SCAN_TIME;
  int freq;
  int *wave = malloc(N * 2 * sizeof(int));

  for (int i = 0; i < N; i++) {
    wave[i] =  analogRead(0); //convert to volts
    delay(1 / SCAN_SPEED);
  }

  for (int i = 0; i < N; i++) {
    wave[i + N] = 0;                       //zero padding
    wave[i] *= (0.54 + (double)0.46 * cos(2 * PI * i / (N - 1))); //hanning window
  }

  freq = findFrequency(wave, N);

  //display frequency onto the lcd screen
  lcd.clear();
  double fA = 440.0; // Base frequency
  double n = 12 * ((log10(freq) / log10(2.0)) - (log10(fA) / log10(2.0))); //formula to find the number of half steps from A4
  int hStep = n > 0 ? round(n) % 12 : (12 + (round(n) % 12)) == 12 ? 0 : (12 + (round(n) % 12));
  int sharp = 0, flat = 0, perf = 0;
  double distance = n - round(n);
  //variable declarations

  lcd.setCursor(0, 0);
  lcd.write(byte(5));
  lcd.write(byte(5));
  lcd.print ("NOTE TRANSCRIBER");
  lcd.write(byte(5));
  lcd.write(byte(5));
  //Print out first line
  lcd.setCursor(0, 1);
  lcd.print("Note: ");

  switch (hStep) {
    case 0:
      lcd.print("A");
      break;
    case 1:
      lcd.print("A#/Bb");
      break;
    case 2:
      lcd.print("B");
      break;
    case 3:
      lcd.print("C");
      break;
    case 4:
      lcd.print("C#/Db");
      break;
    case 5:
      lcd.print("D");
      break;
    case 6:
      lcd.print("D#/Eb");
      break;
    case 7:
      lcd.print("E");
      break;
    case 8:
      lcd.print("F");
      break;
    case 9:
      lcd.print("F#/Gb");
      break;
    case 10:
      lcd.print("G");
      break;
    case 11:
      lcd.print("G#/Ab");
      break;
  }
  //displays note played
  lcd.setCursor(0, 2);
  lcd.print("Freq: ");
  lcd.print(freq);
  lcd.print ("Hz");
  //displays frequency

  if (abs(distance) < 0.05)
    perf = 1;
  else if (distance > 0)
    sharp = 1;
  else
    flat = 1;
  lcd.setCursor(9, 3);
  lcd.print (sharp);
  lcd.print (flat);

  if (perf) {
    /*
      lcd.setCursor(10,2);
      lcd.write(byte(0));

      lcd.setCursor(11,2);
      lcd.write(byte(1));
      lcd.setCursor(9,2);
      lcd.write(byte(3));
      delay (300);
      lcd.setCursor(11,2);
      lcd.write(byte(2));
      lcd.setCursor(9,2);
      lcd.write(byte(4));
      delay (300);
    */
    lcd.setCursor(9, 3);
    lcd.write(byte(0));
    lcd.write(byte(0));
  } else if (sharp) {
    lcd.setCursor(9, 3);
    lcd.print("01");
    lcd.setCursor(11, 3);
    for (int i = 0; i < round(distance * 18); i++) {
      lcd.print ('#');
    }
  } else {
    lcd.setCursor(9, 3);
    lcd.print("10");
    lcd.setCursor(8, 3);
    lcd.rightToLeft();
    for (int i = 0; i < round(distance * -18); i++) {
      lcd.print ('b');
    }
    lcd.leftToRight();
  }
}
