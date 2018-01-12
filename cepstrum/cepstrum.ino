/** Cepstrum Pitch Detection
  The cepstrum analysis involves DFT and iDFT. It is suitable for complex audio
  inputs such as human voices, but musical notes may not work as well. This
  program may not be accurate, due to lack of memory on the Arduino. Ideally,
  there should be 1024-4096 data samples.
*/

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

#define SCAN_SPEED 8//scans per millisecond, maximum is 10
#define SCAN_TIME 20

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
//create a ariety of special characters

double** DFT(double* frequencies, double* wave, int N) { //scanTime is in millis
  double real = 0;
  double imag = 0;

  for (int k = 0; k < N; k++) {
    for (int i = 0; i < N; i++) {
      real += wave[i] * cos(-2 * PI * k * i / N);
      imag += wave[i] * sin(-2 * PI * k * i / N);
    }

    frequencies[k] = real * real + imag * imag;
    real = 0;
    imag = 0;
  }
  //return frequencies;
}

double* iDFT(double* frequencies, double *wave, int N) {
  double real = 0;
  double imag = 0;

  for (int k = 0; k < N; k++) {
    for (int i = 0; i < N; i++) {
      real += wave[i] * cos(2 * PI * k * i / N);
      imag += wave[i] * sin(2 * PI * k * i / N);

    }
    real /= N;
    imag /= N;
    frequencies[k] = sqrt(real * real + imag * imag);
    real = 0;
    imag = 0;
  }
}

int findFrequency(double* wave, int N) {
  double temp[N];

  //step 1 of cepstrum
  DFT(temp, wave, N);
  //from here on, only half of the values should be used.

  //step 2 of cepstrum
  for (int i = 0; i < N / 2; i++) {
    temp[i] = log10(temp[i]);
  }

  //step 3 of cepstrum
  iDFT(wave, temp, N / 2);

  //peak finding
  //start at 10 for high time liftering
  double peak = 0;
  int peakIndex = 9;

  for (int i = 10; i < (N / 2) - 1; i++) {
    if (peak <= wave[i] && (wave[i] > wave[i - 1]) && (wave[i] > wave[i + 1])) {
      peakIndex = i;
      peak = wave[i];
    }
  }

  free(temp);
  return round(1000 * SCAN_SPEED / peakIndex);
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

  lcd.createChar(0, smiley);
  lcd.createChar(1, waveR1);
  lcd.createChar(2, waveR2);
  lcd.createChar(3, waveL1);
  lcd.createChar(4, waveL2);
  lcd.createChar(5, note);
}

void loop() {
  // put your main code here, to run repeatedly:
  int freq;
  int N = SCAN_SPEED * SCAN_TIME;
  double *wave = malloc(N * 2 * sizeof(double));

  for (int i = 0; i < N; i++) {
    wave[i] = (double)5 * analogRead(0) / 1024;
    delay(1 / SCAN_SPEED);
  }

  //hanning window function
  for (int i = 0; i < N; i++) {
    wave[i] *= (0.54 + (double)0.46 * cos(2 * PI * i / (N - 1)));
  }

  freq = findFrequency(wave, N);
  free(wave);

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
