#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define MAX_HEIGHT 1000
#define SCAN_SPEED 32 //ticks per millisecond (mHZ)
#define PI 3.1415926
#define MAX_PEAKS 6
#define MIN_FREQ 100
#define MAX_FREQ 1000

int randMax(int max){
    return rand()%max;
}

//generate a wave from a combination of sine waves of different frequencies
//frequencies is in Hz
//scanTime is the total amount of time the wave lasts for, in milliseconds
int generateWave(int n, int* frequencies, int* amplitudes, int totalTime, double* wave, int totalLength){
    int totalSamples = totalTime*SCAN_SPEED;
    int lengths[n]; //length of each wave in terms of array size
    for (int i = 0; i < n; i++){
        lengths[i] = totalLength*1000/frequencies[i]/totalTime;
    }

    for (int i = 0; i < totalSamples; i++){
        wave[i] = 0;
    }

    //generate waves
    for (int i = 0; i < totalLength; i++){
        for (int j = 0; j < n; j++){
            wave[i]+= amplitudes[j]*sin(2*PI*(i%lengths[j])/(lengths[j]+0.0));
        }
    }

    return 0;
}

//scanTime is in milliseconds
double** DFT(double frequencies[][2], double* wave, int N){
    double a = 0; //a is real part, b is imaginary part
    double b = 0;

    //DFT algorithm
    for (int k = 0; k < N; k++){
        for (int i = 0; i < N; i++){
            a += wave[i]*cos(-2*PI*k*i/N);
            b += wave[i]*sin(-2*PI*k*i/N);
        }

        frequencies[k][0] = a;
        frequencies[k][1] = b;
        a = 0;
        b = 0;
    }
    return frequencies;
}

//iterative version
double* iDFT(double frequencies[][2], int timeDomain, int N){
    double *samples = malloc(sizeof(double)*(timeDomain));
    double a = 0; //a is real, b is imaginary
    double b = 0;

    for (int k = 0; k < timeDomain; k++){
        for (int i = 0; i < N; i++){
            a = a + frequencies[i][0]*cos(2*PI*k*i/N) - frequencies[i][1]*sin(2*PI*k*i/N);
            b = b + frequencies[i][1]*cos(2*PI*k*i/N) + frequencies[i][0]*sin(2*PI*k*i/N);
        }
        a /= N;
        b /= N;
        samples[k] = a;
        a = 0;
        b = 0;
    }
    return samples;
}

//instead of making a new wave array every time, we use leap and offset to simulate this
//frequencies is a N by 2 array
//assumes real only input
int FFT(double frequencies[][2], double* wave, int N, int leap, int offset){
    //base case
    if (N == 1){
        frequencies[0][0] = wave[offset];
        frequencies[0][1] = 0;
        return 0;
    }else{
        double* eFreq = (double*)malloc(N*sizeof(double));
        double* oFreq = (double*)malloc(N*sizeof(double));
        double temp1;
        double temp2;

        //take fft of even and odd indexed members
        FFT(eFreq, wave, N/2, leap*2, offset);
        FFT(oFreq, wave, N/2, leap*2, offset+leap);

        //recombine
        for (int j = 0; j < N/2; j++){
            temp1 = cos(-2*PI*j/N)**(oFreq+j*2) - sin(-2*PI*j/N)**(oFreq+j*2+1);
            temp2 = cos(-2*PI*j/N)**(oFreq+j*2+1) + sin(-2*PI*j/N)**(oFreq+j*2);

            frequencies[j][0] = *(eFreq+j*2) + temp1;
            frequencies[j+N/2][0] = *(eFreq+j*2) - temp1;
            frequencies[j][1] = *(eFreq+j*2+1) + temp2;
            frequencies[j+N/2][1] = *(eFreq+j*2+1) - temp2;
        }
        free(eFreq);
        free(oFreq);
    }
    return 0;
}

//instead of making a new wave array every time, we use leap and offset to simulate this
//frequencies is a N by 2 array of input
//answer is a N by 2 array for output
//assumes input with real and imaginary parts
int iFFT(double frequencies[][2], double answer[][2], int N, int leap, int offset){
    //base case
    if (N == 1){
        answer[0][0] = frequencies[offset][0];
        answer[0][1] = frequencies[offset][1];
        return 0;
    }else{
        double* eAns = (double*)malloc(N*sizeof(double));
        double* oAns = (double*)malloc(N*sizeof(double));
        double temp1;
        double temp2;

        //take ifft of even and odd indexed members
        iFFT(frequencies, eAns, N/2, leap*2, offset);
        iFFT(frequencies, oAns, N/2, leap*2, offset+leap);

        //recombine
        for (int j = 0; j < N/2; j++){
            temp1 = cos(2*PI*j/N)**(oAns+j*2) - sin(2*PI*j/N)**(oAns+j*2+1);
            temp2 = cos(2*PI*j/N)**(oAns+j*2+1) + sin(2*PI*j/N)**(oAns+j*2);

            answer[j][0] = (*(eAns+j*2) + temp1);
            answer[j+N/2][0] = (*(eAns+j*2) - temp1);
            answer[j][1] = (*(eAns+j*2+1) + temp2);
            answer[j+N/2][1] = (*(eAns+j*2+1) - temp2);

            //divide final answers by N
            if (leap == 1){
                answer[j][0] /= N;
                answer[j+N/2][0] /= N;
                answer[j][1] /= N;
                answer[j+N/2][1] /= N;
            }
        }
        free(eAns);
        free(oAns);
    }
}

//tester function - compares DFT and FFT output
int FFTCompare(double* wave, int N){
    double freq1[N][2];
    double freq2[N][2];
    DFT(freq1,wave,N);
    FFT(freq2,wave,N,1,0);

    for (int i = 0; i < N; i++){
        printf("%f %f\n",freq1[i][0], freq2[i][0]);
    }
}

//tester function - compares original wave and iFFT(FFT(wave))
int testFFT(double* wave, int N){
    double freq[N][2];

    FFT(freq,wave,N,1,0);
    double samples [N][2];
    iFFT(freq, samples, N, 1,0);

    for (int i = 0; i < N; i++){
        printf("%f %f\n",wave[i],samples[i][0]);
    }
    free(samples);
}

double findFrequency(double *wave, int N){
    double* freq = malloc(N*4*sizeof(double));

    //Apply FFT, set first element to 0
    FFT(freq,wave,N*2,1,0);
    *freq = 0;
    *(freq+1) = 0;

    //Take magnitudes of all elements
    for (int i = 0; i < N*2; i++){
        *(freq+i*2) = *(freq+i*2)**(freq+i*2)+*(freq+i*2+1)**(freq+i*2+1);
        *(freq+i*2+1) = 0;
    }

    //Apply iFFT
    double samples[N*2][2];
    iFFT(freq, samples, N*2, 1, 0);
    for (int i = 1; i < N; i++){
        samples[i][0] = samples[i][0]/samples[0][0];
    }

    //Find suitable peak - distance from highest peak to 0 is used to calculate frequency
    double max = 0;
    int ind = 0;

    for (int i = 1; i < N-1; i++){
        if(max <= samples[i][0] && (samples[i][0] > samples[i-1][0]) && (samples[i][0] > samples[i+1][0])){
            ind = i;
            max = samples[i][0];
        }
    }

    free(samples);
    return SCAN_SPEED*1000/ind;
}

int main(){
    int scanTime = 16; //in milliseconds
    int n = 1;
    int frequencies[3] = {320,200,400};
    int amplitudes[3] = {30,30,40};
    int max=1;
    int N = scanTime*SCAN_SPEED;

    // constant coefficients
    double a = 0.54;
    double b = 1-a;

    double *wave = malloc(N*2*sizeof(double));
    generateWave(n,frequencies,amplitudes,scanTime,wave,N);


    for (int i = 0; i < N; i++){
        wave[i+N] = 0;                         //zero pad
        wave[i] *= (a + b*cos(2*PI*i/(N-1)));  //hanning window
    }
    printf("Frequency: %f\n",findFrequency(wave,N));

    free(wave);
    return 0;
}




