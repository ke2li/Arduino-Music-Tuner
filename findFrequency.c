#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define MAX_HEIGHT 1000
#define TICK_SPEED 10 //ticks per millisecond
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
int generateWave2(int n, int* frequencies, int* amplitudes, int totalTime, double* wave, int totalLength){
    int totalSamples = totalTime*TICK_SPEED;
    int lengths[n];
    for (int i = 0; i < n; i++){
        lengths[i] = totalLength*1000/frequencies[i]/totalTime;
    }

    for (int i = 0; i < totalSamples; i++){
        wave[i] = 0;
    }

    for (int i = 0; i < totalLength; i++){
        for (int j = 0; j < n; j++){
            wave[i]+= amplitudes[j]*sin(2*PI*(i%lengths[j])/(lengths[j]+0.0));
            //printf("%f\n",wave[i]);
        }
    }

    return 0;
}

int* generateWave(int frequency, int variation, int length, int *vals){ //frequency in Hz, length in milliseconds of scan time
    srand(time(NULL));
    int numVals = 1000*TICK_SPEED/frequency;
    int peaks = 1+randMax(MAX_PEAKS-1); //number of troughs is same as number of peaks

    int peakVals[peaks];
    int peakTimes[peaks];
    int troughVals[peaks];
    int troughTimes[peaks];

    int frontEnd = 0;
    for (int i = 0; i < peaks; i++){ //generate times at which peaks and troughs occur
        peakTimes[i] = frontEnd + randMax(numVals-1-frontEnd)/(peaks);
        frontEnd = peakTimes[i];
        troughTimes[i] = frontEnd + randMax(numVals-1-frontEnd)/(peaks);
        frontEnd = troughTimes[i];
    }

    int lastVal = 0;
    int temp;
    for (int i = 0; i < peaks-1; i++){ //generate peaks and troughs
        peakVals[i] = lastVal+randMax(MAX_HEIGHT-lastVal);
        while(temp=randMax(variation)>MAX_HEIGHT-peakVals[i]);
        peakVals[i]+=i;
        lastVal = peakVals[i];

        troughVals[i] = randMax(lastVal);
        while(temp=randMax(variation)>lastVal-troughVals[i]);
        troughVals[i]+=temp;
        lastVal = troughVals[i];
    }
    peakVals[peaks-1] = lastVal+randMax(MAX_HEIGHT-lastVal);
    lastVal = peakVals[peaks-1];
    troughVals[peaks-1] = randMax((lastVal<peakVals[0])?lastVal:peakVals[0]);

    /**
    printf("Peaks and troughs: ");
    for (int i = 0; i < peaks; i++){
        printf("%d %d ",peakVals[i],troughVals[i]);
    }
    printf("\n");
    **/
    double amplitude, axis, period;
    for (int i = 0; i < peaks; i++){
        //generate values between peak and next trough, then trough and next peak
        //use cos function to bridge gap between peak and trough

        //peak and next trough
        amplitude = (peakVals[i] - troughVals[i])/2.0;
        axis = (peakVals[i] + troughVals[i])/2.0;
        period = (troughTimes[i] - peakTimes[i])*2.0;
        //printf("Function for peak %d: amp: %.2f, axis: %.2f, period: %.2f\n",i,amplitude,axis,period);

        vals[peakTimes[i]] = peakVals[i];
        //printf("%d ",vals[peakTimes[i]]);
        for (int j = peakTimes[i]+1; j < troughTimes[i]; j++){ //generate from peak to trough
            vals[j] = round(amplitude*cos(2*PI*(j-peakTimes[i])/period) + axis);
            //printf("%d ",vals[j]);
        }
        //printf("\n\n");

        //trough and next peak
        amplitude = (peakVals[(i+1)%peaks]-troughVals[i])/-2.0;
        axis = (peakVals[(i+1)%peaks]+troughVals[i])/2.0;
        period = (i == peaks-1)? (numVals-troughTimes[i]+peakTimes[0])*2.0 : (peakTimes[i+1] - troughTimes[i])*2.0;
        //printf("Function for trough %d: amp: %.2f, axis: %.2f, period: %.2f\n",i,amplitude,axis,period);

        vals[troughTimes[i]] = troughVals[i];
        //printf("%d ",vals[troughTimes[i]]);
        for (int j = troughTimes[i]+1; j < troughTimes[i]+period/2; j++){ //generate from trough to peak
            vals[numVals%j] = round(amplitude*cos(2*PI*(j-troughTimes[i])/period)+axis);
            //printf("%d ",vals[numVals%j]);
        }
        //printf("\n\n");
    }

    //repeat wave
    int i,j;
    for (i = 0, j = numVals; j < length*TICK_SPEED; i++, j++){
        vals[j] = vals[i%numVals];
    }
    return vals;
}
//1000/(n/TICK_SPEED)
// F = 200
int getFrequency (double *wave, int size){
    int n = -1;
    double closest = 100000;
    double total1 = 0;
    double total2 = 0;
    double total3 = 0;
    int minN = 1000*TICK_SPEED/MAX_FREQ;
    double average;
    double tempClosest;

    for (int i = 0; i < minN-1; i++){
        total1 += wave[i];
        total2 += wave[minN-1+i];
        total3 += wave[2*(minN-1)+i];
    }
    for (int i = minN-1; i <= size/3; i++){
        total1 += wave[i];
        total2 -= wave[i];
        total2 += wave[i*2];
        total2 += wave[i*2+1];
        total3 -= wave[i*2];
        total3 -= wave[i*2+1];
        total3 += wave[i*3];
        total3 += wave[i*3+1];
        total3 += wave[i*3+2];
        average = (total1+total2+total3)/3;
        tempClosest = sqrt((pow(total1-average,2)+pow(total2-average,2)+pow(total3-average,2))/3);
        //printf("%f\n",tempClosest);
        if (tempClosest < closest){
            closest = tempClosest;
            n = i;
        }
    }
    return 1000*TICK_SPEED/n;
}

int main(){
    int freq = 400;
    int scanTime = 100; //in milliseconds
    double wave[scanTime*TICK_SPEED];
    int variance = 30;
    generateWave(freq,variance,scanTime, wave);
    //double *wave = malloc(TICK_SPEED*scanTime*sizeof(double));
    int n = 3;
    int frequencies[] = {100,200,300};
    int amplitudes[] = {30,45,73};
    //generateWave2(n,frequencies,amplitudes,scanTime,wave,TICK_SPEED*scanTime);
    printf("%d\n",getFrequency(wave,scanTime*TICK_SPEED));
}








