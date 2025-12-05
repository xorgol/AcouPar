#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <windows.h>
#include "strength.h"
#include "functions.h"

int main(int argc, char** argv)
{
    MLS* lpAmp;
    int* stereoMode;
    char* filename;
    char* txtname;
    // FILE* FP;
    float fAP[NPAR * NBANDS * 2]; // parametri acustici
    char  cAP[NPAR * NBANDS * 2 + NBANDS]; // validit� o meno dei parametri acustici
    char ok = 1;
    
    lpAmp = (MLS*)malloc((size_t)sizeof(MLS)); //Structure containing all data
    filename = (char*)malloc(sizeof(char) * 1024);
    txtname = (char*)malloc(sizeof(char) * 1024);
    
    if (filename == NULL)
    {
        printf("Error in mallocn");
        exit(1);
    }
    
    // Parse input arguments
    if (argc != 3) {
        printf("Two input arguments are expected.\n");
        exit(1);
    }
    
    if (strcmp(argv[1], "--omni") == 0) {
        stereoMode = 0;
        printf("Omni elaboration applied.\n");
    } else if (strcmp(argv[1], "--wy") == 0) {
        stereoMode = 1;
        printf("Soundfield WY elaboration applied.\n");
    } else if (strcmp(argv[1], "--pu") == 0){
        stereoMode = 2;
        printf("Omni/FigureOfEight PU Pressure-Velocity elaboration applied.\n");
    } else if (strcmp(argv[1], "--pp") == 0){
        stereoMode = 3;
        printf("P-P probe elaboration applied.\n");
    } else if (strcmp(argv[1], "--bin") == 0){
        stereoMode = 4;
        printf("Binaural elaboration applied.\n");
    } else {
        printf("Unknown option: %s\n", argv[1]);
        printf("The stereo mode option must be one from "
               "--omni, --wy, --pu, --pp or --bin.\n");
        exit(1);
    }
    //0=2 mono mikes, 1=Soundfield WY, 2=Omni/Eight PU, 3=p-p probe, 4=binaural

    // get file name tobe analyzed from command line
    strcpy(filename, argv[2]);

    printf("\nTry opening %s \n",filename);
        drwav wav;
        if (!drwav_init_file(&wav, filename, NULL)) {
            return -1;
        }
        printf("N. of channels = %i \n", wav.channels);
        int IRL = (int)wav.totalPCMFrameCount;
        float* pSampleData = (float*)malloc((size_t)IRL * wav.channels * sizeof(float));
        float* fIR = (float*)malloc((size_t)IRL * sizeof(float));
        float* fIR2 = (float*)malloc((size_t)IRL * sizeof(float));

        drwav_read_pcm_frames_f32(&wav, wav.totalPCMFrameCount, pSampleData);
        // At this point pSampleData contains every decoded sample as float-32, one channel after the other<s.

        //FP = fopen("out.txt", "w");
        //fprintf(FP, "SAMPLES:\t%i\n",IRL);
        //fprintf(FP, "BITSPERSAMPLE:\t32\n");
        //fprintf(FP, "CHANNELS:\t%i\n",wav.channels);
        //fprintf(FP, "SAMPLERATE:\t%i\n",wav.sampleRate);
        //fprintf(FP, "NORMALIZED:\tTRUE\n");

        if (wav.channels==2)
        {
            for (int i = 0; i < IRL; i++)
            {
                // printf("i= %i, ",i);
            fIR[i] = pSampleData[2*i];
            fIR2[i] = pSampleData[2*i+1];
            //fprintf(FP,"%f\t", fIR[i]);
            //fprintf(FP,"%f\n", fIR2[i]);
            }
        }
        else
        {
            for (int i = 0; i < IRL; i++)
            {
                fIR[i] = pSampleData[i];
                fIR2[i] = 0.0;
             }
        }
        //fclose(FP);

        // the left and right channels are now available in fIR and fIR2 in float 32-bits format
        // I can use them for calculating Acoustical Parameters

        // allocate additional buffers fIRFilt, fIRFilt2, LeftBuf, RightBuf
        float* fIRFilt = (float*)malloc((size_t)IRL * sizeof(float) * 2);
        float* fIRFilt2 = (float*)malloc((size_t)IRL * sizeof(float) * 2);
        int NbytesBuf = (int)(1000.0 * 4.0 * 12.0 * (double)IRL / (double)wav.sampleRate); //4 bytes, 12 bands, 1000 ms
        float* LeftBuf = (float*)malloc((size_t)NbytesBuf);
        float* RightBuf = (float*)malloc((size_t)NbytesBuf);

        // set lpAmp variables
        lpAmp->ref1 = 69.0;     // Liv. riferimento sorgente a 10m - 31.5 Hz
        lpAmp->ref2 = 69.0;
        lpAmp->ref3 = 69.0;
        lpAmp->ref4 = 69.0;
        lpAmp->ref5 = 69.0;
        lpAmp->ref6 = 69.0;
        lpAmp->ref7 = 69.0;
        lpAmp->ref8 = 69.0;
        lpAmp->ref9 = 69.0;
        lpAmp->ref10 = 69.0;   // Liv. riferimento sorgente a 10m - 16 kHz
        lpAmp->ref11 = 77.0;  // Liv. riferimento sorgente a 10m - A
        lpAmp->ref12 = 79.0;  // Liv. riferimento sorgente a 10m - LIN

        lpAmp->Distance     = 10.0;    //m  - reference distance for G
        lpAmp->Space        = 120;       //mm - p-p probe mic spacing
        lpAmp->Speed        = 340.0;      //m/s - speed of sound
        lpAmp->fThreshold   = -20.0; //dB - Threshold for detecting the peak of direct sound
        lpAmp->fRTUdBstart  = -5.0; //dB - start point for RT-user
        lpAmp->fRTUdBend    = -15.0;  //dB - end point for RT-user
        lpAmp->StereoMode   = stereoMode; // 0=2 mono mikes, 1=Soundfield WY, 2=Omni/Eight PU, 3=p-p probe, 4=binaural
        lpAmp->IACCmode     = 0;   // 0=Early, 1=Late, 2=Whole
        lpAmp->StageMode    = 0;  // 0=Normal, 1=Stage (compute ST instead of ts,D)
        lpAmp->AVGmode      = 0;	// 0=Normal, 1=Average Mode , 250 to 2k
        lpAmp->AppendMode   = 1; // 0=not append, 1=append
        lpAmp->White2Pink   = 0; // 0=leave unchanged, 1=correct
        lpAmp->fFS          = 120.0;  // Full Scale in dB-SPL
        lpAmp->lSampRate    = wav.sampleRate;
        lpAmp->lN           = IRL;
        lpAmp->wChannels    = wav.channels;
        strcpy(lpAmp->wTitle, filename);

        // Ricerca valori MAX picco
        float MaxL = 0.0;
        float MaxR = 0.0;
        for (int i = 0; i < IRL; i++)
        {
            if (fabs(fIR[i]) > MaxL) MaxL = (float)fabs(fIR[i]);
            if (fabs(fIR2[i]) > MaxR) MaxR = (float)fabs(fIR2[i]);
        }
        // Max = dB(Max * Max) + (double)lpAmp->fFS;// MAX SPL istantaneo
        lpAmp->MaxL = (float)MaxL; // Massimo canale Left
        lpAmp->MaxR = (float)MaxR; // Massimo canale Right
        printf("MaxL = %f - MaxR = %f \n", MaxL, MaxR);

        // compute acoustical parameters
        ok = CalculateAcoustics(fIR, fIRFilt, fIR2, fIRFilt2, fAP, cAP, lpAmp, LeftBuf, RightBuf, &wav);

        // now I save the results in a TXT file
        strcpy(txtname, "acoupar");
        switch(lpAmp->StereoMode)
        {
            case 0: strcat(txtname, "_omni.txt"); //2 omni
                break;
            case 1: strcat(txtname, "_WY.txt"); //Soundfield, WY
                break;
            case 2: strcat(txtname, "_pu.txt"); //p-u probe, Ambix WY
                break;
            case 3: strcat(txtname, "_pp.txt"); //p-p probe
                break;
            case 4: strcat(txtname, "_BIN.txt"); //binaural
                break;
            default:strcat(txtname, ".txt"); //unknown stereo mode
        }
 
        // save results, averaging L&R if required
        int avgmode = 0;
        if ((lpAmp->StereoMode == 0) | (lpAmp->StereoMode == 4))
        {
            avgmode = 1;
        }
        SaveResults(&wav,  txtname, filename,cAP,fAP,avgmode,lpAmp);
        drwav_uninit(&wav);
        return 0;
}
