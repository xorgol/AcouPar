//#include "dr_wav.h"
//#include <stdio.h>
//#include <stdint.h>
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
//#include <time.h>

double NoiseCorrection(float* fIRFilt, long IRL)
{
	double NC = 0.0;
	double sample = 0.0;
	long lIRL16 = IRL / 16;
	long i;

	// calcolo la media delle energie dell'ultimo sedicesimo di IR
	for (i = (IRL - lIRL16); i < IRL; i++)
	{
		sample = fIRFilt[i];
		sample *= sample;
		NC += sample;
	}
	NC = NC / (double)lIRL16; // è il quadrato del valore RMS del rumore, di fatto è già un valore di Leq

	return NC;
}

void CalculateSchInt(float* fIRFilt, long IRL, long FAT, double NC, double* tP2, drwav* wav)
{
	// calcola l'integrale di Schroeder memorizzandolo nel vettore d'ingresso
	// e, contemporaneamente, calcola anche il termine short(t*p^2) per la determinazione di TS
	float fSamp = (float)wav->sampleRate;
	long i, i2 = IRL - FAT;
	double TP2 = 0.0;
	double Tot = 0.0;
	char NoiseC = (NC != 0.0);
	double sample;
	//	if (NoiseC) MessageBox(NULL,"Noise Correction enabled","Warning",MB_ICONSTOP|MB_OK);
	if (NoiseC)
	{
		for (i = IRL - 1; i >= 0; i--)
		{
			// quadrato dell'ampiezza
			sample = fIRFilt[i] * fIRFilt[i] - NC;
			// se non ho raggiunto il FAT aggiorno TP2 e i2
			if (i2) TP2 += sample * (i2--);
			// aggiorno il totale dei quadrati
			Tot += sample;
			// memorizzo il valore attuale dell'integrale
			fIRFilt[i] = (float)Tot;
		}
	}
	else
	{
		for (i = IRL - 1; i >= 0; i--)
		{
			// quadrato dell'ampiezza
			sample = fIRFilt[i] * fIRFilt[i];
			// se non ho raggiunto il FAT aggiorno TP2 e i2
			if (i2) TP2 += sample * (i2--);
			// aggiorno il totale dei quadrati
			Tot += sample;
			// memorizzo il valore attuale dell'integrale
			fIRFilt[i] = (float)Tot;
		}
	}
	// trasformo TP2 da campioni a secondi
	TP2 /= fSamp;
	// memorizzo TP2 in *tP2
	*tP2 = TP2;
}

void Schroeder(float* fIRFilt, long IRL, double NC)
{
	// calcola l'integrale di Schroeder memorizzandolo nel vettore d'ingresso
	long i;
	double Tot = 0.0;
	char NoiseC = (NC != 0.0);
	double sample;
	//	if (NoiseC) MessageBox(NULL,"Noise Correction enabled","Warning",MB_ICONSTOP|MB_OK);
	if (NoiseC)
	{
		for (i = IRL - 1; i >= 0; i--)
		{
			sample = fIRFilt[i] * fIRFilt[i] - NC;
			// aggiorno il totale dei quadrati
			Tot += sample;
			// memorizzo il valore attuale dell'integrale
			fIRFilt[i] = (float)Tot;
		}
	}
	else
	{
		for (i = IRL - 1; i >= 0; i--)
		{
			sample = fIRFilt[i] * fIRFilt[i];
			// aggiorno il totale dei quadrati
			Tot += sample;
			// memorizzo il valore attuale dell'integrale
			fIRFilt[i] = (float)Tot;
		}
	}
}

double dB(double a)
{
	double db;
	double zero = 0.0;

	if (a > zero) db = 10.0 * log10(a); //10.0
	else db = -120.0; // valore minimo assumibile, in dB, da un float

	return (db);
}

void reverb(float* fIRFilt, long IRL, long FAT, char NoiseC, char cEDT2,
	float dBstart, float dBdelta, double fSamp, float* T)
{
	long i1; // istante subito successivo all'arrivo dell'onda diretta
	long x = 0, n = 0;
	long iValid;
	double dB1, dB0, dBend, y, B = 0.0, r = 0.0, t = 0.0, denom = 0.0, denom1 = 0.0; //a;
	double SumX = 0.0, SumY = 0.0, SumXY = 0.0, SumX2 = 0.0, SumY2 = 0.0;
	long ok = 1;

	// istante d'arrivo dell'onda diretta
	i1 = FAT;
	dB1 = fIRFilt[i1]; // valore in dB del primo campione
	dB0 = dB1 - dBstart;   // valore in dB dell'inizio del tempo di riverbero
	dBend = -dBstart - dBdelta; // valore in dB della fine del tempo di riverbero

	// cerco l'istante di inizio decadimento
	while ((dB0 < fIRFilt[i1]) && (i1 < IRL)) i1++;
	if (i1 != IRL)
	{
		double Y1=0, Y2=0;
		// inizializzo x
		x = 1;
		Y1 = fIRFilt[i1] - dB1;
		// ciclo di ricerca e somma
		do
		{
			y = fIRFilt[i1] - dB1;
			if ((y >= dBend) && (i1 < IRL))
			{
				SumX += (double)x;
				SumY += y;
				SumXY += y * (double)x;
				SumX2 += (double)x * (double)x;
				SumY2 += y * y;
				x++;
				Y2 = y;
			}
			i1++;
		} while ((y >= dBend) && (i1 < IRL));

		// limite di validità dei parametri, minore nel caso
		// ci sia stata Noise Correction
		if (NoiseC) iValid = (IRL - IRL / 16);
		else iValid = (IRL - IRL / 64);

		// calcolo della regressione
		if (i1 <= iValid)
		{
			n = x - 1;
			denom = ((double)n * SumX2 - (SumX * SumX));
			if (denom != 0.0)
				B = ((double)n * SumXY - SumX * SumY) / denom;
			else { ok = ok + 10; t = 0.0; }
			if (cEDT2) B = (Y2 - Y1) / (n - 1);
			if (B != 0.0)
				t = -60.0 / B / fSamp; // tempo di riverbero alla Sabine
			else { ok = ok + 20; t = 0.0; }
			denom1 = ((double)n * SumX2 - SumX * SumX) * ((double)n * SumY2 - SumY * SumY);
			if (denom1 != 0.0)
				r = fabs(((double)n * SumXY - SumX * SumY) * ((double)n * SumXY - SumX * SumY) / denom1);
			else { ok = ok + 40; }
			if ((r < 0.85) && (dBstart > 0.5))
			{
				//coeff.correl. insufficiente, ma non per EDT
				ok = ok + 100;
				t = 0.0; r = 0.0;
			}
		}
		else ok = 0;
	}
	else ok = 0;
	// controllo finale
	if (ok != 1)
	{
		t = t - t;
		r = (double)ok;
	} // segnalo la non validità dei parametri

	// ricopio i valori
	*T = (float)t;
}

void CalculateParameters(float* fIRFilt, float* fIRFilt2, float* fAP, char* cAP, MLS* lpAmp, long FAT, long FAT2, double NC, double NC2, drwav* wav, long Band)
{
	long IRL = (long)wav->totalPCMFrameCount; // lunghezza imp. res.  
	double tP, tP2;  // parametro per il calcolo del tempo baricentrico
	long i;
	long Nchannels = wav->channels;
	float fSamp = (float)wav->sampleRate;
	float FS = (float)lpAmp->fFS;
	// printf("CalcParam FS = %f \n", FS);
	double delay = 0.0;
	// printf("Nchannels = %ld, LpAmp->Stereomode = %ld \n", Nchannels, lpAmp->StereoMode);
	
	// definisco il ritardo di gruppo del filtro d'ottava impiegato
	if (Band < 11)
	{
		delay = 0.05; // 50 ms a 31.5 Hz
		for (i = 1; i < Band; i++) delay /= 2.0; // continua a dimezzare per ciascuna ottava
	}

	// Calcolo Peakiness se il segnale è mono oppure se è stereo ma con 2 micr. indipendenti
	if ((Nchannels == 1) | (lpAmp->StereoMode == 0))
	{
		float xmax, x, sumx;
		xmax = 0.0; sumx = 0.0;
		for (i = 0; i < IRL; i++)
		{
			x = fIRFilt[i] * fIRFilt[i];
			if (x > xmax) xmax = x;
			sumx += x;
		}
		sumx = sumx / IRL; // Valore Medio Quadrato dei campioni
		// 11: Peakiness, 
		cAP[11] = 'y';
		fAP[11] = (float)(10.0 * log10(xmax / sumx));

		if (Nchannels == 2)
		{
			xmax = 0.0; sumx = 0.0;
			for (i = 0; i < IRL; i++)
			{
				x = fIRFilt2[i] * fIRFilt2[i];
				if (x > xmax) xmax = x;
				sumx += x;
			}
			// 11: Peakiness, 
			sumx = sumx / IRL; // Valore Medio Quadrato dei campioni
			cAP[11 + NPAR * NBANDS] = 'y';
			fAP[11 + NPAR * NBANDS] = (float)(10.0 * log10(xmax / sumx));
		}
	}

	// Calcolo LFC=12
	// 12: LFC, occorre che la IR, a partire da FAT, sia lunga almeno 80ms
	if ((((IRL - FAT) / fSamp) > (0.08f + delay)) && (Nchannels == 2) && (lpAmp->StereoMode > 0) && (lpAmp->StereoMode < 4))
	{
		double sumOmni2 = 0.0;
		double sumOmniOtto = 0.0;
		for (i = FAT + (long)((0.005 + delay) * fSamp); i < FAT + (long)((0.080 + delay) * fSamp); i++)
			sumOmniOtto += fabs(fIRFilt2[i] * fIRFilt[i]);
		for (i = FAT; i < FAT + (long)((0.080 + delay) * fSamp); i++)
			sumOmni2 += fIRFilt[i] * fIRFilt[i];
		fAP[12] = (float)(sumOmniOtto / sumOmni2);
		if (lpAmp->StereoMode == 1) fAP[12] /= (float)sqrt(2.0); // E' un Soundfield, debbo raddoppiare l'energia del W
		cAP[12] = 'y';
		// duplico sul Right
		fAP[12 + NPAR * NBANDS] = fAP[12];
		cAP[12 + NPAR * NBANDS] = cAP[12];
		// printf("Band= %ld - Sum8 = %f - SumOmni = %f \n", Band, sumOmniOtto , sumOmni2);
		// printf("Band= %ld - Jlfc = %f \n", Band, fAP[12]);
	}

	// Calcolo IACC=11, Tau IACC=12, w IACC=13
	// 14, etc., occorre che la IR, a partire da FAT, sia lunga almeno 80ms
	if ((((IRL - FAT) / fSamp) > (0.08f + delay)) && (Nchannels == 2) && (lpAmp->StereoMode == 4))
	{
		long i1, i2, j, jmax, j1, j2;
		long L1ms;
		double Den = 0.0;
		double SumA = 0.0;
		double SumB = 0.0;
		double SumAB = 0.0;
		double IACC = 0.0;
		double TauIACC = 0.0;
		double wIACC = 0.0;
		double Limit;
		// char textbuf[64];

		L1ms = (long)(0.001 * fSamp); // Numero di campioni pari ad 1ms

		// Definisco i limiti di integrazione i1 e i2
		switch (lpAmp->IACCmode)
		{
		case 0: // 0-80 ms
			i1 = FAT;
			i2 = FAT + (long)((0.080 + delay) * fSamp);
			break;
		case 1: // 80-inf ms
			i1 = FAT + (long)((0.080 + delay) * fSamp);
			i2 = IRL - L1ms;
			break;
		case 2: // 0-inf ms
			i1 = FAT;
			i2 = IRL - L1ms;
			break;
		default: // errore, forzo 0
			lpAmp->IACCmode = 0;
			i1 = FAT;
			i2 = FAT + (long)((0.080 + delay) * fSamp);
			break;
		}
		if ((i1 - L1ms) < 0) i1 = L1ms; //evita che l'indice diventi negativo

		// Calcolo il denominatore
		for (i = i1; i < i2; i++)
			SumA += fIRFilt[i] * fIRFilt[i];
		for (i = i1; i < i2; i++)
			SumB += fIRFilt2[i] * fIRFilt2[i];
		Den = sqrt(SumA * SumB); // Questo è il denominatore

		// Ora calcolo la CC e cerco il suo massimo
		IACC = 0.0; jmax = 0;
		for (j = 0; j <= L1ms; j++)
		{
			SumAB = 0.0;
			for (i = i1; i < i2; i++)
				SumAB += fIRFilt[i] * fIRFilt2[i - j];
			if (SumAB > IACC)
			{
				IACC = SumAB;
				jmax = -j;     // Posizione del massimo, lato sinistro
			}
			SumAB = 0.0;
			for (i = i1; i < i2; i++)
				SumAB += fIRFilt[i] * fIRFilt2[i + j];
			if (SumAB > IACC)
			{
				IACC = SumAB;
				jmax = j;     // Posizione del massimo, lato destro
			}
		}
	
		// Ora cerco i valori di J in cui CC scende a 0.9*Max
		Limit = 0.9 * IACC;
		j1 = 0; j2 = 0;
		for (j = 0; j <= L1ms; j++)
		{
			SumAB = 0.0;
			for (i = i1; i < i2; i++) SumAB += (double)(fIRFilt[i] * fIRFilt2[i - j + jmax]);
			if (SumAB <= Limit)
			{
				j1 = -j;
				break;     // ho trovato il fianco sinistro, esco
			}
		}
		for (j = 0; j <= L1ms; j++)
		{
			SumAB = 0.0;
			for (i = i1; i < i2; i++) SumAB += (double)(fIRFilt[i] * fIRFilt2[i + j + jmax]);
			if (SumAB <= Limit)
			{
				j2 = j;
				break;     // ho trovato il fianco destro, esco
			}
		}

		IACC /= Den;					// Questo è il valore della IACC
		TauIACC = jmax * 1000.0 / fSamp;  // Questo è il tempo di ritardo del picco del max., in ms
		wIACC = (j2 - j1) * 1000.0 / fSamp; // Questa è la larghezza del picco del max., in ms
		fAP[11] = (float)IACC; cAP[11] = 'y';
		fAP[12] = (float)TauIACC; cAP[12] = 'y';
		fAP[13] = (float)wIACC; cAP[13] = 'y';
		fAP[11 + NPAR * NBANDS] = (float)IACC; cAP[11 + NPAR * NBANDS] = 'y';
		fAP[12 + NPAR * NBANDS] = (float)TauIACC; cAP[12 + NPAR * NBANDS] = 'y';
		fAP[13 + NPAR * NBANDS] = (float)wIACC; cAP[13 + NPAR * NBANDS] = 'y';
	}

	// calcolo l'integrale di Schroeder
	// MessageBox(NULL,"Calculate Schroeder","Warning",MB_ICONSTOP|MB_OK);
	CalculateSchInt(fIRFilt, IRL, FAT, NC, &tP, wav);
	if (Nchannels == 2) CalculateSchInt(fIRFilt2, IRL, FAT2, NC2, &tP2, wav);

	// printf("Band= %ld - Sch-L = %f - Sch-R = %f \n", Band, fIRFilt[0], fIRFilt2[0]);

	// a questo punto
	// procedo al calcolo effettivo di alcuni parametri

	// 3: C50, occorre che la IR, a partire da FAT, sia lunga almeno 50ms
	// 5: D50, occorre che la IR, a partire da FAT, sia lunga almeno 50ms
	if (((IRL - FAT) / fSamp) < (0.05f + delay)) { cAP[3] = 'n'; cAP[5] = 'n'; }
	else
	{   // C50
		fAP[3] = (float)(10.0 * log10((fIRFilt[FAT] - fIRFilt[FAT + (long)((0.05 + delay) * fSamp)]) / fIRFilt[FAT + (long)((0.05 + delay) * fSamp)]));
		cAP[3] = 'y';
		// D50  
		fAP[5] = ((fIRFilt[FAT] - fIRFilt[FAT + (long)((0.05 + delay) * fSamp)]) / fIRFilt[FAT]);
		cAP[5] = 'y';
		// esprimo D50 in percentuale 
		fAP[5] *= 100.0f;
	}
	if (Nchannels == 2)
	{
		if (((IRL - FAT2) / fSamp) < (0.05f + delay)) { cAP[3 + NPAR * NBANDS] = 'n'; cAP[5 + NPAR * NBANDS] = 'n'; }
		else
		{
			fAP[3 + NPAR * NBANDS] = (float)(10.0 * log10((fIRFilt2[FAT2] - fIRFilt2[FAT2 + (long)((0.05 + delay) * fSamp)]) / fIRFilt2[FAT2 + (long)((0.05 + delay) * fSamp)]));
			cAP[3 + NPAR * NBANDS] = 'y';
			// D50  
			fAP[5 + NPAR * NBANDS] = ((fIRFilt2[FAT2] - fIRFilt2[FAT2 + (long)((0.05 + delay) * fSamp)]) / fIRFilt2[FAT2]);
			cAP[5 + NPAR * NBANDS] = 'y';
			// esprimo D50 in percentuale 
			fAP[5 + NPAR * NBANDS] *= 100.0f;
		}
	}

	// 4: C80, occorre che la IR, a partire da FAT, sia lunga almeno 80ms
	if (((IRL - FAT) / fSamp) < (0.08f + delay)) cAP[4] = 'n';
	else
	{
		fAP[4] = (float)(10.0 * log10((fIRFilt[FAT] - fIRFilt[FAT + (long)((0.08 + delay) * fSamp)]) / fIRFilt[FAT + (long)((0.08 + delay) * fSamp)]));
		cAP[4] = 'y';
	}
	if (Nchannels == 2)
	{
		if (((IRL - FAT2) / fSamp) < (0.08f + delay)) cAP[4 + NPAR * NBANDS] = 'n';
		else
		{
			fAP[4 + NPAR * NBANDS] = (float)(10.0 * log10((fIRFilt2[FAT2] - fIRFilt2[FAT2 + (long)((0.08 + delay) * fSamp)]) / fIRFilt2[FAT2 + (long)((0.08 + delay) * fSamp)]));
			cAP[4 + NPAR * NBANDS] = 'y';
		}
	}

	// Se non ho lo StageMode, calcolo ts, D50 e EDT. Altrimenti, li sostitusco con ST1, ST2, STlate
	if (lpAmp->StageMode)
	{
		// 5: St1, occorre che la IR, a partire da FAT, sia lunga almeno 100ms
		if (((IRL - FAT) / fSamp) < (0.100f + delay)) { cAP[5] = 'n'; }
		else
		{	// St1
			fAP[5] = (float)(10.0 * log10((fIRFilt[FAT + (long)((0.020 + delay) * fSamp)] - fIRFilt[FAT + (long)((0.100 + delay) * fSamp)]) / (fIRFilt[FAT] - fIRFilt[FAT + (long)((0.010 + delay) * fSamp)])));
			cAP[5] = 'y';
		}
		if (Nchannels == 2)
		{
			// 5: St1, occorre che la IR, a partire da FAT, sia lunga almeno 100ms
			if (((IRL - FAT2) / fSamp) < (0.100f + delay)) { cAP[5 + NPAR * NBANDS] = 'n'; }
			else
			{	// St1
				fAP[5 + NPAR * NBANDS] = (float)(10.0 * log10((fIRFilt2[FAT2 + (long)((0.020 + delay) * fSamp)] - fIRFilt2[FAT2 + (long)((0.100 + delay) * fSamp)]) / (fIRFilt2[FAT2] - fIRFilt2[FAT2 + (long)((0.010 + delay) * fSamp)])));
				cAP[5 + NPAR * NBANDS] = 'y';
			}
		}

		// 6: St2, occorre che la IR, a partire da FAT, sia lunga almeno 200ms
		if (((IRL - FAT) / fSamp) < (0.200f + delay)) { cAP[6] = 'n'; }
		else
		{	// St2
			fAP[6] = (float)(10.0 * log10((fIRFilt[FAT + (long)((0.020 + delay) * fSamp)] - fIRFilt[FAT + (long)((0.200 + delay) * fSamp)]) / (fIRFilt[FAT] - fIRFilt[FAT + (long)((0.010 + delay) * fSamp)])));
			cAP[6] = 'y';
		}
		if (Nchannels == 2)
		{
			// 6: St2, occorre che la IR, a partire da FAT, sia lunga almeno 200ms
			if (((IRL - FAT2) / fSamp) < (0.200f + delay)) { cAP[6 + NPAR * NBANDS] = 'n'; }
			else
			{	// St2
				fAP[6 + NPAR * NBANDS] = (float)(10.0 * log10((fIRFilt2[FAT2 + (long)((0.020 + delay) * fSamp)] - fIRFilt2[FAT2 + (long)((0.200 + delay) * fSamp)]) / (fIRFilt2[FAT2] - fIRFilt2[FAT2 + (long)((0.010 + delay) * fSamp)])));
				cAP[6 + NPAR * NBANDS] = 'y';
			}
		}

		// 7: Stlate, occorre che la IR, a partire da FAT, sia lunga almeno 1000ms
		if (((IRL - FAT) / fSamp) < (1.000f + delay)) { cAP[7] = 'n'; }
		else
		{	// St2
			fAP[7] = (float)(10.0 * log10((fIRFilt[FAT + (long)((0.100 + delay) * fSamp)] - fIRFilt[FAT + (long)((1.000 + delay) * fSamp)]) / (fIRFilt[FAT] - fIRFilt[FAT + (long)((0.010 + delay) * fSamp)])));
			cAP[7] = 'y';
		}
		if (Nchannels == 2)
		{
			// 7: Stlate, occorre che la IR, a partire da FAT, sia lunga almeno 1000ms
			if (((IRL - FAT2) / fSamp) < (1.000f + delay)) { cAP[7 + NPAR * NBANDS] = 'n'; }
			else
			{	// St2
				fAP[7 + NPAR * NBANDS] = (float)(10.0 * log10((fIRFilt2[FAT2 + (long)((0.100 + delay) * fSamp)] - fIRFilt2[FAT2 + (long)((1.000 + delay) * fSamp)]) / (fIRFilt2[FAT2] - fIRFilt2[FAT2 + (long)((0.010 + delay) * fSamp)])));
				cAP[7 + NPAR * NBANDS] = 'y';
			}
		}
	}
	else
	{
		// 6: TS, si può calcolare sempre (N.B. è in ms e non in s!!!)
		fAP[6] = (float)(tP * 1000.0 / fIRFilt[FAT]); cAP[6] = 'y';
		if (Nchannels == 2) { fAP[6 + NPAR * NBANDS] = (float)(tP * 1000.0 / fIRFilt2[FAT2]); cAP[6 + NPAR * NBANDS] = 'y'; }
	}


	// Calcolo LF=11
	// 11: LF, occorre che la IR, a partire da FAT, sia lunga almeno 80ms
	if ((((IRL - FAT) / fSamp) > (0.08f + delay)) && (Nchannels == 2) && (lpAmp->StereoMode > 0) && (lpAmp->StereoMode < 4))
	{
		double sumOmni2 = 0.0;
		double sumOmniOtto = 0.0;
		fAP[11] = (float)((fIRFilt2[FAT2 + (long)((0.005 + delay) * fSamp)] - fIRFilt2[FAT2 + (long)((0.080 + delay) * fSamp)]) / (fIRFilt[FAT] - fIRFilt[FAT + (long)((0.080 + delay) * fSamp)]));
		if (lpAmp->StereoMode == 1) fAP[11] *= 0.5f; // E' un Soundfield, debbo raddoppiare l'energia del W
		cAP[11] = 'y';
		fAP[11 + NPAR * NBANDS] = fAP[11];
		cAP[11 + NPAR * NBANDS] = cAP[11];
	}


	// trasformo l'integrale di Schroeder in dB 
	// per il calcolo dei tempi di riverbero
	for (i = 0; i < IRL; i++) fIRFilt[i] = (float)dB(fIRFilt[i]);
	if (Nchannels == 2) for (i = 0; i < IRL; i++) fIRFilt2[i] = (float)dB(fIRFilt2[i]);

	// Calcolo Lj=13
	// 13: Lj, occorre che la IR, a partire da FAT, sia lunga almeno 80ms
	if ((((IRL - FAT) / fSamp) > (0.08f + delay)) && (Nchannels == 2) && (lpAmp->StereoMode > 0) && (lpAmp->StereoMode < 4))
	{
		float L80ms = (float)(fIRFilt2[FAT + (long)((0.08 + delay) * fSamp)] + lpAmp->fFS - 10 * log10(fSamp / 100.0));
		// è data dal livello delle Schroeder a 80ms del can 2 (otto), cui va sottratto il valore di riferimento a 10m (che è funz. della frequenza)
		switch (Band)
		{
		case 1:		fAP[13] = L80ms - lpAmp->ref1; break;
		case 2:		fAP[13] = L80ms - lpAmp->ref2; break;
		case 3:		fAP[13] = L80ms - lpAmp->ref3; break;
		case 4:		fAP[13] = L80ms - lpAmp->ref4; break;
		case 5:		fAP[13] = L80ms - lpAmp->ref5; break;
		case 6:		fAP[13] = L80ms - lpAmp->ref6; break;
		case 7:		fAP[13] = L80ms - lpAmp->ref7; break;
		case 8:		fAP[13] = L80ms - lpAmp->ref8; break;
		case 9:		fAP[13] = L80ms - lpAmp->ref9; break;
		case 10:	fAP[13] = L80ms - lpAmp->ref10; break;
		case 11:	fAP[13] = L80ms - lpAmp->ref11; break;
		case 12:	fAP[13] = L80ms - lpAmp->ref12; break;
		}
		cAP[13] = 'y';
		// duplico sul Right
		fAP[13 + NPAR * NBANDS] = fAP[13];
		cAP[13 + NPAR * NBANDS] = cAP[13];
	}

	//MessageBox(NULL,"Compute EDT","Warning",MB_ICONSTOP|MB_OK);

		// 7: EDT-10dB
	if (!(lpAmp->StageMode))
	{
		reverb(fIRFilt, IRL, FAT, (char)(NC != 0.0), lpAmp->cEDT2, 0.1f, 10.0f, fSamp, (fAP + 7));
		if (fAP[7] != 0.0f) cAP[7] = 'y';
		else cAP[7] = 'n';
		if (Nchannels == 2)
		{
			reverb(fIRFilt2, IRL, FAT2, (char)(NC2 != 0.0), lpAmp->cEDT2, 0.1f, 10.0f, fSamp, (fAP + 7 + NPAR * NBANDS));
			if (fAP[7 + NPAR * NBANDS] != 0.0f) cAP[7 + NPAR * NBANDS] = 'y';
			else cAP[7 + NPAR * NBANDS] = 'n';
		}
	}

		//8: RT-USER
	reverb(fIRFilt, IRL, FAT, (char)(NC != 0.0), lpAmp->cEDT2, (-(lpAmp->fRTUdBstart)), (lpAmp->fRTUdBstart - lpAmp->fRTUdBend), fSamp, (fAP + 8));
	if (fAP[8] != 0.0f) cAP[8] = 'y';
	else cAP[8] = 'n';
	if (Nchannels == 2)
	{
		reverb(fIRFilt2, IRL, FAT2, (char)(NC != 0.0), lpAmp->cEDT2, (-(lpAmp->fRTUdBstart)), (lpAmp->fRTUdBstart - lpAmp->fRTUdBend), fSamp, (fAP + 8 + NPAR * NBANDS));
		if (fAP[8 + NPAR * NBANDS] != 0.0f) cAP[8 + NPAR * NBANDS] = 'y';
		else cAP[8 + NPAR * NBANDS] = 'n';
	}

	// 9: RT-20dB
	reverb(fIRFilt, IRL, FAT, (char)(NC != 0.0), lpAmp->cEDT2, 5.0f, 20.0f, fSamp, (fAP + 9));
	if (fAP[9] != 0.0f) cAP[9] = 'y'; else cAP[9] = 'n';
	if (Nchannels == 2)
	{
		reverb(fIRFilt2, IRL, FAT2, (char)(NC2 != 0.0), lpAmp->cEDT2, 5.0f, 20.0f, fSamp, (fAP + 9 + NPAR * NBANDS));
		if (fAP[9 + NPAR * NBANDS] != 0.0f) cAP[9 + NPAR * NBANDS] = 'y';
		else cAP[9 + NPAR * NBANDS] = 'n';
	}

	// 10, RT-30dB
	reverb(fIRFilt, IRL, FAT, (char)(NC != 0.0), lpAmp->cEDT2, 5.0f, 30.0f, fSamp, (fAP + 10));
	if (fAP[10] != 0.0f) cAP[10] = 'y';
	else cAP[10] = 'n';
	if (Nchannels == 2)
	{
		reverb(fIRFilt2, IRL, FAT2, (char)(NC2 != 0.0), lpAmp->cEDT2, 5.0f, 30.0f, fSamp, (fAP + 10 + NPAR * NBANDS));
		if (fAP[10 + NPAR * NBANDS] != 0.0f) cAP[10 + NPAR * NBANDS] = 'y';
		else cAP[10 + NPAR * NBANDS] = 'n';
	}

	// Signal è l'inizio dello Schroeder, corretto per il calcolo di INR sulla base di T20

	// 0 Signal
	fAP[0] = (float)(fIRFilt[0] + FS - 10 * log10(fSamp / 100.0)); //scalato come il grafico
	cAP[0] = 'y';
	if (Nchannels == 2)
	{
		fAP[0 + NPAR * NBANDS] = (float)(fIRFilt2[0] + FS - 10 * log10(fSamp / 100.0)); //scalato come il grafico
		cAP[0 + NPAR * NBANDS] = 'y';
	}
	// printf("Band= %ld - Sch-L = %f - Sch-R = %f - Ref1 = %f \n", Band, fIRFilt[0], fIRFilt2[0],lpAmp->ref1);
	// printf("Band= %ld - Signal-L = %f - Signal-R = %f \n", Band, fAP[0], fAP[0 + NPAR * NBANDS]);

	// 1 Noise
	if (NC > 0.0f)
	{
		// NC è un Leq, per cui va moltiplicato per il n. di campioni e poi diluito su 1s
		// fAP[1]=(float)(dB(NC*IRL)+FS-10*log10(fSamp)); 
		// NC è un Leq, va già bene cosi'
		fAP[1] = (float)(dB(NC) + FS);
		cAP[1] = 'y';
	}
	else
	{
		fAP[1] = 0.0f;
		cAP[1] = 'n';
	}
	if (Nchannels == 2)
	{
		if (NC2 > 0.0f)
		{
			// NC è un Leq, per cui va moltiplicato per il n. di campioni e poi diluito su 1s
			// fAP[1+NPAR*NBANDS]=(float)(dB(NC2*IRL)+FS-10*log10(fSamp)); 
			// NC è un Leq, per cui va già bene così
			fAP[1 + NPAR * NBANDS] = (float)(dB(NC2) + FS);
			cAP[1 + NPAR * NBANDS] = 'y';
		}
		else
		{
			fAP[1 + NPAR * NBANDS] = 0.0f;
			cAP[1 + NPAR * NBANDS] = 'n';
		}
	}
	// 2 strenGht;
	// è data dal Signal, cui va sottratto il valore di riferimento a 10m (che è funz. della frequenza)
	switch (Band)
	{
	case 1:		fAP[2] = fAP[0] - lpAmp->ref1; break;
	case 2:		fAP[2] = fAP[0] - lpAmp->ref2; break;
	case 3:		fAP[2] = fAP[0] - lpAmp->ref3; break;
	case 4:		fAP[2] = fAP[0] - lpAmp->ref4; break;
	case 5:		fAP[2] = fAP[0] - lpAmp->ref5; break;
	case 6:		fAP[2] = fAP[0] - lpAmp->ref6; break;
	case 7:		fAP[2] = fAP[0] - lpAmp->ref7; break;
	case 8:		fAP[2] = fAP[0] - lpAmp->ref8; break;
	case 9:		fAP[2] = fAP[0] - lpAmp->ref9; break;
	case 10:	fAP[2] = fAP[0] - lpAmp->ref10; break;
	case 11:	fAP[2] = fAP[0] - lpAmp->ref11; break;
	case 12:	fAP[2] = fAP[0] - lpAmp->ref12; break;
	}
	cAP[2] = 'y';
	if (Nchannels == 2)
	{
		switch (Band)
		{
		case 1:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref1; break;
		case 2:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref2; break;
		case 3:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref3; break;
		case 4:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref4; break;
		case 5:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref5; break;
		case 6:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref6; break;
		case 7:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref7; break;
		case 8:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref8; break;
		case 9:		fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref9; break;
		case 10:	fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref10; break;
		case 11:	fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref11; break;
		case 12:	fAP[2 + NPAR * NBANDS] = fAP[0 + NPAR * NBANDS] - lpAmp->ref12; break;
		}
		cAP[2 + NPAR * NBANDS] = 'y';
	}
}

void IIR(float* fIRFilt, long IRL, double Alfa, double Beta, double Gamma, double Mu, double Sigma)
{
	long n = 0;
	double Xn;
	double Xnm1 = 0.0;
	double Xnm2 = 0.0;
	double Ynm1 = 0.0;
	double Ynm2 = 0.0;

	Xn = fIRFilt[n];
	// X(n) = 2 * (Alfa * (Xn + Mu * Xnm1 + Sigma * Xnm2) + Gamma * Ynm1 - Beta * Ynm2) - Basic
	fIRFilt[n] = (float)(2.0 * (Alfa * (Xn + Mu * Xnm1 + Sigma * Xnm2) + Gamma * Ynm1 - Beta * Ynm2));
	Xnm2 = Xnm1;
	Xnm1 = Xn;

	n = 1;
	Xn = fIRFilt[n];
	Ynm1 = fIRFilt[n - 1];
	fIRFilt[n] = (float)(2.0 * (Alfa * (Xn + Mu * Xnm1 + Sigma * Xnm2) + Gamma * Ynm1 - Beta * Ynm2));
	Xnm2 = Xnm1;
	Xnm1 = Xn;

	for (n = 2; n < IRL; n++)
	{
		Xn = fIRFilt[n];
		Ynm1 = fIRFilt[n - 1];
		Ynm2 = fIRFilt[n - 2];
		fIRFilt[n] = (float)(2.0 * (Alfa * (Xn + Mu * Xnm1 + Sigma * Xnm2) + Gamma * Ynm1 - Beta * Ynm2));
		Xnm2 = Xnm1;
		Xnm1 = Xn;
	}
}

void AFilter(float* X, float fSamp, long IRL)
{
	// filtro IIR a 2 poli
	double PI = 3.141592653589793;
	//  filtro 1 - passa alto 20.6 Hz
	double coef1 = 1.0 - exp(-2.0 * PI * (20.6) / fSamp);
	double coef1m = 1.0 - coef1;
	//  filtro 2 - passa alto 20.6 Hz
	double coef2 = 1.0 - exp(-2.0 * PI * (20.6) / fSamp);
	double coef2m = 1.0 - coef2;
	//  filtro 3 - passa alto 107.7 Hz
	double coef3 = 1.0 - exp(-2.0 * PI * (107.7) / fSamp);
	double coef3m = 1.0 - coef3;
	//  filtro 4 - passa.0 alto 737.9 Hz
	double coef4 = 1.0 - exp(-2.0 * PI * (737.9) / fSamp);
	double coef4m = 1.0 - coef4;
	//  filtro 5 - passa basso 12200 Hz
	double coef5 = 1.0 - exp(-2.0 * PI * (12200.0) / fSamp);
	double coef5m = 1.0 - coef5;
	//  filtro 6 - passa basso 12200 Hz
	double coef6 = 1.0 - exp(-2.0 * PI * (12200.0) / fSamp);
	double coef6m = 1.0 - coef6;
	long k;
	double ritardo1 = 0.0, ritardo2 = 0.0, ritardo3 = 0.0, ritardo4 = 0.0;
	double A, B, C, D, E, F;

	k = 0;

	// elaboro separatamente il primo campione
//	primi 4 filtri passa alto in cascata
	A = X[k] - ritardo1;
	ritardo1 = X[k] * coef1 + ritardo1 * coef1m;
	B = A - ritardo2;
	ritardo2 = A * coef2 + ritardo2 * coef2m;
	C = B - ritardo3;
	ritardo3 = B * coef3 + ritardo3 * coef3m;
	D = C - ritardo4;
	ritardo4 = C * coef4 + ritardo4 * coef4m;

	E = D;
	F = D;
	//	2 filtri passa basso conclusivi
	E = D * coef5 + E * coef5m;
	F = E * coef6 + F * coef6m;
	X[k] = (float)(F / 0.858);

	for (k = 1; k < IRL; k++)
	{
		//		primi 4 filtri passa alto in cascata
		A = X[k] - ritardo1;
		ritardo1 = X[k] * coef1 + ritardo1 * coef1m;
		B = A - ritardo2;
		ritardo2 = A * coef2 + ritardo2 * coef2m;
		C = B - ritardo3;
		ritardo3 = B * coef3 + ritardo3 * coef3m;
		D = C - ritardo4;
		ritardo4 = C * coef4 + ritardo4 * coef4m;

		//		2 filtri passa basso conclusivi
		E = D * coef5 + E * coef5m;
		F = E * coef6 + F * coef6m;
		X[k] = (float)(F / 0.858); // correzione per rendere il guadagno ad 1kHz = 0 dB
	}

	/*
	for (k=1; k<IRL;k++)
	{
//		faccio il solo HP a 737.9 Hz
		ritardo4 = X[k] * coef4 + ritardo4 * coef4m;
		D = X[k] - ritardo4;
		X[k]=(float)(D);
	}
	*/

}

void BandFilter(float f0, float fQ, float* fIRFilt, float fSamp, long IRL)
{
	// filtro IIR a 2 poli
	double PI = 3.141592653589793;
	double Teta0 = 2.0 * PI * f0 / fSamp;
	double d = 2.0 * tan(Teta0 / 2.0 / fQ) / sin(Teta0);
	double Beta = 1.0 / 2.0 * (1.0 - 0.5 * d * sin(Teta0)) / (1.0 + 0.5 * d * sin(Teta0));
	double Gamma = (1.0 / 2.0 + Beta) * cos(Teta0);
	//  BP
	double Mu = 0.0;
	double Alfa = (1.0 / 2.0 - Beta) / 2.0;
	double Sigma = -1.0;

	IIR(fIRFilt, IRL, Alfa, Beta, Gamma, Mu, Sigma);
}
void LPFilter(float f0, float fQ, float* fIRFilt, float fSamp, long IRL)
{
	// filtro IIR a 2 poli
	double PI = 3.141592653589793;
	double Teta0 = 2.0 * PI * f0 / fSamp;
	double d = 2.0 * tan(Teta0 / 2.0 / fQ) / sin(Teta0);
	double Beta = 1.0 / 2.0 * (1.0 - 0.5 * d * sin(Teta0)) / (1.0 + 0.5 * d * sin(Teta0));
	double Gamma = (1.0 / 2.0 + Beta) * cos(Teta0);
	//  LP
	double Mu = 2.0;
	double Alfa = (0.5 + Beta - Gamma) / 4.0;
	double Sigma = 1.0;

	IIR(fIRFilt, IRL, Alfa, Beta, Gamma, Mu, Sigma);
}
void HPFilter(float f0, float fQ, float* fIRFilt, float fSamp, long IRL)
{
	// filtro IIR a 2 poli
	double PI = 3.141592653589793;
	double Teta0 = 2.0 * PI * f0 / fSamp;
	double d = 2.0 * tan(Teta0 / 2.0 / fQ) / sin(Teta0);
	double Beta = 1.0 / 2.0 * (1.0 - 0.5 * d * sin(Teta0)) / (1.0 + 0.5 * d * sin(Teta0));
	double Gamma = (1.0 / 2.0 + Beta) * cos(Teta0);
	//  HP
	double Mu = -2.0;
	double Alfa = (0.5 + Beta + Gamma) / 4.0;
	double Sigma = 1.0;

	IIR(fIRFilt, IRL, Alfa, Beta, Gamma, Mu, Sigma);
}
void NotchFilter(float f0, float fQ, float* fIRFilt, float fSamp, long IRL)
{
	// filtro IIR a 2 poli
	double PI = 3.141592653589793;
	double Teta0 = 2.0 * PI * f0 / fSamp;
	double d = 2.0 * tan(Teta0 / 2.0 / fQ) / sin(Teta0);
	double Beta = 1.0 / 2.0 * (1.0 - 0.5 * d * sin(Teta0)) / (1.0 + 0.5 * d * sin(Teta0));
	double Gamma = (1.0 / 2.0 + Beta) * cos(Teta0);
	//  BR
	double Mu = -2 * cos(Teta0);
	double Alfa = (0.5 + Beta) / 2.0;
	double Sigma = 1.0;

	IIR(fIRFilt, IRL, Alfa, Beta, Gamma, Mu, Sigma);
}

void LFilter(float* X, float fSamp, long IRL)
{
	float f0, fQ;
	//  filtro 1 - passa alto 10 Hz
	f0 = 10.0f; // Hz
	fQ = 0.707f;
	HPFilter(f0, fQ, X, fSamp, IRL);
	//  filtro 2 - passa alto 14 Hz
	f0 = 14.0f; // Hz
	fQ = 0.707f;
	HPFilter(f0, fQ, X, fSamp, IRL);
	//  filtro 3 - passa alto 20 Hz
	f0 = 20.0f; // Hz
	fQ = 0.707f;
	HPFilter(f0, fQ, X, fSamp, IRL);
	/*
//  filtro 4 - passa basso 20 kHz
		f0=20000.0f; // Hz
		fQ=0.707f;
		LPFilter(f0,fQ,X,fSamp,IRL);
//  filtro 5 - passa basso 20.4 kHz
		f0=20400.0f; // Hz
		fQ=0.707f;
		LPFilter(f0,fQ,X,fSamp,IRL);
//  filtro 6 - passa basso 20.9 kHz
		f0=20900.0f; // Hz
		fQ=0.707f;
		LPFilter(f0,fQ,X,fSamp,IRL);
	*/
}

void OctaveFilter(float* fIRFilt, float fFc, float fSamp, long IRL)
{
	float f0, fQ;
	// faccio il filtraggio in banda con un IIR a 6 poli
	// primo filtraggio (2poli)
	f0 = 0.755f * fFc;
	fQ = 6.0f;
	BandFilter(f0, fQ, fIRFilt, fSamp, IRL);
	// secondo filtraggio (2poli)
	f0 = fFc;
	fQ = 3.0f;
	BandFilter(f0, fQ, fIRFilt, fSamp, IRL);
	// terzo filtraggio (2poli)
	f0 = 1.33f * fFc;
	fQ = 6.0f;
	BandFilter(f0, fQ, fIRFilt, fSamp, IRL);
	return;
}

void Demean(float* fIRFilt, long IRL)
{
	double mean = 0.0;
	long i;
	for (i = 0; i < IRL; i++) 	mean += (double)fIRFilt[i];
	mean = mean / (double)IRL;
	for (i = 0; i < IRL; i++) 	fIRFilt[i] -= (float)mean;
	return;
}

char CalculateAcoustics(float* fIR, float* fIRFilt, float* fIR2, float* fIRFilt2, float* fAP, char* cAP, MLS* lpAmp, float* LBuf, float* RBuf, drwav* wav)
{
	float fFrq[] = { 31.5f,63.0f,125.0f,250.0f,500.0f,1000.0f,2000.0f,4000.0f,8000.0f,16000.0f };
	long IRL = (int)wav->totalPCMFrameCount; // lunghezza imp. res.
	float FS = lpAmp->fFS;
	char cNC =1;                    // flag per la noise correction
	double NC, NC2;                 // noise correction
	long i;
	long FAT, FAT2, Nbuf, ibuf;
	long lIRL16 = IRL / 16;
	float fIRL16Samp = ((float)lIRL16) / wav->sampleRate;
	float fSamp = (float)wav->sampleRate;   // parametri utili al filtraggio in bande
	float fSamp2 = ((float)wav->sampleRate) / 2.1f; // un po' meno della f. di Nyquist
	double fTh = (double)(lpAmp->fThreshold);
	long Band;   // indice della banda 
	float fFc;  // frequenza di centro banda
	char ok = TRUE;
	float gain;
	long Nchannels = wav->channels;
	float ScaleFactor, X;
	//printf("CalcAcous FS = %f - dBstart = %f - dBend = %f \n", FS,lpAmp->fRTUdBstart, lpAmp->fRTUdBend);

	gain = 12.8825f; // valore da verificare..... = 22.2 dB ???

	// inizializzo tutti i cAP a 'n'
	for (i = 0; i < (NPAR * NBANDS * Nchannels); i++) cAP[i] = 'n';

	// Elaborazione sonda intensimetrica p-p
	if ((Nchannels == 2) && (lpAmp->StereoMode == 3)) // PP probe
	{
		double p, v = 0.0;
		double d = lpAmp->Space / 1000.0;  // era in mm, cosi' ora è in m
		double c = lpAmp->Speed;         // è già in m/s
		double dt = 1.0 / fSamp;		   // intervallo di campionamento in s
		double Factor = c * dt / d;		   // adimensionale

		// Ora applico il filtraggio passa-banda 20-20k Hz
		// LFilter(fIR,fSamp,IRL);LFilter(fIR2,fSamp,IRL);

		// Calcolo p e v
		for (i = 0; i < IRL; i++)
		{
			p = (fIR[i] + fIR2[i]) / 2.0;		// p=(p1+p2)/2
			v += (fIR2[i] - fIR[i]) * Factor;	// v=Integral((p2-p1)/d*c*dt) - Legge di Eulero
			fIR[i] = (float)p;			// left =p
			fIR2[i] = (float)v;			// right=v
		}
	}

	// cerco il first arrival time, lo cerco anzitutto sul canale Left

	// cerco anzitutto il valore max della risposta all'impulso
	X = lpAmp->MaxL;
	// Scendo al valore di threshold sotto il max
	X *= (float)pow(10, (fTh / 20));
	FAT = 0;
	while ((FAT < IRL) && (X > fabs(fIR[FAT]))) FAT++;
	FAT2 = FAT;

	if (Nchannels == 2)
	{	// guardo anche il Right
		// cerco anzitutto il valore max della risposta all'impulso
		X = lpAmp->MaxR;
		// Scendo al valore di threshold sotto il max
		X *= (float)pow(10, (fTh / 20));
		FAT2 = 0;
		while ((FAT2 < IRL) && (X > fabs(fIR2[FAT2]))) FAT2++;
		if (lpAmp->StereoMode > 0)
		{
			// PP probe o Soundfield o Binaural, quindi assumo un unico FAT, pari al minore dei 2
			if (FAT2 < FAT) FAT = FAT2;
			else FAT2 = FAT;
		}
	}

	// se arrivo in fondo al vettore -> non ho trovato FAT -> la soglia è troppo alta -> esco
	if ((FAT == IRL) || (FAT2 == IRL))
	{
		printf( "Can not locate First Arrival Time.\nThreshold is too high.\nTry again with a lower Threshold.\n");
		return (FALSE);
	}

	// Parto 4 campioni prima del trigger point vero e proprio
	if (FAT > 4) FAT -= 4; else FAT = 0;
	if (FAT2 > 4) FAT2 -= 4; else FAT2 = 0;

	// creo il progress meter
	// ProgressCreate(ci, "Calculating Acoustical Parameters...", NULL);

	// Qui inizia la parte modificata nel 2006 (v. 4.2)

	// calcolo noise correction, parametri acustici per le bande 1-10 e riempio Lbuf ed Rbuf
	Band = 0;
	while ((ok) && (Band < 10) && ((fFc = fFrq[Band]) <= fSamp2))
	{
		float fgain = gain;
		if (lpAmp->White2Pink)
			fgain *= (float)(pow(10.0, ((5.0 - (double)Band) * 3.0102999566398119521373889472449) / 20));
		// copio i valori da fIR a fIRFilt 
		for (i = 0; i < IRL; i++) fIRFilt[i] = fIR[i] * fgain;
		if (Nchannels == 2) for (i = 0; i < IRL; i++) fIRFilt2[i] = fIR2[i] * fgain;

		// Applico il filtro d'ottava
		//MessageBox(NULL,"Octave filter","Warning",MB_ICONSTOP|MB_OK);
		OctaveFilter(fIRFilt, fFc, fSamp, IRL);
		if (Nchannels == 2) OctaveFilter(fIRFilt2, fFc, fSamp, IRL);

		// sottraggo il valore medio (demean)
		//MessageBox(NULL,"Demean","Warning",MB_ICONSTOP|MB_OK);
		Demean(fIRFilt, IRL);
		if (Nchannels == 2) Demean(fIRFilt2, IRL);

		// se è il caso calcolo la Noise Correction per la banda corrente
		// altrimenti la pongo = 0.0 per segnalare che non è stata applicata
		NC = 0.0;
		NC2 = 0.0;
		if (cNC)
		{
			NC = NoiseCorrection(fIRFilt, IRL);
			if (Nchannels == 2) NC2 = NoiseCorrection(fIRFilt2, IRL);
		}

		// Svuoto il contenuto dei Buffers
		Nbuf = (long)(1000.0 * (double)IRL / (double)wav->sampleRate);
		for (ibuf = 0; ibuf < Nbuf; ibuf++) LBuf[ibuf + Band * Nbuf] = 0.0f;
		for (ibuf = 0; ibuf < Nbuf; ibuf++) RBuf[ibuf + Band * Nbuf] = 0.0f;

		// Memorizzo i dati filtrati in Buf, con risoluzione di 1ms
		for (i = FAT; i < IRL; i++)
		{
			ibuf = (long)(1000.0 * i / wav->sampleRate);
			LBuf[ibuf + Band * Nbuf] += fIRFilt[i] * fIRFilt[i];
		}
		if (Nchannels == 2)
			for (i = FAT2; i < IRL; i++)
			{
				ibuf = (long)(1000.0 * i / wav->sampleRate);
				RBuf[ibuf + Band * Nbuf] += fIRFilt2[i] * fIRFilt2[i];
			}

		// Ora riscalo i valori diLBuf in modo che i valori in dB siano poi corretti
		//MessageBox(NULL,"Rescale Buf","Warning",MB_ICONSTOP|MB_OK);
		ScaleFactor = (float)(pow(10.0, FS / 10) / (double)wav->sampleRate);
		for (ibuf = 0; ibuf < Nbuf; ibuf++) LBuf[ibuf + Band * Nbuf] *= ScaleFactor;
		if (Nchannels == 2) for (ibuf = 0; ibuf < Nbuf; ibuf++) RBuf[ibuf + Band * Nbuf] *= ScaleFactor;

		// calcolo i parametri acustici per la banda corrente (ricorda che la banda 11 è la A e la 12 è la LIN)
		//MessageBox(NULL,"Acoustical Params","Warning",MB_ICONSTOP|MB_OK);
		CalculateParameters(fIRFilt, fIRFilt2, (fAP + (NPAR * (Band))), (cAP + (NPAR * (Band))), lpAmp, FAT, FAT2, NC, NC2, wav, (short)(Band + 1));

		// Se il segnale è mono, o stereo ma proveniente da due microfoni indipendenti, calcolo la Millisecondness (par=12)
		// e la Impulsiveness (par=13)
		if ((Nchannels == 1) | (lpAmp->StereoMode == 0))
		{
			float xmax, x, sumx, impx, impxmax;
			long i2;
			xmax = 0.0; sumx = 0.0;  impxmax = 0.0;
			for (i = 0; i < (Nbuf - 35); i++)
			{
				x = LBuf[i + Band * Nbuf];
				if (x > xmax) xmax = x;
				sumx += x;
				impx = 0.0;
				for (i2 = 0; i2 < 35; i2++) impx += LBuf[i + Band * Nbuf + i2]; //Impx è il valore medio su 35 ms, quindi è il valore Impulse
				if (impx > impxmax) impxmax = impx;
			}
			sumx = sumx / (Nbuf - 35); // Valore Medio Quadrato dei buffer millisecondici
			// 12: Millisecondness 
			cAP[12 + Band * NPAR] = 'y';
			fAP[12 + Band * NPAR] = (float)(10.0 * log10(xmax / sumx));
			// 13: Impulsiveness 
			impxmax = impxmax / 35; // Valore max della media Impulse su 35 ms;
			cAP[13 + Band * NPAR] = 'y';
			fAP[13 + Band * NPAR] = (float)(10.0 * log10(impxmax / sumx));

			if (Nchannels == 2)
			{
				xmax = 0.0; sumx = 0.0; impxmax = 0.0;
				for (i = 0; i < (Nbuf - 35); i++)
				{
					x = RBuf[i + Band * Nbuf];
					if (x > xmax) xmax = x;
					sumx += x;
					impx = 0.0;
					for (i2 = 0; i2 < 35; i2++) impx += RBuf[i + Band * Nbuf + i2]; //Impx è il valore medio su 35 ms, quindi è il valore Impulse
					if (impx > impxmax) impxmax = impx;
				}
				// 12: Millisecondness 
				sumx = sumx / (Nbuf - 35); // Valore Medio Quadrato dei buffer millisecondici
				cAP[12 + Band * NPAR + NPAR * NBANDS] = 'y';
				fAP[12 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(xmax / sumx));
				// 13: Impulsiveness
				impxmax = impxmax / 35; // Valore max della media Impulse su 35 ms;
				cAP[13 + Band * NPAR + NPAR * NBANDS] = 'y';
				fAP[13 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(impxmax / sumx));
			}
		}

		// aggiorno il progress meter
		// ProgressMeter(ci, (Band + 1), 12);
		// if (*ci->lpProgressCanceled) ok = FALSE;

		// aggiorno la banda corrente
		Band++;
	}

	Band = 10;
	// Calcolo banda A, è la undicesima - se AVGmode è selezionato, anzichè la Banda A si calcola la media da 500 a 1000 di ciascun parametro
	if (lpAmp->AVGmode)
	{
		// Media parametri fra 500 e 1000 Hz
		for (i = 0; i < NPAR; i++)
			if ((cAP[i + NPAR * 4] == 'y') && (cAP[i + NPAR * 5] == 'y'))
			{
				if (i > 4) // media aritmetica, sono percentuali
					fAP[i + NPAR * Band] = (float)((fAP[i + NPAR * 4] + fAP[i + NPAR * 5]) / 2.0);
				else      // media in dB, sono livelli
					fAP[i + NPAR * Band] = (float)dB((pow(10.0, fAP[i + NPAR * 4] / 10.0) + pow(10.0, fAP[i + NPAR * 5] / 10.0)) / 2.0);
				cAP[i + NPAR * Band] = 'y';
			}
			else
				cAP[i + NPAR * Band] = 'n';
		if (Nchannels == 2)
			for (i = 0; i < NPAR; i++)
				if ((cAP[i + NPAR * 4 + NPAR * NBANDS] == 'y') && (cAP[i + NPAR * 5 + NPAR * NBANDS] == 'y'))
				{
					if (i > 4) // media aritmetica, sono percentuali
						fAP[i + NPAR * Band + NPAR * NBANDS] = (float)((fAP[i + NPAR * 4 + NPAR * NBANDS] + fAP[i + NPAR * 5 + NPAR * NBANDS]) / 2.0);
					else      // media in dB, sono livelli
						fAP[i + NPAR * Band + NPAR * NBANDS] = (float)dB((pow(10.0, fAP[i + NPAR * 4 + NPAR * NBANDS] / 10.0) + pow(10.0, fAP[i + NPAR * 5 + NPAR * NBANDS] / 10.0)) / 2.0);
					cAP[i + NPAR * Band + NPAR * NBANDS] = 'y';
				}
				else
					cAP[i + NPAR * Band + NPAR * NBANDS] = 'n';

		// se sto calcolando LF, etc., allora la media va da 125 a 1000
		if ((Nchannels == 2) && (lpAmp->StereoMode > 0) && (lpAmp->StereoMode < 4))
			// Media parametri fra 125 e 1000 Hz
			for (i = 11; i < NPAR; i++)
				if ((cAP[i + NPAR * 2] == 'y') && (cAP[i + NPAR * 3] == 'y') && (cAP[i + NPAR * 4] == 'y') && (cAP[i + NPAR * 5] == 'y'))
				{
					if (i < 13) // media aritmetica sono percentuali
						fAP[i + NPAR * Band] = (float)((fAP[i + NPAR * 2] + fAP[i + NPAR * 3] + fAP[i + NPAR * 4] + fAP[i + NPAR * 5]) / 4.0);
					else      // media in dB per Lj
						fAP[i + NPAR * Band] = (float)dB((pow(10.0, fAP[i + NPAR * 2] / 10.0) + pow(10.0, fAP[i + NPAR * 3] / 10.0) + pow(10.0, fAP[i + NPAR * 4] / 10.0) + pow(10.0, fAP[i + NPAR * 5] / 10.0)) / 4.0);
					cAP[i + NPAR * Band] = 'y';

					// duplico sul Right
					fAP[i + NPAR * Band + NPAR * NBANDS] = fAP[i + NPAR * Band];
					cAP[i + NPAR * Band + NPAR * NBANDS] = cAP[i + NPAR * Band];

				}
				else
				{
					cAP[i + NPAR * Band] = 'n';
					cAP[i + NPAR * Band + NPAR * NBANDS] = cAP[i + NPAR * Band];
				}
	}
	else
	{
		// Faccio il vero filtraggio "A" a norma IEC
		// copio i valori da fIR a fIRFilt 
		for (i = 0; i < IRL; i++) fIRFilt[i] = fIR[i];
		if (Nchannels == 2) for (i = 0; i < IRL; i++) fIRFilt2[i] = fIR2[i];
		// Applico la ponderazione "A"
		AFilter(fIRFilt, fSamp, IRL);
		if (Nchannels == 2) AFilter(fIRFilt2, fSamp, IRL);
		// sottraggo il valore medio (demean)
		Demean(fIRFilt, IRL);
		if (Nchannels == 2) Demean(fIRFilt2, IRL);
		// se è il caso calcolo la Noise Correction per la banda A
		// altrimenti la pongo = 0.0 per segnalare che non è stata applicata
		NC = 0.0;
		NC2 = 0.0;
		if (cNC)
		{
			NC = NoiseCorrection(fIRFilt, IRL);
			if (Nchannels == 2) NC2 = NoiseCorrection(fIRFilt2, IRL);
		}
		// Svuoto il contenuto di Buf
		// MessageBox(NULL,"Clear Buf","Warning",MB_ICONSTOP|MB_OK);
		Nbuf = (long)(1000.0 * (double)IRL / (double)wav->sampleRate);
		for (ibuf = 0; ibuf < Nbuf; ibuf++) LBuf[ibuf + Band * Nbuf] = 0.0f;
		if (Nchannels == 2) for (ibuf = 0; ibuf < Nbuf; ibuf++) RBuf[ibuf + Band * Nbuf] = 0.0f;
		// Memorizzo i dati filtrati in Buf
		for (i = FAT; i < IRL; i++)
		{
			ibuf = (long)(1000.0 * i / wav->sampleRate);
			LBuf[ibuf + Band * Nbuf] += fIRFilt[i] * fIRFilt[i];
		}
		if (Nchannels == 2)
			for (i = FAT2; i < IRL; i++)
			{
				ibuf = (long)(1000.0 * i / wav->sampleRate);
				RBuf[ibuf + Band * Nbuf] += fIRFilt2[i] * fIRFilt2[i];
			}
		// Ora riscalo i valori diLBuf in modo che i valori in dB siano poi corretti
		//MessageBox(NULL,"Rescale Buf","Warning",MB_ICONSTOP|MB_OK);
		ScaleFactor = (float)(pow(10.0, FS / 10) / (double)wav->sampleRate);
		for (ibuf = 0; ibuf < Nbuf; ibuf++) LBuf[ibuf + Band * Nbuf] *= ScaleFactor;
		if (Nchannels == 2) for (ibuf = 0; ibuf < Nbuf; ibuf++) RBuf[ibuf + Band * Nbuf] *= ScaleFactor;

		// calcolo i parametri acustici per la A band
		CalculateParameters(fIRFilt, fIRFilt2, (fAP + (NPAR * (Band))), (cAP + (NPAR * (Band))), lpAmp, FAT, FAT2, NC, NC2, wav, (short)(Band + 1));

		// Se il segnale è mono, o stereo ma proveniente da due microfoni indipendenti, calcolo la Millisecondness (par=12)
		// e la Impulsiveness (par=13)
		if ((Nchannels == 1) | (lpAmp->StereoMode == 0))
		{
			float xmax, x, sumx, impx, impxmax;
			long i2;
			xmax = 0.0; sumx = 0.0;  impxmax = 0.0;
			for (i = 0; i < (Nbuf - 35); i++)
			{
				x = LBuf[i + Band * Nbuf];
				if (x > xmax) xmax = x;
				sumx += x;
				impx = 0.0;
				for (i2 = 0; i2 < 35; i2++) impx += LBuf[i + Band * Nbuf + i2]; //Impx è il valore medio su 35 ms, quindi è il valore Impulse
				if (impx > impxmax) impxmax = impx;
			}
			sumx = sumx / (Nbuf - 35); // Valore Medio Quadrato dei buffer millisecondici
			// 12: Millisecondness 
			cAP[12 + Band * NPAR] = 'y';
			fAP[12 + Band * NPAR] = (float)(10.0 * log10(xmax / sumx));
			// 13: Impulsiveness 
			impxmax = impxmax / 35; // Valore max della media Impulse su 35 ms;
			cAP[13 + Band * NPAR] = 'y';
			fAP[13 + Band * NPAR] = (float)(10.0 * log10(impxmax / sumx));

			if (Nchannels == 2)
			{
				xmax = 0.0; sumx = 0.0; impxmax = 0.0;
				for (i = 0; i < (Nbuf - 35); i++)
				{
					x = RBuf[i + Band * Nbuf];
					if (x > xmax) xmax = x;
					sumx += x;
					impx = 0.0;
					for (i2 = 0; i2 < 35; i2++) impx += RBuf[i + Band * Nbuf + i2]; //Impx è il valore medio su 35 ms, quindi è il valore Impulse
					if (impx > impxmax) impxmax = impx;
				}
				// 12: Millisecondness 
				sumx = sumx / (Nbuf - 35); // Valore Medio Quadrato dei buffer millisecondici
				cAP[12 + Band * NPAR + NPAR * NBANDS] = 'y';
				fAP[12 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(xmax / sumx));
				// 13: Impulsiveness
				impxmax = impxmax / 35; // Valore max della media Impulse su 35 ms;
				cAP[13 + Band * NPAR + NPAR * NBANDS] = 'y';
				fAP[13 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(impxmax / sumx));
			}
		}

	}	// Fine dell'IF sul calcolo della banda A

	// aggiorno il progress meter (è 11 su 12 perchè ho 10 bande + la Wide + A = 12)
	// ProgressMeter(ci, 11, 12);
	// if (*ci->lpProgressCanceled) ok = FALSE;

	Band = 11;
	// Calcolo banda LIN, è la n. 12 di 12
	// copio i valori da fIR a fIRFilt 
	for (i = 0; i < IRL; i++) fIRFilt[i] = fIR[i];
	if (Nchannels == 2) for (i = 0; i < IRL; i++) fIRFilt2[i] = fIR2[i];
	// Applico la ponderazione "Lin"
	LFilter(fIRFilt, fSamp, IRL);
	if (Nchannels == 2) LFilter(fIRFilt2, fSamp, IRL);
	// sottraggo il valore medio (demean)
	// sottraggo il valore medio (demean)
	Demean(fIRFilt, IRL);
	if (Nchannels == 2) Demean(fIRFilt2, IRL);
	// se è il caso calcolo la Noise Correction per la banda LIN
	// altrimenti la pongo = 0.0 per segnalare che non è stata applicata
	if (cNC)
	{
		NC = NoiseCorrection(fIRFilt, IRL);
		if (Nchannels == 2) NC2 = NoiseCorrection(fIRFilt2, IRL);
	}
	else
	{
		NC = 0.0;
		NC2 = 0.0;
	}
	// Svuoto il contenuto di Buf
	//MessageBox(NULL,"Clear Buf","Warning",MB_ICONSTOP|MB_OK);
	Nbuf = (long)(1000.0 * (double)IRL / (double)wav->sampleRate);
	for (ibuf = 0; ibuf < Nbuf; ibuf++) LBuf[ibuf + Band * Nbuf] = 0.0f;
	if (Nchannels == 2) for (ibuf = 0; ibuf < Nbuf; ibuf++) RBuf[ibuf + Band * Nbuf] = 0.0f;
	// Memorizzo i dati filtrati in Buf
	for (i = FAT; i < IRL; i++)
	{
		ibuf = (long)(1000.0 * i / wav->sampleRate);
		LBuf[ibuf + Band * Nbuf] += fIRFilt[i] * fIRFilt[i];
	}
	if (Nchannels == 2)
		for (i = FAT2; i < IRL; i++)
		{
			ibuf = (long)(1000.0 * i / wav->sampleRate);
			RBuf[ibuf + Band * Nbuf] += fIRFilt2[i] * fIRFilt2[i];
		}
	// Ora riscalo i valori diLBuf in modo che i valori in dB siano poi corretti
	//MessageBox(NULL,"Rescale Buf","Warning",MB_ICONSTOP|MB_OK);
	ScaleFactor = (float)(pow(10.0, FS / 10) / (double)wav->sampleRate);
	for (ibuf = 0; ibuf < Nbuf; ibuf++) LBuf[ibuf + Band * Nbuf] *= ScaleFactor;
	if (Nchannels == 2) for (ibuf = 0; ibuf < Nbuf; ibuf++) RBuf[ibuf + Band * Nbuf] *= ScaleFactor;
	// calcolo i parametri acustici per la banda LIN
	CalculateParameters(fIRFilt, fIRFilt2, (fAP + (NPAR * (Band))), (cAP + (NPAR * (Band))), lpAmp, FAT, FAT2, NC, NC2, wav, (short)(Band + 1));

	// Se il segnale è mono, o stereo ma proveniente da due microfoni indipendenti, calcolo la Millisecondness (par=12)
	// e la Impulsiveness (par=13)
	if ((Nchannels == 1) | (lpAmp->StereoMode == 0))
	{
		float xmax, x, sumx, impx, impxmax;
		long i2;
		xmax = 0.0; sumx = 0.0;  impxmax = 0.0;
		for (i = 0; i < (Nbuf - 35); i++)
		{
			x = LBuf[i + Band * Nbuf];
			if (x > xmax) xmax = x;
			sumx += x;
			impx = 0.0;
			for (i2 = 0; i2 < 35; i2++) impx += LBuf[i + Band * Nbuf + i2]; //Impx è il valore medio su 35 ms, quindi è il valore Impulse
			if (impx > impxmax) impxmax = impx;
		}
		sumx = sumx / (Nbuf - 35); // Valore Medio Quadrato dei buffer millisecondici
		// 12: Millisecondness 
		cAP[12 + Band * NPAR] = 'y';
		fAP[12 + Band * NPAR] = (float)(10.0 * log10(xmax / sumx));
		// 13: Impulsiveness 
		impxmax = impxmax / 35; // Valore max della media Impulse su 35 ms;
		cAP[13 + Band * NPAR] = 'y';
		fAP[13 + Band * NPAR] = (float)(10.0 * log10(impxmax / sumx));

		if (Nchannels == 2)
		{
			xmax = 0.0; sumx = 0.0; impxmax = 0.0;
			for (i = 0; i < (Nbuf - 35); i++)
			{
				x = RBuf[i + Band * Nbuf];
				if (x > xmax) xmax = x;
				sumx += x;
				impx = 0.0;
				for (i2 = 0; i2 < 35; i2++) impx += RBuf[i + Band * Nbuf + i2]; //Impx è il valore medio su 35 ms, quindi è il valore Impulse
				if (impx > impxmax) impxmax = impx;
			}
			// 12: Millisecondness 
			sumx = sumx / (Nbuf - 35); // Valore Medio Quadrato dei buffer millisecondici
			cAP[12 + Band * NPAR + NPAR * NBANDS] = 'y';
			fAP[12 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(xmax / sumx));
			// 13: Impulsiveness
			impxmax = impxmax / 35; // Valore max della media Impulse su 35 ms;
			cAP[13 + Band * NPAR + NPAR * NBANDS] = 'y';
			fAP[13 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(impxmax / sumx));
		}
	}

	// aggiorno il progress meter (è 12 su 12 perchè ho 10 bande + la Wide + A = 12)
	// ProgressMeter(ci, 12, 12);
	// if (*ci->lpProgressCanceled) ok = FALSE;

	// memorizzo il SEL Left (S della banda LIN)
	lpAmp->SEL = fAP[NPAR * 10];

	// Ora posso fare la correzione White2Pink, senza spaciugare il calcolo della Strenght
	if (lpAmp->White2Pink)
	{
		double SA = 0.0, SL = 0.0, NA = 0.0, NL = 0.0;
		double AW[10];
		AW[0] = -39.4;
		AW[1] = -26.2;
		AW[2] = -16.1;
		AW[3] = -8.6;
		AW[4] = -3.2;
		AW[5] = 0.0;
		AW[6] = 1.2;
		AW[7] = 1.0;
		AW[8] = -1.1;
		AW[9] = -6.6;

		for (Band = 0; Band < 10; Band++)
		{
			SA += pow(10.0, (fAP[0 + Band * NPAR] + AW[Band]) / 10.0);
			NA += pow(10.0, (fAP[1 + Band * NPAR] + AW[Band]) / 10.0);
			SL += pow(10.0, fAP[0 + Band * NPAR] / 10.0);
			NL += pow(10.0, fAP[1 + Band * NPAR] / 10.0);
		}
		Band = 10; // A
		fAP[0 + Band * NPAR] = (float)(10.0 * log10(SA)); // S A
		fAP[1 + Band * NPAR] = (float)(10.0 * log10(NA)); // N A 
		Band = 11; // LIN
		fAP[0 + Band * NPAR] = (float)(10.0 * log10(SL)); // S L
		fAP[1 + Band * NPAR] = (float)(10.0 * log10(NL)); // N L 
		if (Nchannels == 2)
		{
			SA = 0.0; SL = 0.0; NA = 0.0; NL = 0.0;
			for (Band = 0; Band < 10; Band++)
			{
				SA += pow(10.0, (fAP[0 + Band * NPAR + NPAR * NBANDS] + AW[Band]) / 10.0);
				NA += pow(10.0, (fAP[1 + Band * NPAR + NPAR * NBANDS] + AW[Band]) / 10.0);
				SL += pow(10.0, fAP[0 + Band * NPAR + NPAR * NBANDS] / 10.0);
				NL += pow(10.0, fAP[1 + Band * NPAR + NPAR * NBANDS] / 10.0);
			}
			Band = 10; // A
			fAP[0 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(SA)); // S A
			fAP[1 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(NA)); // N A 
			Band = 11; // LIN
			fAP[0 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(SL)); // S L
			fAP[1 + Band * NPAR + NPAR * NBANDS] = (float)(10.0 * log10(NL)); // N L 
		}
	}


	// distruggo il progress meter, ho finito!
	// ProgressDestroy(ci);

	// N.B.: in uscita fIRFilt contiene lo Schroeder Plot non filtrato, utile per il display.

	// ritorno ok
	return (ok);
}


void SaveResults(drwav* wav,  char* TXTname, char* filename, char* cAP, float* fAP, int avgmode, MLS* lpAmp)
{
	CRLFDEF;
	VTABDEF;
	long			counter = 0;
	long			i, i2;
	long            Error;
	long            lpWritten = 0;
	char            ErrorString[256];
	char			Clipdata[8192];
	long			N = 0;
	long			nChan = wav->channels;
	char			ExistFlag = FALSE;
	char			Lin[] = "--";
	FILE *			hResults;
	char			bResult = TRUE;
	char* Bands[] = { "31.5 ", // 0
				   "63   ", // 1
				   "125  ", // 2
				   "250  ", // 3
				   "500  ", // 4
				   "1000 ", // 5
				   "2000 ", // 6
				   "4000 ", // 7
				   "8000 ", // 8
				   "16000", // 9
				   "A    ", //10
				   "Lin  " };//11

	char* Pars[] = { "Lsg", // 0
					"Lno", // 1
					"StG", // 2
					"C50", // 3
					"C80", // 4 
					"D50", // 5
					"ts ", // 6
					"EDT", // 7
					"Tus", // 8
					"T20", // 9
					"T30", //10
					"Jlf", //11
					"Jlfc", //12
					"Lj " };//13

	char* ParsMono[] = { "Lsg", // 0
				"Lno", // 1
				"StG", // 2
				"C50", // 3
				"C80", // 4 
				"D50", // 5
				"ts ", // 6
				"EDT", // 7
				"Tus", // 8
				"T20", // 9
				"T30", //10
				"Pks", //11
				"Mscn", //12
				"Imp" };//13

	char* ParsBIN[] = { "Lsg", // 0
			"Lno", // 1
			"StG", // 2
			"C50", // 3
			"C80", // 4 
			"D50", // 5
			"ts ", // 6
			"EDT", // 7
			"Tus", // 8
			"T20", // 9
			"T30", //10
			"IACCe", //11
			"tIACC", //12
			"wIACC" };//13


	// A questo punto, debbo scoprire se il file appendFile esiste già o meno.
	// ipotizzo quindi che esista già, e cerco di aprirlo in append, se poi non ci riesco lo creero'

	// verifico se il file esiste
    hResults= fopen(TXTname, "r");
	if (hResults == NULL)
	{
		// Il file non esiste, debbo crearlo ex-novo
		hResults = fopen(TXTname, "w"); 
		strcpy(ErrorString, TXTname);
		strcat(ErrorString, " did not exist\nA new file with this name was created!\n");
		printf(ErrorString);

		// Ora debbo scrivere l'intestazione del file (1a riga)
		strcpy(Clipdata, "Acoupar.exe - ISO3382 Acoustical Parameter File"); CRLF;
		fputs(Clipdata, hResults);

		// Ora debbo scrivere l'intestazione delle colonne (2a riga)
		// parametri: nome + parametro banda per banda
		strcpy(Clipdata, "Filename");
		for (i = 0; i < NPAR; i++)
			for (i2 = 0; i2 < NBANDS; i2++)
			{
				VTAB;
				if((lpAmp->wChannels == 1) | (lpAmp->StereoMode == 0))
				{
					// mono parameters
					strcat(Clipdata, ParsMono[i]);
				}
				else if ((lpAmp->StereoMode == 1) | (lpAmp->StereoMode == 2) | (lpAmp->StereoMode == 3))
				{
					// WY parameters
					strcat(Clipdata, Pars[i]);
				}
				else if (lpAmp->StereoMode == 4)
				{
					// BIN parameters
					strcat(Clipdata, ParsBIN[i]);
				}
				strcat(Clipdata, "_");
				strcat(Clipdata, Bands[i2]);
			}
		CRLF;
		fputs(Clipdata, hResults);
	}
	else
	{
		//The file was existing, and it has been opened for reading, so I close it
		fclose(hResults);
		// and I reopen it for append
		hResults = fopen(TXTname, "a");
	}

	// ora aggiungo i risultati dell'elaborazione appena conclusa
	// copio i parametri nella riga di testo Clipdata
	strcpy(Clipdata, filename); //copy filename at beginning of each row
	for (int par = 0; par < NPAR; par++)
		for (int band = 0; band < NBANDS; band++)
		{
			float tmp = 0.0;
			char Lin[] = "--";
			char textbuf[20];
			long counter = 0;
			long c = 0;

			VTAB;
			// average data among channels
			if (!avgmode) 
			{
				nChan = 1;
			}
			for (c = 0; c < nChan; c++)
				if (cAP[NBANDS * NPAR * c + (band * NPAR) + par] == 'y')
				{
					tmp += fAP[NBANDS * NPAR * c + (band * NPAR) + par];
					counter += 1;
				}
			if (counter > 0)
			{
				tmp = (float)(tmp / counter);
				sprintf(textbuf, "%.3f", tmp);
				strcat(Clipdata, textbuf);
			}
			else strcat(Clipdata, Lin);
		} //next band
	CRLF;

	// Ora scrivo la stringa su file.
	fputs(Clipdata, hResults);

	strcpy(ErrorString, "Data written to file ");
	strcat(ErrorString, TXTname);
	strcat(ErrorString, "\n");
	printf(ErrorString);

	fclose(hResults);
}
