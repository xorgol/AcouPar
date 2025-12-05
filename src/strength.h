/************************************************/
/*                Strength.h                    */
/************************************************/

// Use an Internal Structure for your function, and access it through a handle
typedef struct
{
	float ref1;    // Liv. riferimento sorgente a 10m - 31.5 Hz
	float ref2;
	float ref3;
	float ref4;
	float ref5;
	float ref6;
	float ref7;
	float ref8;
	float ref9;
	float ref10;   // Liv. riferimento sorgente a 10m - 16 kHz
	float ref11;   // Liv. riferimento sorgente a 10m - Lin
	float ref12;   // Liv. riferimento sorgente a 10m - A
	float fFS;     // Full Scale in dB-SPL
	float Distance;
	float fRTUdBstart;
	float fRTUdBend;
	float Space;
	float Speed;
	float fThreshold;
	int cNoiseC;
	int cEDT2;
	int cBand;
	int cChannel;
	long fAP;
	long cAP;
	long LeftBuf;
	long RightBuf;
	int FormatNum;
	int StereoMode; // 0=2 mikes, 1=Soundfield, 2=Omni/Eight, 3=p-p probe, 4=binaural
	int IACCmode;   // 0=Early, 1=Late, 2=Whole
	int StageMode;  // 0=Normal, 1=Stage (compute ST instead of ts,D)
	int AVGmode;	 // 0=Normal, 1=Average Mode , 250 to 2k
	int AppendMode; // 0=not append, 1=append
	int White2Pink; // 0=leave unchanged, 1=correct
	long lSampRate;
	long lN;
	float MaxL;
	float MaxR;
	float SEL;
	int wChannels;
	//long fIR;		// puntatore alla IR corrente (mono o left)
	//long fIRFilt;   // puntatore alla IR corrente filtrata (mono o left)
	//long fIR2;		// puntatore alla IR corrente (right)
	//long fIRFilt2;  // puntatore alla IR corrente filtrata (right)
	//char LastIRDir[256];
	//char cName[256];
	//char IRFilePath[256];
	char wTitle[256];
	//char appendFile[256];
}
MLS;

#define NPAR 14		// numero parametri: S,N,G,C50,C80,D,Ts,EDT,TU,T20,T30,LF,LFC,LG
#define NBANDS 12	// numero bande: 31.5 ... 16 kHz (10 bande) + A (11^a) + LIN (12^a)           
#define PARLEN 4    
#define CRLFDEF char crlf[]={0x0d,0x0a,0}
#define CRLF strcat(Clipdata,crlf)
#define VTABDEF char vtab[]={0x09,0}
#define VTAB strcat(Clipdata,vtab)
#define TRUE 1
#define FALSE 0

