
#include<stdio.h>
#include<stdlib.h>
#include<DSPF_sp_minval.h>
#include<DSPF_sp_cfftr2_dit.h>
#include<DSPF_sp_icfftr2_dif.h>
#include<DSPF_sp_bitrev_cplx.h>
#include<math.h>
#include<string.h>
void function_carre_256_float_to_float (float *);
#define NBECH  19559
#define LENGTH_BUFFER 256
#define N 256
#define M_PI 3.14159265358979323846
const float Tuckey_Windows[256] = { 0.0, 0.015101531982495253,
    0.059493902857107761, 0.13049554138967034, 0.2238175135197471,
    0.33382260026017013, 0.45386582026834893, 0.57669582743934256,
    0.69489293664633955, 0.80131731818962815, 0.88954028726283518,
    0.95423263590976182, 0.991486549841951, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    0.99148654984195117, 0.95423263590976126, 0.88954028726283463,
    0.80131731818962781, 0.69489293664633955, 0.576695827439343,
    0.45386582026834982, 0.33382260026017141, 0.22381751351974855,
    0.13049554138966946, 0.059493902857107317, 0.015101531982495142, 0.0 };
float buffer_in[N*2];
float buffer_out[N*2];
float w[N], w_save[N], Echantillon_Precedent[256], Bruit_temp[256], E_Bruit[256], PS_Bruit[256], E_argument[256], temp_float256[256], temp_phase_float256[256], PS_Ech_Deb[256], stop[256], Echantillon[256], buffer_out_256[256], output[256];
float s, alpha0,alpha,Beta, SNR_trame, rms_Signal, rms_Bruit;
float Tab_ech[80][256];

typedef struct complexe {
   float re, img;
} complexe;

void bit_rev(float* x, int n)
{
     int i, j, k;
     float rtemp, itemp;

     j = 0;
     for(i=1; i < (n-1); i++)
        {
        k = n >> 1;
        while(k <= j)
             {
             j -= k;
             k >>= 1;
             }
        j += k;
        if(i < j)
          {
           rtemp    = x[j*2];
           x[j*2]   = x[i*2];
           x[i*2]   = rtemp;
           itemp    = x[j*2+1];
           x[j*2+1] = x[i*2+1];
           x[i*2+1] = itemp;
          }
        }
}

int gen_twiddle(float *w, int n)
{
    float delta = 2 * M_PI / n;
    int i;
    for(i = 0; i < n/2; i++)
    {
        w[2 * i + 1] = sin(i * delta);
        w[2 * i] = cos(i * delta);
    }

 return n;
}

void bitrev_index(short *index, int n)
{
    int   i, j, k, radix = 2;
    short nbits, nbot, ntop, ndiff, n2, raddiv2;

    nbits = 0;
    i = n;
    while (i > 1)
    {
        i = i >> 1;
        nbits++;
    }

    raddiv2 = radix >> 1;
    nbot    = nbits >> raddiv2;
    nbot    = nbot << raddiv2 - 1;
    ndiff   = nbits & raddiv2;
    ntop    = nbot + ndiff;
    n2      = 1 << ntop;

    index[0] = 0;
    for ( i = 1, j = n2/radix + 1; i < n2 - 1; i++)
    {
        index[i] = j - 1;

        for (k = n2/radix; k*(radix-1) < j; k /= radix)
            j -= k*(radix-1);

        j += k;
    }
    index[n2 - 1] = n2 - 1;
}

void abs_256_complexe_to_256_float(float *x, float *y)
{
  int k;
  float a;
  for (k = 0; k < 512; k+=2) {
    a = x[k]*x[k] + x[k+1]*x[k+1];
    y[k/2] = (float)sqrt(a);
  }
}

/*
void function_carre_256_float_to_float(float *x){
	int k;
	for (k = 0; k < 256; k++) {
	    x[k] = x[k] * x[k];
	  }
}
*/
void argument_256_complexe_to_256_float(float *x, float *y)
{
  int k;
  for (k = 0; k < 512; k+=2) {
    y[k/2] = atan2(x[k+1], x[k]);
  }
}

float moyenne_256_float(float x[256])
{
  float y;
  int k;
  y = x[0];
  for (k = 0; k < 255; k++) {
    y += x[k + 1];
  }

  y /= 256.0F;
  return y;
}

float rms(const float x[256])
{
  float y;
  int i;
  float b_x[256];
  for (i = 0; i < 256; i++) {
    b_x[i] = x[i] * x[i];
  }

  y = b_x[0];
  for (i = 0; i < 255; i++) {
    y += b_x[i + 1];
  }

  return (float)sqrt(y / 256.0F);
}
/*
void cast_log(float *x)
{
  *x = (float)log(*x);
}

void cast_sqrt(float *x)
{
	*x = (float)sqrt(*x);
}


void fft(float *x, float *w, short n){
	//int size;
	int size = gen_twiddle(w, N);
	bit_rev(w, N>>1);
	DSPF_sp_cfftr2_dit(x, w, N);
}*/

int main(void)
{
FILE 	*fp, *fp2;

// Modifier selon votre chemin
char 	name[] = "C:\\signal_bruite.txt";
int 	i, j,k, nb_blocs, nb_blocs_noise, counter;
short	table[16];


// Definition des limites de blocs
//nb_blocs = NBECH/LENGTH_BUFFER;
nb_blocs = ((NBECH-256)/242)+1; // Taille du signal - taille premiere trame / taille des trames avec overlapping +1 (pour la premiere trame)
nb_blocs_noise = 7;
counter=0;

//Variables du projet
Beta = 0.05;

printf("%s \n", name);
if( (fp = fopen(name, "r+")) == NULL)
{
    printf("No file\n");
    exit(1);
}

// Initialisation des variables
for(i=0; i<256; i++){
	Bruit_temp[i] = 0;
	E_Bruit[i] = 0;
	PS_Bruit[i] = 0;
	E_argument[i] = 0;
	temp_float256[i]=0;
	PS_Ech_Deb[i]=0;
	buffer_in[i]=0;
}
for(i=0; i<512; i++){
	buffer_in[i]=0;
	buffer_out[i]=0;
}

//Generation du twiddle table
bitrev_index(table,N);
gen_twiddle(w, N);
bit_rev(w,N>>1);

//Sauvegarde de la twiddle table
for(i=0; i<256; i++){
	w_save[i] = w[i];
}

////////////////////////////////
// Creation du bruit moyen
////////////////////////////////
for(i=0; i<nb_blocs_noise; i++){
	for(j=0; j<512; j++){ 	//initialisation des buffers
		buffer_in[j]=0;
		buffer_out[j]=0;
	}
	for(j=0; j<256; j++){	//recuperation de w
		w[j] = w_save[j];
	}
	for(j=0; j<512; j+=2){	//lecture du bloc
		fscanf(fp,"%f", &buffer_in[j]);
		buffer_in[j+1] = 0;
	}
	DSPF_sp_cfftr2_dit(buffer_in, w, N);
	DSPF_sp_bitrev_cplx((double*)buffer_in,table,N);
	// Ici buffer_in contient les valeurs exploitables de la FFT : equivalent au code matlab :fft(buffer_in)
	abs_256_complexe_to_256_float(buffer_in,temp_float256); 	//calcul valeurs absolues des complexes
	function_carre_256_float_to_float(temp_float256);			//calcul des valeur puissance 2

	//Addition des valeurs du bruit de chaque bloc
	for(j=0; j<LENGTH_BUFFER; j++)
		Bruit_temp[j] = Bruit_temp[j] + temp_float256[j];
}

for(j=0; j<LENGTH_BUFFER; j++){
	PS_Bruit[j]= Bruit_temp[j] / i;
}

rewind(fp); //Remet le pointeur de lecture du fichier au debut du fichier


////////////////////////////////
//Calcul alpha0 et s
////////////////////////////////

alpha0 = 0;

for (i=1; i<8 ;i++){

		for(j=0; j<256; j++){	//recuperation de w
			w[j] = w_save[j];
		}
		for(j=0; j<512; j+=2){	//lecture du bloc
			fscanf(fp,"%f", &buffer_in[j]);
			buffer_in[j+1] = 0;
		}
		DSPF_sp_cfftr2_dit(buffer_in, w, N);
		DSPF_sp_bitrev_cplx((double*)buffer_in,table,N);
		abs_256_complexe_to_256_float(buffer_in,temp_float256); 	//calcul valeurs absolues des complexes
	    function_carre_256_float_to_float(temp_float256);			//calcul des valeur puissance 2
	    for (alpha=1; alpha < 10; alpha++){
	    		for (k=0; k<LENGTH_BUFFER;k++){    				    // PS_Ech_Deb = temp_float256 - Alpha PS_bruit
	    			PS_Ech_Deb[k] = temp_float256[k]- alpha*PS_Bruit[k];
	    			if (PS_Ech_Deb[k] < Beta*PS_Bruit[k]){
	    				PS_Ech_Deb[k] = Beta*PS_Bruit[k];
	    			}
	    			stop[k] = Beta*PS_Bruit[k];
	    		}
	    		if (moyenne_256_float(PS_Ech_Deb)==moyenne_256_float(stop)){
	    			alpha0 += alpha;
	    			break;
	    		}
	    }


}

alpha0 = alpha0/(i-1);

//Calcul de s

s = 1/((1-alpha0)/20);

rewind(fp); //Remet le pointeur de lecture du fichier au debut du fichier

////////////////////////////////
// Debruitage du premier bloc
////////////////////////////////
for(j=0; j<1024; j++){	//lecture du bloc
	fscanf(fp,"%f", &Echantillon[1]);
}
//fseek (fp, 1024, SEEK_SET); //4*256=1024
for(j=0; j<256; j++){	//lecture du bloc
	fscanf(fp,"%f", &Echantillon[j]);
}
rewind(fp);
rms_Bruit = rms(Echantillon);
rms_Bruit = rms_Bruit * rms_Bruit;

//Calcul de la valeur du signal pour le SNR
for(j=0; j<256; j++){	//lecture du bloc
	fscanf(fp,"%f", &Echantillon[j]);
	Echantillon_Precedent[j] = Echantillon[j];
	counter++;
}
//rewind(fp);
rms_Signal = rms(Echantillon);
rms_Signal = rms_Signal * rms_Signal;

//Calcul du SNR
SNR_trame = 10 * log(rms_Signal/rms_Bruit);
if (SNR_trame < -5){
	SNR_trame = -5;
}
if (SNR_trame > 20){
	SNR_trame = 20;
}
alpha = alpha0 - SNR_trame/s;


//Debruitage
for(j=0; j<512; j+=2){	//lecture du bloc
	buffer_in[j] = Echantillon [j/2];
	buffer_in[j+1] = 0;
}
DSPF_sp_cfftr2_dit(buffer_in, w, N);
DSPF_sp_bitrev_cplx((double*)buffer_in,table,N);
abs_256_complexe_to_256_float(buffer_in,temp_float256); 	//calcul valeurs absolues des complexes
function_carre_256_float_to_float(temp_float256);	//PS_Echantillon
argument_256_complexe_to_256_float(buffer_in, temp_phase_float256);

for(j=0; j<256; j++){
	PS_Ech_Deb[j] = temp_float256[j] - alpha * PS_Bruit[j];
	if (PS_Ech_Deb[j] < Beta*PS_Bruit[j]){
		PS_Ech_Deb[j] = Beta*PS_Bruit[j];
	}
	PS_Ech_Deb[j] = sqrt(PS_Ech_Deb[j]); //Spectre_E_Debruite
}
//Reconstruction complexe
for(j=0; j<512; j+=2){	//lecture du bloc
	buffer_out[j] = PS_Ech_Deb[j/2]*cos(temp_phase_float256[j/2]);
	buffer_out[j+1] = PS_Ech_Deb[j/2]*sin(temp_phase_float256[j/2]);
}
// ifft
DSPF_sp_bitrev_cplx((double*)buffer_out,table,N);
DSPF_sp_icfftr2_dif(buffer_out, w, N);
for(j=0; j<512; j++){
	buffer_out[j]=buffer_out[j]/N;
}
// Tuckey_windows
for(j=0; j<512; j+=2){
	Tab_ech[0][j/2]=buffer_out[j]*Tuckey_Windows[j/2];
}

fp2 = fopen("C:\\output.txt", "w+");
/////////////////////////////
// Recuperation des echantillons de 256 overlappes
/////////////////////////////
for(i=1; i<nb_blocs; i++)
{

	for(j=0; j<14;j++){
		Tab_ech[i][j]=Echantillon_Precedent[j+(256-14)];
	}
	//fseek (fp, (i*256)-14, SEEK_SET); //On place le pointeur de lecture au 2eme ech - overlapping
	for(j=14; j<256; j++){	//lecture du bloc overlappe
		fscanf(fp,"%f", &Tab_ech[i][j]);
		counter++;
	}
	for(j=0; j<256; j++){	//lecture du bloc overlappe
		Echantillon_Precedent[j] = Tab_ech[i][j];
	}
	//for(j=14; j<256; j++){	//lecture du bloc
		//fscanf(fp,"%f", &Tab_ech[i][j]);
	//}
	//rewind(fp);

	////////////////////
	//Script Débruitage
	////////////////////

	//Calcul de la valeur du signal pour le SNR
	rms_Signal = rms(Tab_ech[i]);
	rms_Signal = rms_Signal * rms_Signal;

	//Calcul du SNR
	SNR_trame = 10 * log(rms_Signal/rms_Bruit);
	if (SNR_trame < -5){
		SNR_trame = -5;
	}
	if (SNR_trame > 20){
		SNR_trame = 20;
	}
	alpha = alpha0 - SNR_trame/s;

	//Debruitage
	for(j=0; j<512; j+=2){	//lecture du bloc
		buffer_in[j] = Tab_ech[i][j/2];
		buffer_in[j+1] = 0;
	}
	DSPF_sp_cfftr2_dit(buffer_in, w, N);
	DSPF_sp_bitrev_cplx((double*)buffer_in,table,N);
	abs_256_complexe_to_256_float(buffer_in,temp_float256); 	//calcul valeurs absolues des complexes
	function_carre_256_float_to_float(temp_float256);	//PS_Echantillon
	argument_256_complexe_to_256_float(buffer_in, temp_phase_float256);

	for(j=0; j<256; j++){
		PS_Ech_Deb[j] = temp_float256[j] - alpha * PS_Bruit[j];
		if (PS_Ech_Deb[j] < Beta*PS_Bruit[j]){
			PS_Ech_Deb[j] = Beta*PS_Bruit[j];
		}
		PS_Ech_Deb[j] = sqrt(PS_Ech_Deb[j]); //Spectre_E_Debruite
	}
	//Reconstruction complexe
	for(j=0; j<512; j+=2){	//lecture du bloc
		buffer_out[j] = PS_Ech_Deb[j/2]*cos(temp_phase_float256[j/2]);
		buffer_out[j+1] = PS_Ech_Deb[j/2]*sin(temp_phase_float256[j/2]);
	}
	// ifft
	DSPF_sp_bitrev_cplx((double*)buffer_out,table,N);
	DSPF_sp_icfftr2_dif(buffer_out, w, N);
	for(j=0; j<512; j++){
		buffer_out[j]=buffer_out[j]/N;
	}
	// Tuckey_windows
	for(j=0; j<512; j+=2){
		Tab_ech[i][j/2]=buffer_out[j]*Tuckey_Windows[j/2];
	}
	// Prise en compte overlapping
	for(j=242; j<256;j++){
		Tab_ech[i-1][j] += Tab_ech[i][j-242];
	}

	// Reconstruction de la sortie
	for(j=14; j<LENGTH_BUFFER; j++){
		output[j] = Tab_ech[i-1][j];
		fprintf(fp2, "%f \n", output[j]);
	}
}

fclose(fp2);

printf("Last sample = %d\n",nb_blocs*j);

if(feof(fp))
	printf("End of file\n");
else
	printf("Can't read all file\n");

fclose(fp);
return 0;
}


//FFT Reference : https://www.dsprelated.com/showthread/c6x/3086-1.php


