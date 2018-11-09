#include "FFT.h"
#include<ESP8266WiFi.h>
#include "string.h"
#include <Ticker.h>

/* Identification Wifi  */

#define C_Conversion 3.3/1023 // 10 bits de conversion ADC nodeMCU

#define  ESID_WIFI "iPhone de Yassine"
#define  MDP_WIFI "yassine124"
/* Creation object Ticker& client MAIL */
Ticker ticker;



/* Global variables  */
const uint16_t samples = 256; //This value MUST ALWAYS be a power of 2
const double maxTimeAquire =0.5;/*Temps de mesure en seconde*/
const double samplingFrequency = samples/maxTimeAquire;

double vReal[samples]={0.0};
double vImag[samples]={0.0};
double majorPeak=0;


uint32_t indexAquire =0;
bool isFull =false; // Pour debug mode

void mesureFFT();
//void sendDataSocket(double* Vreal,double* Vimag,uint16_t samlples );
void AquireData();




void setup()
{
  WiFi.begin(ESID_WIFI,MDP_WIFI); /* Init Connection */
  Serial.begin(115200); /* Init Serial communication protocol*/
  while(WiFi.status()!= WL_CONNECTED){ /* wait until wifi connected status */
    delay(200);
    Serial.print('.');
  }

  Serial.print(F("Wifi Connected to : "));
  Serial.print(ESID_WIFI);
  Serial.print(F("Adresse IP :  "));
  Serial.println(WiFi.localIP());

  ticker.attach_ms((maxTimeAquire/samples)*1000, AquireData);

}

void loop()
{

}

/*
    Name : AquireData
    Description :Cette fonction mesure les (samples ) points, est il calcule le fréquence et le l'amplitude de peak Major
    Return :  Null
*/
void AquireData(){
  if(indexAquire == samples){
    if(!isFull){
      mesureFFT();

      isFull =true; // Debug mode (lancer mesure une seule fois )
      Serial.println(F("[DEBUG] Fin mesure"));
      //sendDataSocket(vReal, vImag, samples);
    }else {
      //Serial.println(F("[DEBUG] nothing "));

    }

  }else {
   vReal[indexAquire] = analogRead(A0)*C_Conversion;
    //system_adc_read_fast((uint16_t*)vReal,samples,8);
    //for(int i=0;i<samples;i++) Serial.println(vReal[i]);
   //Serial.print(F("[DEBUG] Analog value : "));
   //Serial.println(vReal[indexAquire]);
   indexAquire ++;
  }
}
/*
    Name : mesure FFT
    Description :Mesurer la fréquence du signal ansi l'amlitude de peak max
    Return :Null
*/
void mesureFFT(){

  FFT fft = FFT();
  double amplitude=0;
  fft.FFT_PACK(vReal,samples,FFT_WIN_TYP_HAMMING,&majorPeak,&amplitude);
  Serial.print(F("Major Peak equal : "));
  Serial.println(majorPeak);
  Serial.print(F("Major Peak amplitude equal : "));
  Serial.println(amplitude);
}

/*
void sendDataSocket(double* Vreal,double* Vimag,uint16_t samlples ){
  WiFiClient client;
  if(!client.connect("172.20.10.3", 1555)){
    Serial.println(F("Connection Failed"));
    delay(100);
    return ;
  }
  Serial.println(F("[DEBUG] connection to server successful"));

  for(int i=0;i<samlples;i++){

    //strcpy(dataTosend,Vreal_char + "," +Vimag_char);
    client.print(Vreal[i]);
    //client.write(",");
    //client.write(Vimag[i]);
    //client.write("\n");
  }

  client.println("END");
  client.stop();

}*/
