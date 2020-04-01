                        PROGRAM ALTIMETR;
{        PROGRAM FOR CALULATION DIFFERENT FUNCTION

 1) INOSPHERIC CORRECTION    5 FEBRUARY 2020
 2) DRY TROPOSPHIRIC CORRECTION 7 FEBRUARY 2020
 3) WET TPOPOSHERIC CORRECTION 7 FEBRUARY 2020
 4) CORRECTION OF THE ZEMLI 13 FEBRUARY 2020
 5) speckle-noise  21 february 2020
 6) INFLUENCE OF ZEMLI OF RCS 22/02/2020
 7) CALCULATION OF EMB - 25/2/2020
 8) CALCULATION OF WIND SPEED - 26/02/2020 + CHAPRON
 9) CALCULATION OF IMPULS - AMORACH - 9/3/2020

}

  USES
       GREAD_V,USEFUL_V;

  VAR
     FALT,TEC0,TECK,HTEC,TEC,HCOR,ATM,FIW_LAT,WET,WET0,T : EXTENDED;
     RZ,TECM, K1,X,RZEMLI,SWH,G,ASMALL,BSMALL,DH1 : EXTENDED;
     A0,A1,A2,A0_1,A0_2,A0_3,A0_4,A0_5,A1_1,A1_2,A1_3,A1_4,A1_5 : EXTENDED;
     A2_1,A2_2,A2_3,A2_4,A2_5,DH2,DH3,DH4,DH5,DH6,DH7 : EXTENDED;
     A,B,AA_1,B1,B2,AI_0,AI_1,AI_2,AI_3,AI_4,AI_5 : EXTENDED;
     BA_1,AA_2,BA_2,U10_1,U10_2,U10_3,B1_1,U10_4 : EXTENDED;
     TAU0,TAU_H,SIGMA_P,TET0,HZEM,RSPUTNIK,HSMALL,DX,GAMMA : EXTENDED;
     SIGMA_H,C,SIGMA_C,DELTA,TEK,BETA,US,ALFA,VS,A_BIG,IMP1 : EXTENDED;
     IMP2,IMP3,ALFA1,ALFA2,US1,VS1,VS2,US2,IMP4 : EXTENDED;

     NTEC, I,I0,NM, N,NI,ITAU : INTEGER;

     F : TEXT;
     FNAME, STR : STRING;

     LOG_ION, LOG_TROP,LOG_IMP : BOOLEAN;

{
 --------------------------------------------------------------------------
}
         FUNCTION  IONOSPHERA(F0, NTC : EXTENDED) : EXTENDED;
{             CALCULATION OF CORRECTION  5/02-20
}
    VAR
        A,Y : EXTENDED;
    BEGIN
        A:= 40.3;

        Y:= A*NTC/SQR(F0);

        IONOSPHERA:= Y;
    END;

{
 --------------------------------------------------------------------------
}
         FUNCTION  TROPOSPHERA_D(PA, FI : EXTENDED) : EXTENDED;
{             CALCULATION OF DRY CORRECTION  7/02-20
        INPUT PA V MM RTUTNOGO STOLBA
}
    VAR
        A,Y,KBAR : EXTENDED;
    BEGIN
        KBAR:= 750.06156;
        A:= 2.277;

        Y:= A*(PA/KBAR)*( 1.0 + 0.0026*COS(2.0*FI) );

        TROPOSPHERA_D:= Y;
    END;

{
 --------------------------------------------------------------------------
}
         FUNCTION  TROPOSPHERA_W(WT : EXTENDED) : EXTENDED;
{             CALCULATION OF WET CORRECTION  12/02-20
        WET VVODITSYA V G/SM**2
        OTVET POLUCHAEM V SM
}
    VAR
        Y : EXTENDED;
    BEGIN

        Y:= 6.36*WT;

        TROPOSPHERA_W:= Y;
    END;

{
 --------------------------------------------------------------------------
}
         FUNCTION  TROPOSPHERA_WT(WT, TEMP : EXTENDED) : EXTENDED;
{             CALCULATION OF WET CORRECTION  12/02-20
        WET VVODITSYA V KG/M**2
        OTVET POLUCHAEM V M,  S UCHETOM TEMPERATU
}
    VAR
        Y : EXTENDED;
    BEGIN

        Y:= 1.723*WT/TEMP;

        TROPOSPHERA_WT:= Y;
    END;

{
 --------------------------------------------------------------------------
}
         FUNCTION ZEMLYA(XT, RZ : EXTENDED) : EXTENDED;
{             CALCULATION OF ZEMLYA CORRECTION  13/02-20
}
    VAR
        Y : EXTENDED;
    BEGIN

        Y:= SQR(XT)  / (2.0*RZ);

        ZEMLYA:= Y;
    END;

{
 --------------------------------------------------------------------------
}
         FUNCTION FORM(XT,KT,TI : EXTENDED; N0,NIMP : INTEGER) : EXTENDED;
{             CALCULATION OF SPECKL  20/02-20
}
    VAR
        Y,I0,FIT,A0, SUM : EXTENDED;
        I,J : INTEGER;

    BEGIN
        IF XT <= TI
        THEN A0:= XT*KT
        ELSE A0:= TI*KT;


        RANDOMIZE;
        SUM:= 0.0;
        FOR J:= 1 TO NIMP
        DO BEGIN

              Y:= 0.0;

              FOR I:= 1 TO N0
              DO BEGIN
                   FIT:= RANDOM(256);

                   FIT:= (FIT/255.0)*(2.0*PI);

                   Y:= Y + A0*SIN(FIT);

{WRITELN(' FIT ',FIT*180/PI:6:2,' I ',I,' Y ',Y:12);}

                END;

              Y:= SQR(Y/N0);

              SUM:= SUM + Y;
           END;


        FORM:= SUM/NIMP;
    END;


{
 --------------------------------------------------------------------------
}
               FUNCTION WIND_F1(SIGMA0,SWH : EXTENDED) : EXTENDED;
{               model CHAPRONA F1
}
   VAR
      W2_1,W2_2,B2,MKU,BKU,MSWH,BSWH,P_1,P_2: EXTENDED;
      W1_1,W1_2,W02_1,W02_2,B1_1,B1_2,SIGMAKU : EXTENDED;
      W1P_1,W1P_2,X_1,X_2,W2X,U,BWIND,MWIND : EXTENDED;

   BEGIN
      W02_1:= 0.5401201;
      W02_2:= 10.4048140;

      B2:= -2.2838729;

      W1_1:= -33.9506170;
      W1_2:= -11.0339400;
      W2_1:= -3.9342847;
      W2_2:= -0.0583444;

      B1_1:= 18.0637810;
      B1_2:= -0.3722814;

      SIGMAKU:= SIGMA0;
      MKU:= 0.0690915;
      BKU:= -0.3433598;
      MSWH:= 0.0637450;
      BSWH:= 0.0872510;

      P_1:= MKU*SIGMAKU + BKU;
      P_2:= MSWH*SWH + BSWH;

      W1P_1:= W1_1*P_1 + W1_2*P_2;
      W1P_2:= W2_1*P_1 + W2_2*P_2;

      X_1:= 1.0 / (1.0 + EXPO(-(W1P_1+B1_1)));
      X_2:= 1.0 / (1.0 + EXPO(-(W1P_2+B1_2)));

      W2X:= W02_1*X_1 + W02_2*X_2;
      U:= 1.0 / (1.0+ EXPO(-(W2X+B2)));

      BWIND:= 0.1;
      MWIND:= 0.0284394;

      WIND_F1:= (U - BWIND)/MWIND;

   END;


    FUNCTION  ERF_1(X : REAL) : REAL;
{ FUNCTION ERF 22/06-2015 }
      VAR
        Y,A1,A2,A3,P,T,SIGN : REAL;
      BEGIN
         A1:= 0.3480242;
         A2:= -0.0958798;
         A3:= 0.7478556;
         P:= 0.47047;
         SIGN:= 1;

         IF X < 0
         THEN BEGIN
                 X:= -X;
                 SIGN:= -1.0;
              END;
         T:= 1.0/(1.0+P*X);

         Y:= 1.0 - (A1*T + A2*T*T + A3*T*T*T)*EXP(-X*X);

         ERF_1:= Y*SIGN;
      END;


    FUNCTION  IMPULS_1(TEK,SIGMA_P,TET0,HZEM,RSPUTNIK,DX,SWH,
                        C,A0 : EXTENDED) : EXTENDED;
{        FUNCTION FORMA IMPULSA 11/03-2020
            VSE UGLY V RADIANAH }

      VAR
         HSMALL, GAMMA,SIGMA_H,SIGMA_C,G,DELTA,BETA,ALFA : EXTENDED;
         US,VS,A_BIG,Y : EXTENDED;

      BEGIN

           HSMALL:= RSPUTNIK*(1.0 + RSPUTNIK/HZEM);

           GAMMA:= 2.0*SQR( SIN(DX/2.0) )/LN(2.0);

           SIGMA_H:= SWH/4.0;

           SIGMA_C:= SQRT( SQR(SIGMA_P) + SQR(2.0*SIGMA_H/C) );

           G:= 9.81;

           DELTA:= 4.0*C*COS(2.0*TET0)/(GAMMA*HSMALL);
           BETA:= 4.0*SQRT(C/HSMALL)*SIN(2.0*TET0)/GAMMA;

           ALFA:= DELTA - SQR(BETA)/4.0;

           US:= (TEK-ALFA*SQR(SIGMA_C))/( SQRT(2.0)*SIGMA_C);
           VS:= ALFA*( TEK-ALFA*SQR(SIGMA_C)/2.0 );

           A_BIG:= A0*EXPO( -4.0*SQR(SIN(TET0))/GAMMA );

           Y:= A_BIG*EXPO(-VS)*( 1.0 + ERF_1(US) );

           IMPULS_1:= Y;


      END;


      FUNCTION  IMPULS_2(TEK,SIGMA_P,TET0,HZEM,RSPUTNIK,DX,SWH,
                        C,A0 : EXTENDED) : EXTENDED;
{        FUNCTION FORMA IMPULSA 11/03-2020
            VSE UGLY V RADIANAH - FORMULA AMORACH  }

      VAR
         HSMALL,GAMMA,SIGMA_H,SIGMA_C,G,DELTA,BETA,ALFA1,ALFA2 : EXTENDED;
         US1,US2,VS1,VS2,A_BIG,Y : EXTENDED;

      BEGIN

           HSMALL:= RSPUTNIK*(1.0 + RSPUTNIK/HZEM);

           GAMMA:= 2.0*SQR( SIN(DX/2.0) )/LN(2.0);

           SIGMA_H:= SWH/4.0;

           SIGMA_C:= SQRT( SQR(SIGMA_P) + SQR(2.0*SIGMA_H/C) );

           G:= 9.81;

           DELTA:= 4.0*C*COS(2.0*TET0)/(GAMMA*HSMALL);
           BETA:= 4.0*SQRT(C/HSMALL)*SIN(2.0*TET0)/GAMMA;

           ALFA1:= DELTA - SQR(BETA)/8.0;
           ALFA2:= DELTA;

           US1:= (TEK-ALFA1*SQR(SIGMA_C))/( SQRT(2.0)*SIGMA_C );
           US2:= (TEK-ALFA2*SQR(SIGMA_C))/( SQRT(2.0)*SIGMA_C );

           VS1:= ALFA1*( TEK-ALFA1*SQR(SIGMA_C)/2.0 );
           VS2:= ALFA2*( TEK-ALFA2*SQR(SIGMA_C)/2.0 );

           A_BIG:= A0*EXPO( -4.0*SQR(SIN(TET0))/GAMMA );

           Y:=  A_BIG*EXPO(-VS1)*(1.0 + ERF_1(US1) )
                        -0.5*A_BIG*EXPO(-VS2)*(1.0 + ERF_1(US2) );
           Y:= 2.0*Y;

           IMPULS_2:= Y;

      END;


{***************************************************************************}
{                             ƒ‹€‚€Ÿ —€‘’œ                                 }

 BEGIN

WRITELN;
WRITELN('                ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
WRITELN('                ³       Žƒ€ŒŒ€ ALTIMETR.PAS    ³');
WRITELN('                ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
WRITELN;

{repeat
    RANDOMIZE;
    I0:= RANDOM(256);

    IF I0 = 0
    THEN X:= 0
    ELSE X:= I0/255.0;

    WRITELN(' I0 ',I0,' X = ',(2.0*PI*X)*180.0/PI:6:2);
until NOT antwrt(' CONTINUE ? ');  }

   IF ANT(' REWRITTEN FILE ? ')
   THEN BEGIN
           ASSIGN(F,'ALTIMETR.DAT');
           REWRITE(F);
        END
   ELSE HALT;

   LOG_ION:= ANTWRT(' IONOSHERA CORRECTIONS ? ');

   IF LOG_ION
   THEN BEGIN
           FALT:= RREAD(' INPUT FREQUENCY , GHz ');
           FALT:= FALT * 1000000000.0;

           TEC0:= RREAD(' INPUT TEC0 ');
           TECK:= RREAD(' INPUT TECK ');
           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);

           TEC0:= TEC0*10000000000000000.0;
           TECK:= TECK*10000000000000000.0;
           HTEC:= HTEC*10000000000000000.0;

           WRITELN(F,'  TEC  ','   HCOR_m   ','  HCOR_CM  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 TEC:= TEC0 + (I-1)*HTEC;

                 HCOR:= IONOSPHERA(FALT,TEC);

                 TEC:= TEC/10000000000000000.0;
                 WRITELN(F,TEC:5:2,'   ',HCOR:7:5,'   ',HCOR*100.0:6:4);

              END;

        END;

   LOG_TROP:= ANTWRT(' DRY TROPOSHERA CORRECTIONS ? ');

   IF LOG_TROP
   THEN BEGIN
           FIW_LAT:= RREAD(' INPUT LATITUDE ');
           FIW_LAT:= FIW_LAT*PI/180.0;

           TEC0:= RREAD(' INPUT P0 (mm rtutnogo stolba) ');
           TECK:= RREAD(' INPUT PK (mm rtutnogo stolba) ');
           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);

           WRITELN(F,'  P_ATMSP  ','   HCOR_m   ','  HCOR_CM  ',
                     '   P_BAR   ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 ATM:= TEC0 + (I-1)*HTEC;

                 HCOR:= TROPOSPHERA_D(ATM,FIW_LAT);

                 WRITELN(F,ATM:5:2,'   ',HCOR:7:5,'   ',HCOR*100.0:6:4,
                          '    ',ATM/750.06156:6:5);

              END;

        END;

   LOG_TROP:= ANTWRT(' WET TROPOSHERA CORRECTION 1 ? ');

   IF LOG_TROP
   THEN BEGIN
           TEC0:= RREAD(' INPUT WET_0 (GRAM/CM) ');
           TECK:= RREAD(' INPUT WET_  (GRAM/CM) ');
           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);

           WRITELN(F,'  WET  ','   HCOR_m   ','  HCOR_CM  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 WET:= TEC0 + (I-1)*HTEC;

                 HCOR:= TROPOSPHERA_W(WET);


                 WRITELN(F,WET:5:2,'   ',HCOR/100.0:7:5,'   ',HCOR:6:4);

              END;

        END;

   LOG_TROP:= ANTWRT(' WET TROPOSHERA CORRECTION 2 ? ');

   IF LOG_TROP
   THEN BEGIN
           TEC0:= RREAD(' INPUT TEMP_0 (K) ');
           TECK:= RREAD(' INPUT TEMP_K (K) ');
           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);
           WET0:= RREAD(' INPUT WET , KG/M**2  ');

           WRITELN(F,'  WET  ','   HCOR_m   ','  HCOR_CM  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 T:= TEC0 + (I-1)*HTEC;

                 HCOR:= TROPOSPHERA_WT(WET0,T);

                 WRITELN(F,T:5:2,'   ',HCOR:7:5,'   ',HCOR*100.0:6:4);

              END;

        END;

   LOG_TROP:= ANTWRT(' RADIUS ZEMLI ? ');

   IF LOG_TROP
   THEN BEGIN
           TEC0:= RREAD(' INPUT X0, M ');
           TECK:= RREAD(' INPUT XK, M ');
           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);
           RZ:= RREAD(' INPUT R ZEMLI (6371 KM)  ');
           RZ:= RZ * 1000.0;

           WRITELN(F,'  XT  ','   HCOR_m   ','  HCOR_CM  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 T:= TEC0 + (I-1)*HTEC;

                 HCOR:= ZEMLYA(T,RZ);

                 WRITELN(F,T:5:2,'   ',HCOR:7:5,'   ',HCOR*100.0:6:4);

              END;

        END;

   LOG_TROP:= ANTWRT(' SPECKLE NOISE ? ');

   IF LOG_TROP
   THEN BEGIN
           TEC0:= RREAD(' INPUT T0, S ');
           TECM:= RREAD(' IMPUT MIDLE TM, S ');
           TECK:= RREAD(' INPUT TK, M ');
           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);

           K1:= RREAD(' INPUT K1  ');
           NM:= IREAD(' INPUT POINTS NM ');
           NI:= IREAD(' INPUT N IMP ');

           WRITELN(F,'  T  ','   IMP   ','  IM_SUM  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 T:= TEC0 + (I-1)*HTEC;

                 HCOR:= FORM(T,K1,TECM,NM,NI);

                 IF T <= TECM
                 THEN X:= T*K1
                 ELSE X:= TECM*K1;

                 WRITELN(F,T:5:2,'   ',X:12,'     ',SQRT(HCOR));

              END;

        END;

   LOG_TROP:= ANTWRT(' INFLUENCE R_ZEMLI ON RCS ? ');

   IF LOG_TROP
   THEN BEGIN
           TEC0:= RREAD(' INPUT H0, KM  ');
           TECK:= RREAD(' INPUT HK, KM  ');
           RZEMLI:= RREAD(' INPUT R_ZEMLI, KM (6371) ');

           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);

           WRITELN(F,'  H0  ','   SIGMA_0  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 T:= TEC0 + (I-1)*HTEC;

                 HCOR:= 10.0*LOG10(1.0 + T/RZEMLI);

                 WRITELN(F,T:5:2,'   ',HCOR);

              END;

        END;


   LOG_TROP:= ANTWRT(' CALCULATION OF EMB ? ');

   IF LOG_TROP
   THEN BEGIN
           TEC0:= RREAD(' INPUT U10_0, M/S  ');
           TECK:= RREAD(' INPUT U10_K, M/S  ');
           SWH:= RREAD(' INPUT SWH, M ');

           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);
           G:= 9.81;

           ASMALL:= 3.21; {RREAD(' INPUT COEFFICIENT A (3.21) ');}
           BSMALL:= 0.62; {RREAD(' INPUT COEFFICIENT B (0.62) ');}

           A0:= -0.0146;
           A1:= -0.00215;
           A2:= -0.00389;

           A0_1:= 0.00666;
           A1_1:= 0.0015;

           A0_2:= 0.0327;
           A2_2:= -0.0022;

           A0_3:= 0.0245;
           A1_3:= 0.00122;
           A2_3:= -0.0034;

           A1_4:= -0.66;
           A2_4:= -0.15;

           A0_5:= -0.0029;
           A1_5:= -0.0038;
           A2_5:= 0.000155;

           WRITELN(F,'  U10  ','   DH1  ','  DH2  ','   DH3  ',
                     '  DH4  ','  DH5  ','  DH6   ','  DH7  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 T:= TEC0 + (I-1)*HTEC;

                 DH1:= ASMALL*POWER(G*SWH/T,BSMALL);

                 DH2:= -(A0 + A1*T+A2*SWH);
                 DH3:= A0_1 + A1_1*T;
                 DH4:= A0_2 + A2_2*SWH;
                 DH5:= A0_3 + A1_3*T + A2_3*SWH;
                 DH6:= -(A1_4*SWH + A2_4*SWH*T);
                 DH7:= -(A0_5 + A1_5*T + A2_5*T*T);

                 WRITELN(F,T:5:2,'   ',DH1:7:5,'   ',DH2*100.0:7:5,'   ',
                          DH3*100.0:7:5,'   ',DH4*100.0:7:5,'   ',
                          DH5*100.0:7:5,'   ',
                          DH6:7:5,'   ',DH7*100.0:7:5);

              END;

        END;

   LOG_TROP:= ANTWRT(' CALCULATION OF U10 (ONE PARAMETER) ? ');

   IF LOG_TROP
   THEN BEGIN
           TEC0:= RREAD(' INPUT SIGMA_0, dB  ');
           TECK:= RREAD(' INPUT SIGMA_K, dB  ');
           SWH:= RREAD(' INPUT SWH, M ');

           NTEC:= IREAD(' INPUT N ');
           HTEC:= (TECK-TEC0)/(NTEC-1);
           G:= 9.81;

           A:= 0.080074;
           B:= -0.124651;

           A1:= 0.039893;
           B1:= -0.031996;

           A2:= 0.01595;
           B2:= 0.017215;

           AI_0:= -15.383;
           AI_1:= 16.077;
           AI_2:= -2.305;
           AI_3:= 0.09896;
           AI_4:= 0.00018;
           AI_5:= -0.00006414;

           AA_1:= 46.5;
           BA_1:= 3.6;

           AA_2:= 1690;
           BA_2:= 0.5;


           WRITELN(F,'  SIGMA_0  ','   U10_1  ','  U10_2  ','   U10_3  ',
                     '  U10_CHAPR  ');

           FOR I:= 1 TO NTEC
           DO BEGIN
                 T:= TEC0 + (I-1)*HTEC;

                 IF T < 10.12
                 THEN BEGIN
                         U10_1:= EXPO( (POWER(10.0,-(0.21+T/10.0))- B)/A );
                      END
                 ELSE IF T <= 10.9
                      THEN BEGIN
                             U10_1:=EXPO( (POWER(10.0,-(0.21+T/10.0))- B1)
                                                                       /A1 );
                           END
                      ELSE BEGIN
                             U10_1:=EXPO( (POWER(10.0,-(0.21+T/10.0))- B2)
                                                                       /A2 );
                           END;

                 IF (8.0 <= T) AND (T<= 15.0)
                 THEN BEGIN
                         U10_2:= AI_0 + AI_1*T + AI_2*T*T + AI_3*T*T*T
                                      + AI_4*SQR(T*T) + AI_5*SQR(T*T)*T;
                      END
                 ELSE BEGIN
                         WRITELN(' SIGMA OUT SIDE ');
                         READLN;
                         HALT;
                      END;

                 IF (T <= 10.917)
                 THEN BEGIN
                         U10_3:= AA_1 - BA_1*T;
                      END
                 ELSE BEGIN
                         U10_3:= AA_2*EXPO(-BA_2*T);
                      END;


                 U10_4:= WIND_F1(T,SWH);


                 WRITELN(F,T:5:2,'   ',U10_1:7:5,'   ',U10_2:7:5,'   ',
                          U10_3:7:5,'   ',U10_4:7:5);

              END;

        END;

   LOG_IMP:= ANTWRT(' CALCULATION OF IMPULS SHAPE (BROWN) ? ');
   IF LOG_IMP
   THEN BEGIN
           TAU0:= RREAD(' INPUT TAU_0, NS  ');
           TAU0:= TAU0*0.000000001;
           ITAU:= IREAD(' INPUT N POINTS ');
           TAU_H:= RREAD(' INPUT TAU STEP, NS ');
           TAU_H:= TAU_H*0.000000001;

           IF ANTWRT(' SIGMA_P - BROWN ? ')
           THEN BEGIN
                   WRITELN(' COEF = 0.425 ');
                   SIGMA_P:= 0.425*TAU_H;
                END
           ELSE BEGIN
                   WRITELN(' COEF = 0.513 ');
                   SIGMA_P:= 0.513*TAU_H;
                END;

           TET0:= RREAD(' INPUT INCIDENCE ANGLE, GRAD. ');
           TET0:= TET0*PI/180.0;

           HZEM:= RREAD(' INPUT  R ZEMLI, KM (6370) ');
           HZEM:= HZEM*1000.0;

           RSPUTNIK:= RREAD(' INPUT  H SPUTNIK, KM ');
           RSPUTNIK:= RSPUTNIK*1000.0;

           HSMALL:= RSPUTNIK*(1.0 + RSPUTNIK/HZEM);

           DX:= RREAD(' INPUT ANTENNA BEAM, GRAD. ');
           DX:= DX*PI/180.0;
           GAMMA:= 2.0*SQR( SIN(DX/2.0) )/LN(2.0);

           SWH:= RREAD(' INPUT SWH, M ');
           SIGMA_H:= SWH/4.0;

           C:= RREAD(' INPUT C, KM/S ');
           C:= C * 1000.0;

           SIGMA_C:= SQRT( SQR(SIGMA_P) + SQR(2.0*SIGMA_H/C) );
           A0:= RREAD(' INPUT A0 ');


           G:= 9.81;

           WRITELN(F,' TAU  ','  IMP_1  ','  IMP_2 ','   IMP3  ','  IMP4  ',
                     '  I  ');

           FOR I:= 1 TO ITAU
           DO BEGIN
                 TEK:= TAU0 + (I-1)*TAU_H;

                 DELTA:= 4.0*C*COS(2.0*TET0)/(GAMMA*HSMALL);
                 BETA:= 4.0*SQRT(C/HSMALL)*SIN(2.0*TET0)/GAMMA;

                 ALFA:=  DELTA - SQR(BETA)/4.0;
                 ALFA1:= DELTA - SQR(BETA)/8.0;
                 ALFA2:= DELTA;

                 US:=  (TEK- ALFA*SQR(SIGMA_C))/( SQRT(2.0)*SIGMA_C );
                 US1:= (TEK-ALFA1*SQR(SIGMA_C))/( SQRT(2.0)*SIGMA_C );
                 US2:= (TEK-ALFA2*SQR(SIGMA_C))/( SQRT(2.0)*SIGMA_C );

                 VS:=  ALFA* ( TEK- ALFA*SQR(SIGMA_C)/2.0 );
                 VS1:= ALFA1*( TEK-ALFA1*SQR(SIGMA_C)/2.0 );
                 VS2:= ALFA2*( TEK-ALFA2*SQR(SIGMA_C)/2.0 );

                 A_BIG:= A0*EXPO( -4.0*SQR(SIN(TET0))/GAMMA );

                 IMP1:=  A_BIG*EXPO(-VS)*( 1.0 + ERF_1(US) );
                 IMP2:= IMPULS_1(TEK,SIGMA_P,TET0,HZEM,RSPUTNIK,DX,SWH,C,A0);
                 IMP3:=  A_BIG*EXPO(-VS1)*(1.0 + ERF_1(US1) )
                        -0.5*A_BIG*EXPO(-VS2)*(1.0 + ERF_1(US2) );
                 IMP3:= 2.0*IMP3;
                 IMP4:= IMPULS_2(TEK,SIGMA_P,TET0,HZEM,RSPUTNIK,DX,SWH,C,A0);

                 WRITELN(F,TEK:12,'   ',IMP1:12,'   ',IMP2:12,'  ',IMP3:12,
                           '    ',IMP4:12,'   ',I);

              END;

        END;

   WRITELN('    ‘—…’     ‡€ŠŽ—…,            ENTER ');
   READLN;
   CLOSE(F);

 END.
