//
//  RINGsize2.7_C.c
//  
//
//  Created by Miyako Oyadomari on 2019/07/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// global variables from python script
char g_string1[256],g_string2[256],g_string3[256];
int g_InitialCP_X, g_InitialCP_Y, g_Initial_R, g_search_range;
double g_flux_lo, g_flux_up;



void share_string(char* string1, char* string2, char* string3){
    printf("******The C-library is called******\n");
    printf ("-----Module: share_string ------\n");
    strcpy(g_string1,string1);
    strcpy(g_string2,string2);
    strcpy(g_string3,string3);
    printf("shared_string1: %s\n", g_string1);
    printf("shared_string2: %s\n", g_string2);
    printf("shared_string3: %s\n", g_string3);
    printf("*****************END***************\n");
}

void share_int(int P1, int P2, int P3, int P4){
    printf("******The C-library is called******\n");
    printf ("-----Module: share_int ------\n");
    g_InitialCP_X = P1; g_InitialCP_Y = P2;
    g_Initial_R = P3; g_search_range = P4;
    printf("integer g_InitialCP_X: %d\n", P1);
    printf("integer g_InitialCP_Y: %d\n", P2);
    printf("integer g_Initial_R: %d\n", P3);
    printf("integer g_search_range: %d\n", P4);
    printf("*****************END***************\n");
}

void share_float(float P1, float P2){
    printf("******The C-library is called******\n");
    printf ("-----Module: share_float ------\n");
    g_flux_lo = P1; g_flux_up = P2;
    printf("float g_flux_lo: %f\n", P1);
    printf("float g_flux_up: %f\n", P2);
    printf("*****************END***************\n");
}


struct RESULTS{
    int CP_X;
    int CP_Y;
    int R_median;
    int R_lo;
    int R_up;
    float D1;
    float D2;
    float D2D1;
    float D_D2D1;
};

struct CP_His{
    int CP_X;
    int CP_Y;
    int CP_count;
};



void CalculateRadius(float **matrix, int task_no, float flux_total) {
    printf("******The C-library is called ******\n");
    printf ("-----Module: CalculateRadius ------\n");

    int CP_X, CP_Y; // Center position in a searching area
    float D_D2D1;

    //Open the output file for "RESULTS extractD2D1"
    printf("# Open an output file\n");
    FILE *fp;
    if (task_no == 4){
        printf("# file_nema: %s\n", g_string1);
        fp = fopen(g_string1, "a");
    }
    if (task_no == 5){
        printf("# file_nema: %s\n", g_string2);
        fp = fopen(g_string2, "a");
    }
    if (fp == NULL) {
        printf("can not open the file\n");
    }

    
    
    printf ("Iteration: Sweep the searching area\n");
    int min_D1 = g_Initial_R;
    float min_D_D2D1 = 0.5;
    
    struct RESULTS whole[(g_search_range +1)*(g_search_range +1)];
    struct RESULTS extractD1[(g_search_range +1)*(g_search_range +1)];
    struct RESULTS extractD2D1[(g_search_range +1)*(g_search_range +1)];

    // ----- Iteration: Sweep the searching area ------
    int num = 0;
    
    for (int dy=0; dy <= g_search_range; dy++){
        CP_Y = g_InitialCP_Y - (g_search_range/2) + dy;
        for (int dx=0; dx <= g_search_range; dx++){
        CP_X = g_InitialCP_X - (g_search_range/2) + dx;
    
    
    // ----- Slice the matrix ------
    int k=0; int l=0;
    float slice[2*g_Initial_R+1][2*g_Initial_R+1];
    for(int i=(CP_Y - g_Initial_R); i<=(CP_Y + g_Initial_R); i++){
        for(int j=(CP_X - g_Initial_R); j<=(CP_X + g_Initial_R); j++){
            slice[k][l] = matrix[i][j];
            l++;
        }
        l=0;k++;
    }     // ----- End of Slice the matrix ------
    
    
    int R_median, R_lo, R_up;
    float D1, D2, D2D1, rate_max;
    float min_D_rate = 0.5;
    float min_D_rate_lo = 0.5;
    float min_D_rate_up = 0.5;


    // ----- Iteration: R=Initial radius to R=1 ------

    for (int radius=g_Initial_R; radius> 1; radius--){
    
    // ----- Calcurate the flux in the circle ------
    int RR; float flux_IN, flux_rate, D_rate;
    flux_IN=0.0;
    for (int k=0; k<= 2*g_Initial_R; k++){
        for (int l=0; l<= 2*g_Initial_R; l++){
            RR = ((l - g_Initial_R)*(l - g_Initial_R)) + ((k - g_Initial_R)*(k - g_Initial_R));
            if (RR <= (radius*radius)){
                flux_IN = flux_IN + slice[k][l];
            }
        }
    }  // ------ END of Calcuration:flux in the circle ---------
        
    flux_rate = flux_IN/flux_total;
//    printf( "R=%d, flux_rate=%f\n", radius,flux_rate);

    // ----- Flux rate at the R maximun ------
    if (radius == g_Initial_R){
        rate_max = flux_rate;
    }
        
    // ----- Radius median ------
//    D_rate = fabsf(0.5 - flux_rate);
    D_rate = fabs(0.5 - flux_rate);
    if (min_D_rate > D_rate){
        min_D_rate = D_rate;
        R_median = radius;
    }

    // ----- Radius lower limit ------
//    D_rate = fabsf(g_flux_lo - flux_rate);
    D_rate = fabs(g_flux_lo - flux_rate);
    if (min_D_rate_lo > D_rate){
        min_D_rate_lo = D_rate;
        R_lo = radius;
    }
    // ----- Radius upper limit ------
//    D_rate = fabsf(g_flux_up - flux_rate);
    D_rate = fabs(g_flux_up - flux_rate);
    if (min_D_rate_up > D_rate){
        min_D_rate_up = D_rate;
        R_up = radius;
    }
        
    } // ----- END of Iteration: R=Initial radius to R=1 ------


    // ----- judge the effectiveness of R _up -----
    // If the rate_max does not reach the g_flux_up,
    // the results at the Centet Position is ignored.
    if (rate_max < g_flux_up){
        continue;
    }
            
    D1 = R_up - R_lo; D2 = R_median - R_lo;
    D2D1 = D2 / D1;
//    D_D2D1 = fabsf(0.5 - D2D1);
    D_D2D1 = fabs(0.5 - D2D1);

    // ----- D1 minimum ------
    if (min_D1 > D1){
        min_D1 = D1;
    }

    // ----- struct RESULTS wole initialize ----
    whole[num].CP_X = CP_X;
    whole[num].CP_Y = CP_Y;
    whole[num].R_median = R_median;
    whole[num].R_lo = R_lo;
    whole[num].R_up = R_up;
    whole[num].D1 = D1;
    whole[num].D2 = D2;
    whole[num].D2D1 = D2D1;
    whole[num].D_D2D1 = D_D2D1;
    num++;

    } // --- END of Iteration: Sweep the searching area (loop for dx)
    CP_X = g_InitialCP_X - (g_search_range/2); // reset the CP_X
    } // --- END of Iteration: Sweep the searching area (loop for dy)
    printf("--- \n");

    for (int i=0; i < num; i++){
        printf("%d,%d,%d,%d,%d,%f,%f,%f,%f\n",whole[i].CP_X,whole[i].CP_Y,whole[i].R_median,whole[i].R_lo,whole[i].R_up,whole[i].D1, whole[i].D2,whole[i].D2D1,whole[i].D_D2D1);
    }
    printf("min_D1=%d\n", min_D1);
    printf("NUM=%d\n", num);

    // --- Extract the D1 minimu from RESULTS whole ---
    int id1 = 0;
    for (int i=0; i < num; i++){
        if (whole[i].D1 == min_D1){
        extractD1[id1].CP_X = whole[i].CP_X;
        extractD1[id1].CP_Y = whole[i].CP_Y;
        extractD1[id1].R_median = whole[i].R_median;
        extractD1[id1].R_lo = whole[i].R_lo;
        extractD1[id1].R_up = whole[i].R_up;
        extractD1[id1].D1 = whole[i].D1;
        extractD1[id1].D2 = whole[i].D2;
        extractD1[id1].D2D1 = whole[i].D2D1;
        extractD1[id1].D_D2D1 = whole[i].D_D2D1;
            
            // ----- D2D1 ~ 0.5 ------
            if (min_D_D2D1 > whole[i].D_D2D1){
                min_D_D2D1 = whole[i].D_D2D1;
            }

        id1++;
        }
    } // --- END of Extract the D1 minimu from RESULTS whole ---
    
    for (int i=0; i < id1; i++){
        printf("%d,%d,%d,%d,%d,%f,%f,%f,%f\n",extractD1[i].CP_X,extractD1[i].CP_Y,extractD1[i].R_median,extractD1[i].R_lo,extractD1[i].R_up,extractD1[i].D1, extractD1[i].D2,extractD1[i].D2D1,extractD1[i].D_D2D1);
    }
    printf("id1=%d\n",id1);
    printf("min_D_D2D1=%f\n",min_D_D2D1);

    // --- Extract the D2D1~0.5 from RESULTS extractD1 ---
    int id2 = 0;
    for (int i=0; i < id1; i++){
        if (extractD1[i].D_D2D1 == min_D_D2D1){
            extractD2D1[id2].CP_X = extractD1[i].CP_X;
            extractD2D1[id2].CP_Y = extractD1[i].CP_Y;
            extractD2D1[id2].R_median = extractD1[i].R_median;
            extractD2D1[id2].R_lo = extractD1[i].R_lo;
            extractD2D1[id2].R_up = extractD1[i].R_up;
            extractD2D1[id2].D1 = extractD1[i].D1;
            extractD2D1[id2].D2 = extractD1[i].D2;
            extractD2D1[id2].D2D1 = extractD1[i].D2D1;
            extractD2D1[id2].D_D2D1 = extractD1[i].D_D2D1;
            id2++;
        }
    } // Extract the D2D1~0.5 from RESULTS extractD1  ---

    
    
    printf("id2=%d\n",id2);
    for (int i=0; i < id2; i++){
    printf("%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n",i,extractD2D1[i].CP_X,extractD2D1[i].CP_Y,extractD2D1[i].R_median,extractD2D1[i].R_lo,extractD2D1[i].R_up,extractD2D1[i].D1, extractD2D1[i].D2,extractD2D1[i].D2D1,extractD2D1[i].D_D2D1);
        
    fprintf(fp,"%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n",i,extractD2D1[i].CP_X,extractD2D1[i].CP_Y,extractD2D1[i].R_median,extractD2D1[i].R_lo,extractD2D1[i].R_up,extractD2D1[i].D1, extractD2D1[i].D2,extractD2D1[i].D2D1,extractD2D1[i].D_D2D1);
    }
    printf("id2=%d\n",id2);



    fclose(fp);

    printf("*****************END***************\n");

    
} // ----- END of Module: CalculateRadius -----


void HistCal() {
    printf("******The C-library is called ******\n");
    printf ("-----Module: HistCal ------\n");
    
    
    //Open the output file for "RESULTS extractD2D1"
    FILE *fp;
    char line[256];
    int line_num = 0;
    
    // count the lines in the result file
    printf("# file_nema: %s\n", g_string2);
    fp = fopen(g_string2, "r");
    if (fp == NULL) {
        printf("can not open the file\n");
    }

    while(fgets(line, 256, fp) != NULL) {
        line_num++;
    }
    fclose(fp);

    
    printf("line number in tmep_CenterPosition.txt: %d\n", line_num);
    struct RESULTS CenterPosition[line_num];

    
    // Make a file of Center position
    printf("# file_nema: %s\n", g_string2);
    fp = fopen(g_string2, "r");
    if (fp == NULL) {
        printf("can not open the file\n");
    }
    
    float data[10]; int ret; int num = 0;
    while( (ret=fscanf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &data[0], &data[1], &data[2], &data[3], &data[4], &data[5], &data[6], &data[7], &data[8], &data[9])) != EOF){
        CenterPosition[num].CP_X = (int)data[1];
        CenterPosition[num].CP_Y = (int)data[2];
        num++;
    }
    
    fclose(fp);

    // Minimum & Maximum X-Position
    // Minimum & Maximum Y-Position
    int CP_X_max = CenterPosition[0].CP_X;
    int CP_X_min = CenterPosition[0].CP_X;
    int CP_Y_max = CenterPosition[0].CP_Y;
    int CP_Y_min = CenterPosition[0].CP_Y;
    for (int i=0; i < line_num; i++){
        //printf("%d %d %d\n", i, CenterPosition[i].CP_X,CenterPosition[i].CP_Y);
        if (CP_X_max < CenterPosition[i].CP_X){
            CP_X_max = CenterPosition[i].CP_X;
        }
        if (CP_X_min > CenterPosition[i].CP_X){
            CP_X_min = CenterPosition[i].CP_X;
        }
        if (CP_Y_max < CenterPosition[i].CP_Y){
            CP_Y_max = CenterPosition[i].CP_Y;
        }
        if (CP_Y_min > CenterPosition[i].CP_Y){
            CP_Y_min = CenterPosition[i].CP_Y;
        }
    }
    
//    printf("CP_X_max=%d\n",CP_X_max);
//    printf("CP_X_min=%d\n",CP_X_min);
//    printf("CP_Y_max=%d\n",CP_Y_max);
//    printf("CP_Y_min=%d\n",CP_Y_min);
    

    // Make CP_histgram data
    struct CP_His CP_histgram[(CP_X_max-CP_X_min+1)*(CP_Y_max-CP_Y_min+1)];

    int sum_i = 0; int hist_num = 0;
    for(int x=CP_X_min; x <= CP_X_max; x++){
        for (int y=CP_Y_min; y <= CP_Y_max; y++){
            for (int i=0; i < line_num; i++){
                if ((CenterPosition[i].CP_X == x) && (CenterPosition[i].CP_Y == y) ){
                    hist_num++;
                }
            }
            CP_histgram[sum_i].CP_X = x;
            CP_histgram[sum_i].CP_Y = y;
            CP_histgram[sum_i].CP_count = hist_num;
        //printf("%d,%d,%d\n",CP_histgram[sum_i].CP_X,CP_histgram[sum_i].CP_Y, CP_histgram[sum_i].CP_count);
            hist_num = 0;
            sum_i++;
        }
    }// -- End of Make CP_histgram data
    
    // Make a file for CP_histgram data
    printf("# file_nema: %s\n", g_string3);
    fp = fopen(g_string3, "w");
    if (fp == NULL) {
        printf("can not open the file\n");
    }
    int area_sum;
    int count_left, count_right, count_up, count_lo;
    for (int i=0; i < sum_i; i++){
        // count_left
        int count_left_X = CP_histgram[i].CP_X - 1;
        int count_left_Y = CP_histgram[i].CP_Y;
        if (count_left_X < CP_X_min){
            count_left = 0;
        }
        else{
            for (int j=0; j < sum_i; j++){
                if ((CP_histgram[j].CP_X == count_left_X) && (CP_histgram[j].CP_Y == count_left_Y)){
                    count_left = CP_histgram[j].CP_count;
                }
            }
        }

        // count_right
        int count_right_X = CP_histgram[i].CP_X + 1;
        int count_right_Y = CP_histgram[i].CP_Y;
        if (count_right_X > CP_X_max){
            count_right = 0;
        }
        else{
            for (int j=0; j < sum_i; j++){
                if ((CP_histgram[j].CP_X == count_right_X) && (CP_histgram[j].CP_Y == count_right_Y)){
                    count_right = CP_histgram[j].CP_count;
                }
            }
        }

        // count_up
        int count_up_X = CP_histgram[i].CP_X;
        int count_up_Y = CP_histgram[i].CP_Y - 1;
        if (count_up_Y < CP_Y_min){
            count_up = 0;
        }
        else{
            for (int j=0; j < sum_i; j++){
                if ((CP_histgram[j].CP_X == count_up_X) && (CP_histgram[j].CP_Y == count_up_Y)){
                    count_up = CP_histgram[j].CP_count;
                }
            }
        }

        // count_lo
        int count_lo_X = CP_histgram[i].CP_X;
        int count_lo_Y = CP_histgram[i].CP_Y + 1;
        if (count_lo_Y > CP_Y_max){
            count_lo = 0;
        }
        else{
            for (int j=0; j < sum_i; j++){
                if ((CP_histgram[j].CP_X == count_lo_X) && (CP_histgram[j].CP_Y == count_lo_Y)){
                    count_lo = CP_histgram[j].CP_count;
                }
            }
        }

        area_sum = CP_histgram[i].CP_count + count_left + count_right + count_up + count_lo;
        
     // printf("%d,%d,%d,%d,left(%d,%d,%d),right(%d,%d,%d),up(%d,%d,%d),lo(%d,%d,%d)\n",CP_histgram[i].CP_X,CP_histgram[i].CP_Y, CP_histgram[i].CP_count, area_sum, count_left_X,count_left_Y,count_left,count_right_X,count_right_Y,count_right,count_up_X,count_up_Y,count_up,count_lo_X,count_lo_Y,count_lo);
        fprintf(fp,"%d,%d,%d,%d\n",CP_histgram[i].CP_X,CP_histgram[i].CP_Y, CP_histgram[i].CP_count, area_sum);
    }
    fclose(fp);
    printf("*****************END***************\n");

}

