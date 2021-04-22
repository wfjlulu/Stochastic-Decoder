//shuffle partially_MTFM for (2048,1723) (32,6) LDPC code
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"mex.h"

#define CN_DEGREE 32
#define VN_DEGREE 6
#define CN 384
#define VN 2048

void d3_subunit_mex(int in1, int in2, int in3,int *IM_REG, int rand01, int *b);
void d2_subunit_mex( int b1, int b2, int * r, int *s );

void Shuffled_MTFM_mex ( double *LLR, double N0, double * H , double * codeword )
{
    const int It_Max = 500;
    const double BETA = 1.0/16;
    const double alpha = 3.0;
    const double Y = 6.0;
    double LLR_NDS[VN];
    double Pn[VN];
    double rand01[VN][2];
    int rand_IM[VN][8];
    int ChannelBits[VN];
    int Hc[CN][VN] = {0};
    int ColumIndex[VN_DEGREE][VN];
    int RowIndex[CN][CN_DEGREE];
    int CN2VN[CN][VN] = {0};
    int VN2CN[CN][VN] = {0};
    int IM_REG [VN][16] = {0};
    double Pt_TFM[VN];
    int NG = 128;
    int G = VN/NG;

        // find Hc from H
    for ( int i = 0 ;i < CN ; i++)  //  loop
    {
        for ( int j = 0; j < VN; j++)  //  loop
        {
            Hc[i][j] = (int)H[j*CN+i];  // Hc is defined as int 
        }
    }

    // get index for cn and vn update
    // different from the matlab , index of c begins from 0, be careful!
    for ( int i = 0; i < VN; i++ )
    {   
        int count = 0;
        for ( int j = 0; j < CN; j++ )
        {   
            if ( Hc[j][i] == 1)
            {
                ColumIndex[count][i] = j;
                count = count + 1;
            }
        }
    }

    for ( int i = 0; i< CN; i++)
    {
        int count = 0;
        for ( int j = 0; j < VN; j++)
        {
            if ( Hc[i][j] == 1)
            {
                RowIndex[i][count] = j;
                count = count + 1;
            }
        }
    }

    //NDS and get bit stream
    for ( int i = 0; i< VN ; i++)
    {
        LLR_NDS[i] = alpha*N0/Y*LLR[i];
        Pn[i] = 1.0 / ( 1.0 + exp(LLR_NDS[i]) );
        Pt_TFM[i] = Pn[i];
    }
    
    // generate channel bit and random number
    srand((unsigned int)time(NULL));    //reset seed every fram
    for (int iter = 0; iter < It_Max; iter++) //iteration loop
    {
        for ( int i = 0;i < VN ;i++)
        {
            rand01[i][0] = rand() / ( RAND_MAX + 1.0 ); //for channel bit
            rand01[i][1] = rand() / ( RAND_MAX + 1.0 ); //for TFM
            for (int k = 0; k < 8; k++ )               //for 8 im_reg
            {
                rand_IM[i][k] = (rand() % 2) ;
            }
            ChannelBits[i] = ( Pn[i] > rand01[i][0] );
        }

        // VN2CN initial by channel
        if ( iter == 0)
        {
            for ( int i = 0;i < CN; i++)
            {
                for ( int j = 0; j < VN ; j++)
                {
                    VN2CN[i][j] = ChannelBits[j];
                }
            }
        }

        for ( int g = 0 ; g < G; g++ )   // group loop 
        {
            // horizontical step  ( cn update )
            for ( int n = g*NG ; n < ( g +1 )* NG ; n++ )   // current group to be update
            {
                for ( int m = 0 ; m < VN_DEGREE ; m++ )  // the cn number to be updated
                {
                    int flag ;
                    int cn_index = ColumIndex[m][n];
                    flag = VN2CN [cn_index] [ RowIndex[cn_index][0] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][1] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][2] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][3] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][4] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][5] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][6] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][7] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][8] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][9] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][10] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][11] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][12] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][13] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][14] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][15] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][16] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][17] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][18] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][19] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][20] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][21] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][22] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][23] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][24] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][25] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][26] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][27] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][28] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][29] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][30] ]^
                           VN2CN [cn_index] [ RowIndex[cn_index][31] ];
                    for ( int current_vn = 0 ; current_vn < CN_DEGREE ; current_vn++ )
                    {
                        if ( RowIndex[cn_index][current_vn] >= g*NG && RowIndex[cn_index][current_vn] < (g+1)*NG )
                        {
                            CN2VN [cn_index] [ RowIndex[cn_index][current_vn] ] = flag ^ VN2CN [cn_index] [ RowIndex[cn_index][current_vn] ];
                        }

                    }


                }

            }
            
            // vertical step ( vn update )
            for ( int i = g*NG ; i < ( g + 1 )*NG ; i++ )
            {
                int in1,in2,in3,in4,in5,in6;

                in1 = CN2VN[ColumIndex[0][i]][i];
                in2 = CN2VN[ColumIndex[1][i]][i];
                in3 = CN2VN[ColumIndex[2][i]][i];
                in4 = CN2VN[ColumIndex[3][i]][i];
                in5 = CN2VN[ColumIndex[4][i]][i];
                in6 = CN2VN[ColumIndex[5][i]][i];
    
                int b[8];
                int r[6];
                int s[6];
    
                d3_subunit_mex(ChannelBits[i],in2,in3,&(IM_REG[i][0]),rand_IM[i][0],&(b[0]));
                d3_subunit_mex(ChannelBits[i],in1,in3,&(IM_REG[i][2]),rand_IM[i][1],&(b[1]));
                d3_subunit_mex(ChannelBits[i],in1,in2,&(IM_REG[i][4]),rand_IM[i][2],&(b[2]));
                d3_subunit_mex(ChannelBits[i],in5,in6,&(IM_REG[i][6]),rand_IM[i][3],&(b[3]));
                d3_subunit_mex(ChannelBits[i],in4,in6,&(IM_REG[i][8]),rand_IM[i][4],&(b[4]));
                d3_subunit_mex(ChannelBits[i],in4,in5,&(IM_REG[i][10]),rand_IM[i][5],&(b[5]));
                d3_subunit_mex(in4,in5,in6,&(IM_REG[i][12]),rand_IM[i][6],&(b[6]));
                d3_subunit_mex(in1,in2,in3,&(IM_REG[i][14]),rand_IM[i][7],&(b[7]));
    
                d2_subunit_mex(b[0],b[6],&r[0],&s[0]);
                d2_subunit_mex(b[1],b[6],&r[1],&s[1]);
                d2_subunit_mex(b[2],b[6],&r[2],&s[2]);
                d2_subunit_mex(b[3],b[7],&r[3],&s[3]);
                d2_subunit_mex(b[4],b[7],&r[4],&s[4]);
                d2_subunit_mex(b[5],b[7],&r[5],&s[5]);
    
                int St;
                int U;
                St = s[0] + s[1] +s[2] + s[3] +s[4] + s[5] ;  // edges updated
                U = St;
                //printf("U = %d\n",U);
    
                int Xt;
                Xt = r[0]*s[0] + r[1]*s[1] +r[2]*s[2] +r[3]*s[3] +r[4]*s[4] +r[5]*s[5] ; // edges 1 updated
                int rr;
                if ( Xt > 3 )
                {
                    rr = 1;
                }
                else
                {
                    rr = 0;
                }
    
                int TFM_OUT;
                TFM_OUT = ( Pt_TFM[i] > rand01[i][1]);
    
                if( U == 6 )
                {
                    Pt_TFM[i] = (1-BETA) * Pt_TFM[i] + rr * BETA;
                    //printf("run here\n");
                }
                //printf("Pt_TFM[1] = %f\n",Pt_TFM[1]);           
                for (int cnt = 0; cnt < VN_DEGREE; cnt++)
                {    
                    if ( s[cnt] == 1)
                    {
                        VN2CN[ColumIndex[cnt][i]][i] = r[cnt];
                    }
                    else
                    {
                        VN2CN[ColumIndex[cnt][i]][i] = TFM_OUT; 
                    }
                }
                int insum = in1 + in2 + in3 +in4 + in5 + in6; // decision
                codeword[i] = (( ChannelBits[i] + insum ) > 3 );
    
            }
        } // end group loop

        int check[CN] = {0};
        int sum = 0;

        if( iter == It_Max )
        {
            break;
        }
        else
        {
            for ( int kcn = 0; kcn < CN; kcn++)
            {
                for ( int kvn = 0; kvn < VN; kvn++)
                {
                    check[kcn] = check[kcn] + Hc[kcn][kvn] * (int)codeword[kvn];
                }
            }
            for(int kcn = 0; kcn < CN ; kcn++)
            {
                check[kcn] = ( check[kcn] % 2 );
                sum = sum + check[kcn];
            }
            if ( sum == 0 )
            {
                break;
            }
        }

    } // end  iteration loop

}  // end function

void d3_subunit_mex(int ChannelBit, int in1,int in2, int *IM_REG, int rand_im , int * b )
{
    if( (ChannelBit == in1) && ( in1 == in2 ) )
    {
        *b = in1;
        *(IM_REG + 1) = *IM_REG;
        *IM_REG = in1;
    }
    else
    {
        *b = *( IM_REG + rand_im );
    }
}

//d2_subunit_mex(b[0],b[6],&r[0],&s[0]);
void d2_subunit_mex(int b1, int b2, int *r, int *s)
{
    if( b1 == b2)
    {
        *r = b1;
        *s = 1;
    }
    else
    {
        *r = 0;
        *s = 0;
    }
}

// mexfunction interface
void mexFunction( int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    double * LLR;
    double N0;
    double * H;
    double * codeword;

    LLR = mxGetPr(prhs[0]);
    N0 = mxGetScalar(prhs[1]);
    H = mxGetPr(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(VN,1,mxREAL);
    codeword = mxGetPr(plhs[0]);

    Shuffled_MTFM_mex(LLR,N0,H,codeword);
}