
/* Copyright 2016 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Scott Munro, Matthew Shin
                  
                                                                     */

#include "syrk2_blk_var1.h"
#include "../UnblockedVar1/syrk2_un_var1.h"
#include "FLAME.h"

#define UPLO FLA_UPPER_TRIANGULAR
#define TRANS FLA_NO_TRANSPOSE

int syrk2_blk_var1( FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  FLA_Obj CTL,   CTR,      C00, C01, C02, 
          CBL,   CBR,      C10, C11, C12,
                           C20, C21, C22;

  int b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_TOP );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    b = min( FLA_Obj_length( AB ), nb_alg );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                        /* ** */            /* ** */
                                              &B1, 
                           BB,                &B2,        b, FLA_BOTTOM );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00, /**/ &C01, &C02,
                        /* ************* */   /* ******************** */
                                                &C10, /**/ &C11, &C12,
                           CBL, /**/ CBR,       &C20, /**/ &C21, &C22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    // unblocked version
    /*
    //c01 = A0*transpose(b1t) + B0*transpose(a1t) + c01;
    FLA_Gemv(TRANS, FLA_ONE, A0, b1t, FLA_ONE, c01);
    FLA_Gemv(TRANS, FLA_ONE, B0, a1t, FLA_ONE, c01);

    //gamma11 = a1t*transpose(b1t) +  b1t*transpose(a1t) + gamma11;
    FLA_Dots(FLA_ONE, a1t, b1t, FLA_ONE, gamma11);
    FLA_Dots(FLA_ONE, a1t, b1t, FLA_ONE, gamma11);
    */

    // blocked
    // C01 = A0 * transpose(B1) + B0 * transpose(A1) + C01;
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A0, B1, FLA_ONE, C01);
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, B0, A1, FLA_ONE, C01);

    // C11 = A1 * transpose(B1) + B1 * transpose(A1) + C11;
    syrk2_un_unb_var1( A1, B1, C11 );



    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  A1, 
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                                                  B1, 
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00, C01, /**/ C02,
                                                     C10, C11, /**/ C12,
                            /* ************** */  /* ****************** */
                              &CBL, /**/ &CBR,       C20, C21, /**/ C22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}
