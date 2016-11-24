#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"

#include "cblas.h"
#include "UnblockedVar1/syrk2_un_var1.h"
#include "BlockedVar1/syrk2_blk_var1.h"

#define UPLO FLA_UPPER_TRIANGULAR
#define TRANS FLA_NO_TRANSPOSE

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0
#define BLOCKSIZE 4

int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i, ii, jj, irep,
    nrepeats;

  double
    dtime, dtime_best, 
    diff;

  FLA_Obj
    Aobj, Bobj, Cobj, Cold, Cref;
  
  /* Initialize FLAME.  Notice, we will only use the FLAME routines to time the routines  */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  for ( n=nfirst; n<= nlast; n+=ninc ){

    /* Allocate space for the matrices */
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Aobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cold );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cref );

    /* Generate random matrices A, B, and C */
    FLA_Random_matrix( Aobj );
    FLA_Random_matrix( Bobj );
    FLA_Random_matrix( Cold );



    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time reference implementation (from libflame) */
      FLA_Copy( Cold, Cref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute Cref = A * B' + B * A' + Cref or Cref = A' * B + B' * A + Cref,
	 where C is symmetric stored in the
	 indicated triangular part of array C, by calling FLA_Syr2k.  The
	 result ends up in Cref, which we will consider to be the correct
	 result. */
      FLA_Syr2k( UPLO, TRANS, FLA_ONE, Aobj, Bobj, FLA_ONE, Cref );
      
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n,
	    dtime_best );
    fflush( stdout );

    /* Time your unblocked Variant 1 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( Cold, Cobj );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      syrk2_un_unb_var1( Aobj, Bobj, Cobj );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( Cobj, Cref );
  
    printf( "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );

    /* Time your blocked Variant 1 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( Cold, Cobj );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      syrk2_blk_var1( Aobj, Bobj, Cobj, BLOCKSIZE );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( Cobj, Cref );
  
    printf( "data_blk_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );

    FLA_Obj_free( &Aobj );
    FLA_Obj_free( &Bobj );
    FLA_Obj_free( &Cobj );
    FLA_Obj_free( &Cref );
    FLA_Obj_free( &Cold );

    i++;
  }
  FLA_Finalize( );

  exit( 0 );
}

