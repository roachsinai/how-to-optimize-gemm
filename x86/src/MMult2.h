#include <stdio.h>

#define A(i,j) a[ (i)*lda + (j) ]
#define B(i,j) b[ (i)*ldb + (j) ]
#define C(i,j) c[ (i)*ldc + (j) ]

/* Routine for computing C = A * B + C */

/* Create macro to let X( i ) equal the ith element of x */

#define Y(i) y[ (i)*incx ]

void AddDot( int k, float *x, int incx,  float *y, float *gamma )
{
  /* compute gamma := x' * y + gamma with vectors x and y of length n.
     Here x starts at location x with increment (stride) incx and y starts at location y and has (implicit) stride of 1.
  */

  int p;

  for ( p=0; p<k; p++ ){
    *gamma += x[p] * Y(p);
  }
}

// 相对于 MY_MMult1，虽然 MY_MMult2 的想法是一次计算结果矩阵的 4 个元素
// 但是由于 AddDot 是完整的计算了两个向量的结果之后才结束
// 导致 MY_MMult2 中矩阵 B 的空间局部性其实和 MY_MMult1 是一样的
void MY_MMult2( int m, int n, int k, float *a, int lda,
                                    float *b, int ldb,
                                    float *c, int ldc )
{
  int i, j;

  for ( j=0; j<n; j+=4 ){        /* Loop over the columns of C, unrolled by 4 */
    for ( i=0; i<m; i+=1 ){        /* Loop over the rows of C */
      /* Update the C( i,j ) with the inner product of the ith row of A
	 and the jth column of B */

      AddDot( k, &A( i,0 ), lda, &B( 0,j ), &C( i,j ) );

      /* Update the C( i,j+1 ) with the inner product of the ith row of A
	 and the (j+1)th column of B */

      AddDot( k, &A( i,0 ), lda, &B( 0,j+1 ), &C( i,j+1 ) );

      /* Update the C( i,j+2 ) with the inner product of the ith row of A
	 and the (j+2)th column of B */

      AddDot( k, &A( i,0 ), lda, &B( 0,j+2 ), &C( i,j+2 ) );

      /* Update the C( i,j+3 ) with the inner product of the ith row of A
	 and the (j+1)th column of B */

      AddDot( k, &A( i,0 ), lda, &B( 0,j+3 ), &C( i,j+3 ) );
    }
  }
}

