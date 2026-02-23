/*
 * Copyright (c) 2010 Adina Ciomaga <ciomaga@cmla.ens-cachan.fr>
 *                    Marco Mondelli <m.mondelli@sssup.it>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @mainpage Finite Difference Scheme for the Mean Curvature Motion
 *
 * README.txt:

 * @verbinclude README.txt

 */

/**
 * @file FDS_MCM.c
 *
 * @brief Finite Difference Scheme for the Mean Curvature Motion
 *
 * @version 1.3
 * @author Adina Ciomaga;       <adina.ciomaga@cmla.ens-cachan.fr> ;
 * @author Marco Mondelli;      <m.mondelli@sssup.it>
 *
 * requirements:
 * - ANSI C compiler
 * - libtiff
 *
 * compilation:
 * - type "gcc FDS_MCM.c io_tiff_all.c -ltiff -o mcm.out" to obtain the
 * executable file
 * - type "mcm.out" followed by the required arguments to run the executable
 * file
 *
 * usage: this programs requires 3 arguments written in the following order,
 * - argv[1] contains the normalized scale you want to reach
 * - argv[2] contains the name of the input file
 * - argv[3] contains the name of the output file
 *
 * The normalized scale represents the scale @f$ R @f$ at which a circle with
 * radius @f$ R @f$ disappears.
 * The number of iterations needed is linked to normalized scale and time step
 * by the following formula:
 * @f[
 * n_{iter}= \frac{R^2}{2 \cdot dt}
 * @f]
 * where @f$ R @f$ is the normalized scale, @f$ dt @f$ the time step and
 * @f$ n_{iter} @f$ the number of iterations.
 *
 */


/* This line was added by Jose-Luis Lisani to prevent compilation warning */
#include <string.h>
/*                                                                        */

/*
 * Pre processing
 */

#include <stdlib.h>
/* 
#include <tiffio.h>
#include "io_tiff_routine.h"
#include <time.h>
*/
#include <math.h>

/*
 * Print error messages
 */


/** print a message and abort the execution */
#define FATAL(MSG)\
        do {\
          fprintf(stderr, MSG "\n");\
        abort();\
        } while (0);

/** print an info message */
#define INFO(MSG)\
    do {\
        fprintf(stderr, MSG "\n");\
    } while (0);

/** print a debug message with detailed information */
#define DEBUG(MSG)\
    do {\
        if (debug_flag)\
            fprintf(stderr, __FILE__ ":%03i " MSG "\n", __LINE__);\
    } while (0);


/*
 * Compute one iteration of Heat Equation
 */

/**
 * @brief Computes one iteration of Heat Equation.
 *
 * This function applies the following Finite Difference Scheme for Heat 
 * Equation
 * @f[
 * u_{n+1} = u_n + dt \cdot \Delta u_n
 * @f]
 * where the discrete laplacian @f$ \Delta u_n @f$ is given by
 * @f[
 * \Delta u_n =  u_n(i+1,j) + u_n(i-1,j)
 *     + u_n(i,j+1) + u_n(i,j-1)
 *     - 4u_n(i,j)
 * @f]
 *
 *
 * @param[out] ptr_out  pointer to the output array
 * @param[in]  ptr_in   pointer to the input array
 * @param[in]  n_col    number of columns
 * @param      i        current row
 * @param      i_prev   previous row
 * @param      i_next   next row
 * @param      j        current column
 * @param      j_prev   previous column
 * @param      j_next   next column
 * @param[in]  dt       time step
 *
 *
 */

static void heat (float *ptr_out, float *ptr_in, size_t n_col, int i, 
                  int i_prev, int i_next, int j, int j_prev, int j_next, 
                  float dt)
{
    float laplacian;
    laplacian = -4.0 * (*(ptr_in + i*n_col + j)) +
                *(ptr_in + i_prev*n_col + j) + *(ptr_in + i*n_col + j_prev) +
                *(ptr_in + i*n_col + j_next) + *(ptr_in + i_next*n_col + j);

    *(ptr_out+i*n_col + j)= *(ptr_in+i*n_col+j) + 0.5 * dt * laplacian;
}

/*
 * Compute one iteration of the Mean Curvature Motion
 */

/**
 * @brief Computes one iteration of the Mean Curvature Motion
 *
 * This function applies a Finite Difference Scheme for the Mean Curvature 
 * Motion
 * @f[
 * u_t=|Du|curv(u)
 * @f]
 * When @f$ |Du| \neq 0 @f$, we denote by @f$ \xi @f$ the direction orthogonal 
 * to @f$ Du @f$ and by @f$ \theta @f$ the angle between the @f$ x @f$ positive
 * semiaxis and the gradient.
 * Since
 * @f[
 * u_{\xi\xi}=D^2u(\xi,\xi)=D^2u\Bigl(\frac{Du^{\perp}}{|Du|},\frac{Du^{\perp}}
 * {|Du|}\Bigr)= \frac{1}{|Du|^2}D^2u(Du^{\perp},Du^{\perp})= |Du|curv(u)
 * @f]
 * the continuous equation can be translated in the followig discrete version:
 * @f[
 * u_{n+1}(i,j)=u_n(i,j)+\Delta t (u_{\xi\xi})_n(i,j)
 * @f]
 * For the numerical evaluation of @f$ (u_{\xi \xi})_n(i,j) @f$ we adopt a 
 * linear finite difference scheme based on a @f$ 3 \times 3 @f$ stencil.
 * @f[
 * (u_{\xi \xi})_n(i,j)=
 *      \frac{1}{\Delta x^2}( \hspace{0.3em} -4\lambda_0 \cdot u_n(i,j)
 * \hspace{0.3em} + \hspace{0.3em} \lambda_1 \cdot (u_n(i,j+1) + u_n(i,j-1))
 * \hspace{0.3em} + \hspace{0.3em} \lambda_2 \cdot (u_n(i+1,j) + u_n(i-1,j))
 * \hspace{0.3em} + \hspace{0.3em} \lambda_3 \cdot (u_n(i-1,j-1) +u_n(i+1,j+1))
 * \hspace{0.3em} + \hspace{0.3em} \lambda_4 \cdot (u_n(i-1,j+1) +u_n(i+1,j-1))
 * \hspace{0.3em})
 * @f]
 * The later can be rewritten in a more synthetic form, using a discrete 
 * convolution with a variable kernel @f$ A @f$
 * @f[
 * u_{\xi\xi_n}=\frac{1}{\Delta x^2}\cdot (A \star u_n) \quad \mbox{where} 
 * \quad
 * A=\Biggl(\begin{array}{ccc}
 * \lambda_3(\theta) & \lambda_2(\theta)   & \lambda_4(\theta) \\
 * \lambda_1(\theta) & -4\lambda_0(\theta) & \lambda_1(\theta) \\
 * \lambda_4(\theta) & \lambda_2(\theta)   & \lambda_3(\theta)
 * \end{array}
 * \Biggr)
 * @f]
 * We talk about @f$ \emph{variable kernel} @f$ since the @f$ \lambda_i @f$
 * depend on the direction of @f$ |Du| @f$ and can be calculated either from
 * @f$ \cos\theta @f$ and @f$ \sin\theta\ @f$ or (as we do in the code to 
 * optimize the number of variables allocated) directly from the FDS of the 
 * partial derivatives @f$ u_x @f$ and @f$ u_y @f$.
 * @f[
 * \begin{array}{l}
 * u_x= \displaystyle \frac{2\cdot (u(i,j+1) - u(i,j-1)) +
 *              u(i-1,j+1) - u(i-1,j-1) +
 *              u(i+1,j+1) - u(i+1,j-1)} {8\Delta x}=
 *      \displaystyle \frac{s_x}{8\Delta x} \\
 *                                          \\
 * u_y= \displaystyle \frac{2\cdot (u(i+1,j) - u(i-1,j)) +
 *              u(i+1,j+1) - u(i-1,j+1) +
 *              u(i+1,j-1) - u(i-1,j-1)} {8\Delta x}=
 *      \displaystyle \frac{s_y}{8\Delta x} \\
 * \end{array}
 * @f]
 * where we have denoted respectively by @f$ s_x @f$ and @f$ s_y @f$ the
 * algebric sums at the numerator of @f$ u_x @f$ and @f$ u_y @f$.
 *
 * With the choice @f$\lambda_0 = 0.5 - (\cos\theta \sin\theta)^2 @f$ the
 * following formulas hold for the @f$ \lambda_i@f$:
 * @f[
 * \begin{array}{l}
 * \lambda_0= \displaystyle \frac{1}{2} - \frac{s_x^2 \cdot s_y^2}
 * {(s_x^2+s_y^2)^2} \\
 * \\
 * \lambda_1= 2\lambda_0 - \displaystyle \frac{s_x^2}{s_x^2+s_y^2} \\
 * \\
 * \lambda_2= 2\lambda_0 - \displaystyle \frac{s_y^2}{s_x^2+s_y^2} \\
 * \\
 * \lambda_3= -\lambda_0 +\displaystyle \frac{1}{2} \cdot 
 *                          \Bigl(1-\frac{s_x \cdot s_y} {s_x^2+s_y^2}\Bigr) \\
 * \\
 * \lambda_4= -\lambda_0 +\displaystyle \frac{1}{2} \cdot 
 *                          \Bigl(1+\frac{s_x \cdot s_y} {s_x^2+s_y^2}\Bigr)
 * \end{array}
 * @f]
 *
 * Remember that the array representing the image is examinated as a matrix, 
 * i.e.
 * @f$ u(i,j)=pixel(i,j) @f$: @f$ i @f$ is the current row number and
 * @f$ j @f$ is the current column number. The infinitezimal increment 
 * @f$\Delta x @f$ is the pixel
 * width and is set to @f$ 1 @f$.
 *
 * At the borders, we symmetrize the image with respect to the axes 
 * @f$ x=0 @f$ and @f$
 * y=0 @f$, i.e.
 * @f[

 * \begin{array}{l}
 * u(-1,j)    = u(0,j)       \\
 * \\
 * u(n_{row},j) = u(n_{row}-1,j) \\
 * \\
 * u(i,-1)    = u(i,0)       \\
 * \\
 * u(i,n_{col}) = u(i,n_{col}-1)
 * \end{array}
 * @f]
 *
 * All these speculations are based on the hypothesis that @f$ |Du|\neq 0 @f$.
 * In fact @f$ curv(u)=\displaystyle\frac{1}{|Du|^3} D^2u(Du^{\perp},
 * Du^{\perp})@f$, which means that, if @f$ Du = 0@f$, the Mean
 * Curvature Motion is undefined and has no sense talking about @f$ \xi @f$ or
 * @f$ \theta @f$. Numerically, even if the gradient is @f$\neq 0 @f$, but
 * small in modulus, its direction becomes substantially random, because it 
 * will be driven by rounding and approximation errors.
 *
 * These considerations bring us to replace @f$u_{\xi\xi} @f$ by half the
 * laplacian, if @f$ |Du| @f$ is below a given threshold (experimentally set to
 * @f$ 4 @f$).
 *
 * @param[out] ptr_out  pointer to the output array
 * @param[in]  ptr_in   pointer to the input array
 * @param[in]  n_row    number of rows
 * @param[in]  n_col    number of columns
 * @param[in]  dt               time step
 * @param[in]  t_g      threshold for the gradient
 *
 */

static void mcm(float *ptr_out, float *ptr_in, size_t n_row, size_t n_col,
                float dt, float t_g)
{
    unsigned int  i, j, i_next, i_prev, j_next, j_prev;
    float         grad, s_x, s_y, lambda0, lambda1, lambda2, lambda3, lambda4;

    /* iterate on j and i, following the array order */
    for (i=0; i<n_row; i++ )  {

        /* symmetry for borders */
        i_next = (i<n_row-1?i+1:i);
        i_prev = (i>0?i-1:i);

        for (j=0; j<n_col; j++)  {

            /* symmetry for borders */
            j_next = (j<n_col-1?j+1:j);
            j_prev = (j>0?j-1:j);

            /* computation of s_x, s_y and gradient */
            s_x= 2 * ( *(ptr_in+i*n_col+j_next) - *(ptr_in+i*n_col+j_prev) ) +
                *(ptr_in+i_prev*n_col+j_next) - *(ptr_in+i_prev*n_col+j_prev) +
                 *(ptr_in+i_next*n_col+j_next) - *(ptr_in+i_next*n_col+j_prev);
            s_y= 2 * ( *(ptr_in+i_next*n_col+j) - *(ptr_in+i_prev*n_col+j) ) +
                *(ptr_in+i_next*n_col+j_next) - *(ptr_in+i_prev*n_col+j_next) +
                *(ptr_in+i_next*n_col+j_prev) - *(ptr_in+i_prev*n_col+j_prev);
            grad = 0.125 * sqrt( (s_x*s_x) + (s_y*s_y) );

            if (grad>t_g) {
                /* MCM diffusion */
                grad*=8;   /* i.e. grad=\sqrt{s_x^2 + s_y^2} */
                lambda0 = 0.5 - ( (s_x*s_x*s_y*s_y) / ( pow(grad,4) ) );
                lambda1= 2 * lambda0 - ( (s_x*s_x) / (grad*grad) );
                lambda2= 2 * lambda0 - ( (s_y*s_y) / (grad*grad) );
                lambda3= -lambda0 + 0.5 * ( 1 - ( (s_x*s_y)/(grad*grad) ) );
                lambda4= -lambda0 + 0.5 * ( 1 + ( (s_x*s_y)/(grad*grad) ) );

                *(ptr_out+i*n_col + j)= *(ptr_in+i*n_col+j) + dt * (
                                         -4 * lambda0 * (*(ptr_in+i*n_col+j)) +
                                         lambda1 * (*(ptr_in+i*n_col+j_next)  +
                                                    *(ptr_in+i*n_col+j_prev)) +
                                         lambda2 * (*(ptr_in+i_next*n_col+j)  +
                                                    *(ptr_in+i_prev*n_col+j)) +
                                         lambda3 * (*(ptr_in+i_prev*n_col
                                                                     +j_prev) +
                                                    *(ptr_in+i_next*n_col
                                                                     +j_next))+
                                         lambda4 * (*(ptr_in+i_prev*n_col
                                                                     +j_next) +
                                                    *(ptr_in+i_next*n_col
                                                                    +j_prev)));
            }

            else
                /* Heat Equation */
                heat (ptr_out, ptr_in, n_col, i, i_prev, i_next, 
                      j, j_prev, j_next, dt);
        }
    }
}

/* José Luis */
/* Output: in */
void mcm_main(float *in, float *aux, int w, int h, float scale, float gradth)
{
  float *out;
  int n, m;

  int          n_iter=0;                    /* number of iterations */
  float        R;                           /* normalization scale */
  float        t_g=gradth;                       /* gradient threshold */
  /*float        t_g=4;*/                       /* gradient threshold */
  /*float        dt=0.1;*/
  float        dt=0.5;
  float        n_iter_f;                    /* final number of iterations */

  R=scale;
  /* number of iterations needed */
  n_iter_f= (R * R) / (2 * dt);

  /* round n_iter_f (must be integer) */
  n_iter= (int) (n_iter_f+0.5);

  if (aux) out=aux;
  else out=(float *) malloc(w*h*sizeof(float));

  for (m=0; m<n_iter; m++)   {
    /* apply mean curvature motion */
    mcm(out, in, h, w, dt, t_g);
    memcpy(in, out, w*h*sizeof(float));
  }
  for (n=0; n < w*h; n++) {
     /* saturate values to [0,255] */
     if (in[n] < 0) in[n]=0;
     if (in[n] > 255) in[n]=255;
  }

  if (!aux) free(out);
}

/*
 * Main function section
 */

/**
 * @brief Main function call
 *
 * The program reads an image.tiff given as input,
 * applies a Finite Difference Scheme for the Mean Curvature Motion at 
 * renormalized scale @f$ R @f$,
 * writes the result as a new image.tiff and returns it as output.
 *
 */

/* Removed by Jose Luis */
/* int main(int argc, char *argv[]) */

