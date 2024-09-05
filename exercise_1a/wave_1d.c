#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>


// Option to change numerical precision.
typedef int64_t int_t;
typedef double real_t;

// Simulation parameters: size, step count, and how often to save the state.
const int_t
    N = 1024,
    max_iteration = 4000,
    snapshot_freq = 10;

// Wave equation parameters, time step is derived from the space step.
const real_t
    c  = 1.0,
    dx  = 1.0;
real_t
    dt;

// Buffers for three time steps, indexed with 2 ghost points for the boundary.
real_t
    *buffers[3] = { NULL, NULL, NULL };


#define U_prv(i) buffers[0][(i)+1]
#define U(i)     buffers[1][(i)+1]
#define U_nxt(i) buffers[2][(i)+1]


// Save the present time step in a numbered file under 'data/'.
void domain_save ( int_t step )
{
    char filename[256];
    sprintf ( filename, "data/%.5ld.dat", step );
    FILE *out = fopen ( filename, "wb" );
    fwrite ( &U(0), sizeof(real_t), N, out );
    fclose ( out );
}


// TASK: T1
// Set up our three buffers, fill two with an initial cosine wave,
// and set the time step.
void domain_initialize ( void )
{
// BEGIN: T1
    // 1. Allocate memory for three sets of spatial domains
    for (int i = 0; i < 3; i++) {
        buffers[i] = (real_t*) malloc((N + 2) * sizeof(real_t));
        if (!buffers[i]) {
            printf("Failed to allocate memory for buffer %d", i);
            exit(EXIT_FAILURE);
        }
    }

    // 2. Apply the initial condition from Equation 4.
    for (int i = 0; i < N; i++) {
        real_t x_value = (real_t)i / (N - 1);
        U(i) = cos(M_PI * x_value);
        U_prv(i) = cos(M_PI * x_value);
    }

    // 3. Set the time step, dt, based on Equation 5.
    dt = dx / c;
    ;
// END: T1
}


// TASK T2:
// Return the memory to the OS.
// BEGIN: T2
void domain_finalize ( void )
{
    for (int i = 0; i < 3; i++) {
        free(buffers[i]);
    }
    ;
}
// END: T2


// TASK: T3
// Rotate the time step buffers.
// BEGIN: T3
void domain_shift_window( void )
{
    real_t* outdated_buffer = buffers[0];
    buffers[0] = buffers[1];
    buffers[1] = buffers[2];
    free(outdated_buffer);

    buffers[2] = (real_t*) malloc((N + 2) * sizeof(real_t));
    if (!buffers[2]) {
        printf("Failed to allocate memory for buffer 2");
        exit(EXIT_FAILURE);
    }

}
    ;
// END: T3


// TASK: T4
// Derive step t+1 from steps t and t-1.
// BEGIN: T4
    ;
void derive_time_step( void ) 
{
    real_t stability_term = ((pow(dt, 2) * pow(c, 2)) / pow(dx, 2));
    for (int i = 0; i < N; i++) {
        real_t x_value = 
            -U_prv(i) 
            + (2 * U(i)) 
            + stability_term * (U(i - 1) + U(i + 1) - (2 * U(i)));
        U_nxt(i) = x_value;
    }
}
// END: T4


// TASK: T5
// Neumann (reflective) boundary condition.
// BEGIN: T5
void set_ghost_values( void )
{
    U(-1) = U(0);
    U(N) = U(N - 1);
}
    ;
// END: T5


// TASK: T6
// Main time integration.
void simulate( void )
{
// BEGIN: T6
    int_t iteration=0;
    for (int i = 0; i < max_iteration; i++) {
        set_ghost_values();
        derive_time_step();
        domain_shift_window();
        if (iteration % snapshot_freq == 0) {
            domain_save ( iteration / snapshot_freq );
        }
        ++iteration;
    }
// END: T6
}


int main ( void )
{
    domain_initialize();

    simulate();

    domain_finalize();
    exit ( EXIT_SUCCESS );
}
