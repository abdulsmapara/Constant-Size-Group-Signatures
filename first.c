/*
	author: @abdulsmapara
	Paper: CSGS
*/
#include <pbc/pbc.h>
#include <pbc/pbc_test.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>


// generates a random prime number of n bits
void random_prime_bits(mpz_t result, mpz_t n) {

	mpz_t random;
	mpz_init(random);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	mpz_init(random);	
	gmp_randseed_ui(state, (rand()+1)*(rand()+1));
	if (mpz_cmp_ui(n,1) <= 0) {
		printf("NO PRIME EXISTS\n");
	} else {
		while (1){
			mpz_urandomb(random, state ,mpz_get_ui(n));
			if (mpz_probab_prime_p(random,mpz_get_ui(n))) {
				mpz_set(result,random);
				break;
			}
		}			
	}

}


// setup operation Core Construction 5.1
void setup(mpz_t security_parameter) {


	// choose k and m = O(security_parameter)
	mpz_t k,m;
	
	mpz_init(k);
	mpz_init(m);

	mpz_set(k, security_parameter);
	mpz_set(m, security_parameter);

	// choose p, q - random prime numbers
	mpz_t p_bits, q_bits;
	
	mpz_init(p_bits);
	mpz_init(q_bits);

	mpz_add_ui(p_bits, k, 0);
	mpz_add_ui(q_bits, k, 0);

	mpz_t p,q;

	mpz_init(p);
	mpz_init(q);

	random_prime_bits(p, p_bits);
	random_prime_bits(q, q_bits);

	// gmp_printf("Chosen p,q = %Zd\n%Zd\n", p,q);

	// calculate n

	mpz_t n;

	mpz_init(n);
	
	mpz_mul(n, p, q);

	gmp_printf("N = %Zd\n", n);
	gmp_printf("P = %Zd\n", p);
	gmp_printf("Q = %Zd\n", q);


	// Cyclic Bilinear group G of order n with symmetric pairing
	
	pairing_t pairing;
	
	pbc_param_t par;

	pbc_param_init_a1_gen(par, n);
	pairing_init_pbc_param(pairing, par);
	element_t g1, gt1, identity, h, add, temp, mk0;

	element_init_G1(g1, pairing);
	element_init_G2(g1, pairing);
	element_init_GT(gt1, pairing);

	element_init_G1(add, pairing);
	element_init_G1(temp, pairing);
	element_init_G1(h, pairing);
	element_init_G1(identity, pairing);
	element_init_G2(mk0, pairing);
	element_set0(identity);

	do {
		element_random(g1);	
		element_pow_mpz(h, g1, p);
	} while(element_cmp(h, identity) == 0);
	
	element_printf("h = %B\n", h);

	/* ******************************
		Verify that h is a generator 
		Run the below code only for small values
	   ******************************
	*/
	// if (element_cmp(h, identity) == 0) {
	// 	printf("h is not a generator of group Gq\n");
	// 	return;
	// }


	// element_set(add, h);
	// mpz_t count;
	// mpz_init(count);
	// mpz_set_ui(count, 1);
	// while(element_cmp(h, identity) != 0) {
	// 	element_add(h, h, add);
	// 	element_printf("h generates %B\n", h);
	// 	mpz_add_ui(count, count, 1);
	// }
	// if (mpz_cmp(count, q) == 0) {
	// 	printf("h is a generator of Gq\n");
	// } else {
	// 	printf("h is not a generator of Gq\n");
	// 	return;
	// }
	
	/*
		*****************************
		Verification Done for H is a generator
		*****************************
	*/
	mpz_t required, gen;
	mpz_init(required);
	mpz_init(gen);
	mpz_add_ui(required, m,3);
	element_t* generators = (element_t*)malloc(sizeof(element_t)*(3+mpz_get_ui(m)));
	
	for (unsigned long int i = 0; i < (3+mpz_get_ui(m)); i++) {
		element_init_G1(generators[i], pairing);
	}
	unsigned long long int index = 0;
	// Choosing m+3 (one generator g and m+2 others) generator for G
	do {
		element_random(g1);
		element_pow_mpz(temp, g1, n);

		if (element_cmp(temp, identity) == 0) {
			element_pow_mpz(temp, g1, p);
			if (element_cmp(temp, identity)) {
				element_pow_mpz(temp, g1, q);
				if (element_cmp(temp, identity)) {
					mpz_add_ui(gen, gen, 1);
					element_set(generators[index], g1);
					element_printf("Generator = %B\n", generators[index]);
					index++;
				}
			}
		}
	} while(mpz_cmp(gen, required));


	// Pick random exponents alpha, omega from Zn
	mpz_t alpha, omega;
	mpz_init(alpha);
	mpz_init(omega);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, (rand()+1)*(rand()+1));
	mpz_urandomm(alpha, state, n);
	mpz_urandomm(omega, state, n);
	gmp_printf("Alpha = %Zd\nOmega = %Zd\n", alpha, omega);

	element_pairing(gt1, g1, g1);
	element_pow_mpz(gt1, gt1, alpha);
	element_printf("A = %B\n", gt1);

	element_pow_mpz(g1, g1, omega);
	element_printf("B_Omega = %B\n", g1);

	printf("PUBLIC INFORMATION\n");

	printf("MASTER KEY: ");

	element_pow_mpz(mk0, generators[0], alpha);
	element_printf("(%B,", mk0);
	gmp_printf("%Zd)\n", omega);

	printf("TRACING KEY: ");
	gmp_printf("%Zd\n", q);
}


int main (int argc, char **argv) {
	srand(time(NULL));
	
	mpz_t security_parameter;
	mpz_init(security_parameter);
	mpz_set_ui(security_parameter, 4);
	setup(security_parameter);

}