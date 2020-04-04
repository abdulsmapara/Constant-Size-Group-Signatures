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
#include <assert.h>

typedef struct setup_output {
	element_t mk0;
	mpz_t mk1;
	mpz_t tk;
	mpz_t n;
	element_t G;
	pairing_t pairing;
	// e can be computed by using element_pairing function
	element_t gen_subgroup_q;
	element_t* generators;
	element_t B_Omega;
	element_t A_val;
	mpz_t alpha;
	element_t identity;	// identity of G
	pbc_param_t p;
	mpz_t omega;
	mpz_t pval, qval;
}setup_result;
setup_result ret_setup;


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
		mpz_t lower_limit;
		mpz_init(lower_limit);
		mpz_ui_pow_ui(lower_limit, 2, mpz_get_ui(n)-1);
		while (1){
			mpz_urandomb(random, state ,mpz_get_ui(n));

			if (mpz_cmp(random, lower_limit) > 0 && mpz_probab_prime_p(random,mpz_get_ui(n))) {
				mpz_set(result,random);
				break;
			}
		}			
	}

}


// setup operation Core Construction 5.1
void setup(setup_result* retval, mpz_t security_parameter) {

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

	mpz_add_ui(p_bits, k, 1);
	mpz_add_ui(q_bits, k, 2);

	mpz_t p,q;

	mpz_init(p);
	mpz_init(q);
	mpz_init(retval->pval);
	mpz_init(retval->qval);
	
	random_prime_bits(p, p_bits);
	do{
		random_prime_bits(q, q_bits);
	}while(mpz_cmp(p,q) == 0);

	// gmp_printf("Chosen p,q = %Zd\n%Zd\n", p,q);

	// calculate n

	mpz_t n;

	mpz_init(n);
	
	mpz_mul(n, p, q);

	mpz_init(retval->n);
	mpz_set(retval->n, n);

	gmp_printf("N = %Zd\n", n);
	gmp_printf("P = %Zd\n", p);
	gmp_printf("Q = %Zd\n", q);

	mpz_set(retval->pval, p);
	mpz_set(retval->qval, q);
	
	// Cyclic Bilinear group G of order n with symmetric pairing
	
	pairing_t pairing;
	
	pbc_param_t par;

	pbc_param_init_a1_gen(par, n);
	pbc_param_init_a1_gen(retval->p, n);
	pairing_init_pbc_param(pairing, par);
	element_t g1, gt1, identity, h, add, temp, mk0;
	pairing_init_pbc_param(retval->pairing, par);
	element_init_G1(g1, pairing);
	element_init_G2(g1, pairing);
	element_init_GT(gt1, pairing);

	element_init_G1(add, pairing);
	element_init_G1(temp, pairing);
	element_init_G1(h, pairing);
	element_init_G1(identity, pairing);
	element_init_G1(retval->G, pairing);
	element_init_GT(retval->A_val, pairing);
	element_init_G2(mk0, pairing);
	element_set0(identity);
	element_init_G1(retval->identity, pairing);
	element_set0(retval->identity);

	element_set(retval->G, g1);
	
	do {
		element_random(g1);	
		element_pow_mpz(h, g1, p);
	} while(element_cmp(h, identity) == 0);
	
	element_init_G1(retval->gen_subgroup_q, pairing);
	element_set(retval->gen_subgroup_q, h);

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
					index++;
				}
			}
		}
	} while(mpz_cmp(gen, required));

	retval->generators = generators;

	// Pick random exponents alpha, omega from Zn
	mpz_t alpha, omega;
	mpz_init(alpha);
	mpz_init(omega);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, (rand()+1)*(rand()+1));
	mpz_urandomm(alpha, state, n);
	mpz_urandomm(omega, state, n);
	mpz_init(retval->omega);
	mpz_set(retval->omega, omega);
	gmp_printf("Alpha = %Zd\nOmega = %Zd\n", alpha, omega);
	mpz_init(retval->alpha);
	mpz_set(retval->alpha, alpha);

	element_pow_mpz(g1, g1, omega);
	element_init_G1(retval->B_Omega, pairing);
	element_set(retval->B_Omega, g1);

	element_pow_mpz(mk0, generators[0], alpha);
	element_init_G1(retval->mk0, pairing);
	element_set(retval->mk0, mk0);

	mpz_init(retval->mk1);
	mpz_set(retval->mk1, omega);

	mpz_init(retval->tk);
	mpz_set(retval->tk, q);
}

/* 
	Input: PP, MK, ID
	0 <= ID < 2^k < p
*/
static int num_user = 0;
int add_user_val;
void enroll(mpz_t userID, element_t k1, element_t k2, element_t k3) {
	// Choose unique SID
	/*
		Using the following property -
		Zn = {0, 1, 2, 3, ...., p, p+1, p+2,...., q, q+1, q+2,.. n}
		Zn* = {1, 2, ..., min(p,q) > 2^k, ..... }
		       <------Say F-------> 
		       Remaining = (p-1)*(q-1) - min(p,q)
		       say q is min of p,q
		       therefore, rem = pq + 1 - p - q - q
						= pq + 1 - p - 2*q
						  if q > 4 (can safely assume)
						  pq > 4p
						  and q < p
						  therefore, rem > 4p - p - 2p = p i.e there are sufficient values left
	*/
	mpz_t gcd, sid;
	mpz_init(gcd);
	mpz_init(sid);
	mpz_gcd(gcd, userID, ret_setup.n);
	assert(mpz_cmp_ui(gcd,1) == 0);
	mpz_t val;
	mpz_init(val);
	mpz_set(val, ret_setup.pval);
	if (mpz_cmp(ret_setup.pval,ret_setup.qval) > 0) {
		mpz_set(val, ret_setup.qval);
	}
	
	while (1){
		mpz_add(sid, userID, ret_setup.omega);
		mpz_mod(sid, sid, ret_setup.n);
		mpz_gcd(gcd, sid, ret_setup.n);
		if (mpz_cmp_ui(gcd, 1) == 0) {
			break;
		}
		num_user+=1;
		mpz_add_ui(userID, val, num_user);
	}
	mpz_t inverse;
	mpz_init(inverse);
	assert(mpz_invert(inverse, sid, ret_setup.n));
	gmp_printf("sid (not to be disclosed to user): %Zd\n", sid);
	// k1, k2, k3 = ((g^alpha)^(1/w+sid), g^sid, u^sid)

	element_init_G1(k1, ret_setup.pairing);
	element_init_G1(k2, ret_setup.pairing);
	element_init_G1(k3, ret_setup.pairing);
	element_pow_mpz(k1, ret_setup.generators[0], ret_setup.alpha);
	element_pow_mpz(k1, k1, inverse);

	element_pow_mpz(k2, ret_setup.generators[0], sid);
	element_pow_mpz(k3, ret_setup.generators[1], sid);

}
int main (int argc, char **argv) {
	srand(time(NULL));
	mpz_t security_parameter;
	mpz_init(security_parameter);
	mpz_set_ui(security_parameter, 256);
	

	/*
		SETUP starts
	*/
	setup(&ret_setup, security_parameter);
	printf("PUBLIC INFORMATION:\n");
	gmp_printf("N: %Zd\n", ret_setup.n);
	printf("G: ");
	element_t grp, add;
	element_init_G1(grp, ret_setup.pairing);
	element_init_G1(add, ret_setup.pairing);
	element_set(grp, ret_setup.generators[0]);
	element_set(add, grp);
	/*
		For printing group, uncomment following lines
		Try for only small values of security_parameter
	*/
	// do {
	// 	element_printf("%B, ",grp);
	// 	element_add(grp, grp, add);
	// }while(element_cmp(grp, ret_setup.identity));
	// printf("\n");
	/*
		End for printing security_parameter
	*/
	printf("FOR e and GT, the pairing parameters are as follows:\n");
	pbc_param_out_str(stdout, ret_setup.p);

	// Generators
	// printf("m+3 GENERATORS of GROUP G: ");
	// for (unsigned long int i = 0; i < mpz_get_ui(security_parameter)+3; i++) {
	// 	element_printf("%B,",ret_setup.generators[i]);
	// }
	// printf("\n");

	// Value of A
	element_t gt1;
	element_init_GT(gt1, ret_setup.pairing);
	element_pairing(gt1, ret_setup.generators[0], ret_setup.generators[0]);
	element_pow_mpz(gt1, gt1, ret_setup.alpha);
	printf("A: ");
	element_printf("%B\n", gt1);

	printf("OMEGA: ");
	element_printf("%B\n", ret_setup.B_Omega);
	printf("MASTER KEY: ");
	element_printf("%B, ",ret_setup.mk0);
	gmp_printf("%Zd\n", ret_setup.mk1);
	printf("TRACING KEY: ");
	element_printf("%B\n", ret_setup.gen_subgroup_q);
	
	/*
		SETUP ends
	*/
	add_user_val = rand()%6+1;
	

	/*
		ENROLL starts
	*/
	//	Generate some user ID
	mpz_t random, userID;
	mpz_init(random);
	mpz_init(userID);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, (rand()+1)*(rand()+1));

	mpz_t count, limit;
	mpz_init(count);
	mpz_init(limit);
	mpz_set_ui(limit, 1);
	while(mpz_cmp(count, security_parameter)) {
		mpz_mul_ui(limit, limit, 2);
		mpz_add_ui(count, count, 1);
	}
	do {
		mpz_urandomm(userID, state, limit);	
	} while(mpz_cmp_ui(userID, 0) == 0);
	element_t k1, k2, k3;
	enroll(userID, k1, k2, k3);
	element_printf("K1: %B\nK2: %B\nK3: %B\n", k1, k2, k3);

}