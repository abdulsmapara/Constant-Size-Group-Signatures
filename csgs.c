/*
	author: @abdulsmapara
	Course: Number Theory and Modern Cryptography
	Paper Implemented: CSGS
*/
#include <pbc/pbc.h>
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
	// e = element_pairing() function
	element_t gen_subgroup_q;
	element_t* generators;
	element_t B_Omega;
	element_t A_val;
	mpz_t alpha;
	element_t identity;	// identity of G
	pbc_param_t p;
	mpz_t omega;
	mpz_t pval, qval;
	mpz_t mval;
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

	mpz_init(retval->mval);
	mpz_set(retval->mval, m);
	// choose p, q - random prime numbers
	mpz_t p_bits, q_bits;
	
	mpz_init(p_bits);
	mpz_init(q_bits);

	mpz_add_ui(p_bits, k, 2);
	mpz_add_ui(q_bits, k, 3);

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
	#ifdef DEBUG
		element_printf("h = %B\n", h);
	#endif
	element_set(retval->gen_subgroup_q, h);

	/* ******************************
		Verify that h is a generator 
		Uncomment the below code only for small values
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

	// element_printf("h = %B\n", h);
	// while(element_cmp(h, identity) != 0) {
	// 	element_add(h, h, add);
	// 	element_printf("h generates %B\n", h);
	// 	mpz_add_ui(count, count, 1);
	// }
	// if (mpz_cmp(count, q) == 0) {
	// 	printf("----------------------------h is a generator of Gq--------------------------\n");
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
	mpz_t gcd_alpha_p, gcd_alpha_q;
	mpz_init(gcd_alpha_p);
	mpz_init(gcd_alpha_q);
	do{
		mpz_urandomm(alpha, state, n);
		mpz_gcd(gcd_alpha_p, alpha, p);
		mpz_gcd(gcd_alpha_q, alpha, q);
	} while(mpz_cmp_ui(gcd_alpha_p, 1) != 0 || mpz_cmp_ui(gcd_alpha_q, 1) != 0);
	mpz_urandomm(omega, state, n);
	mpz_init(retval->omega);
	mpz_set(retval->omega, omega);
	gmp_printf("Alpha = %Zd\nOmega = %Zd\n", alpha, omega);
	mpz_init(retval->alpha);
	mpz_set(retval->alpha, alpha);

	element_pow_mpz(g1, generators[0], omega);
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
static int num_user;
int add_user_val;
void enroll(mpz_t final_sid, mpz_t userID, element_t k1, element_t k2, element_t k3) {
	// Choose unique SID
	

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
		num_user+=1;
		if (mpz_cmp_ui(gcd, 1) == 0) {
			break;
		}
		
		mpz_add_ui(userID, val, num_user);
	}
	mpz_t inverse;
	mpz_init(inverse);
	assert(mpz_invert(inverse, sid, ret_setup.n));
	gmp_printf("sid (not to be disclosed to user): %Zd\n", userID);
	mpz_init(final_sid);
	mpz_set(final_sid, userID);
	// gmp_printf("DEBUG: inverse of sid+omega: %Zd\n", inverse);
	
	// k1, k2, k3 = ((g^alpha)^(1/w+sid), g^sid, u^sid)

	element_init_G1(k1, ret_setup.pairing);
	element_init_G1(k2, ret_setup.pairing);
	element_init_G1(k3, ret_setup.pairing);
	element_pow_mpz(k1, ret_setup.generators[0], ret_setup.alpha);
	element_pow_mpz(k1, k1, inverse);

	element_pow_mpz(k2, ret_setup.generators[0], userID);
	element_pow_mpz(k3, ret_setup.generators[1], userID);

}
void sign(element_t pi1, element_t pi2,element_t sigma1,element_t sigma2,element_t sigma3,element_t sigma4, element_t k1, element_t k2, element_t k3, char* message) {
	// CHOOSE A RANDOM s in Zn
	mpz_t s, t1, t2, t3, t4, t_temp;
	
	mpz_init(s);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(t3);
	mpz_init(t4);
	mpz_init(t_temp);

	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, (rand()+1)*(rand()+1));
	mpz_urandomm(s, state, ret_setup.n);
	// 01 = k1
	// 02 = k2
	// 03 -
	element_t theta3, val1;
	element_init_G1(theta3, ret_setup.pairing);
	element_init_G1(val1, ret_setup.pairing);
	element_set(theta3, ret_setup.generators[2]);
	element_t temp_pw;
	element_init_G1(temp_pw, ret_setup.pairing);
	mpz_t msg_val;
	mpz_init(msg_val);
	for (int i = 3; i < (3+mpz_get_ui(ret_setup.mval)); i++) {
		if (message[i-3] == '1') {
			mpz_set_ui(msg_val,1);
		} else {
			mpz_set_ui(msg_val, 0);
		}
		// gmp_printf("%Zd\n", msg_val);
		element_pow_mpz(temp_pw, ret_setup.generators[i], msg_val);
		element_mul(theta3, theta3, temp_pw);
	}
	element_set(val1, theta3);
	element_pow_mpz(theta3, theta3, s);
	element_mul(theta3, theta3, k3);
	// element_printf("theta3: %B\n", theta3);
	// 04 -
	element_t theta4;
	element_init_G1(theta4, ret_setup.pairing);
	element_set(theta4, ret_setup.generators[0]);
	element_pow_mpz(theta4, theta4, s);
	element_invert(theta4, theta4);
	// element_printf("theta4: %B\n", theta4);
	// Initial signature formed

	// Verify
	element_t ver1, ver2, ver3;
	element_init_GT(ver1, ret_setup.pairing);
	element_init_GT(ver2, ret_setup.pairing);
	element_init_GT(ver3, ret_setup.pairing);
	element_pairing(ver1, theta3, ret_setup.generators[0]);
	// element_printf("ver1: %B\nval1:%B\n", ver1, val1);
	element_pairing(ver2, theta4, val1);
	// element_printf("ver2: %B\n", ver2);
	element_mul(ver1, ver1, ver2);
	// element_printf("ver1: %B\tver2: %B\n",ver1, ver2);
	element_pairing(ver3, k2, ret_setup.generators[1]);
	// The other verification check is same as 1st verification check of key, so to unnecessary increase time complexity, it is not done again
	// element_printf("ver1-final: %B\nver3-final: %B\n",ver1, ver3);

	if (element_cmp(ver1, ver3) == 0) {
		printf("-----------------Verification of initial signature successful------------------\n");
	} else {
		printf("Verification of initial signature UNSUCCESSFUL\n");
	}
	mpz_urandomm(t1, state, ret_setup.n);
	mpz_urandomm(t2, state, ret_setup.n);
	mpz_urandomm(t3, state, ret_setup.n);
	mpz_urandomm(t4, state, ret_setup.n);

	#ifdef DEBUG
		gmp_printf("t1, t2, t3, t4, s: %Zd, %Zd, %Zd, %Zd, %Zd\n", t1, t2, t3, t4, s);
	#endif

	element_t h, temp;
	element_init_G1(sigma1, ret_setup.pairing);
	element_init_G1(sigma2, ret_setup.pairing);
	element_init_G1(sigma3, ret_setup.pairing);
	element_init_G1(sigma4, ret_setup.pairing);
	element_init_G1(h, ret_setup.pairing);
	element_init_G1(temp, ret_setup.pairing);

	element_set(h, ret_setup.gen_subgroup_q);

	element_set(sigma1, k1);
	element_set(sigma2, k2);
	element_set(sigma3, theta3);
	element_set(sigma4, theta4);

	element_pow_mpz(temp, h, t1);
	element_mul(sigma1, temp, sigma1);

	element_pow_mpz(temp, h, t2);
	element_mul(sigma2, temp, sigma2);
	
	element_pow_mpz(temp, h, t3);
	element_mul(sigma3, temp, sigma3);
	
	element_pow_mpz(temp, h, t4);
	element_mul(sigma4, temp, sigma4);

	// gmp_printf("T = %Zd, %Zd, %Zd, %Zd\n", t1, t2, t3, t4);

	element_t theta1, theta2, u;
	element_init_G1(pi1, ret_setup.pairing);
	element_init_G1(pi2, ret_setup.pairing);
	element_init_G1(theta1, ret_setup.pairing);
	element_init_G1(theta2, ret_setup.pairing);
	element_init_G1(u, ret_setup.pairing);
	element_set(theta1, k1);
	element_set(theta2, k2);
	element_set(u, ret_setup.generators[1]);
	mpz_mul(t_temp, t1, t2);
	element_pow_mpz(h, h, t_temp);
	element_pow_mpz(theta1, theta1,t2);
	element_mul(theta2, theta2, ret_setup.B_Omega);
	element_pow_mpz(theta2, theta2, t1);
	element_mul(pi1, theta1, h);
	element_mul(pi1, pi1, theta2);

	element_pow_mpz(pi2, u, t2);
	element_pow_mpz(u, ret_setup.generators[0], t3);
	element_invert(u, u);
	element_mul(pi2, pi2, u);
	element_pow_mpz(val1, val1, t4);
	element_invert(val1, val1);
	element_mul(pi2, pi2, val1);


	
}
void verify(element_t Aval, element_t sigma1, element_t sigma2, element_t sigma3, element_t sigma4, element_t pi1, element_t pi2, char* message) {

	element_t Acopy, val1, val2, pairing_result1, pairing_result2, T1, T2;
	element_init_G1(val2, ret_setup.pairing);
	element_init_G1(val1, ret_setup.pairing);
	element_init_GT(Acopy, ret_setup.pairing);
	element_init_GT(pairing_result1, ret_setup.pairing);
	element_init_GT(pairing_result2, ret_setup.pairing);
	element_init_GT(T1, ret_setup.pairing);
	element_init_GT(T2, ret_setup.pairing);

	element_set(Acopy, Aval);
	element_invert(Acopy, Acopy);
	element_mul(val1, sigma2, ret_setup.B_Omega);
	element_pairing(T1, sigma1, val1);
	element_mul(T1, T1, Acopy);

	
	
	element_set(val2, ret_setup.generators[2]);
	element_t temp_pw;
	element_init_G1(temp_pw, ret_setup.pairing);
	mpz_t msg_val;
	mpz_init(msg_val);
	for (int i = 3; i < (3+mpz_get_ui(ret_setup.mval)); i++) {
		if (message[i-3] == '1') {
			mpz_set_ui(msg_val,1);
		} else {
			mpz_set_ui(msg_val, 0);
		}
		
		element_pow_mpz(temp_pw, ret_setup.generators[i], msg_val);
		element_mul(val2, val2, temp_pw);
	}
	#ifdef DEBUG
		element_printf("sigma4 %B\n", sigma4);
	#endif
	element_pairing(T2, ret_setup.generators[1],sigma2);
	element_pairing(pairing_result2, ret_setup.generators[0], sigma3);
	
	element_invert(pairing_result2, pairing_result2);
	element_mul(T2, T2, pairing_result2);
	element_pairing(pairing_result2, sigma4, val2);
	element_invert(pairing_result2, pairing_result2);	
	element_mul(T2, T2, pairing_result2);

	element_t T1_verify, T2_verify;
	element_init_GT(T1_verify, ret_setup.pairing);
	element_pairing(T1_verify, ret_setup.gen_subgroup_q, pi1);
	element_init_GT(T2_verify, ret_setup.pairing);
	element_pairing(T2_verify, ret_setup.gen_subgroup_q, pi2);

	#ifdef DEBUG
		element_printf("%B\n%B\n",T2, T2_verify);
	#endif
	
	if (element_cmp(T1, T1_verify) == 0 && element_cmp(T2, T2_verify) == 0) {
		printf("---------------------VERIFICATION VALID-------------------------\n");
	} else {
		printf("VERIFICATION INVALID\n");
	}



}
void trace(element_t sigma2, mpz_t* sids, unsigned long num_of_users) {
	element_t sigma2_copy, val;
	element_init_G1(sigma2_copy, ret_setup.pairing);
	element_init_G1(val, ret_setup.pairing);
	element_set(sigma2_copy, sigma2);
	element_pow_mpz(sigma2_copy, sigma2_copy, ret_setup.tk);
	for (unsigned long i = 0; i < (num_of_users); i++) {
		element_pow_mpz(val, ret_setup.generators[0], sids[i]);
		element_pow_mpz(val, val, ret_setup.tk);
		if (element_cmp(val, sigma2_copy) == 0) {
			gmp_printf("----------------------------SID %Zd traced successfully--------------------------\n", sids[i]);
		} else {
			gmp_printf("----------------------------SID %Zd unsuccessful in passing trace test--------------------------\n", sids[i]);
		}
	}
}
int main (int argc, char **argv) {
// for(int j=0; j<5000;j++){
	num_user = 0;
	srand(time(NULL));
	mpz_t security_parameter;
	mpz_init(security_parameter);

	/*
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		PLEASE READ THIS BEFORE UPDATING security_parameter
		For small values, verification may fail. (as pairing returns 0,0 sometimes for small values of security_parameter like 3,4)
		You may keep values 16, ... 128, 256, 512
		Number of bits in p & q is greater than the security_parameter
	*/

	mpz_set_ui(security_parameter, 512);
	

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
		End for printing group
	*/
	printf("FOR e and GT, the pairing parameters are as follows:\n");
	pbc_param_out_str(stdout, ret_setup.p);

	// Print Generators, uncomment the below code for printing the 1+m+2 generators
	// printf("m+3 GENERATORS of GROUP G: ");
	// for (unsigned long int i = 0; i < mpz_get_ui(security_parameter)+3; i++) {
	// 	element_printf("%B,",ret_setup.generators[i]);
	// }
	// printf("\n");

	// Value of A
	element_t Aval;
	element_init_GT(Aval, ret_setup.pairing);
	element_pairing(Aval, ret_setup.generators[0], ret_setup.generators[0]);
	element_pow_mpz(Aval, Aval, ret_setup.alpha);
	printf("A: ");

	element_printf("%B\n", Aval);

	printf("OMEGA: ");
	element_printf("%B\n", ret_setup.B_Omega);
	printf("MASTER KEY: ");
	element_printf("%B, ",ret_setup.mk0);
	gmp_printf("%Zd\n", ret_setup.mk1);
	printf("TRACING KEY: ");
	gmp_printf("%Zd\n", ret_setup.qval);
	
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

	// Generate users
	unsigned long num_of_users = 10;
	mpz_t* sids = (mpz_t*)malloc(sizeof(mpz_t)*num_of_users);
	for (unsigned long i = 0; i < num_of_users;i++) {
		element_t k1, k2, k3;
		mpz_init(sids[i]);
		do {
			mpz_urandomm(userID, state, limit);	
		} while(mpz_cmp_ui(userID, 0) == 0 || mpz_cmp(userID, ret_setup.pval) == 0 || mpz_cmp(userID, ret_setup.qval) == 0);
		enroll(sids[i], userID, k1, k2, k3);
	}

	
	
	element_t k1, k2, k3;
	gmp_printf("------------------Demonstrating for USERID - %Zd-----------------\n", userID);
	mpz_t sid;
	mpz_init(sid);
	enroll(sid, userID, k1, k2, k3);
	
	element_printf("K1: %B\nK2: %B\nK3: %B\n", k1, k2, k3);

	element_t val1, val2, val3, val4;
	element_init_G1(val4, ret_setup.pairing);
	element_init_GT(val1, ret_setup.pairing);
	element_init_GT(val2, ret_setup.pairing);
	element_init_GT(val3, ret_setup.pairing);
	element_pairing(val1, k2, ret_setup.generators[1]);
	element_pairing(val2, k3, ret_setup.generators[0]);
	element_mul(val4, k2, ret_setup.B_Omega);
	element_pairing(val3, k1, val4);
	#ifdef DEBUG
		if (element_cmp(val1, val2) == 0) {
			printf("SUCCESS - 1\n");
		} 
		if (element_cmp(val3, Aval) == 0) {
			printf("SUCCESS - 2\n");
		}
	#endif
	if (element_cmp(val1, val2) == 0 && element_cmp(val3, Aval) == 0) {
		printf("--------------------Verification of key successful - key well formed by enroll()----------------------------\n");
	} else {
		printf("Verification of key UNSUCCESSFUL - key not well formed by enroll()\n");
	}
	/*
		ENROLL ends
	*/
	/*
		SIGN STARTS
	*/
	char* msg = (char*)malloc(sizeof(char)*(1+mpz_get_ui(ret_setup.mval)));

	for (unsigned int i = 0; i < mpz_get_ui(ret_setup.mval); i++) {
		int b = rand()%2;
		msg[i] = '1';
		if (b == 0) {
			msg[i] = '0';
		}
	}
	msg[mpz_get_ui(ret_setup.mval)] = '\0';
	printf("Message (m bits): %s\n", msg);
	element_t sigma1, sigma2, sigma3, sigma4, pi1, pi2;
	sign(pi1, pi2,sigma1, sigma2, sigma3, sigma4, k1,k2,k3, msg);


	element_printf("Final Signature:\nSigma1: %B\nSigma2: %B\nSigma3: %B\nSigma4: %B\nPi1: %B\nPi2: %B\n",sigma1, sigma2, sigma3, sigma4, pi1, pi2);
	/*
		SIGN ENDS
	*/

	/*
		VERIFY STARTS
	*/

	verify(Aval, sigma1, sigma2, sigma3, sigma4, pi1, pi2, msg);

	/*
		VERIFY ENDS
	*/

	/*
		TRACE STARTS
	*/
	trace(sigma2, sids, num_of_users);
	/*
		TRACE ENDS
	*/

	return 0;

}