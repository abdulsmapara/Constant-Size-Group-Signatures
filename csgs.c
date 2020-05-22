/*
			ASSIGNMENT 1

		Number Theory & Modern Cryptography
				VNIT, Nagpur

Submitted by:			Guided by:
Abdul Sattar Mapara		Dr. Syed Taqi Ali Sir
BT16CSE053

ABOUT:
- Guide:  Dr. Syed Taqi Ali Sir
- Author: Abdul Sattar Mapara
- Course: Number Theory and Modern Cryptography
- Paper Implemented: CSGS (https://crypto.stanford.edu/pbc/)

Short SUMMARY of the Paper:
- The aim of the scheme is to allow any member of a certain group to sign a message on behalf of the group, but the signer remains anonymous within the group. 
- However, in certain situations, an authority should have the ability to evoke the anonymity of a signer and trace the signature.
- Use in Anonymous attestation, which has practical applications such as in building Trusted Platform Modules (TPMs)

CONCEPTS explained:
- Pairing based cryptography
	* Cryptography with use of a pairing between elements of two cryptographic groups to a third group with a mapping e: G1*G2-> GT.
	* If G1=G2, then the pairing is called SYMMETRIC PAIRING, which is used in this program
	* If G1 != G2, then the pairing is called ASYMMETRIC PAIRING
- BILINEAR MAP:
	* The map from two groups G1, G2 to a third group GT is the bilinear map.
	(In pbc library, the G1 and G2 groups associated with pairings, are groups of points on an ellitpic curve; GT group is currently implemented as a subgroup of a finite field)
	* Denoted (as observed at many places, including the research paper being implemented) with e
- Properties of bilinear pairing:
	* Bilinearity: e(g^a, h^b) = e(g^b, h^a) = e(g,h)^(ab), where g,h are generators of group G1 & G2 respectively.
	* Non-degenerate: If g,h are generator of G,H then e(g,h) is generator of GT

The comments/description should be read according to the following: 
- Read the comments before the struct declared just after importing the libraries.
- Read the comments just after the struct referred above
- Jump to the main function and start reading the comments
- Whenever a function call is made, jump to the function body and read comments written just above & inside the function body
- References are mentioned in the format [REF-<Number>]
- IMPORTANT NOTE: Whenever a FUNCTION/DATA-TYPE is encountered for the first time according to the flow mentioned above, it is explained. For the next time, the Duplicate explanation is avoided. If an explanation to FUNCTION/DATA-TYPE is required next time it is being read, please search it in index_of_references.txt (OR simply look at the end of this file), note the reference number & search the explanation of function/data-type by searching +[REF-<Number>].
  For example, if I want to search for details on a function whose reference number is [REF-2], then I will search <Plus symbol '+'>[REF-2] in this file.
[It took more effort to avoid duplicate explanation, than it would have taken by writing explanation again. Since, it was a good practice to avoid duplicate explanation, I have used such a technique.]
- Preferably, set word wrap to ON in the text editor to avoid horizontal scrolling at certain places.

Following is the WORKING CODE (with comments & descriptions at suitable places) in C programming language for the Section "5.2 Core Construction" of the paper. The output for the code is present in the file output_large.txt (for large values) and output_small.txt (for small values).
The code can be compiled using - 'gcc csgs.c -lpbc -lgmp' and can be executed using './a.out'
*/

// Importing the required libraries
#include <pbc/pbc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <assert.h>

/*
	struct setup_output OR setup_result is a structure.
	The ret_setup is a variable of this type (type - setup_output) that holds values of the variables set by setup function to be used by other functions.
	The data type of each variable is explained here. Use of each variable is explained later.
*/
typedef struct setup_output {
	/*
	+[REF-1] element_t
	It is a data type (offered by pbc library), that stores the elements of groups, rings and fields.
	*/
	element_t mk0;
	/*
	+[REF-2] mpz_t
	Data type offered by GMP library for storing large integers.
	*/
	mpz_t mk1;
	mpz_t tk;
	mpz_t n;
	
	// element_t explained in [REF-1]
	element_t G;

	/*
	+[REF-3] pairing_t
	pairing_t is a data type offered by pbc library, that stores the bilinear pairings.
	*/
	pairing_t pairing;
	
	// e = element_pairing() function

	// element_t explained in [REF-1]
	element_t gen_subgroup_q;
	element_t* generators;
	element_t B_Omega;
	element_t Aval;
	
	// mpz_t explained in [REF-2]
	mpz_t alpha;

	// element_t explained in [REF-1]
	element_t identity;	// identity of G

	/*
	+[REF-4] pbc_param_t
	Data type offered by pbc library to generate/store pairing parameters
	*/
	pbc_param_t p;

	// mpz_t explained in [REF-2]
	mpz_t omega;
	mpz_t pval, qval;
	mpz_t mval;
}setup_result;
/*
	+[REF-5] ret_setup is a variable of type setup_result (which is a structure). It stores/handles/manipulates the following values -

		a. element_t mk0;
		   mk0 is the first part of the master enrollment key, generated by the function setup. It is of type element_t, which is explained in[REF-1].
		b. mpz_t mk1;
			mk1 is the second part of the master enrollment key, generated by the function setup. It is of type mpz_t, which is explained in [REF-2].
		c. mpz_t tk;
			tk is the group manager’s tracing key, which is generated (assigned value) by the function setup. It is of type mpz_t, which is explained in [REF-2]
		d. mpz_t n;
			The setup function chooses n = p*q, where p and q are random primes of bit size ceil(log p (base 2)), ceil(log q (base 2)) = Θ(λ) > k. n is the variable storing the value of p*q. It is of type mpz_t, which is explained in [REF-2]
		e. element_t G;
			The function setup builds a bilinear group G of order n. 
			The bilinear group is stored in the variable G. It is of type element_t, which is explained in [REF-1].
		f. pairing_t pairing;
			The A1 type pairing generated is required by many other functions later as well. So, it is stored in the variable pairing of ret_setup. It is of the type pairing_t, which is explained in [REF-3]	
		g. element_t gen_subgroup_q;
			gen_subgroup_q is the generator of Gq (cyclic subgroup of G of order q).
			It is of type element_t, which is explained in [REF-1].
		h. element_t* generators;
			generators is an array which stores the required generators of the group G.
			Every value of generators is of type element_t, which is explained in [REF-1].
		i. element_t B_Omega;
			B_Omega is a variable that stores Ω = g^(ω) ∈ G, where ω (omega, explained in n.) is a random exponent ∈ Zn.
			It is of type element_t, explained in [REF-1].
		j. element_t Aval;
			g is a generator of G, stored at 0th index in array generators and e is a bilinear map (details explained while using the map in setup function). Aval stores the value of A = e(g, g)^(alpha), where alpha (explained in k.) is a random exponent ∈ Zn
		k. mpz_t alpha;
			alpha is a random exponent ∈ Zn
			It is of type mpz_t, explained in [REF-2]
		l. element_t identity;
			identity is the identity of group G
			It is of type element_t explained in [REF-1].
		m. pbc_param_t p;
			p generates/stores the pairing parameters, used to initialize a pairing and later display the pairing parameters on the console.
			It is of type pbc_param_t, explained in [REF-4].
		n. mpz_t omega;
			omega is a random exponent ∈ Zn.
			It is of type mpz_t, explained in [REF-2].
		o. mpz_t pval, qval;
			Recall that n = p*q, where p and q are primes.
			pval stores the value of p, while qval stores the value of q.
			pval and qval are of type mpz_t, explained in [REF-2]
		p. mpz_t mval;
			mval stores the value of m used in setup (Refer the paper).
			It is of type mpz_t,explained in [REF-2].

		The bilinear map e, is not explicitly stored as it can be used by calling the function element_pairing provided by the pbc library. 
*/
setup_result ret_setup;


// generates a random prime number of n bits and stores it in the parameter param.
void random_prime_bits(mpz_t result, mpz_t n) {

	// Declare variable random of type mpz_t (Refer [REF-2] for details on mpz_t)
	mpz_t random;
	// Initialize the variable random (Refer [REF-7] for details on mpz_init function)
	mpz_init(random);

	/*
	+[REF-11] gmp_randstate_t
		Random state means an algorithm selection and current state data. The data type (provided by gmp library) for such objects is gmp_randstate_t. For example: gmp_randstate_t rstate;
	*/
	// Declare state
	gmp_randstate_t state;
	
	/*
	+[REF-12] void gmp_randinit_default (gmp randstate t state)
		Function provided by gmp library that initializes the state variable with a default algorithm.
	*/
	// Initialize state
	gmp_randinit_default(state);

	/*
	+[REF-13] void gmp_randseed_ui (gmp randstate t state, unsigned long int seed)
		Function provided by gmp library to set an initial seed value into state.
	*/
	// Set initial seed value to state. We pass seed as a random number defined as (random_number + 1)*(another-random-number) + 1. 1 was added to avoid seed = 0.
	gmp_randseed_ui(state, (rand()+1)*(rand()+1));

	/*
	+[REF-14] int mpz_cmp_ui (const mpz t op1, unsigned long int op2)
		Function provided by GMP library to compare op1 and op2. 
		It returns a positive value if op1 > op2, zero if op1 == op2, or a negative value if op1 < op2.
	
	Pseudo code for self-implementation:
		int retval = 0
		if (op1 > op2) {
			retval = +1;
		}
		else if (op1 < op2) {
			retval = -1;
		}
		return retval;
	*/
	// If number of bits <= 1, then no prime exists for that number of bits
	if (mpz_cmp_ui(n,1) <= 0) {
		printf("NO PRIME EXISTS\n");
	} else {
		// The prime generated should be greater than 2^(n-1) [As it should be n-bits]
		// Declare a variable to store the lower limit specified above
		mpz_t lower_limit;

		// Initialize the variable declared of type mpz_t
		mpz_init(lower_limit);
		/*
		+[REF-15] void mpz_ui_pow_ui (mpz t rop, unsigned long int base, unsigned long int exp)
			Function provided by gmp library to set rop to base^exp. According to the manual, 0^0 returns 1.
		Pseudo code for self-implementation:
			def power_helper(base, exp):
				if base == 1:
					return 1
				else if base == 0:
					return 0
				else if exp == 0:
					return 1
				else:
					save_num = base
					p = power_helper(base, exp/2)
					p = p*p
					if exp%2 == 1:
						p = p*save_num
					return p
			def power(base, exp):
				int rop
				if base == 0 and exp == 0:
					rop = 1
				else:
					rop = power_helper(base, exp)
		*/
		/*
		+[REF-30] unsigned long int mpz_get_ui (const mpz t op)
			Provided by gmp library, it returns the value of op as an unsigned long.
		*/
		mpz_ui_pow_ui(lower_limit, 2, mpz_get_ui(n)-1);
		// Loop till we find a prime number of n bits
		while (1){
			/*
			+[REF-16] void mpz_urandomb (mpz t rop, gmp randstate t state, mp bitcnt t n)
				Function provided by gmp library to generate a uniformly distributed random integer in the range 0 to 2^n − 1, inclusive.
			*/
			// Store a random value from 0 to (2^n)-1 in the variable random
			mpz_urandomb(random, state ,mpz_get_ui(n));

			/*
			+[REF-17] int mpz_cmp (const mpz t op1, const mpz t op2)
				Function provided by GMP library to compare op1 and op2. 
				It returns a positive value if op1 > op2, zero if op1 == op2, or a negative value if op1 < op2.
	
			Pseudo code for self-implementation:
				int retval = 0
				if (op1 > op2) {
					retval = +1;
				}
				else if (op1 < op2) {
					retval = -1;
				}
				return retval;
			*/
			/*
			+[REF-18] int mpz_probab_prime_p (const mpz t n, int reps)
				Function provided by gmp library to check if n is prime. Returns 2 if n is definitely prime, returns 1 if n is probably prime (without being certain), or return 0 if n is definitely composite.
				About the argument reps: It controls how many such tests are done. Larger value of reps will reduce the chances of a composite being returned as probably prime.

				The function (using Rabin-Miller primality testing algorithm) is also implemented on own (without using mpz_probab_prime_p) as follows -

			*/
			// typedef enum{FALSE,TRUE} boolean;
			// boolean isPrimeUtil(mpz_t num, mpz_t base_mpz){

			// 	boolean ret = FALSE;
			// 	//to check is num is prime or not
				
			// 	mpz_t temp,res,count;
			// 	mpz_t num_minus_1;
			// 	mpz_t T;

			// 	mpz_init(T);
			// 	mpz_init(temp);
			// 	mpz_init(res);
			// 	mpz_init(count);
			// 	mpz_init(num_minus_1);

				

			// 	mpz_set_ui(count,0);
			// 	mpz_sub_ui(temp,num,1);
			// 	mpz_sub_ui(num_minus_1,num,1);
				
			// 	mpz_mod_ui(res,temp,2);
			// 	while(mpz_cmp_ui(res,0) == 0){
			// 		mpz_add_ui(count,count,1);
			// 		//temp /= 2
			// 		mpz_div_ui(temp,temp,2);
			// 		mpz_mod_ui(res,temp,2);
			// 	}
			// 	// mpz_out_str(stdout,10,temp);
			// 	// printf("\n");
			// 	// mpz_out_str(stdout,10,count);
			// 	// printf("\n");

			// 	//	m = temp and k = count
			// 	mpz_powm(T,base_mpz,temp,num);
			// 	//	T = (base_mpz ^ temp ) % NUM

			// 	if((mpz_cmp_ui(T,1) == 0) || (mpz_cmp(T,num_minus_1) == 0)){
			// 		ret = TRUE;
			// 	}else{
			// 		/*int i = 1;
			// 		while(i <= count-1){
			// 			T = (T*T)%num;
			// 			if(T == 1){
			// 				ret = FALSE;
			// 				break;
			// 			}else if (T == -1){
			// 				ret = TRUE;
			// 				break;
			// 			}else{
			// 				i = i+1;	
			// 			}
			// 		}*/
			// 		mpz_t i;
			// 		mpz_init(i);
			// 		mpz_set_ui(i,1);
			// 		boolean break_loop = FALSE;
			// 		while((!break_loop) && mpz_cmp(i,count) < 0){
			// 			mpz_powm_ui(T,T,2,num);
			// 			if(mpz_cmp_ui(T,1) == 0){
			// 				ret = FALSE;
			// 				break_loop = TRUE;
			// 			}else if(mpz_cmp(T,num_minus_1) == 0){
			// 				ret = TRUE;
			// 				break_loop = TRUE;
			// 			}else{
			// 				mpz_add_ui(i,i,1);
			// 			}
			// 		}

			// 	}
			// 	return ret;
			// }
			// boolean isPrime(mpz_t num){
			// 	if(mpz_cmp_ui(num,2) == 0 || mpz_cmp_ui(num,3) == 0){
			// 		return TRUE;
			// 	}else if(mpz_cmp_ui(num,2) < 0){
			// 		return FALSE;
			// 	}else{
			// 		mpz_t modulo;
			// 		mpz_init(modulo);
			// 		mpz_mod_ui(modulo,num,2);
			// 		if(mpz_cmp_ui(modulo,0) == 0){
			// 			return FALSE;
			// 		}
			// 	}
			// 	mpz_t random_num;
			// 	mpz_init(random_num);
				
			// 	mpz_set_ui(random_num,2);

			// 	if(!isPrimeUtil(num,random_num)){
			// 		// first checked using base 2
			// 		return FALSE;
			// 	}

			// 	gmp_randstate_t state;
			// 	gmp_randinit_mt(state);

			// 	mpz_urandomb(random_num,state,2048);

			// 	mpz_t num_minus_3;
			// 	mpz_init(num_minus_3);
			// 	mpz_sub_ui(num_minus_3,num,3);
			// 	mpz_mod(random_num,random_num,num_minus_3);

			// 	mpz_add_ui(random_num,random_num,2);
			// 	return isPrimeUtil(num,random_num);
			// }
			
			// We are sure that random is smaller than (2^n)-1. To confirm that random is of n-bits, we check if random is greater than the lower limit.
			// Then check whether random number is prime we use the probabilistic function
			if (mpz_cmp(random, lower_limit) > 0 && mpz_probab_prime_p(random,mpz_get_ui(n))) {
				// If random is of n-bits and is prime (probably or for sure), then set result to the random number generated and return
				mpz_set(result,random);	// mpz_set is explained in [REF-9]
				break;
			}
		}			
	}

}

void helper_setup() {
	// The public information consists of the bilinear group, (n, G, GT , e), and a few public values.
	// Display those to the user.
	printf("PUBLIC INFORMATION:\n");
	// Print the value of n to console
	gmp_printf("N: %Zd\n", ret_setup.n);
	/*
		We print the group in the following way (using generator)-
		- Pick up the generator and store it in grp
		- Set add = generator
		- Print grp
		- Compute grp = grp <operator> add // Refer [REF-33] for explanation on how grp <operator> add is computed (using element_add).
		- If grp is identity, print and break
		- Else, print and continue 
	*/
	printf("G: ");
	// Declare, initialize and set some variables used for printing group G
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
	// 	element_add(grp, grp, add); // element_add explained in [REF-33]
	// }while(element_cmp(grp, ret_setup.identity));
	// printf("\n");
	/*
		End for printing group
	*/
	// Print the pairing parameters for e and GT.
	// e is the bilinear map, which can be accessed using the function element_pairing()
	printf("FOR e and GT, the pairing parameters are as follows:\n");
	/*
	+[REF-34] void pbc_param_out_str(FILE *stream, pbc_param_t p)
		Provided by pbc library, it is used to write pairing parameters to 'stream' in a text format. To print on console, we provide stream = stdout.
	*/
	pbc_param_out_str(stdout, ret_setup.p);
	// Can print GT using one of its generators also. Pseudo code for printing GT -
	/*
		// g is the generator of G
		// u is another generator of G
		// From non-degenerate property, e(g,u) is the generator of GT
		//So, 
		gen_GT = e(g,u)
		store_gen = gen_GT
		while (gen_GT != identity-of-GT) {
			print(gen_GT)
			gen_GT = element_add(gen_GT, gen_GT, store_gen) // equivalent to gen_GT = gen_GT <operator> store_gen
		}
	*/
	/* C CODE USING ABOVE EXPLANATION to print elements of group GT. UNCOMMENT BELOW CODE FOR PRINTING GT*/
	// element_t gt, identity_gt, store_gen;
	// element_init_GT(identity_gt, ret_setup.pairing);
	// element_init_GT(store_gen, ret_setup.pairing);
	// element_init_GT(gt, ret_setup.pairing);
	// element_pairing(gt, ret_setup.generators[0], ret_setup.generators[1]);
	// element_set0(identity_gt);
	// element_set(store_gen, gt);
	// printf("GT-");
	// while(element_cmp(gt, identity_gt)) {
	// 	element_printf("%B, ", gt);
	// 	element_mul(gt, store_gen, gt);
	// }

	
	/* 
		Print Generators of G using generators array
		Uncomment the below code for printing the 1+m+2 generators
	*/
	// printf("m+3 GENERATORS of GROUP G: ");
	// for (unsigned long int i = 0; i < mpz_get_ui(security_parameter)+3; i++) {
	// 	element_printf("%B,",ret_setup.generators[i]);
	// }
	// printf("\n");
	/*
		Done printing of generators of G
	*/

	// Next, we compute A = e(g, g)^α ∈ GT
	// Declare variable for storing value of A and initialize with pairing (type-A1-gen computed and store earlier)
	element_t Aval;
	element_init_GT(Aval, ret_setup.pairing);
	// To compute e(g,g) i.e to map g,g from G1,G1 respectively to an element in GT, we use the following function element_pairing
	/*
	+[REF-35] void element_pairing(element_t out, element_t in1, element_t in2)
		Provided by pbc library, it computes a pairing: out = e(in1, in2), where in1, in2, out must be in the groups G1, G2, GT. (In our case, G1=G2)
	*/
	element_pairing(Aval, ret_setup.generators[0], ret_setup.generators[0]);// Currently, Aval = e(g,g)
	element_pow_mpz(Aval, Aval, ret_setup.alpha); // Set Aval = e(g,g) ^ alpha
	printf("A: ");

	// Set value of A computed in ret_setup
	element_set(ret_setup.Aval, Aval);

	// Print value of A (part of public values) to console (element_printf explained in [REF-29])
	element_printf("%B\n", Aval);
	element_printf("%B\n", ret_setup.Aval);

	// Print remaining values - B_Omega (part of public values), master key and tracing key
	printf("OMEGA: ");
	element_printf("%B\n", ret_setup.B_Omega); // element_printf explained in [REF-29]
	printf("MASTER KEY: ");
	element_printf("%B, ",ret_setup.mk0); // element_printf explained in [REF-29]
	gmp_printf("%Zd\n", ret_setup.mk1);
	printf("TRACING KEY: ");
	gmp_printf("%Zd\n", ret_setup.qval);// element_printf explained in [REF-29]
	
	// Hence, displayed the Public Information and the public values	to the console
	// Also, displayed the Master enrollment Key, MK and the group manager's tracing key to the console.
}

/* 
Setup algorithm from Core Construction 5.2
The input is a security parameter in unary, 1^λ (here in code just λ i.e security_parameter in the code below represents the value of λ)
*/
void setup(setup_result* retval, mpz_t security_parameter) {

	// NOTE: ret_setup defined in main function is referred in this function by the name retval. The address of ret_setup is passed in retval, so that we can access the original ret_setup defined in main function. 

	// start choosing k and m = O(λ)
	mpz_t k,m; // (Refer [REF-2] for explanation on mpz_t data type)
	
	// Initialize k, m (Refer [REF-7] for explanation on mpz_init)
	mpz_init(k);
	mpz_init(m);

	/*
	+[REF-9] void mpz_set(mpz_t a, mpz_t b)
		The function (provided by GMP library) is used to set the value of a equal to the value of b.

	Pseudo code for self-implementation (i.e without using predefined function):
		b = a;
	*/
	// set the value of k, m equal to value of λ
	mpz_set(k, security_parameter);
	mpz_set(m, security_parameter);

	// Initialize the value of m to be stored in ret_setup  (Refer [REF-7] for explanation on mpz_init)
	mpz_init(retval->mval);

	// Store the value of m in the ret_setup variable, so that it can be used later by other functions (Refer [REF-9] for explanation on mpz_set)
	mpz_set(retval->mval, m);

	// start choosing p, q - random prime numbers
	// Declare no. of bits in p and q as p_bits, q_bits
	mpz_t p_bits, q_bits; // (Refer [REF-2] for explanation on mpz_t data type)
	
	// Initialize the number of bits in p and q. (Refer [REF-7] for explanation on mpz_init)
	mpz_init(p_bits);
	mpz_init(q_bits);

	/*
	+[REF-10] void mpz_add_ui (mpz t rop, const mpz t op1, unsigned long int op2) 
		Sets rop equal to addition of op1 and op2
	
	Pseudo code for self-implementation (without using pre-defined function):
		rop = op1 + op2;
	*/
	// add some number to p_bits and q_bits to achieve p_bits, q_bits = Θ(λ) > k
	// Since, k is order(λ), p_bits and q_bits will be order(λ) and by the following operation, we ensure that p_bits and q_bits > k.
	mpz_add_ui(p_bits, k, 2);
	mpz_add_ui(q_bits, k, 3);

	// Declare actual p and q.
	mpz_t p,q; // (Refer [REF-2] for explanation on mpz_t data type)

	// Initialize p and q (Refer [REF-7] for explanation on mpz_init)
	mpz_init(p);
	mpz_init(q);

	// Initialize values of p and q in the ret_setup, that is used by other functions
	// (Refer [REF-7] for explanation on mpz_init)
	mpz_init(retval->pval);
	mpz_init(retval->qval);
	
	// Generate the value of p, with no. of bits = p_bits
	/*
	void random_prime_bits(mpz_t result, mpz_t n) 
		A custom function (written on own), that generates a prime number of n bits and stores in the parameter result.
		Function details explained in the function body.
	*/
	random_prime_bits(p, p_bits);

	// Generate the value of q, with no. of bits = q_bits. Ensure that p is not equal to q.
	do{
		random_prime_bits(q, q_bits);
	}while(mpz_cmp(p,q) == 0);	// mpz_cmp is explained in [REF-17]

	// gmp_printf("Chosen p,q = %Zd\n%Zd\n", p,q);

	// Start calculating the value of n
	// Declare n (large integer using gmp library)
	mpz_t n; // mpz_t explained in [REF-2]

	// Initialize n
	mpz_init(n);	// mpz_init explained in [REF-7]
	
	/*
	+[REF-19] void mpz_mul (mpz t rop, const mpz t op1, const mpz t op2)
		Function provided by gmp library to set rop = op1*op2
	
	Pseudo code:
		rop = op1*op2;
	*/
	// Set n = p*q (n is composite)
	mpz_mul(n, p, q);

	// Initialize value of n in ret_setup variable (that is kind of returned by setup function)
	mpz_init(retval->n);
	// Set value of n in ret_setup variable
	mpz_set(retval->n, n);

	/*
	+[REF-20] int gmp_printf (const char *fmt, . . .)
		Function provided by gmp library that prints the variable of a data type provided by gmp library
	*/
	// Print values of n, p and q to console
	gmp_printf("N = %Zd\n", n);
	gmp_printf("P = %Zd\n", p);
	gmp_printf("Q = %Zd\n", q);

	// Set value of p and q in ret_setup variable
	mpz_set(retval->pval, p);
	mpz_set(retval->qval, q);
	
	// Start generating cyclic Bilinear group G of order n with SYMMETRIC PAIRING
	
	// Declare pairing to store bilinear pairings
	pairing_t pairing; // pairing_t is explained in [REF-3]
	
	// Declare par to store pairing parameters
	pbc_param_t par; // pbc_param_t explained in [REF-4]

	/*
	+[REF-21] void pbc_param_init_a1_gen(pbc_param_t param, mpz_t n)
		Function provided by pbc library.
		According to pbc manual, it generate type A1 pairing parameters and store them in param. 
		The order of bilinear group (that will be generated using param) will be n. 
		Example provided in manual:  n a product of two primes. (Exactly what we needed to generate the group).
	*/
	// Generate type A1 pairing parameter (to be used to initialize pairings, which will be used to generete groups, etc.) and store in par variable
	pbc_param_init_a1_gen(par, n);
	// Generate type A1 pairing parameter and store in ret_setup variable
	pbc_param_init_a1_gen(retval->p, n);

	/*
	+[REF-22] void pairing_init_pbc_param(pairing, pbc_param_t p)
		Function provided by pbc library.
		According to PBC manual, it initializes a pairing (to be directly used for generating groups etc.) with pairing parameters p.
	*/
	// Initialize pairings with parameter par generated above (type- A1 pairing parameter)
	pairing_init_pbc_param(pairing, par);

	/*
		bilinear map: G1*G2 -> GT
		We use G1 (or g1) and G2 to denote the input groups to the pairing, and GT (or gt1) for the output group. 
		All have order n.
		Here, we use symmetric pairing, so G2=G1 (say = G)
	*/
	// Declare variables to store elements of group
	element_t g1, gt1, identity, h, add, temp, mk0; // Refer [REF-1] for details on element_t
	
	// Initialize pairings in ret_setup with parameter par generated above (type- A1 pairing parameter)
	pairing_init_pbc_param(retval->pairing, par); // Refer [REF-22] for details on pairing_init_pbc_param
	/*
	+[REF-23] void element_init_G1(element_t e, pairing_t pairing)
		The function initializes e to be an element of the group G1 of the pairing (initialized above in the code using type A1 pairing parameters).
		When an element is initialized it is associated with an algebraic structure, such as a particular finite field or elliptic curve group.
	*/
	element_init_G1(g1, pairing);
	// Here, G1 = G2, so element_init_G1 and element_init_G2 will do exactly the same operations.
	element_init_G2(g1, pairing);

	/*
	+[REF-24] void element_init_GT(element_t e, pairing_t pairing)
		The function initializes e to be an element of the group GT of the pairing (initialized above in the code using type A1 pairing parameters).
	*/
	element_init_GT(gt1, pairing);

	/*
	NOTE: Please refer [REF-23] for details on element_init_G1() and [REF-24] for details on element_init_GT(). 
	Since, G1 = G2, element_init_G2() does exactly what element_init_G1() does.
	*/
	element_init_G1(add, pairing);
	element_init_G1(temp, pairing);
	element_init_G1(h, pairing);
	element_init_G1(identity, pairing);
	element_init_G1(retval->G, pairing);
	element_init_GT(retval->Aval, pairing);
	element_init_G2(mk0, pairing);

	/*
	+[REF-25] void element_set0(element_t e)
		For groups of points on an ellitpic curve, such as the G1 and G2 groups associated with pairings, both 0 and 1 represent the identity element. 
		Hence, setting element to 0 means setting element to identity element of G1/G2.
	*/
	// set identity to identity of group G (of order n)
	element_set0(identity);
	// Initialize identity of group G stored in ret_setup with pairing
	element_init_G1(retval->identity, pairing);

	// Set identity in ret_setup to the identity of group G (of order n)
	element_set0(retval->identity);

	/*
	+[REF-26] void element_set(element_t e, element_t a)
		Provided by pbc library, it sets element e to a
	*/
	// Set G in ret_setup with g1 used here.
	element_set(retval->G, g1);
	// Below, we select a generator h of cyclic subgroup of G of order q
	/*
	How to select generator h of cyclic subgroup of G of order q.
	We use the following for selecting h -
	From https://crypto.stanford.edu/pbc/notes/numbertheory/cyclic.html -
	- Corollary: Euler’s Theorem (and Fermat’s Theorem). Any a belongs to Zn* (The * in Zn* stresses that we are only considering multiplication and forgetting about addition) generates a cyclic subgroup {a, a^2, ... a^d = 1} thus d|Phi(n) and hence a^Phi(n) = 1
	- Theorem:  All subgroups of a cyclic group are cyclic. If G=<g> is a cyclic group of order n  then for each divisor d of n there exists exactly one subgroup of order d and it can be generated by a^(n/d)
	So, we can choose some ‘a’ that belongs to G and generate the subgroup by generators- a^(n/p) = a^q and similarly a^p
	So, we choose a random element from group G and calculate (that element)^(p) to get h (generator of subgroup of G of order q)
	*/
	do {
		/*
		+[REF-27] void element_random(element_t e)
			Function provided by pbc library, it assigns a uniformly random element to e (as e lies in a group).
		*/
		// Choose some random element from group G
		element_random(g1);	
		
		/*
		+[REF-28] void element_pow_mpz(element_t rop, element_t e, mpz_t a)
			Provided by pbc library
			Computes e^a and assigns to rop. Here, e is element of group (or any other structure i.e of type element_t) and a is of type mpz_t.
		*/

		// Compute h as g1^ p
		element_pow_mpz(h, g1, p);

		// The loop makes sure that h is not identity element of group G and some other element
	} while(element_cmp(h, identity) == 0);
	
	// The below code sets the value of generator of subgroup of G, of order q in ret_setup to be used by other functions later on
	element_init_G1(retval->gen_subgroup_q, pairing); /*Initialize that value in ret_setup with pairing*/
	#ifdef DEBUG
		/*
		+[REF-29] int element_printf(const char *format, …)
			Provided by pbc library, it prints the element to console.
		*/
		element_printf("h = %B\n", h);
	#endif
	// Set h in ret_setup
	element_set(retval->gen_subgroup_q, h); // Details on element_set available in [REF-26]

	/*
	+[REF-33] void element_add(element_t n, element_t a, element_t b)
		Provided by pbc library, it sets n = a + b
		NOTE: For groups of points on an ellitpic curve, such as the G1 and G2 groups associated with pairings, the operation - addition and multiplication represent the group operation.
		Hence, we can use element_add to compute a <operation> b, where a and b are elements of group
	*/
	/* ******************************
		Verify that h is a generator (for small values)
		Uncomment the below code only for small values
	   ******************************
	   Here, we start with h, and compute h = h operator h, till h becomes identity element.
	   We keep the count of number of steps required for h to become identity element i.e we store the order of h. If order of h == q, then h is a generator of subgroup of G order q.
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
	// 	element_add(h, h, add);	// element_add explained in [REF-33]
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

	// Below we choose m+3 generators of group G - one is g and other m+2 generators.
	mpz_t required, gen;
	mpz_init(required);
	mpz_init(gen);
	/*Set required number of generators to m + 3 */
	mpz_add_ui(required, m,3); // For details on mpz_add_ui, refer [REF-10]

	// Declare array of size m+3 to store all generators of group G
	element_t* generators = (element_t*)malloc(sizeof(element_t)*(3+mpz_get_ui(m))); // mpz_get_ui explained in [REF-30]
	
	// Initialize all generators with pairing (as all generators are also elements of a group)
	for (unsigned long int i = 0; i < (3+mpz_get_ui(m)); i++) { // mpz_get_ui explained in [REF-30]
		element_init_G1(generators[i], pairing);
	}

	unsigned long long int index = 0;
	// Choosing m+3 (one generator g and m+2 others) generators for G by the following method -
	/*	To select a generator g from G, we use the following trick -
		- Select a random element, say ‘r’ from G
		- Check if r^n is identity element or not
		- If not, choose random element again and continue
		- If r^(n) is identity element, then proceed
		- If r^p is not identity element and r^q is not identity element, then r is the generator
		- Else, not. We repeat the process from selecting the random element
		This trick will work-
		- Since p and q are primes, and n = p*q, n has only 2 factors p and q.
		- We need to show that the order of r is n i.e from r^1 to r^n, only r^n is an identity element.If r^n is an identity element, only then r can be a generator. So, let r^n be an identity element. We know that things get repeated after generating an identity element. So, if there is some x, such that r^x = identity element, then n should be a multiple of x. Since, p and q are the only factors of n, and we check r^p is not an identity element and r^q is not an identity element, we can argue that there is no x such that r^x = identity element and x != n.
	*/
	do {
		// Choose a random element from G
		element_random(g1);	// element_random explained in [REF-27]
		// Compute random element power n
		element_pow_mpz(temp, g1, n); // element_pow_mpz explained in [REF-28]

		// Compare random element power n with identity, if matched then proceed, else discard
		if (element_cmp(temp, identity) == 0) {
			// Compute random element power p
			element_pow_mpz(temp, g1, p);
			// Check if random element power p equals identity. If equals, discard, else proceed
			if (element_cmp(temp, identity)) {
				// Compute random element power q
				element_pow_mpz(temp, g1, q);
				// Check if random element power q equals identity. If equals, discard, else proceed
				if (element_cmp(temp, identity)) {
					// Got a generator of G. Add to the count.
					mpz_add_ui(gen, gen, 1);
					// Put the generator at the next index in the array for storing generators.
					element_set(generators[index], g1);
					index++;
				}
			}
		}
	// Loop till number of generators generated < number of generators required to be generated
	} while(mpz_cmp(gen, required)); // mpz_cmp explained in detail in [REF-17]

	retval->generators = generators; // Provide address of array to ret_setup, so that the generators of G can be used by other functions

	// Below, we pick random exponents alpha, omega from Zn (a step from implementation descried in paper)
	// declare and initialize alpha and omega
	mpz_t alpha, omega;
	mpz_init(alpha);
	mpz_init(omega);

	// Using gmp library to pick a random number (in the same way as picking random no. in the function random_prime_bits)
	gmp_randstate_t state; // gmp_randstate_t described in [REF-11]
	gmp_randinit_default(state); // gmp_randinit_default described in [REF-12]
	gmp_randseed_ui(state, (rand()+1)*(rand()+1)); // gmp_randseed_ui described in [REF-13]
	
	// Below, we make sure that alpha's gcd with p and q both is 1 i.e. alpha is relatively prime with p and q.
	// Declare and initialize variables for storing gcd of p,q and alpha
	mpz_t gcd_alpha_p, gcd_alpha_q;
	mpz_init(gcd_alpha_p);
	mpz_init(gcd_alpha_q);

	do{
		/*
		+[REF-31] void mpz_urandomm (mpz t rop, gmp randstate t state, const mpz t n)
			Provided by gmp library, it generates a uniform random integer in the range 0 to n − 1, both included.
		*/
		// Choose alpha randomly between & including 0 and n-1 i.e Zn
		mpz_urandomm(alpha, state, n);
		/*
		+[REF-32] unsigned long int mpz_gcd (mpz_t rop, const mpz_t op1, const mpz_t op2)
			Provided by gmp library, it computes the greatest common divisor of op1 and op2 and stores the result obtained after computation to rop
			Pseudo-code/Code for gcd implemented without function mpz_gcd is as follows-
			void gcd_self(mpz_t gcd, mpz_t a, mpz_t b) {
				mpz_t zero;
				mpz_init(zero);
				while (mpz_cmp(b,zero) > 0 ) {
					mpz_t r;
					mpz_init(r);
					mpz_mod(r, a, b);
					mpz_set(a, b);
					mpz_set(b, r);
				}
				mpz_set(gcd, a);
				
			}			
		*/

		// Compute gcd of alpha, p and alpha, q
		mpz_gcd(gcd_alpha_p, alpha, p);
		mpz_gcd(gcd_alpha_q, alpha, q);

	// Loop till gcd(alpha, p) and gcd(alpha, q) is not 1.
	} while(mpz_cmp_ui(gcd_alpha_p, 1) != 0 || mpz_cmp_ui(gcd_alpha_q, 1) != 0);

	// Choose omega randomly between & including 0 and n-1 i.e Zn
	mpz_urandomm(omega, state, n); // mpz_urandomm explained in [REF-31]

	// Below, we set values of omega and alpha in ret_setup
	mpz_init(retval->omega); // Initialize value of omega in ret_setup
	mpz_set(retval->omega, omega);	// Set value of omega in ret_setup
	gmp_printf("Alpha = %Zd\nOmega = %Zd\n", alpha, omega); // Print values of alpha and omega to console
	mpz_init(retval->alpha);	// Initialize value of alpha in ret_setup
	mpz_set(retval->alpha, alpha);	// Set value of alpha in ret_setup

	// Next, we compute Ω = g^ω ∈ G (as directed by research paper) and set to variable B_Omega in ret_setup
	element_pow_mpz(g1, generators[0], omega);	// element_pow_mpz explained in [REF-28] 
	// Initialize B_Omega in ret_setup with pairing
	element_init_G1(retval->B_Omega, pairing);
	// Set B_Omega with g^w computed (& stored in g1)
	element_set(retval->B_Omega, g1);

	// Next, we compute the first part of master enrollment key, MK below
	// Compute g^alpha ∈ G and store in mk0
	element_pow_mpz(mk0, generators[0], alpha);
	// Initialize variable mk0 that stores the first part of MK in ret_setup with pairing
	element_init_G1(retval->mk0, pairing);
	// Set mk0 in ret_setup with the first part of master enrollment key, MK computed above
	element_set(retval->mk0, mk0);

	// Now, we calculate second part of master enrollment key, MK below
	// Intialize value of second part of MK in ret_setup
	mpz_init(retval->mk1);
	// Second part of MK = omega, so set the second part of MK to omega
	mpz_set(retval->mk1, omega);

	// Next, we compute the group manager’s tracing key, TK
	// Initialize the value of TK in ret_setup
	mpz_init(retval->tk);
	// TK = q, so set the value of TK in ret_setup to q
	mpz_set(retval->tk, q);

	// Remaining operations of setup are done by the function helper_setup (written on own)
	helper_setup();

	/* 
	The setup is complete and it stores all required values in ret_setup (so that it can be used later by other functions)
	*/
}

/* 
	enroll function
	Input: PP, MK, ID
	0 <= ID < 2^k < p (& q)
	Sets SID (secret unique ID) and signing key - K = (k1,k2,k3)
*/
static int num_user; // used to get a unique number for every user
void enroll(mpz_t final_sid, mpz_t userID, element_t k1, element_t k2, element_t k3) {
	// Choose unique SID for a given USER ID
	// SID is later used for tracing purposes
	/*
		Conditions for generating SID - 
		1. Should be unique for given user
		2. omega + sid must lie in Zn* (multiplicative group modulo n)
		
		SID generated follows the above conditions as we perform the following steps -
		1. We first assign sid = user ID
		2. If sid + ω belongs to Zn* (i.e. if gcd(sid+ω, n) is 1), we are done
		3. If not, we set sid = min(p,q) + num_user + 1(where num_user is initially 0) and check if sid + ω belongs to Zn*
		4. If not, we increment, num_user and check if sid + ω belongs to Zn*
		5. We repeat till we get such an sid
		6. Since num_user will be different for each user (as we use static variable), sid is unique
		7. The maximum number of users that can be enrolled = 2^ k (as userID has upper bound 2^k < p).
		   We know that 2^k < p (and q also, as no. of bits in q are greater than k).
		   Hence sid assigned in 3. is always different than any sid assigned in 1. (sid in 1. < p & q, while sid in 3. > p & q)
		   Moreover, num_user is unique for every user, so sid is uniquely generated in 3.
	*/
	// Declare & Initialize variables - sid to store sid, gcd to store gcd(sid + ω, n) 
	mpz_t gcd, sid;
	mpz_init(gcd);
	mpz_init(sid);

	// Compute gcd of userID and n
	mpz_gcd(gcd, userID, ret_setup.n);
	// Ensure that gcd(userID, n) is always 1 (as userID < p & q)
	assert(mpz_cmp_ui(gcd,1) == 0);

	// val stores the minimum of p and q
	mpz_t val;
	mpz_init(val);
	mpz_set(val, ret_setup.pval);
	if (mpz_cmp(ret_setup.pval,ret_setup.qval) > 0) {
		mpz_set(val, ret_setup.qval);
	}

	// Loop till we find sid such that gcd(sid + ω, n) == 1
	while (1){
		/*
		+[REF-37] void mpz_add (mpz_t rop, const mpz_t op1, const mpz_t op2)
			Provided by gmp library
			Sets rop = op1 + op2
		*/
		// Compute (sid + ω) % n, then take gcd of sid + ω with n, then compare with 1
		mpz_add(sid, userID, ret_setup.omega);
		mpz_mod(sid, sid, ret_setup.n);
		mpz_gcd(gcd, sid, ret_setup.n);

		// Increment value of num_user to make it unique every time
		num_user+=1;
		// If gcd is 1, sid is found, we can break
		if (mpz_cmp_ui(gcd, 1) == 0) {
			break;
		}
		// Otherwise set userID = min(p,q) + num_user (1 already added above to num_user)
		// Note that here userID stores the temporary value of sid.
		mpz_add_ui(userID, val, num_user);
	}

	// Now,we found the suitable sid (as per conditions)
	// NOTE: variable sid stores the value sid + ω and real sid is stored in variable userID
	// Next, we need to generate the signing key - k = (k1,k2,k3)

	// First step is to take inverse of sid + ω wrt modulo n (inverse exists as sid + ω is in Zn*)
	mpz_t inverse;
	mpz_init(inverse);
	/*
	+[REF-38] int mpz_invert (mpz_t rop, const mpz_t op1, const mpz_t op2)
		It computes the inverse of op1 modulo op2 and puts the result in rop. 
		If the inverse exists, the return value is non-zero and rop will satisfy 0 < rop < |op2|. If an inverse doesn’t exist then the return value is zero and rop is undefined.

		Pseudo code/code for mpz_invert (without library function mpz_invert):
		void ExtendedEucledian(mpz_t rop,mpz_t a, mpz_t b) {
			// rop = inverse(a) wrt modulo b
			mpz_t r,r1,r2,q,s,t, s1, s2, t1,t2;
			mpz_init(r);
			mpz_init(t1);
			mpz_init(t2);
			mpz_init(s1);
			mpz_init(s2);
			mpz_init(r1);
			mpz_init(r2);
			mpz_init(q);
			mpz_init(s);
			mpz_init(t);
			mpz_set(r1, a);
			mpz_set(r2, b);
			mpz_set_ui(t1, 0);
			mpz_set_ui(s1, 1);
			mpz_set_ui(t2 , 1);
			mpz_set_ui(s2, 0);
			while (mpz_cmp_ui(r2, 0) > 0) {
				
				mpz_tdiv_q(q,r1,r2);
				mpz_t mul;
				mpz_init(mul);
				mpz_mul(mul, q, r2);
				mpz_sub(r, r1, mul);
				mpz_set(r1, r2);
				mpz_set(r2 , r);
				mpz_mul(mul, q, s2);
				mpz_sub(s, s1, mul);
				mpz_set(s1, s2);
				mpz_set(s2, s);
				mpz_mul(mul, q, t2);
				mpz_sub(t, t1, mul);
				mpz_set(t1, t2);
				mpz_set(t2, t);
			}
			//gmp_printf("GCD = %Zd\n",r1);
			//gmp_printf("%Zd\n",s1);
			//gmp_printf("%Zd\n",t1);
			
			mpz_set(rop, s1);
			//gmp_printf("%Zd\n",rop);
		}

	*/
	// store inverse of sid + ω in variable inverse
	assert(mpz_invert(inverse, sid, ret_setup.n));

	// Print SID to console
	gmp_printf("sid (not to be disclosed to user): %Zd\n", userID);
	
	// Store SID in the parameter final_sid to be used by other functions [As SID is stored in variable with name userID]
	mpz_init(final_sid);
	mpz_set(final_sid, userID);

	// gmp_printf("DEBUG: inverse of sid+omega: %Zd\n", inverse);
	
	// k1, k2 and k3 are defined as follows -
	// k1, k2, k3 = ((g^alpha)^(1/w+sid), g^sid, u^sid)
	// NOTE: u is the next generator stored in the array generators

	// Initialize k1,k2,k3 to be elements of group G1 (Remember the bilinear map: G1*G1->GT)
	element_init_G1(k1, ret_setup.pairing);
	element_init_G1(k2, ret_setup.pairing);
	element_init_G1(k3, ret_setup.pairing);
	// compute k1 = g^alpha
	element_pow_mpz(k1, ret_setup.generators[0], ret_setup.alpha);
	// compute k1 = (current_k1)^(inverse of sid + ω)
	element_pow_mpz(k1, k1, inverse); // k1 is ready

	// compute k2 = g^sid [sid is stored in variable userID]
	element_pow_mpz(k2, ret_setup.generators[0], userID);

	// compute k3 = u^sid [sid is stored in variable userID].
	// NOTE: u is another generator stored in the array generators
	element_pow_mpz(k3, ret_setup.generators[1], userID);

	// ENROLL sets the parameters final_sid with sid assigned to user with userID - userID
	// It sets k1, k2 and k3 so that signing key is ready for use by other functions.
}
/*
	sign function
	For a given user with Kid = k1, k2, k3 and a message m (of m-bits), it returns a signature-
	Signature = (sigma1, sigma2, sigma3, sigma4, pi1, pi2)
*/
void sign(element_t pi1, element_t pi2,element_t sigma1,element_t sigma2,element_t sigma3,element_t sigma4, element_t k1, element_t k2, element_t k3, char* message) {
	// First, Kid is used to create a two-level hybrid signature with the message at the second level. 

	// Declare & Initialize integers (using gmp library)
	mpz_t s, t1, t2, t3, t4, t_temp;
	mpz_init(s);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(t3);
	mpz_init(t4);
	mpz_init(t_temp);

	// Below, we generate a random number s ∈ Zn
	// Declare state
	gmp_randstate_t state; //For details on gmp_randstate, refer [REF-11] 
	// Initialize state
	gmp_randinit_default(state); // For details on gmp_randinit_default, refer [REF-12]
	// Set initial seed value to state. We pass seed as a random number defined as (random_number + 1)*(another-random-number) + 1. 1 was added to avoid seed = 0.
	gmp_randseed_ui(state, (rand()+1)*(rand()+1)); // For details on gmp_randseed_ui, refer [REF-13]
	
	// Generate random number s in the range [0,n-1]
	mpz_urandomm(s, state, ret_setup.n);
	
	// First step in signature generation is to compute the (randomized but unblinded) hybrid signature - 01, 02, 03, 04 [Read 0 as theta]

	// 01 = k1 (straightforward)
	// 02 = k2 (straightforward)
	
	// Declare and initialize elements of groups (G1, GT suitably) with pairing generated before to be used for above tasks
	element_t theta3, val1; 
	element_init_G1(theta3, ret_setup.pairing); // theta3 belongs to G1
	element_init_G1(val1, ret_setup.pairing); // val1 belongs to G1

	// Compute 03 = v dash, the 3rd generator of group G (order n) stored in array generators
	element_set(theta3, ret_setup.generators[2]);

	// Declare & initialize another element of group G1
	element_t temp_pw;
	element_init_G1(temp_pw, ret_setup.pairing);

	// msg_val stores the current message bit in the loop below
	mpz_t msg_val;
	mpz_init(msg_val);

	// We need to multiply generators from 4th generator to (m+3)-th generator, each raised to power = message bit value. This value multiplied with v dash is stored in val1
	for (int i = 3; i < (3+mpz_get_ui(ret_setup.mval)); i++) {
		/*
		+[REF-40] void mpz_set_ui(mpz_t rop, unsigned long int a)
			Provided by gmp library, sets rop to value of a
		*/
		// Store the current message bit
		if (message[i-3] == '1') {
			mpz_set_ui(msg_val,1);
		} else {
			mpz_set_ui(msg_val, 0);
		}
		// gmp_printf("%Zd\n", msg_val);
		// Compute temp_pw = (current-generator) raised to power current message bit value
		element_pow_mpz(temp_pw, ret_setup.generators[i], msg_val);
		// Compute 03=03*(temp_pw) [Note that 03 was initialized with v dash (3rd generator of G)]
		element_mul(theta3, theta3, temp_pw);
	}
	// Set val1 = 03 so that current computed value of theta3 can be used later
	element_set(val1, theta3); //[REF-VAL-1]
	// Compute 03 = 03 raised to power s [s is the random number generated at the start of the function]
	element_pow_mpz(theta3, theta3, s);

	// Compute theta3 = theta3*k3
	element_mul(theta3, theta3, k3);
	// theta3 has been computed successfully

	// element_printf("theta3: %B\n", theta3);
	
	// Start computing theta4 (04) -
	// Declare and initialize theta4 - element belong to group G1 (with the pairing generated earlier)
	element_t theta4;
	element_init_G1(theta4, ret_setup.pairing);
	// Initialize 04 = g (1st generator of G)
	element_set(theta4, ret_setup.generators[0]);
	// Compute 04 = 04 raised to power s
	element_pow_mpz(theta4, theta4, s);
	/*
	+[REF-41] void element_invert(element_t a, element_t b)
		Provided by pbc library, it assigns inverse of element b to element a
	*/
	// Compute 04 = inverse(current value of theta4)
	element_invert(theta4, theta4);
	// Theta4 computed successfully

	// element_printf("theta4: %B\n", theta4);
	// ----Initial signature formed---

	/*
		----------TEST 2----------
		Verify the initial signature using the regular verification equations (refer research paper)
	*/
	// Declare & Initialize elements of groups with pairing
	element_t ver1, ver2, ver3;
	element_init_GT(ver1, ret_setup.pairing);	// ver1 is an element of group GT
	element_init_GT(ver2, ret_setup.pairing);	// ver2 is an element of group GT
	element_init_GT(ver3, ret_setup.pairing);	// ver3 is an element of group GT

	// Set ver1 = e(03, g)
	element_pairing(ver1, theta3, ret_setup.generators[0]);
	// element_printf("ver1: %B\nval1:%B\n", ver1, val1);
	
	// Set ver2 = e(04, val1) [val1 is stored at [REF-VAL-1]]
	element_pairing(ver2, theta4, val1);
	// element_printf("ver2: %B\n", ver2);
	
	// compute ver1 = current value of ver1 * current value of ver2
	element_mul(ver1, ver1, ver2);
	// element_printf("ver1: %B\tver2: %B\n",ver1, ver2);

	// Compute ver3 = e(02 = k2, u)
	element_pairing(ver3, k2, ret_setup.generators[1]);

	// The other verification check here is EXACTLY THE SAME as 1st verification check of key (in TEST 1)
	// It is done in main function just after enroll function returns
	// If it passes there, it will pass here. If it fails there, we will know at that point only (as TEST 1 is conducted before call to sign function). 
	// As this code is for demonstration, to unnecessary increase time complexity of the program, it is not done again

	// element_printf("ver1-final: %B\nver3-final: %B\n",ver1, ver3);

	// ver1 and ver3 should match for success of the verification
	if (element_cmp(ver1, ver3) == 0) {
		// Verification successful
		printf("-----------------Verification of initial signature successful------------------\n");
	} else {
		// Verification unsuccessful
		printf("Verification of initial signature UNSUCCESSFUL\n");
	}
	// Generate random numbers t1, t2, t3, t4 belonging to Zn
	mpz_urandomm(t1, state, ret_setup.n);
	mpz_urandomm(t2, state, ret_setup.n);
	mpz_urandomm(t3, state, ret_setup.n);
	mpz_urandomm(t4, state, ret_setup.n);

	// Print t1, t2, t3, t4, s in debug mode
	#ifdef DEBUG
		gmp_printf("t1, t2, t3, t4, s: %Zd, %Zd, %Zd, %Zd, %Zd\n", t1, t2, t3, t4, s);
	#endif

	// Declare, initialize some elements of group (G1/GT)
	element_t h, temp;
	// sigma1, sigma2, sigma3, sigma4 belong to group G1
	element_init_G1(sigma1, ret_setup.pairing);
	element_init_G1(sigma2, ret_setup.pairing);
	element_init_G1(sigma3, ret_setup.pairing);
	element_init_G1(sigma4, ret_setup.pairing);

	// h is the generator of cyclic subgroup of G of order q.
	// Initialize it with pairing as an element of group G1
	element_init_G1(h, ret_setup.pairing);

	// Initialize temp variable with pairing
	element_init_G1(temp, ret_setup.pairing);

	// Set h to generator of cyclic subgroup of G of order q using ret_setup
	element_set(h, ret_setup.gen_subgroup_q);

	// Initialize sigma1 = k1 = theta1
	element_set(sigma1, k1);
	// Initialize sigma2 = k2 = theta2
	element_set(sigma2, k2);
	// Initialize sigma3 = theta3
	element_set(sigma3, theta3);
	// Initialize sigma4 = theta4
	element_set(sigma4, theta4);

	// Compute h^t1 & store in temp
	element_pow_mpz(temp, h, t1);
	// Set sigma1 = temp*sigma1 = (h^t1) * current value of sigma1
	element_mul(sigma1, temp, sigma1);

	// Compute h^t2 and store in temp
	element_pow_mpz(temp, h, t2);
	// Set sigma2 = temp*sigma2 = (h^t2)*current value of sigma2
	element_mul(sigma2, temp, sigma2);
	
	// Compute h^t3 and store in temp
	element_pow_mpz(temp, h, t3);
	// Set sigma3 = temp*sigma3 = (h^t3)*current value of sigma3
	element_mul(sigma3, temp, sigma3);
	
	// Compute h^t4 and store in temp
	element_pow_mpz(temp, h, t4);
	// Set sigma4 = temp*sigma4 = (h^t4)*current value of sigma4
	element_mul(sigma4, temp, sigma4);

	// Additionally, we need to compute two group elements pi1, pi2 according to paper.
	// We do so in the following part of the code
	// Declare & Initialize some elements of group (G1/GT) using the pairing generated earlier
	element_t theta1, theta2, u;
	// pi1, pi2, theta1, theta2, u belong to group G1 (theta1, theta2 are computed earlier)
	element_init_G1(pi1, ret_setup.pairing);
	element_init_G1(pi2, ret_setup.pairing);
	element_init_G1(theta1, ret_setup.pairing);
	element_init_G1(theta2, ret_setup.pairing);
	element_init_G1(u, ret_setup.pairing);
	
	// Assign value k1 to theta1 as theta1 = k1
	element_set(theta1, k1);
	// Assign value k2 to theta2 as theta2 = k2
	element_set(theta2, k2);

	// u is the second generator of bilinear group G of order n
	element_set(u, ret_setup.generators[1]);

	// Compute & set t_temp = t1*t2
	mpz_mul(t_temp, t1, t2);
	// Compute h^t_temp = h^(t1*t2) & store in h
	element_pow_mpz(h, h, t_temp);

	// Compute theta1^t2 & store in theta1 only
	element_pow_mpz(theta1, theta1,t2);
	// Compute theta2*Ω & store in theta2 only
	element_mul(theta2, theta2, ret_setup.B_Omega);
	// Compute (original theta2*Ω)^t1 i.e. (current value of variable theta2)^t1
	element_pow_mpz(theta2, theta2, t1);

	// Compute pi1 = (current value of variable theta1)*(current value of variable h) [Refer above comments to know current values of the variables]
	element_mul(pi1, theta1, h);
	// Compute pi1 = (current value of pi1)*(current value of variable theta2)
	element_mul(pi1, pi1, theta2);
	// pi1 is ready. It is equal to 
	// [(original h, generator of subgroup of G of order q)^(t1*t2)]*[(original theta1)^t2]*[(original theta2*Ω)^t1]

	// Now, we compute pi2
	// Set pi2 = u^t2
	element_pow_mpz(pi2, u, t2);
	// Assign u = g^t3 (g is first generator of group G)
	element_pow_mpz(u, ret_setup.generators[0], t3);
	// Compute inverse of current value of u (=g^t3) and store in variable u
	element_invert(u, u);
	// Compute pi2 = 
	// [(original value of u i.e. the second generator of bilinear group G of order n)^t2]*(inverse of g^t3)
	element_mul(pi2, pi2, u);

	// Details on val1 - Multiply generators from 4th generator to (m+3)-th generator, each raised to power = message bit value. This value multiplied with v dash (3rd generator of bilinear group G) is stored in val1
	// Compute val1 = val1^t4. 
	element_pow_mpz(val1, val1, t4);
	// Inverse current value of val1 and store in val1
	element_invert(val1, val1);
	// Compute pi2 = current value of pi2 * current (modified) value of val1
	element_mul(pi2, pi2, val1);
	// pi2 is computed as [(original value of u i.e. the second generator of bilinear group G of order n)^t2]*(inverse of g^t3)*[inverse(original value of val1 as defined by [REF-VAL-1] ^ t4)]
	// Refer the paper for clear equation of pi1 & pi2

	// The sign function has computed the signature = (sigma1, sigma2, sigma3, sigma4, pi1, pi2) & set the parameters so that this signature can be used by other functions
}
/*
	verify function
	It validates a given group signature σ on a given message M
	If the group signature is validated successfully, the verifier outputs valid; otherwise, it outputs invalid.
*/
void verify(element_t Aval, element_t sigma1, element_t sigma2, element_t sigma3, element_t sigma4, element_t pi1, element_t pi2, char* message) {

	// Declare & intialize group elements using the pairing generated earlier
	element_t Acopy, val1, val2, pairing_result1, pairing_result2, T1, T2;
	element_init_G1(val2, ret_setup.pairing); // val2 belongs to group G1
	element_init_G1(val1, ret_setup.pairing); // val1 belongs to group G1
	element_init_GT(Acopy, ret_setup.pairing); // Acopy belongs to group GT
	element_init_GT(pairing_result1, ret_setup.pairing);	// pairing_result1 belongs to group GT
	element_init_GT(pairing_result2, ret_setup.pairing);	// pairing_result2 belongs to group GT
	element_init_GT(T1, ret_setup.pairing);	// T1 belongs to group GT
	element_init_GT(T2, ret_setup.pairing);	// T2 belongs to group GT

	// Set Acopy to the value of A (computed by setup, stored in ret_setup also, passed to verify by main function)
	element_set(Acopy, Aval);
	// Compute inverse of A (element of group G1) and assign to Acopy
	element_invert(Acopy, Acopy);

	// Set val1 = sigma2*Ω
	element_mul(val1, sigma2, ret_setup.B_Omega);
	// Compute e(sigma1, val1 = sigma2* Ω), where e is bilinear map using element_pairing function (explained in [REF-35])
	element_pairing(T1, sigma1, val1);
	// Set T1 = Value of T1 set above * current value of Acopy
	element_mul(T1, T1, Acopy); // Finally, T1 = (inverse of A)* e(sigma1, sigma2*Ω)

	// set/initialize val2 = v dash (3rd generator of cyclic bilinear group G) stored in generators array	
	element_set(val2, ret_setup.generators[2]);

	// Declare & Initialize temp_pw (which belongs to G1) using the pairing generated earlier
	element_t temp_pw;
	element_init_G1(temp_pw, ret_setup.pairing);

	// msg_val stores the current message bit in the loop below
	mpz_t msg_val;
	mpz_init(msg_val);

	// We need to multiply generators from 4th generator to (m+3)-th generator, each raised to power = message bit value. This value multiplied with v dash is stored in val2
	for (int i = 3; i < (3+mpz_get_ui(ret_setup.mval)); i++) {
		// Store the current message bit
		if (message[i-3] == '1') {
			mpz_set_ui(msg_val,1);
		} else {
			mpz_set_ui(msg_val, 0);
		}
		// Compute temp_pw = (current-generator) raised to power of current message bit value
		element_pow_mpz(temp_pw, ret_setup.generators[i], msg_val);
		// Compute val2=val2*(temp_pw) [Note that val2 was initialized with v dash (3rd generator of G)]
		element_mul(val2, val2, temp_pw);
	}

	// e - bilinear map
	// Compute e(sigma2, u) using pbc function element_pairing
	element_pairing(T2, ret_setup.generators[1],sigma2);
	// Compute e(sigma3, g) using pbc function element_pairing
	element_pairing(pairing_result2, ret_setup.generators[0], sigma3);
	
	// Compute inverse of e(sigma3, g) stored in variable pairing_result2 & store the inverse in the same variable pairing_result2
	element_invert(pairing_result2, pairing_result2);
	// Compute T2 = e(sigma2, u)*(current value of pairing_result2) = (current value of T2)*(current value of pairing_result2)
	element_mul(T2, T2, pairing_result2);

	// Compute e(sigma4, val2) & overwrite pairing_result2 with this value
	element_pairing(pairing_result2, sigma4, val2);
	
	// Compute inverse of pairing_result2's current value & store back in the same variable pairing_result2
	element_invert(pairing_result2, pairing_result2);
	// Compute T2 = (current value of T2)*(current value of pairing_result2)
	element_mul(T2, T2, pairing_result2); // Final value of T2 = [e(sigma2, u)*inverse(e(sigma3,g))*inverse(e(sigma4, original value of val2))]

	/*
		-----------TEST 3-------------
		We verify whether the signature is valid or not
	*/
	// Declare & Initialize GT group elements for usage in verification test
	element_t T1_verify, T2_verify;
	element_init_GT(T1_verify, ret_setup.pairing);
	// Compute T1_verify as e((generator of subgroup of G of order q), pi1), where e is the bilinear map
	element_pairing(T1_verify, ret_setup.gen_subgroup_q, pi1);
	element_init_GT(T2_verify, ret_setup.pairing);
	// Compute T2_verify as e((generator of subgroup of G of order q), pi2), where e is the bilinear map
	element_pairing(T2_verify, ret_setup.gen_subgroup_q, pi2);

	#ifdef DEBUG
		element_printf("%B\n%B\n",T2, T2_verify);
	#endif
	
	// Compare T1 & T1_verify and T2 & T2_verify
	if (element_cmp(T1, T1_verify) == 0 && element_cmp(T2, T2_verify) == 0) {
		// If both pairs match, verification is valid
		printf("---------------------VERIFICATION VALID-------------------------\n");
	} else {
		// Otherwise, it is invalid
		printf("VERIFICATION INVALID\n");
	}
	/*
		Verify function does the task of validating a given group signature σ on a given message
	*/
}
/*
	trace function
	It recovers the sid of the user, whose message's signature component (sigma2) is passed here, from a list of sids of all users in the group.
	It requires the use of tracing key TK
*/
void trace(element_t sigma2, mpz_t* sids, unsigned long num_of_users) {
	// Declare & Initialize elements of group G1 using the pairing generated earlier (available via ret_setup)
	element_t sigma2_copy, val;
	element_init_G1(sigma2_copy, ret_setup.pairing);
	element_init_G1(val, ret_setup.pairing);
	
	// Copy the value of sigma2 to sigma2_copy variable
	element_set(sigma2_copy, sigma2);
	// Compute sigma2^(Tracing Key) and store in sigma2_copy variable
	element_pow_mpz(sigma2_copy, sigma2_copy, ret_setup.tk);

	/*
		----------TEST 4----------
		Test to check if a given sid is the identity to be recovered or not
	*/

	// Traverse the list of all SIDs and figure out the identity to be recovered
	for (unsigned long i = 0; i < (num_of_users); i++) {
		// Compute val = g^(current sid)
		element_pow_mpz(val, ret_setup.generators[0], sids[i]);
		// Compute val = (current value of val) ^ TK = (g^(current sid))^TK
		element_pow_mpz(val, val, ret_setup.tk);

		// Check the condition that val = sigma2^(Tracing Key) OR not
		if (element_cmp(val, sigma2_copy) == 0) {
			// If condition satisfied, then the SID used is the recovered identity
			gmp_printf("----------------------------SID %Zd traced successfully--------------------------\n", sids[i]);
		} else {
			// Otherwise, the SID used is not the identity to be recovered
			gmp_printf("----------------------------SID %Zd unsuccessful in passing trace test--------------------------\n", sids[i]);
		}
	}
	/*
		Trace function, thus does the task of tracing the identity of the signer using a component of signature, the tracing key & a list of SIDs of all users (along with public values)
	*/
}
int main (int argc, char **argv) {

	num_user = 0;
	srand(time(NULL));

	/* 
	+[REF-6] security_parameter: 
		Declare the variable for storing the security parameter
	
	The following declaration declares an integer with the name security_parameter.
	*/
	mpz_t security_parameter;
	/*
	+[REF-7] void mpz_init (mpz_t x):
		Implemented by the GMP library, the function initializes the variable passed and sets its value to 0
	
	Pseudo code for self-implementation (instead of using mpz_init(mpz_t x)):
		x = 0;

	The following function call initializes the security_parameter and sets its value to 0.
	*/
	mpz_init(security_parameter);

	/*
	+[REF-8] void mpz_set_ui (mpz_t rop, unsigned long int op):
		Implemented by the GMP library, the function sets the first parameter to the value of second parameter op

	Pseudo code for self-implementation (instead of using mpz_set_ui(mpz_t rop, unsigned long int op)):
		rop = op;

	The following function call sets the value of security_paramter to the value of second paramter (currently 512).
	*/
	mpz_set_ui(security_parameter, 512);
	

	/*
	SETUP(1^lambda) from the research paper starts below.
	Here, 1^lambda, is a security parameter in unary.
	
	We pass the value of lambda = security_parameter, to the below function -setup, written on own (not pre-defined).
	ret_setup is a variable of type structure setup_result (defined near the start of the code) which stores various values returned/processed by the function - setup and used by other functions. The first parameter ret_setup is explained in detail in [REF-5]. 
	(the address of ret_setup is passed to the function to allow setup to set values in ret_setup that can later be used by other functions, when setup returns)
	
	Also, we have explained the second parameter - security_parameter.
	Detailed comments about the function (what it does, how it functions) are present in the function body.
	*/
	setup(&ret_setup, security_parameter);
	/*
		SETUP function ends
	*/

	/*
		Enroll(PP, MK, ID) from the research paper starts below
		PP is the public values generated by setup. MK- Master Enrollment Key
		ID- User ID,  0 ≤ ID < 2^k < p

		Enroll creates a signing key for user ID - ID.

	*/
	/*	
		We generate some random userIDs
	*/
	// declare, initialize variables required to achieve above tasks.
	mpz_t random, userID;
	mpz_init(random);
	mpz_init(userID);

	// Declare state 
	gmp_randstate_t state; // gmp_randstate_t explained in detail in [REF-11]
	// Initialize state
	gmp_randinit_default(state); // gmp_randinit_default explained in detail in [REF-12]
	// Set initial seed value to state. We pass seed as a random number defined as (random_number + 1)*(another-random-number) + 1. 1 was added to avoid seed = 0.
	gmp_randseed_ui(state, (rand()+1)*(rand()+1)); // gmp_randseed_ui explained in detail in [REF-13]

	// Below, we compute the UPPER LIMIT of the user ID , which is 2^k
	mpz_t count, limit;
	mpz_init(count);
	mpz_init(limit);
	// start upper limit with 1
	mpz_set_ui(limit, 1);

	// till k iterations continue
	while(mpz_cmp(count, security_parameter)) {
		/*
		+[REF-36] void mpz_mul_ui(mpz_t rop, mpz_t a, unsigned long int x)
			Provided by gmp library, used for setting rop = a*x
		*/
		mpz_mul_ui(limit, limit, 2); // we do limit = limit*2
		mpz_add_ui(count, count, 1); // add 1 to count so that loop runs for k times. (we have used k = security_parameter/lambda as only k = O(security_parameter) is the required condition)
	}

	// Below, we generate 10 random user IDs to demonstrate the functionality
	unsigned long num_of_users = 10; // num of users are 10

	// Prepare an array (of size num_of_users) for storing secret unique value SID that will be generated by the function enroll
	mpz_t* sids = (mpz_t*)malloc(sizeof(mpz_t)*num_of_users);


	for (unsigned long i = 0; i < num_of_users;i++) {
		element_t k1, k2, k3;

		// Initialize all user IDs
		mpz_init(sids[i]);
		do {
			// Generate a random userID within the UPPER LIMIT stored in variable limit
			mpz_urandomm(userID, state, limit);	
			/* 
			Loop till the userID generated is NOT ZERO and userID generated is not equal to the value of p or q. mpz_cmp_ui is explained in [REF-14], mpz_cmp is explained in [REF-17]
			*/
		} while(mpz_cmp_ui(userID, 0) == 0 || mpz_cmp(userID, ret_setup.pval) == 0 || mpz_cmp(userID, ret_setup.qval) == 0);

		/*
		The function enroll enrolls the user with user ID and generates the signing key in the form of three values - k1, k2, k3. The function also assigns a secret unique value SID
		*/
		enroll(sids[i], userID, k1, k2, k3);
	}

	// WE NOW PICK ONE USER ID THAT WILL BE USED TO TRACE A USER LATER IN THE PROGRAM.
	
	// k1, k2, k3 that will be generated by enroll function for the picked user ID
	element_t k1, k2, k3;

	// Print the User ID that will be used for tracing & demonstration.
	gmp_printf("------------------Demonstrating for USERID - %Zd-----------------\n", userID);
	
	mpz_t sid;
	mpz_init(sid);

	/* 
		Pass userID. 
		Other parameters passed will be assigned suitable values by the enroll function.
		Jump to enroll function body for details on it.
	*/
	enroll(sid, userID, k1, k2, k3);
	/*
		enroll function has assigned sid to user with user ID = userID
		It has also assigned suitable values to k1, k2, k3
	*/
	
	// Print the elements k1, k2, k3 (signing key) given to the user with the user ID picked (for demonstration) - userID
	element_printf("K1: %B\nK2: %B\nK3: %B\n", k1, k2, k3);

	/*
		--------------TEST 1-------------
		Below, we check that the signing key returned by enroll is WELL-FORMED or NOT
	*/
	// Declare and initialize variables with pairing generated earlier in the program
	// Below, val4 belongs to group G1, val1,val2, val3 belongs to group GT
	element_t val1, val2, val3, val4;
	element_init_G1(val4, ret_setup.pairing);
	element_init_GT(val1, ret_setup.pairing);
	element_init_GT(val2, ret_setup.pairing);
	element_init_GT(val3, ret_setup.pairing);

	// Compute val1 = e(k2, u), where e is the bilinear map according to research paper
	// Bilinear map can be called by using the function element_pairing (EXPLAINED in [REF-35])
	element_pairing(val1, k2, ret_setup.generators[1]);

	// Compute val2 = e(k3, g) where e is the bilinear map
	element_pairing(val2, k3, ret_setup.generators[0]);
	
	// Compute val4 = k2*Ω. Multiplication is done using pre-defined function element_mul
	/*
	+[REF-39] void element_mul(element_t rop, element_t p, element_t q)
		Provided by PBC library. Sets rop to p*q.
		The function is used for multiplying elements of groups etc.
	*/
	element_mul(val4, k2, ret_setup.B_Omega);

	// Compute val3 = e(k1, val4 = k2*Ω) where e is the bilinear map
	element_pairing(val3, k1, val4);

	// FOR DEBUG purpose while writing/testing code, avoid reading
	#ifdef DEBUG
		if (element_cmp(val1, val2) == 0) {
			printf("SUCCESS - 1\n");
		} 
		if (element_cmp(val3, ret_setup.Aval) == 0) {
			printf("SUCCESS - 2\n");
		}
	#endif

	/* 
		VERIFY WHETHER 
		1. val1 = e(k2, u) is equal to val2 = e(k3, g)
			AND
		2. val3 = e(k1, k2*Ω) is equal to A computed by setup function
	*/
	if (element_cmp(val1, val2) == 0 && element_cmp(val3, ret_setup.Aval) == 0) {
		// VERIFICATION SUCCEEDS
		printf("--------------------Verification of key successful - key well formed by enroll()----------------------------\n");
	} else {
		// VERIFICATION FAILS
		printf("Verification of key UNSUCCESSFUL - key not well formed by enroll()\n");
	}
	/*
		ENROLL ends
	*/
	/*
		SIGN(PP, ID, KID, M)  STARTS
		where PP = public values generated by setup (availabe in ret_setup)
		ID = userID
		KID = signing key (k1, k2, k3)
		M = message to be signed
	*/
	// Create an array for storing message of size m = O(security_parameter)
	// Message is a binary message with number of bits = m
	char* msg = (char*)malloc(sizeof(char)*(1+mpz_get_ui(ret_setup.mval)));

	// Generate bits in the message randomly
	for (unsigned int i = 0; i < mpz_get_ui(ret_setup.mval); i++) {
		// rand function provided by C standard library produces a random number
		int b = rand()%2;
		// If rand()%2 == 0 message bit is 0, otherwise it is 1
		msg[i] = '1';
		if (b == 0) {
			msg[i] = '0';
		}
	}
	// Set last character of message to '\0' (terminating character of string in C)
	msg[mpz_get_ui(ret_setup.mval)] = '\0'; //mpz_get_ui() is explained in [REF-30]
	
	// Print message to console
	printf("Message (m bits): %s\n", msg);

	// Declare elements of groups to pass to sign function
	element_t sigma1, sigma2, sigma3, sigma4, pi1, pi2;
	/*
		Sign requires as input - PP (available using ret_setup), k1, k2, k3, message and sets pi1, pi2, sigma1, sigma2, sigma3, sigma4 which forms the signature for message msg and user (k1, k2, k3 uniquely represents the user).
		Please jump to function body for further details.
	*/
	sign(pi1, pi2,sigma1, sigma2, sigma3, sigma4, k1,k2,k3, msg);

	/*
		Print the final signature to the console
	*/
	element_printf("Final Signature:\nSigma1: %B\nSigma2: %B\nSigma3: %B\nSigma4: %B\nPi1: %B\nPi2: %B\n",sigma1, sigma2, sigma3, sigma4, pi1, pi2);
	/*
		SIGN ENDS
	*/

	/*
		VERIFY(PP, M, σ) STARTS
		where PP = public values computed by setup & stored in ret_setup, M = message & σ computed by sign function
		Use: To validate a group signature σ on a message M
	*/
	// To verify function,
	// Pass required information out of all public values - A
	// Pass signature σ = (sigma1, sigma2, sigma3, sigma4, pi1, pi2)
	// Pass message msg
	// For more explanation on verify function, please jump to the function body
	verify(ret_setup.Aval, sigma1, sigma2, sigma3, sigma4, pi1, pi2, msg);
	/*
		VERIFY ENDS
	*/


	/*
		TRACE(PP, TK, signature) STARTS, where PP = Public values computed by setup & stored in ret_setup
		TK = group manager’s tracing key (available to TRACE function using ret_setup)
		signature = signature assumed to pass the verification test for some message M, which is not needed here
		The aim is to recover the identity of the signer.
		The tracer outputs the recovered identity of the signer.
	*/
	// Only sigma2 out of the complete signature is required by the trace function
	// We pass all sids and expect the trace function to recover the sid of the user whose message M's signature is passed to the function. We also pass num_of_users to help traverse array of user-IDs.
	// We recall that the use of SID was stated to be tracing purpose only.
	// For more details on trace function, please jump to the function body.
	trace(sigma2, sids, num_of_users);
	/*
		TRACE ENDS
	*/

	return 0;

}
/*
REFERENCES:
- GMP Manual. [Online]. Available: https://gmplib.org/gmp-man-6.0.0a.pdf (accessed May, 2020)
- PBC Library Manual. [Online]. Available: https://crypto.stanford.edu/pbc/manual.html (accessed May, 2020)
*/
/*
	index_of_references
REFERENCE NUMBER	EXPLANATION OF [Function names end with ()]
REF-1					element_t
REF-2					mpz_t
REF-3					pairing_t
REF-4					pbc_param_t
REF-5					ret_setup
REF-6					security_parameter
REF-7					mpz_init()
REF-8					mpz_set_ui()
REF-9					mpz_set()
REF-10					mpz_add_ui()
REF-11					gmp_randstate_t
REF-12					gmp_randinit_default()
REF-13					gmp_randseed_ui()
REF-14					mpz_cmp_ui()
REF-15					mpz_ui_pow_ui()
REF-16					mpz_urandomb()
REF-17					mpz_cmp()
REF-18					mpz_probab_prime_p()
REF-19					mpz_mul()
REF-20					gmp_printf()
REF-21					pbc_param_init_a1_gen()
REF-22					pairing_init_pbc_param()
REF-23					element_init_G1()
REF-24					element_init_GT()
REF-25					element_set0()
REF-26					element_set()
REF-27					element_random()
REF-28					element_pow_mpz()					
REF-29					element_printf()
REF-30					mpz_get_ui()
REF-31					mpz_urandomm()
REF-32					mpz_gcd()
REF-33					element_add()
REF-34					pbc_param_out_str()
REF-35					element_pairing()
REF-36					mpz_mul_ui()
REF-37					mpz_add()
REF-38					mpz_invert()
REF-39					element_mul()
REF-40					mpz_set_ui()
REF-41					element_invert()
*/