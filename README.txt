Thanks for looking at this project.

This is a project I started a few years ago, during the end of 2019 and beginning of 2020.
This is still a work in progress, as I have reorganized the initial project.
It is a working version, you can run the algorithm by running the .exe file in the build folder.
See more in the "running the algorithm" section.


##### Introduction ######

This project is a C implementation of the Quadratic Sieve to factor integers.

Here are some detailed characteristics :
- Small multiplier used
- Multiple polynomials used for sieving
- Single large prime variation used to collect relations
- Batch smoothness and trial division tests available
- 1 CPU or all CPU sieving possible
- Wiedemann algorithm used for the linear algebra step

##### Sources #####

Overall algorithm:
- "Prime numbers, a computational perspective" by Richard Crandall and Carl Pomerance (really good): https://link.springer.com/book/10.1007/0-387-28979-8
- "The Quadratic Sieve Factoring Algorithm" by Eric Landquist: https://www.cs.virginia.edu/crab/QFS_Simple.pdf

Multiple polynomials and double large primes
- "Factoring Integers with Large-Prime Variations of the Quadratic Sieve" by Henk Boender and Herman J. J. te Riele: https://projecteuclid.org/journals/experimental-mathematics/volume-5/issue-4/Factoring-integers-with-large-prime-variations-of-the-quadratic-sieve/em/1047565445.pdf

Batch smoothness test:
- "HOW TO FIND SMOOTH PARTS OF INTEGERS" by DANIEL J. BERNSTEIN: https://cr.yp.to/factorization/smoothparts-20040510.pdf

Gaussian elimination:
- Took the gaussian elimination code from https://github.com/skollmann/PyFactorise/blob/master/factorise.py#L39

Block Lanczos:
- "A Block Lanczos Algorithm for Finding Dependencies over GF(2)" by Peter L. Montgomery: https://scispace.com/pdf/a-block-lanczos-algorithm-for-finding-dependencies-over-gf-2-ezdu2qt0pp.pdf
- "A modified block Lanczos algorithm with fewer vectors" by Emmanuel Thomé (really good): https://eprint.iacr.org/2016/329.pdf

Wiedemann algorithm:
- "SOLVING HOMOGENEOUS LINEAR EQUATIONSOVER GF(2) VIA BLOCK WIEDEMANN ALGORITHM" by Don Coppersmith: https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1192970-7/S0025-5718-1994-1192970-7.pdf

##### running the algortihm #####

You just need to run the build/MPQS_v_1_0_0.exe program in a terminal.
The program will ask you to input the number you want to factor.

##### config parameters #####

Here are the following parameters you can edit in the config file:

- flag_parallel_sieve : If 0, only one CPU is used for sieving. Else, all the available CPUs are used for sieving.

- flag_batch_smooth : If 0, the trial division (naive) test is used. Else, the batch smoothness test is used. 

##### General discussion #####

I tried overall to use as few exernal libraries as possible. While this makes the code 
	way harder to write for me, and understand for you, I viewed it as a good challenge
	to write good code, and get to know hard and foreign concepts to me. It took a lot
	of work, but I got a bug free code, and got much more knowledgable on the subroutines
	thanks to this. Furthermore, this project was a challenge for me to become proficient in C.



Here are some points I consider working on at some point:

- No matter how hard I tried, I am stuck on understanding the block Wiedemann algorithm.
	For now, the best I can do is the scalar one, with some optimizations. Namely, I
	use binary encoding of the blocks of vectors to compute matrix-vector product very
	efficiently. Then, for each scalar sequence, I compute its generator using scalar
	Berlekamp-Massey algorithm, and use it to update the minimal polynomial of the matrix.
	I consider working some way to compute the matrix generator of the matrix sequence
	through some algorithm (matrix Berlekamp-Massey or some other).

- I have to make the linear algebra step more parallel. Even though it is the fastest of the
	two main steps of the algorithm, no one likes to wait if it is possible to go faster.

Here are the next steps:

- I still have to implement many features : choose how many cpu are used for sieving, adding
	gaussian elimination and block Lanczos algorithms for the linear algebra step.

- I have to continue cleaning up the code.

- Finally, implement the GNFS algorithm in C.