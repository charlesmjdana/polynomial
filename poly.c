/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   poly.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: charles <charles.dana@hec.edu>             +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2023/07/15 11:09:45 by charles           #+#    #+#             */
/*   Updated: 2023/07/16 16:30:38 by charles          ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

double complex polynomial[10];
int	deg;

void show(double complex z) {
	printf("(%+.2f + %+.2fi)", creal(z), cimag(z));
}

double complex pot(int n) {
	double complex z;

	z = (double complex)1;
	while (n > 0) {
		z = z * 2;
		n--;
	}
	while (n < 0) {
		z = z / 2;
		n++;
	}
	return z;
}

double complex pos(int n, double complex s) {
	double complex z;

	z = (double complex)1;
	while (n > 0) {
		z = z * s;
		n--;
	}
	while (n < 0 && s != (double complex)0) {
		z = z / s;
		n++;
	}
	return z;
}

double complex construct(long x) {
	double complex	z;
	int	i;

	z = (double complex)0;
	i = 1;
	while (i < 21 && x > 0) {
		if (x % 2 == 1) {
			if (i % 4 == 1) {
				z = z + pot((i - 1) / 4);
			}
			if (i % 4 == 2) {
				z = z + pot(-1 + (2 - i) / 4);
			}
			if (i % 4 == 3) {
				z = z + I * pot((i - 3) / 4);
			}
			if (i % 4 == 0) {
				z = z + I * pot(-i / 4);
			}
		}
		x = (x - x % 2) / 2;
		i++;
	}
	if (x > 0) {
		if (x % 2 == 1)
			z = -z;
		x = (x - x % 2) / 2;
	}
	if (x > 0) {
		if (x % 2 == 1) {
			z = creal(z) - I * cimag(z);
		}
	}
	return z;
}

double complex evaluate(long x) {
	double complex	z;
	double complex	f;
	int	i;

	z = (double complex)0;
	i = 1;
	while (i < 21 && x > 0) {
		if (x % 2 == 1) {
			if (i % 4 == 1) {
				z = z + pot((i - 1) / 4);
			}
			if (i % 4 == 2) {
				z = z + pot(-1 + (2 - i) / 4);
			}
			if (i % 4 == 3) {
				z = z + I * pot((i - 3) / 4);
			}
			if (i % 4 == 0) {
				z = z + I * pot(-i / 4);
			}
		}
		x = (x - x % 2) / 2;
		i++;
	}
	if (x > 0) {
		if (x % 2 == 1)
			z = -z;
		x = (x - x % 2) / 2;
	}
	if (x > 0) {
		if (x % 2 == 1) {
			z = creal(z) - I * cimag(z);
		}
	}
	i = 1;
	f = polynomial[0];
	while (i <= deg) {
		f = f + pos(i, z) * polynomial[i];
		i++;
	}
	return f;
}

double complex evaluatefromC(double complex z) {
	int i;
	double complex f;

	i = 1;
	f = polynomial[0];
	while (i <= deg) {
		f = f + pos(i, z) * polynomial[i];
		i++;
	}
	return f;
}

long Ps[100000];
long Qs[100000];
int	Zs[100000];
int	pnqindex;

long lpot(int n) {
	long z;
	
	z = 1;
	while (n > 0) {
		z = z * 2;
		n--;
	}
	return z;
}

int	compute(long p, long q) {
	int i;
	int z;
	int buffer[24];
	int b;
	long tmp1;
	long tmp2;
	long tmp3;
	long magnitude;

	if (p == q)
		return 1;
	if (p < q)
		return 0;
	i = 0;
	while (i < pnqindex && (Ps[i] != p || Qs[i] != q))
		i++;
	if (i < pnqindex)
		return Zs[i];
	b = 0;
	i = 1;
	tmp1 = p;
	tmp2 = q;
	while (i < 23 && tmp1 > 0) {
		if (tmp1 % 2 == 1 && tmp2 % 2 == 0)
			buffer[b++] = i;
		if (tmp1 % 2 == 0 && tmp2 % 2 == 1)
			return 0;
		i++;
		tmp1 = (tmp1 - tmp1 % 2) / 2;
		tmp2 = (tmp2 - tmp2 % 2) / 2;
	}
	if (b == 0)
		return 1;
	magnitude = lpot(b);
	z = 0;
	tmp1 = 0;
	while (tmp1 < magnitude - 1) {
		tmp2 = q;
		i = 0;
		tmp3 = tmp1;
		while (i < b && tmp3 > 0) {
			if (tmp3 % 2 == 1)
				tmp2 = tmp2 + lpot(buffer[i] - 1);
			tmp3 = (tmp3 - tmp3 % 2) / 2;
			i++;
		}
		z = z + compute(tmp2, q);
		tmp1++;
	}
	if (pnqindex < 99999) {
		Ps[pnqindex] = p;
		Qs[pnqindex] = q;
		Zs[pnqindex] = -z;
		pnqindex++;
	}
	//printf("p: %ld q: %ld magnitude: %ld z:%d\n", p, q, magnitude, -z);
	return (-z);
}

double approximate(long p) {
	int	buffer[100];
	int b;
	int i;
	long tmp1;
	long tmp2;
	long tmp3;
	long q;
	double complex d;
	double complex e;
	
	if (p == 0) {
		d = evaluate(0);
		return -sqrt(creal(d) * creal(d) + cimag(d) * cimag(d));
	}
	tmp1 = p;
	i = 1;
	b = 0;
	while (i < 23 && tmp1 > 0) {
		if (tmp1 % 2 == 1)
			buffer[b++] = i;
		i++;
		tmp1 = (tmp1 - tmp1 % 2) / 2;
	}
	tmp1 = 1;
	i = 0;
	while (i < b) {
		tmp1 = tmp1 * 2;
		i++;
	}
	d = 0;
	tmp2 = 0;
	while (tmp2 < tmp1) {
		q = 0;
		tmp3 = tmp2;
		i = 0;
		while (tmp3 > 0 && i < 23) {
			if (tmp3 % 2 == 1)
				q = q + lpot(buffer[i] - 1);
			tmp3 = (tmp3 - tmp3 % 2) / 2;
			i++;
		}
		if (compute(p, q) != 0)
			e = evaluate(q);
			d = d - compute(p, q) * sqrt(creal(e) * creal(e) + cimag(e) * cimag(e));
		tmp2++;
	}
	return d;
}

long history[100000];
int historyindex;

long generate(int pr) {
	long g;
	int i;

	g = 0;
	while (pr > 0 && rand() % pr > 0) {
		g = (g + 1) % 8388608;
	}
	i = 0;
	while (i < historyindex && history[i] != g)
		i++;
	if (i == historyindex && historyindex < 99999)
		history[historyindex++] = g;
	else
		return generate(pr);
	return g;
}

char* problem[100000];
double weight[100000];
int	problemindex;

int expand(int pr) {
	long g;
	double a;
	int	buffer[100];
	int b;
	int i;
	int j;
	char *clause;

	g = generate(pr);
	a = approximate(g);
	b = 0;
	i = 1;
	while (g > 0) {
		if (g % 2 == 1)
			buffer[b++] = i;
		i++;
		g = (g - g % 2) / 2;
	}
	if (a > 0) {
		clause = malloc(sizeof(char) * (b + 1));
		clause[b] = '\0';
		i = 0;
		while (i < b) {
			clause[i] = -buffer[i];
			i++;
		}
		problem[problemindex] = clause;
		weight[problemindex++] = -a;
		return 1;
	}
	if (a < 0) {
		i = 0;
		while (i < b) {
			clause = malloc(sizeof(char) * (b + 1 - i));
			clause[0] = buffer[i];
			j = 1;
			while (i + j < b) {
				clause[j] = -buffer[i + j];
				j++;
			}
			clause[j] = '\0';
			problem[problemindex] = clause;
			weight[problemindex++] = a;
			i++;
		}
		return b;
	}
	return 0;
}

int shuffle(int permutations) {
	int i;
	int j;
	char *tmpc;
	double tmpw;

	if (problemindex <= 0)
		return 0;
	while (permutations > 0) {
		i = rand() % problemindex;
		j = rand() % problemindex;
		tmpc = problem[i];
		tmpw = weight[i];
		problem[i] = problem[j];
		weight[i] = weight[j];
		problem[j] = tmpc;
		weight[j] = tmpw;
		permutations--;
	}
	return 1;
}

double proxy(long l) {
	double complex z;

	z = evaluate(l);
	return sqrt(creal(z) * creal(z) + cimag(z) * cimag(z));
}

double asymptotic(long l) {
	double complex z;
	
	z = construct(l);
	z = pos(deg, z);
	return sqrt(creal(z) * creal(z) + cimag(z) * cimag(z));
}

long solve(int iterations) {
	int i;
	int j;
	int t;
	int r;
	char solution[23];
	int poss[23];
	int possindex;
	double score;
	long best;
	long candidate;

	best = 0;
	solution[0] = ' ';
	i = 1;
	while (i < 23)
		solution[i++] = '0';
	solution[i] = '\0';
	while (iterations > 0) {
		expand(130);
		shuffle(problemindex);
		i = 0;
		score = 0.0;
		t = 0;
		while (i < problemindex) {
			j = 0;
			possindex = 0;
			while (problem[i][j]) {
				if (problem[i][j] > 0 && solution[(int)problem[i][j]] == '1')
					break;
				if (problem[i][j] < 0 && solution[(int)problem[i][j]] == '2')
					break;
				if (solution[abs((int)problem[i][j])] == '0')
					poss[possindex++] = problem[i][j];
				j++;
			}
			if (problem[i][j]) {
				score = (t * score + weight[i]) / ((double)(t + 1)); // clause secured pos weight
			} else {
				if (possindex > 0) {
					score = (t * score + weight[i]) / ((double)(t + 1)); // clause secured pos weight
					r = rand() % possindex;
					if (poss[r] > 0) {
						solution[r] = '1';
					} else {
						solution[-r] = '2';
					}
				} else if (weight[i] > score) {
					r = rand() % j;
					solution[abs((int)problem[i][r])] = 0;
					score = (t * score) / ((double)(t + 1)); // clause possible weight 0
				} else {
					score = (t * score - weight[i]) / ((double)(t + 1)); // skip clause neg weight
				}
			}
			t++;
			i++;
		}
		//printf("score: %.2f\n", score);
		t = 0;
		while (t < 10000) {
			candidate = 0;
			i = 1;
			while (i < 23) {
				if (solution[i] == '0') {
					candidate = candidate + (rand() % 2) * lpot(i - 1);
				}
				if (solution[i] == '1') {
					candidate = candidate + lpot(i - 1);
				}
				i++;
			}
			if (proxy(best) > proxy(candidate)) {
				//printf("current best: %ld score: %.2f root: ", candidate, proxy(candidate));
				//show(construct(candidate));
				//printf("\n");
				best = candidate;
			}
			t++;
		}
		iterations = iterations - 1;
		if (proxy(best) < 0.1 * asymptotic(best) || proxy(best) < 0.05)
			break;
	}
	return best;
}

int fact(int n) {
	if (n == 0)
		return 1;
	return n * fact(n-1);
}

int reduce(double complex root) {
	int i;
	int j;
	int k;
	double complex newpoly[100];
	double complex eps = 0.01;
	double complex tmp[100];
	double complex tmp2[100];

	printf("(X - ");
	show(root);
	printf(")");
	if (deg <= 1)
		printf("\n");
	if (root == (double complex)0) {
		i = 1;
		while (i <= deg) {
			polynomial[i - 1] = polynomial[i];
			i++;
		}
		deg--;
	} else {
		newpoly[0] = -polynomial[0] / root;
		i = 1;
		while (i <= deg) {
			j = 0;
			while (j <= i) {
				tmp[j] = evaluatefromC(eps * j) / (eps * j - root);
				j++;
			}
			while (j > 1) {
				k = 0;
				while (k < j) {
					tmp2[k] = (tmp[k + 1] - tmp[k]) / eps;
					k++;
				}
				k = 0;
				while (k < j) {
					tmp[k] = tmp2[k];
					k++;
				}
				j--;
			}
			newpoly[i] = tmp[0] / fact(i);
			i++;
		}
		i = 0;
		while (i < deg) {
			polynomial[i] = newpoly[i];
			i++;
		}
		deg--;
	}
	historyindex = 1;
	problemindex = 0;
	expand(130);
	expand(130);
	expand(130);
	expand(130);
	expand(130);
	expand(130);
	expand(130);
	return deg;
}

double proxyfromC(double complex root) {
	double complex f;

	f = evaluatefromC(root);
	return sqrt(creal(f) * creal(f) + cimag(f) * cimag(f));
}

double complex improve(double complex root) {
	int	T;
	double score;
	double complex best;
	double complex candidate;
	double eps = 0.001;
	int	i;
	int j;
	
	best = root;
	score = proxyfromC(root);
	T = 0;
	i = -25;
	while (i < 25) {
		j = -25;
		while (j < 25) {
			candidate = creal(root) * (1.0 + i * eps) + I * cimag(root) * (1.0 + j * eps);
			if (proxyfromC(candidate) < score) {
				T = 1;
				score = proxyfromC(candidate);
				best = candidate;
			}
			j++;
		}
		i++;
	}
	if (T == 1)
		return improve(best);
	return best;
}

int main(int argc, char **argv) {
	int i;
	
	i = 1;
	deg = 0;
	polynomial[0] = (double complex)0;
	while (argc > i + 1) {
		polynomial[deg++] = (double)atoi(argv[i]) + I * ((double)atoi(argv[i + 1]));
		polynomial[deg] = (double complex)0;
		i = i + 2;
	}
	deg--;
	i = 0;
	while (i <= deg) {
		polynomial[i] = polynomial[i] / polynomial[deg];
		i++;
	}
	if (deg <= 0) {
		printf("P(X) = ");
		show(polynomial[0]);
		printf("\n");
		return EXIT_SUCCESS;
	}
	printf("P(X) = ");
	i = 0;
	while (i < deg) {
		show(polynomial[i]);
		if (i > 0) {
			printf(" * X^%d", i);
		}
		if (i < deg)
			printf(" + ");
		i++;
	}
	printf("X^%d\n", deg);
	srand(time(NULL));
	pnqindex = 0;
	historyindex = 1;
	history[0] = 0;
	problemindex = 0;
	expand(130);
	expand(130);
	printf("P(X) ~ ");
	while (deg > 0) {
		deg = reduce(improve(construct(solve(100))));
	}
	return EXIT_SUCCESS;
}
