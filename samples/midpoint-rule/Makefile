all: midpoint-rule

midpoint-rule: midpoint_rule.c
	mpicc -g -Wall -o midpoint_rule midpoint_rule.c -lm -std=c99

clean:
	rm -f midpoint_rule
