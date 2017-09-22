all: midpoint-rule midpoint-rule-ulfm

midpoint-rule-ulfm: midpoint_rule.c
	mpicc -g -Wall -o midpoint_rule midpoint_rule.c -lm -std=c99 -D ULFM_SUPPORT

midpoint-rule: midpoint_rule.c
	mpicc -g -Wall -o midpoint_rule_ulfm midpoint_rule.c -lm -std=c99

clean:
	rm -f midpoint_rule midpoint_rule_ulfm
