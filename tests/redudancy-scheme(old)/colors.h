#ifndef _COLORS_H_
#define _COLORS_H_

/* FOREGROUND */
#define RST  "\x1B[0m"
#define CRED  "\x1B[31m"
#define CGRN  "\x1B[32m"
#define CYEL  "\x1B[33m"
#define CBLU  "\x1B[34m"
#define CMAG  "\x1B[35m"
#define CCYN  "\x1B[36m"
#define CWHT  "\x1B[37m"

#define RED(x) CRED x RST
#define GRN(x) CGRN x RST
#define YEL(x) CYEL x RST
#define BLU(x) CBLU x RST
#define MAG(x) CMAG x RST
#define CYN(x) CCYN x RST
#define WHT(x) CWHT x RST

#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

#endif  /* _COLORS_H_ */

