#if BPI == 16
#define ODD_MASK  0xaaaa
#define EVEN_MASK 0x5555
#else
#define ODD_MASK  0xaaaaaaaa
#define EVEN_MASK 0x55555555
#endif

#define POSITIVE 1
#define NEGATIVE 0

#define PRESENT 1
#define ABSENT  0

#define RAISED 2

/* black_white.c */ void setup_bw(pcover R, pcube c);
/* black_white.c */ void free_bw();
/* black_white.c */ bool black_white();
/* black_white.c */ void split_list(pcover R, int v);
/* black_white.c */ void merge_list();
/* black_white.c */ void print_bw(int size);
/* black_white.c */ void variable_list_alloc(int size);
/* black_white.c */ void variable_list_init(int reduced_c_free_count, int *reduced_c_free_list);
/* black_white.c */ void variable_list_delete(int element);
/* black_white.c */ void variable_list_insert(int element);
/* black_white.c */ bool variable_list_empty();
/* black_white.c */ void get_next_variable(int *pv, int *pphase, pcover R);
/* black_white.c */ void print_variable_list();
/* black_white.c */ extern void reset_black_list();
/* black_white.c */ extern void push_black_list();
/* black_white.c */ extern void pop_black_list();
/* canonical.c */ int is_minterm();
/* canonical.c */ pcover find_canonical_cover(pcover F1, pcover D, pcover R);
/* essentiality.c */ pcover etr_order(pcover F, pcover E, pcover R, pcube c, pcube d);
/* essentiality.c */ void aux_etr_order(pcover F, pcover E, pcover R, pcube c, pcube d);
/* essentiality.c */ pcover get_mins(pcube c);
/* essentiality.c */ int ascending(const void *p1, const void *p2);
/* util_signature.c */ void set_time_limit(int seconds);
/* util_signature.c */ void print_cover(pcover F, char *name);
/* util_signature.c */ int sf_equal(pcover F1, pcover F2);
/* util_signature.c */ int time_usage(char *name);
/* util_signature.c */ void s_totals(long time, int i);
/* util_signature.c */ void s_runtime(long total);
/* sigma.c */ pcube get_sigma(pcover R, register pcube c);
/* sigma.c */ void set_not(pcube c);
/* signature.c */ void cleanup(int sig);
/* signature.c */ pcover signature(pcover F1, pcover D1, pcover R1);
/* signature.c */ pcover generate_primes(pcover F, pcover R);
/* signature_exact.c */ pcover signature_minimize_exact(pcover ESCubes, pcover ESSet);
/* signature_exact.c */ sm_matrix *signature_form_table(pcover ESCubes, pcover ESSet);
