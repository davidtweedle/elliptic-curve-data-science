/*-*- compile-command: "cc -c -o generate-data.gp.o -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -I\"/usr/include/x86_64-linux-gnu\" generate-data.gp.c && cc -o generate-data.gp.so -shared -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -Wl,-shared -Wl,-z,relro generate-data.gp.o -lc -lm -L/usr/lib/x86_64-linux-gnu -lpari"; -*-*/
#include <pari/pari.h>
/*
GP;install("init_generate_data","v","init_generate_data","./generate-data.gp.so");
GP;install("write_aps","vGD0,G,D0,G,sp","write_aps","./generate-data.gp.so");
GP;install("construct_box_curves","vGGGsp","construct_box_curves","./generate-data.gp.so");
GP;install("construct_twists_sqrt_n3","vGGsp","construct_twists_sqrt_n3","./generate-data.gp.so");
GP;install("construct_twists_sqrt_n1","vD0,G,D0,G,D0,G,p","construct_twists_sqrt_n1","./generate-data.gp.so");
GP;install("rescale_aps_sqrt_n1_step","vD0,G,D0,G,D0,G,D0,G,D0,G,D0,G,D0,G,p","rescale_aps_sqrt_n1_step","./generate-data.gp.so");
GP;install("rescale_aps_sqrt_n1","vD0,G,D0,G,D0,G,D0,G,D0,G,D0,G,p","rescale_aps_sqrt_n1","./generate-data.gp.so");
GP;install("construct_twists_sqrt_N","vGGGsp","construct_twists_sqrt_N","./generate-data.gp.so");
GP;install("construct_all_twists","vD0,G,GGGsp","construct_all_twists","./generate-data.gp.so");
GP;install("nth_power_free","GG","nth_power_free","./generate-data.gp.so");
*/
void init_generate_data(void);
void write_aps(GEN P_bd, GEN curves, GEN labels, char *filepath, long prec);
void construct_box_curves(GEN P_bd, GEN A_bd, GEN B_bd, char *filepath, long prec);
void construct_twists_sqrt_n3(GEN P_bd, GEN D_bd, char *filepath, long prec);
void construct_twists_sqrt_n1(GEN P_bd, GEN D_bd, GEN filepath, long prec);
void rescale_aps_sqrt_n1_step(GEN P_bd, GEN start, GEN stepsize, GEN numsteps, GEN N, GEN n, GEN fileprefix, long prec);
void rescale_aps_sqrt_n1(GEN P, GEN D1, GEN D2, GEN N, GEN n, GEN filepath, long prec);
void construct_twists_sqrt_N(GEN P_bd, GEN D_bd, GEN N, char *filepath, long prec);
void construct_all_twists(GEN P, GEN D_bd, GEN A, GEN B, char *filepath, long prec);
GEN nth_power_free(GEN s, GEN n);
/*End of prototype*/

void
init_generate_data(void)	  /* void */
{
  pari_sp ltop = avma;
  avma = ltop;
  return;
}

/* for each E in curves write a row to filepath
* E.label E.conductor E.rank E.ap(2) E.ap(3) E.ap(5) ... E.ap(pn)
* where 2,3,...,pn are the primes between 1 and P_bd
* write the result to the file at filepath
* we assume that the curves and labels have been constructed
* from the LMFDB database
*/
void
write_aps(GEN P_bd, GEN curves, GEN labels, char *filepath, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN P = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN n = gen_0, m = gen_0, file = gen_0;
  if (typ(P_bd) != t_INT)
    pari_err_TYPE("write_aps",P_bd);
  p1 = cgetg(3, t_VEC);
  gel(p1, 1) = gen_1;
  gel(p1, 2) = icopy(P_bd);
  P = primes0(p1);
  n = stoi(glength(P));
  m = stoi(glength(curves));
  file = stoi(gp_fileopen(filepath, "w"));
  if (glength(labels) == 0)
  {
    long l2;
    GEN p3 = gen_0;	  /* vec */
    l2 = glength(curves);
    {
      long i;
      p3 = cgetg(l2+1, t_VEC);
      for (i = 1; i <= l2; ++i)
        gel(p3, i) = stoi(i);
    }
    labels = p3;
  }
  gp_filewrite1(gtos(file), "label discriminant conductor rank");
  {
    pari_sp btop = avma;
    GEN j = gen_0;
    for (j = gen_1; gcmp(j, n) <= 0; j = gaddgs(j, 1))
    {
      gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), gel(P, gtos(j)))));
      if (gc_needed(btop, 1))
        j = gerepilecopy(btop, j);
    }
  }
  gp_filewrite1(gtos(file), "\n");
  {
    pari_sp btop = avma;
    GEN i = gen_0;
    for (i = gen_1; gcmp(i, m) <= 0; i = gaddgs(i, 1))
    {
      {
        GEN p4 = gen_0;	  /* vec */
        p4 = cgetg(8, t_VEC);
        gel(p4, 1) = gcopy(gel(labels, gtos(i)));
        gel(p4, 2) = strtoGENstr(" ");
        gel(p4, 3) = gcopy(member_disc(gel(curves, gtos(i))));
        gel(p4, 4) = strtoGENstr(" ");
        gel(p4, 5) = gcopy(gel(ellglobalred(gel(curves, gtos(i))), 1));
        gel(p4, 6) = strtoGENstr(" ");
        gel(p4, 7) = gcopy(gel(ellanalyticrank_bitprec(gel(curves, gtos(i)), NULL, prec2nbits(prec)), 1));
        gp_filewrite1(gtos(file), GENtostr_unquoted(gconcat1(p4)));
        {
          pari_sp btop = avma;
          GEN j = gen_0;
          for (j = gen_1; gcmp(j, n) <= 0; j = gaddgs(j, 1))
          {
            gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), ellap(gel(curves, gtos(i)), gel(P, gtos(j))))));
            if (gc_needed(btop, 1))
              j = gerepilecopy(btop, j);
          }
        }
        gp_filewrite1(gtos(file), "\n");
      }
      if (gc_needed(btop, 1))
        i = gerepilecopy(btop, i);
    }
  }
  gp_fileclose(gtos(file));
  avma = ltop;
  return;
}

/* P_bd: bound for the set of primes to consider
* A_bd: -A_bd<= A <= A_bd
* B_bd: -B_bd <= B <= B_bd
* filepath: destination to write the results
* for each (A,B) with -A_bd<= A<= A_bd, -B_bd<=B<= B_bd
* construct an elliptic curve E:y^2 = x^3 + Ax + B
* compute the discriminant, conductor and rank
* then write these and the ap(E) for 1<=p<=P_bd, p prime
* to the destination file
*/
void
construct_box_curves(GEN P_bd, GEN A_bd, GEN B_bd, char *filepath, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN P = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN file = gen_0;
  long l2;
  GEN p3 = gen_0;	  /* int */
  GEN E = pol_x(fetch_user_var("E"));
  if (typ(P_bd) != t_INT)
    pari_err_TYPE("construct_box_curves",P_bd);
  if (typ(A_bd) != t_INT)
    pari_err_TYPE("construct_box_curves",A_bd);
  if (typ(B_bd) != t_INT)
    pari_err_TYPE("construct_box_curves",B_bd);
  p1 = cgetg(3, t_VEC);
  gel(p1, 1) = gen_1;
  gel(p1, 2) = icopy(P_bd);
  P = primes0(p1);
  file = stoi(gp_fileopen(filepath, "w"));
  gp_filewrite1(gtos(file), "A,B discriminant conductor rank");
  l2 = glength(P);
  {
    pari_sp btop = avma;
    long j;
    for (j = 1; j <= l2; ++j)
    {
      gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), gel(P, j))));
      avma = btop;
    }
  }
  gp_filewrite1(gtos(file), "\n");
  p3 = negi(A_bd);
  {
    pari_sp btop = avma;
    GEN j = gen_0;
    for (j = p3; gcmp(j, A_bd) <= 0; j = gaddgs(j, 1))
    {
      {
        GEN p4 = gen_0;	  /* int */
        p4 = negi(B_bd);
        {
          pari_sp btop = avma;
          GEN i = gen_0;
          for (i = p4; gcmp(i, B_bd) <= 0; i = gaddgs(i, 1))
          {
            {
              GEN p5 = gen_0;	  /* vec */
              p5 = cgetg(3, t_VEC);
              gel(p5, 1) = gcopy(i);
              gel(p5, 2) = gcopy(j);
              if (!gequal0(E = ellinit(p5, NULL, prec)))
              {
                GEN p6 = gen_0;	  /* vec */
                long l7;
                p6 = cgetg(10, t_VEC);
                gel(p6, 1) = gcopy(i);
                gel(p6, 2) = strtoGENstr(",");
                gel(p6, 3) = gcopy(j);
                gel(p6, 4) = strtoGENstr(" ");
                gel(p6, 5) = gcopy(member_disc(E));
                gel(p6, 6) = strtoGENstr(" ");
                gel(p6, 7) = gcopy(gel(ellglobalred(E), 1));
                gel(p6, 8) = strtoGENstr(" ");
                gel(p6, 9) = gcopy(gel(ellanalyticrank_bitprec(E, NULL, prec2nbits(prec)), 1));
                gp_filewrite1(gtos(file), GENtostr_unquoted(gconcat1(p6)));
                l7 = glength(P);
                {
                  pari_sp btop = avma;
                  long k;
                  for (k = 1; k <= l7; ++k)
                  {
                    gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), ellap(E, gel(P, k)))));
                    avma = btop;
                  }
                }
                gp_filewrite1(gtos(file), "\n");
              }
            }
            if (gc_needed(btop, 1))
              gerepileall(btop, 2, &i, &E);
          }
        }
      }
      if (gc_needed(btop, 1))
        gerepileall(btop, 2, &j, &E);
    }
  }
  gp_fileclose(gtos(file));
  avma = ltop;
  return;
}

static GEN
anon_0(GEN p)
{
  pari_sp ltop = avma;
  p = gerepilecopy(ltop, p);
  return p;
}

static GEN
wrap_anon_0(void * _cargs, GEN p)
{
  (void) _cargs;
  GEN _res = anon_0(p);
  return _res;
}

static long
anon_1(GEN p)	  /* bool */
{
  pari_sp ltop = avma;
  long l1;	  /* bool */
  l1 = kronecker(stoi(-3), p) == 1;
  avma = ltop;
  return l1;
}

static long
wrap_anon_1(void * _cargs, GEN p)
{
  (void) _cargs;
  long _res = anon_1(p);
  return _res;
}

/* given the curve y^2=x^3+1 construct a file
* containing all twists of the curve by D for
* -D_bd <= D <= D_bd
* and the aps for each p<=P_bd, such that p splits in
* QQ(sqrt(-3))
*/
void
construct_twists_sqrt_n3(GEN P_bd, GEN D_bd, char *filepath, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN P = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN file = gen_0;
  long l2;
  GEN p3 = gen_0;	  /* int */
  GEN E = pol_x(fetch_user_var("E"));
  if (typ(P_bd) != t_INT)
    pari_err_TYPE("construct_twists_sqrt_n3",P_bd);
  if (typ(D_bd) != t_INT)
    pari_err_TYPE("construct_twists_sqrt_n3",D_bd);
  p1 = cgetg(3, t_VEC);
  gel(p1, 1) = gen_1;
  gel(p1, 2) = icopy(P_bd);
  P = vecselapply(NULL, wrap_anon_1, NULL, wrap_anon_0, primes0(p1));
  file = stoi(gp_fileopen(filepath, "w"));
  gp_filewrite1(gtos(file), "D discriminant conductor rank");
  l2 = glength(P);
  {
    pari_sp btop = avma;
    long j;
    for (j = 1; j <= l2; ++j)
    {
      gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), gel(P, j))));
      avma = btop;
    }
  }
  gp_filewrite1(gtos(file), "\n");
  p3 = negi(D_bd);
  {
    pari_sp btop = avma;
    GEN j = gen_0;
    for (j = p3; gcmp(j, D_bd) <= 0; j = gaddgs(j, 1))
    {
      {
        GEN p4 = gen_0;	  /* vec */
        p4 = cgetg(3, t_VEC);
        gel(p4, 1) = gen_0;
        gel(p4, 2) = gcopy(j);
        if (!gequal0(E = ellinit(p4, NULL, prec)) && !gequal0(nth_power_free(j, stoi(6))))
        {
          GEN p5 = gen_0;	  /* vec */
          long l6;
          p5 = cgetg(8, t_VEC);
          gel(p5, 1) = gcopy(j);
          gel(p5, 2) = strtoGENstr(" ");
          gel(p5, 3) = gcopy(member_disc(E));
          gel(p5, 4) = strtoGENstr(" ");
          gel(p5, 5) = gcopy(gel(ellglobalred(E), 1));
          gel(p5, 6) = strtoGENstr(" ");
          gel(p5, 7) = gcopy(gel(ellanalyticrank_bitprec(E, NULL, prec2nbits(prec)), 1));
          gp_filewrite1(gtos(file), GENtostr_unquoted(gconcat1(p5)));
          l6 = glength(P);
          {
            pari_sp btop = avma;
            long k;
            for (k = 1; k <= l6; ++k)
            {
              gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), ellap(E, gel(P, k)))));
              avma = btop;
            }
          }
          gp_filewrite1(gtos(file), "\n");
        }
      }
      if (gc_needed(btop, 1))
        gerepileall(btop, 2, &j, &E);
    }
  }
  gp_fileclose(gtos(file));
  avma = ltop;
  return;
}

static GEN
anon_2(GEN p)
{
  pari_sp ltop = avma;
  p = gerepilecopy(ltop, p);
  return p;
}

static GEN
wrap_anon_2(void * _cargs, GEN p)
{
  (void) _cargs;
  GEN _res = anon_2(p);
  return _res;
}

static long
anon_3(GEN p)	  /* bool */
{
  pari_sp ltop = avma;
  long l1;	  /* bool */
  l1 = gequal1(gmodulo(p, stoi(4)));
  avma = ltop;
  return l1;
}

static long
wrap_anon_3(void * _cargs, GEN p)
{
  (void) _cargs;
  long _res = anon_3(p);
  return _res;
}

/* given the curve y^2=x^3+x construct a file
* containing all twists of the curve by D
* for -D_bd <= D <= D_bd and the aps for each 
* p<= P_bd which splits in QQ(sqrt(-1))
*/
void
construct_twists_sqrt_n1(GEN P_bd, GEN D_bd, GEN filepath, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN P = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN file = gen_0;
  long l2;
  GEN p3 = gen_0, E = pol_x(fetch_user_var("E"));
  p1 = cgetg(3, t_VEC);
  gel(p1, 1) = gen_1;
  gel(p1, 2) = gcopy(P_bd);
  P = vecselapply(NULL, wrap_anon_3, NULL, wrap_anon_2, primes0(p1));
  file = stoi(gp_fileopen(GENtostr_unquoted(filepath), "w"));
  gp_filewrite1(gtos(file), "D discriminant conductor rank");
  l2 = glength(P);
  {
    pari_sp btop = avma;
    long j;
    for (j = 1; j <= l2; ++j)
    {
      gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), gel(P, j))));
      avma = btop;
    }
  }
  gp_filewrite1(gtos(file), "\n");
  p3 = gneg(D_bd);
  {
    pari_sp btop = avma;
    GEN j = gen_0;
    for (j = p3; gcmp(j, D_bd) <= 0; j = gaddgs(j, 1))
    {
      /* enumerate over fourth-power free D */
      {
        GEN p4 = gen_0;	  /* vec */
        /* enumerate over fourth-power free D */
        p4 = cgetg(3, t_VEC);
        gel(p4, 1) = gcopy(j);
        gel(p4, 2) = gen_0;
        if (!gequal0(E = ellinit(p4, NULL, prec)) && !gequal0(nth_power_free(j, stoi(4))))
        {
          GEN p5 = gen_0;	  /* vec */
          long l6;
          p5 = cgetg(8, t_VEC);
          gel(p5, 1) = gcopy(j);
          gel(p5, 2) = strtoGENstr(" ");
          gel(p5, 3) = gcopy(member_disc(E));
          gel(p5, 4) = strtoGENstr(" ");
          gel(p5, 5) = gcopy(gel(ellglobalred(E), 1));
          gel(p5, 6) = strtoGENstr(" ");
          gel(p5, 7) = gcopy(gel(ellanalyticrank_bitprec(E, NULL, prec2nbits(prec)), 1));
          gp_filewrite1(gtos(file), GENtostr_unquoted(gconcat1(p5)));
          l6 = glength(P);
          {
            pari_sp btop = avma;
            long k;
            for (k = 1; k <= l6; ++k)
            {
              gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), ellap(E, gel(P, k)))));
              avma = btop;
            }
          }
          gp_filewrite1(gtos(file), "\n");
        }
      }
      if (gc_needed(btop, 1))
        gerepileall(btop, 2, &j, &E);
    }
  }
  gp_fileclose(gtos(file));
  avma = ltop;
  return;
}

static GEN
anon_4(GEN p)
{
  pari_sp ltop = avma;
  p = gerepilecopy(ltop, p);
  return p;
}

static GEN
wrap_anon_4(void * _cargs, GEN p)
{
  (void) _cargs;
  GEN _res = anon_4(p);
  return _res;
}

static long
anon_5(GEN p)	  /* bool */
{
  pari_sp ltop = avma;
  long l1;	  /* bool */
  l1 = gequal1(gmodulo(p, stoi(4)));
  avma = ltop;
  return l1;
}

static long
wrap_anon_5(void * _cargs, GEN p)
{
  (void) _cargs;
  long _res = anon_5(p);
  return _res;
}

void
rescale_aps_sqrt_n1_step(GEN P_bd, GEN start, GEN stepsize, GEN numsteps, GEN N, GEN n, GEN fileprefix, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN P = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN postfix = gen_0;
  p1 = cgetg(3, t_VEC);
  gel(p1, 1) = gen_1;
  gel(p1, 2) = gcopy(P_bd);
  P = vecselapply(NULL, wrap_anon_5, NULL, wrap_anon_4, primes0(p1));
  postfix = strtoGENstr(".data");
  {
    pari_sp btop = avma;
    GEN j = gen_0;
    for (j = gen_1; gcmp(j, numsteps) <= 0; j = gaddgs(j, 1))
    {
      {
        GEN p2 = gen_0;	  /* vec */
        p2 = cgetg(4, t_VEC);
        gel(p2, 1) = gcopy(fileprefix);
        gel(p2, 2) = gcopy(j);
        gel(p2, 3) = gcopy(postfix);
        rescale_aps_sqrt_n1(P, gaddgs(gadd(start, gmul(gsubgs(j, 1), stepsize)), 1), gadd(start, gmul(j, stepsize)), N, n, gconcat1(p2), prec);
      }
      if (gc_needed(btop, 1))
        j = gerepilecopy(btop, j);
    }
  }
  avma = ltop;
  return;
}

void
rescale_aps_sqrt_n1(GEN P, GEN D1, GEN D2, GEN N, GEN n, GEN filepath, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN counts_odd = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN sum_aps_odd = gen_0;
  GEN p2 = gen_0;	  /* vec */
  GEN counts_even = gen_0;
  GEN p3 = gen_0;	  /* vec */
  GEN sum_aps_even = gen_0;
  GEN p4 = gen_0;	  /* vec */
  GEN E = pol_x(fetch_user_var("E")), file = gen_0;
  {
    long l5;
    p1 = cgetg(gtos(n)+1, t_VEC);
    for (l5 = 1; gcmpsg(l5, n) <= 0; ++l5)
      gel(p1, l5) = gen_0;
  }
  counts_odd = p1;
  {
    long l6;
    p2 = cgetg(gtos(n)+1, t_VEC);
    for (l6 = 1; gcmpsg(l6, n) <= 0; ++l6)
      gel(p2, l6) = gen_0;
  }
  sum_aps_odd = p2;
  {
    long l7;
    p3 = cgetg(gtos(n)+1, t_VEC);
    for (l7 = 1; gcmpsg(l7, n) <= 0; ++l7)
      gel(p3, l7) = gen_0;
  }
  counts_even = p3;
  {
    long l8;
    p4 = cgetg(gtos(n)+1, t_VEC);
    for (l8 = 1; gcmpsg(l8, n) <= 0; ++l8)
      gel(p4, l8) = gen_0;
  }
  sum_aps_even = p4;
  {
    pari_sp btop = avma;
    GEN j = gen_0;
    for (j = gcopy(D1); gcmp(j, D2) <= 0; j = gaddgs(j, 1))
    {
      {
        GEN p9 = gen_0;	  /* vec */
        p9 = cgetg(3, t_VEC);
        gel(p9, 1) = gcopy(j);
        gel(p9, 2) = gen_0;
        if (!gequal0(E = ellinit(p9, NULL, prec)) && !gequal0(nth_power_free(j, stoi(4))))
        {
          GEN c = gen_0;
          long l10;
          c = gcopy(gel(ellglobalred(E), 1));
          l10 = glength(P);
          {
            pari_sp btop = avma;
            long i;
            for (i = 1; i <= l10; ++i)
            {
              {
                GEN idx = gen_0;
                if (gcmp(gel(P, i), gmul(N, c)) > 0)
                  break;
                idx = gceil(gdiv(gmul(n, gel(P, i)), gmul(c, N)));
                if (gequal0(gmodulo(gel(ellanalyticrank_bitprec(E, NULL, prec2nbits(prec)), 1), gen_2)))
                {
                  gel(counts_even, gtos(idx)) = gaddgs(gel(counts_even, gtos(idx)), 1);
                  gel(sum_aps_even, gtos(idx)) = gadd(gel(sum_aps_even, gtos(idx)), ellap(E, gel(P, i)));
                }
                else
                {
                  gel(counts_odd, gtos(idx)) = gaddgs(gel(counts_odd, gtos(idx)), 1);
                  gel(sum_aps_odd, gtos(idx)) = gadd(gel(sum_aps_odd, gtos(idx)), ellap(E, gel(P, i)));
                }
              }
              if (gc_needed(btop, 1))
                gerepileall(btop, 4, &counts_odd, &sum_aps_odd, &counts_even, &sum_aps_even);
            }
          }
        }
      }
      if (gc_needed(btop, 1))
        gerepileall(btop, 6, &counts_odd, &sum_aps_odd, &counts_even, &sum_aps_even, &j, &E);
    }
  }
  file = stoi(gp_fileopen(GENtostr_unquoted(filepath), "w"));
  gp_filewrite1(gtos(file), "x sum_aps_even count_even sum_aps_odd count_odd\n");
  {
    pari_sp btop = avma;
    GEN i = gen_0;
    for (i = gen_1; gcmp(i, n) <= 0; i = gaddgs(i, 1))
    {
      {
        GEN p11 = gen_0;	  /* vec */
        p11 = cgetg(11, t_VEC);
        gel(p11, 1) = gdiv(gmul(gmul(real_1(prec), N), i), n);
        gel(p11, 2) = strtoGENstr(" ");
        gel(p11, 3) = gcopy(gel(sum_aps_even, gtos(i)));
        gel(p11, 4) = strtoGENstr(" ");
        gel(p11, 5) = gcopy(gel(counts_even, gtos(i)));
        gel(p11, 6) = strtoGENstr(" ");
        gel(p11, 7) = gcopy(gel(sum_aps_odd, gtos(i)));
        gel(p11, 8) = strtoGENstr(" ");
        gel(p11, 9) = gcopy(gel(counts_odd, gtos(i)));
        gel(p11, 10) = strtoGENstr("\n");
        gp_filewrite1(gtos(file), GENtostr_unquoted(gconcat1(p11)));
      }
      if (gc_needed(btop, 1))
        i = gerepilecopy(btop, i);
    }
  }
  gp_fileclose(gtos(file));
  avma = ltop;
  return;
}

static GEN
anon_6(GEN p)
{
  pari_sp ltop = avma;
  p = gerepilecopy(ltop, p);
  return p;
}

static GEN
wrap_anon_6(void * _cargs, GEN p)
{
  (void) _cargs;
  GEN _res = anon_6(p);
  return _res;
}

static long
anon_7(GEN p, GEN N)	  /* bool */
{
  pari_sp ltop = avma;
  long l1;	  /* bool */
  l1 = kronecker(N, p) == 1;
  avma = ltop;
  return l1;
}

static long
wrap_anon_7(void * _cargs, GEN p)
{
  GEN _args = (GEN) _cargs;
  long _res = anon_7(p, gel(_args,1));
  return _res;
}

/* given a curve with complex multiplication by QQ(sqrt(N))
* construct a file containing all twists of the curve by D
* for -D_bd <= D <= D_bd and the aps for each p <= P_bd
* for which p splits in QQ(sqrt(N))
*/
void
construct_twists_sqrt_N(GEN P_bd, GEN D_bd, GEN N, char *filepath, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN P = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN j_Map = gen_0;
  GEN p2 = gen_0, p3 = gen_0, p4 = gen_0, p5 = gen_0, p6 = gen_0, p7 = gen_0, p8 = gen_0, p9 = gen_0;	  /* vec */
  GEN A = gen_0, B = gen_0, p10 = gen_0;
  if (typ(P_bd) != t_INT)
    pari_err_TYPE("construct_twists_sqrt_N",P_bd);
  if (typ(D_bd) != t_INT)
    pari_err_TYPE("construct_twists_sqrt_N",D_bd);
  if (typ(N) != t_INT)
    pari_err_TYPE("construct_twists_sqrt_N",N);
  p1 = cgetg(3, t_VEC);
  gel(p1, 1) = gen_1;
  gel(p1, 2) = icopy(P_bd);
  P = vecselapply(mkvec(N), wrap_anon_7, NULL, wrap_anon_6, primes0(p1));
  p2 = cgetg(3, t_MAT);
  gel(p2, 1) = cgetg(8, t_COL);
  gel(p2, 2) = cgetg(8, t_COL);
  gcoeff(p2, 1, 1) = gen_m2;
  p3 = cgetg(3, t_VEC);
  gel(p3, 1) = stoi(-270);
  gel(p3, 2) = stoi(1512);
  gcoeff(p2, 1, 2) = p3;
  gcoeff(p2, 2, 1) = stoi(-7);
  p4 = cgetg(3, t_VEC);
  gel(p4, 1) = stoi(-35);
  gel(p4, 2) = stoi(-98);
  gcoeff(p2, 2, 2) = p4;
  gcoeff(p2, 3, 1) = stoi(-11);
  p5 = cgetg(3, t_VEC);
  gel(p5, 1) = stoi(-9504);
  gel(p5, 2) = stoi(365904);
  gcoeff(p2, 3, 2) = p5;
  gcoeff(p2, 4, 1) = stoi(-19);
  p6 = cgetg(3, t_VEC);
  gel(p6, 1) = stoi(-219488);
  gel(p6, 2) = stoi(-39617584);
  gcoeff(p2, 4, 2) = p6;
  gcoeff(p2, 5, 1) = stoi(-43);
  p7 = cgetg(3, t_VEC);
  gel(p7, 1) = stoi(-25442240);
  gel(p7, 2) = negi(readseq("49394836848"));
  gcoeff(p2, 5, 2) = p7;
  gcoeff(p2, 6, 1) = stoi(-67);
  p8 = cgetg(3, t_VEC);
  gel(p8, 1) = stoi(-529342880);
  gel(p8, 2) = negi(readseq("4687634371504"));
  gcoeff(p2, 6, 2) = p8;
  gcoeff(p2, 7, 1) = stoi(-163);
  p9 = cgetg(3, t_VEC);
  gel(p9, 1) = negi(readseq("924354639680"));
  gel(p9, 2) = negi(readseq("34206961763303088"));
  gcoeff(p2, 7, 2) = p9;
  j_Map = listinit(gtomap(p2));
  p10 = mapget(j_Map, N);
  A = gcopy(gel(p10, 1));
  B = gcopy(gel(p10, 2));
  construct_all_twists(P, D_bd, A, B, filepath, prec);
  avma = ltop;
  return;
}

/* construct all twists of the curve y^2 = x^3 + Ax + B
* for -D_bd <= D <= D_bd and write the aps to a file
* for each prime p in P
*/
void
construct_all_twists(GEN P, GEN D_bd, GEN A, GEN B, char *filepath, long prec)	  /* void */
{
  pari_sp ltop = avma;
  GEN file = gen_0;
  long l1;
  GEN p2 = gen_0;	  /* int */
  GEN E = pol_x(fetch_user_var("E"));
  if (typ(D_bd) != t_INT)
    pari_err_TYPE("construct_all_twists",D_bd);
  if (typ(A) != t_INT)
    pari_err_TYPE("construct_all_twists",A);
  if (typ(B) != t_INT)
    pari_err_TYPE("construct_all_twists",B);
  file = stoi(gp_fileopen(filepath, "w"));
  gp_filewrite1(gtos(file), "D discriminant conductor rank");
  l1 = glength(P);
  {
    pari_sp btop = avma;
    long j;
    for (j = 1; j <= l1; ++j)
    {
      gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), gel(P, j))));
      avma = btop;
    }
  }
  gp_filewrite1(gtos(file), "\n");
  p2 = negi(D_bd);
  {
    pari_sp btop = avma;
    GEN j = gen_0;
    for (j = p2; gcmp(j, D_bd) <= 0; j = gaddgs(j, 1))
    {
      /* enumerate over square-free D */
      {
        GEN p3 = gen_0;	  /* vec */
        /* enumerate over square-free D */
        p3 = cgetg(3, t_VEC);
        gel(p3, 1) = gmul(A, gsqr(j));
        gel(p3, 2) = gmul(B, gpowgs(j, 3));
        if (!gequal0(E = ellinit(p3, NULL, prec)) && !gequal0(nth_power_free(j, gen_2)))
        {
          GEN p4 = gen_0;	  /* vec */
          long l5;
          p4 = cgetg(8, t_VEC);
          gel(p4, 1) = gcopy(j);
          gel(p4, 2) = strtoGENstr(" ");
          gel(p4, 3) = gcopy(member_disc(E));
          gel(p4, 4) = strtoGENstr(" ");
          gel(p4, 5) = gcopy(gel(ellglobalred(E), 1));
          gel(p4, 6) = strtoGENstr(" ");
          gel(p4, 7) = gcopy(gel(ellanalyticrank_bitprec(E, NULL, prec2nbits(prec)), 1));
          gp_filewrite1(gtos(file), GENtostr_unquoted(gconcat1(p4)));
          l5 = glength(P);
          {
            pari_sp btop = avma;
            long k;
            for (k = 1; k <= l5; ++k)
            {
              gp_filewrite1(gtos(file), GSTR(gconcat(strtoGENstr(" "), ellap(E, gel(P, k)))));
              avma = btop;
            }
          }
          gp_filewrite1(gtos(file), "\n");
        }
      }
      if (gc_needed(btop, 1))
        gerepileall(btop, 2, &j, &E);
    }
  }
  gp_fileclose(gtos(file));
  avma = ltop;
  return;
}

/* return true if s is not divisible by any n-th power
* return false if p^n divides s for some prime p
*/
GEN
nth_power_free(GEN s, GEN n)
{
  pari_sp ltop = avma;
  GEN e = gen_0, k = gen_0, res = gen_1;
  if (typ(s) != t_INT)
    pari_err_TYPE("nth_power_free",s);
  if (typ(n) != t_INT)
    pari_err_TYPE("nth_power_free",n);
  e = gcopy(gel(Z_factor(s), 2));
  k = stoi(glength(e));
  {
    pari_sp btop = avma;
    GEN i = gen_0;
    for (i = gen_1; gcmp(i, k) <= 0; i = gaddgs(i, 1))
    {
      res = gmulgs(res, gcmp(gel(e, gtos(i)), n) < 0);
      if (gc_needed(btop, 1))
        gerepileall(btop, 2, &res, &i);
    }
  }
  res = gerepilecopy(ltop, res);
  return res;
}

