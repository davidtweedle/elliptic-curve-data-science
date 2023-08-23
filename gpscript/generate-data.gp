/* for each E in curves write a row to filepath
 * E.label E.conductor E.rank E.ap(2) E.ap(3) E.ap(5) ... E.ap(pn)
 * where 2,3,...,pn are the primes between 1 and P_bd
 * write the result to the file at filepath
 * we assume that the curves and labels have been constructed
 * from the LMFDB database
 */
write_aps(P_bd:int,curves,labels,filepath:str)=
{    	
  my (P=primes([1,P_bd]), n=length(P), m=length(curves),\
    file=fileopen(filepath,"w"));
  if (length(labels)==0,labels = vector(length(curves),i,i));
  filewrite1(file,"label discriminant conductor rank");
  for (j=1,n,filewrite1(file,concat(" ",P[j])));
  filewrite1(file,"\n");
  for (i=1,m,
     filewrite1(file,concat(\
        [labels[i]," ",\
         curves[i].disc," ",\
         ellglobalred(curves[i])[1]," ",\
         ellanalyticrank(curves[i])[1]]));
     for (j=1,n,filewrite1(file,concat(" ",ellap(curves[i],P[j]))));
     filewrite1(file,"\n"));
  fileclose(file);
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
construct_box_curves(P_bd:int,A_bd:int,B_bd:int,filepath:str)=
{
  my (P = primes([1,P_bd]));
  my (file=fileopen(filepath,"w"));
  filewrite1(file,"A,B discriminant conductor rank");
  for (j=1,length(P),filewrite1(file,concat(" ",P[j])));
  filewrite1(file,"\n");
  for (j=-A_bd,A_bd,
    for (i=-B_bd,B_bd,
      if (E=ellinit([i,j]),
        filewrite1(file,concat(\
          [i,",",j," ",E.disc," ",ellglobalred(E)[1]," ",\
           ellanalyticrank(E)[1]]));
        for(k=1, length(P),filewrite1(file,concat(" ",ellap(E,P[k]))));
        filewrite1(file,"\n"))));
  fileclose(file);
}
/* given the curve y^2=x^3+1 construct a file
 * containing all twists of the curve by D for
 * -D_bd <= D <= D_bd
 * and the aps for each p<=P_bd, such that p splits in
 * QQ(sqrt(-3))
 */
construct_twists_sqrt_n3(P_bd:int,D_bd:int,filepath:str)=
{
  my (P=[p|p<-primes([1,P_bd]),kronecker(-3,p)==1]);
  my (file=fileopen(filepath,"w"));
  filewrite1(file, "D discriminant conductor rank");
  for (j=1, length(P), filewrite1(file,concat(" ",P[j])));
  filewrite1(file, "\n");
  for (j=-D_bd,D_bd,
    if((E=ellinit([0,j]))&&nth_power_free(j,6),
      filewrite1(file,concat(\
        [j," ",E.disc," ",ellglobalred(E)[1]," ",\
         ellanalyticrank(E)[1]]));
      for (k=1,length(P),filewrite1(file,concat(" ",ellap(E,P[k]))));
      filewrite1(file,"\n")));
  fileclose(file);
}
/* given the curve y^2=x^3+x construct a file
 * containing all twists of the curve by D
 * for -D_bd <= D <= D_bd and the aps for each 
 * p<= P_bd which splits in QQ(sqrt(-1))
 */
construct_twists_sqrt_n1(P_bd,D_bd,filepath)=
{
  my (P=[p|p<-primes([1,P_bd]),Mod(p,4)==1]);
  my (file=fileopen(filepath,"w"));
  filewrite1(file,"D discriminant conductor rank");
  for (j=1, length(P), filewrite1(file,concat(" ",P[j])));
  filewrite1(file,"\n");
  for (j=-D_bd,D_bd, \\ enumerate over fourth-power free D
    if ((E = ellinit([j,0]))&&nth_power_free(j,4),
      filewrite1(file,concat(\
        [j, " ", E.disc," ",ellglobalred(E)[1], " ",\
         ellanalyticrank(E)[1]]));
      for (k=1,length(P),filewrite1(file,concat(" ",ellap(E,P[k]))));
      filewrite1(file,"\n")));
  fileclose(file);
}

/* given a curve with complex multiplication by QQ(sqrt(N))
 * construct a file containing all twists of the curve by D
 * for -D_bd <= D <= D_bd and the aps for each p <= P_bd
 * for which p splits in QQ(sqrt(N))
 */
construct_twists_sqrt_N(P_bd:int,D_bd:int,N:int,filepath:str)=
{
  my (P=[p|p<-primes([1,P_bd]),kronecker(N,p)==1]);
  my (j_Map = Map([-2,[-270,1512];-7,[-35,-98];\
                 -11,[-9504,365904];-19,[-219488,-39617584];\
                 -43,[-25442240,-49394836848];\
                 -67,[-529342880,-4687634371504];\
                 -163,[-924354639680,-34206961763303088]]));
  my (A, B);
  [A,B] = mapget(j_Map,N);
  construct_all_twists(P,D_bd,A,B,filepath);
}
/* construct all twists of the curve y^2 = x^3 + Ax + B
 * for -D_bd <= D <= D_bd and write the aps to a file
 * for each prime p in P
 */
construct_all_twists(P,D_bd:int,A:int,B:int,filepath:str)=
{
  my (file=fileopen(filepath,"w"));
  filewrite1(file,"D discriminant conductor rank");
  for (j=1, length(P), filewrite1(file,concat(" ",P[j])));
  filewrite1(file,"\n");
  for (j=-D_bd,D_bd, \\ enumerate over square-free D
    if ((E = ellinit([A*j^2,B*j^3]))&&nth_power_free(j,2),
      filewrite1(file,concat(\
        [j, " ", E.disc," ",ellglobalred(E)[1], " ",\
         ellanalyticrank(E)[1]]));
      for (k=1,length(P),filewrite1(file,concat(" ",ellap(E,P[k]))));
      filewrite1(file,"\n")));
  fileclose(file);
}

/* return true if s is not divisible by any n-th power
 * return false if p^n divides s for some prime p
 */
nth_power_free(s:int,n:int)=
{
  my (e = factor(s)[,2],k=length(e),res=1);
  for (i=1,k,res *= (e[i]<n));
  res;
}
/* Compute vectors of "rescaled murmurations"
 * for the curves y^2 = x^3+Dx
 * and for primes <= P_bd
 * starting at D = start+1
 * and ending at start+stepsize*numsteps
 * splitting it into numsteps different files
 * N is the upper bound for p/N(E)
 * n is the number of bins
 */
rescale_aps_sqrt_n1_step(P_bd,start,stepsize,numsteps,N,n,fileprefix)=
{ 
  my (P=[p|p<-primes([1,P_bd]),Mod(p,4)==1]);
  my (postfix=".data");
  for (j=1,numsteps,
    rescale_aps_sqrt_n1(P,start+(j-1)*stepsize+1,start+j*stepsize,N,n,concat([fileprefix,j,postfix])));
}
/* Helper function for rescale_aps_sqrt_n1_step
 * for each p in P, and each D=D1,...,D2
 * consider E_D: y^2 = x^3+Dx and
 * bin the value of a_p(E) according to p/N(E_D)
 * where N(E_D) is the conductor fo E_D
 * then print the results to a file
 * N is the maximum of p/N(E_D)
 * n is the number of bins.
 */
rescale_aps_sqrt_n1(P,D1,D2,N,n,filepath)=
{
  my (counts_odd = vector(n));
  my (sum_aps_odd=vector(n));
  my (counts_even = vector(n));
  my (sum_aps_even = vector(n));
  for (j=D1,D2,
    if((E=ellinit([j,0]))&&nth_power_free(j,4),
      my (c = ellglobalred(E)[1]);\
      for (i=1,length(P),
        if(P[i]>N*c,break);
          my(idx = ceil(n*P[i]/(c*N)));
          if(Mod(ellanalyticrank(E)[1],2)==0,
            counts_even[idx]++;
            sum_aps_even[idx]+=ellap(E,P[i]);,
            counts_odd[idx]++;
            sum_aps_odd[idx]+=ellap(E,P[i])))));
  my (file=fileopen(filepath,"w"));
  filewrite1(file,"x sum_aps_even count_even sum_aps_odd count_odd\n");
  for(i=1,n,
    filewrite1(file,concat([1.0*N*i/n," ",sum_aps_even[i]," ",counts_even[i]," ",sum_aps_odd[i]," ",counts_odd[i],"\n"])));
  fileclose(file);
}