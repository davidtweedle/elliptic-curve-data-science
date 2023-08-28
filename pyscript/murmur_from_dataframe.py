import pandas as pd
import numpy as np
import math

def bin_murmuration(ec:pd.DataFrame, primes:pd.Index, x):
  ''' ec is a pandas DataFrame 
      primes is a subset of the columns of ec
      corresponding to the prime columns of ec
      N is the bound on p/c where p is prime, c is hte conductor on E
      n is the number of bins
  '''
  n = len(x); 
  idx = np.digitize(np.outer(np.reciprocal(1.0*ec['conductor']),primes.astype(float)),x);
  counts = np.zeros(shape=x.shape);
  sum_aps = np.zeros(shape=x.shape);
  for i in range(n):
    where = np.equal(idx,i);
    counts[i] = where.sum();
    sum_aps[i] = ec[primes].to_numpy()[where].sum();
  return sum_aps,counts;

def display_murmuration(ec,scale,primes,title,ax,xlim=None,ylim=None):
  ''' ec is a dataframe representing a list of elliptic curves
      rank has been normalized so that ec['rank'] in {-1,1}
      and -1 corresponds to odd rank, 1 to even rank
      primes is a subset of columns of ec corresponding to primes
      title is the title to be displayed
      ax is an array of 3 axes
  '''
  for axis in ax:
    axis.clear();
    if xlim:
      axis.set_xlim(*xlim)
    if ylim:
      axis.set_ylim(*ylim)
  t = primes.astype(int)*scale;
  ax[0].plot(t,ec.loc[ec['rank']==1,primes].mean(),'b.');
  ax[1].plot(t,ec.loc[ec['rank']==-1,primes].mean(),'r.');
  ax[2].plot(t,ec[primes].mul(ec['rank'],axis=0).mean(),'c.');
  ax[0].set_title(title+", even rank");
  ax[1].set_title(title+", odd rank");
  ax[2].set_title(title+", combined ranks");

def display_rescaled_murmuration(x,sum_aps_even,counts_even,sum_aps_odd,counts_odd,title,ax,xlim=None,ylim=None):
  ''' x is an array of x coordinates
      sum_aps_even is the sum_aps as binned by compute_murmur for the curves with even rank
      count_even is the binned counts of the aps as computed by compute_murmur
      sum_aps_odd is the sum_aps as binned by compute_murmur for the curves with odd rank
      count_odd is the binned counts of the aps as computed by compute_murmur for the curves with odd rank
      title is a string which will be the title of each plot
      ax is 3 Axes objects
  '''
    for axis in ax:
    axis.clear();
    if xlim:
      axis.set_xlim(*xlim)
    if ylim:
      axis.set_ylim(*ylim)
  ax[0].plot(x,sum_aps_even/counts_even,'b.');
  ax[1].plot(x,sum_aps_odd/counts_odd,'r.');
  ax[2].plot(x,(sum_aps_even-sum_aps_odd)/(counts_even+counts_odd),'c.');
  ax[0].set_title(title+", even");
  ax[1].set_title(title+", odd");
  ax[2].set_title(title+", combined ranks");
