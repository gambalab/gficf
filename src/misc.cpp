#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <progress.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]


// [[Rcpp::export]]
arma::mat armaColSumFull(const arma::mat& m, int ncores=1, bool verbose=false)
{
  arma::vec s(m.n_cols, arma::fill::zeros);
  int tot = m.n_cols;
  Progress p(tot, verbose);
#pragma omp parallel for num_threads(ncores) shared(s)
  for(int i=0;i<(m.n_cols);i++) {
    s(i) = arma::sum(m.col(i));
  }
  return(s);
}

// [[Rcpp::export]]
arma::mat armaColSumSparse(const arma::sp_mat& m, int ncores=1, bool verbose=false)
{
  arma::vec s(m.n_cols, arma::fill::zeros);
  int tot = m.n_cols;
  Progress p(tot, verbose);
#pragma omp parallel for num_threads(ncores) shared(s)
  for(int i=0;i<(m.n_cols);i++) {
    s(i) = arma::sum(m.col(i));
  }
  return(s);
}

// calculate column mean and variance
// [[Rcpp::export]]
Rcpp::DataFrame colMeanVarS(const arma::sp_mat& m, int ncores=1) {
  
  arma::vec meanV(m.n_cols,arma::fill::zeros); arma::vec varV(m.n_cols,arma::fill::zeros); arma::vec nobsV(m.n_cols,arma::fill::zeros);
#pragma omp parallel for num_threads(ncores) shared(meanV,varV,nobsV)
  for(int i=0;i<(m.n_cols);i++) {
    meanV(i) = arma::mean(m.col(i));
    nobsV(i) = m.col(i).n_nonzero;
    varV(i) = arma::var(m.col(i));
  }
  return Rcpp::DataFrame::create(Named("m")=meanV,Named("v")=varV,Named("nobs",nobsV));
}

// [[Rcpp::export]]
arma::sp_mat armaManhattan(const arma::sp_mat& m, int ncores=1,bool verbose=true, bool full=false, bool diag=true)
{
  typedef arma::sp_mat::const_col_iterator iter;
  arma::sp_mat d(m.n_cols,m.n_cols);
  int tot = (full) ? m.n_cols*m.n_cols - ((diag) ? m.n_cols : 0) : (m.n_cols*m.n_cols + ((diag) ? m.n_cols : 0))/2;
  Progress p(tot, verbose);
#pragma omp parallel for num_threads(ncores) shared(d)
  for(int i=0;i<(m.n_cols);i++) {
    for(int j = diag ? i : i+1;j<m.n_cols;j++) {
      if ( !Progress::check_abort())
      {
        p.increment(); // update progress
        iter i_iter = m.begin_col(i);
        iter j_iter = m.begin_col(j);
        double num=0;
        
        while( (i_iter != m.end_col(i)) && (j_iter != m.end_col(j)) )
        {
          if(i_iter.row() == j_iter.row())
          {
            num+= std::abs((*i_iter)-(*j_iter));
            ++i_iter;
            ++j_iter;
          } else {
            if(i_iter.row() < j_iter.row())
            {
              num+=std::abs(*i_iter);
              ++i_iter;
            } else {
              num+=std::abs(*j_iter);
              ++j_iter;
            }
          }
        }
        for(; i_iter != m.end_col(i); ++i_iter) { num+=std::abs(*i_iter); }
        for(; j_iter != m.end_col(j); ++j_iter){ num+=std::abs(*j_iter); }
        d(j,i) = num;
        if (full) {d(i,j) = d(j,i);}
      }
    }
  }
  
  return(d);
}


// Correlation coefficient between the columns of a dense matrix m
// [[Rcpp::export]]
arma::sp_mat armaCorr(const arma::mat& m, int ncores=1, bool verbose=true, bool full=false, bool diag=true, bool dist=true)
{
  arma::sp_mat d(m.n_cols,m.n_cols);
  d.zeros();
  int tot = (full) ? m.n_cols*m.n_cols - ((diag) ? m.n_cols : 0) : (m.n_cols*m.n_cols + ((diag) ? m.n_cols : 0))/2;
  Progress p(tot, verbose);
#pragma omp parallel for num_threads(ncores) shared(d)
  for(int i=0;i<(m.n_cols);i++) {
    arma::vec li = m.col(i);
    for(int j = diag ? i : i+1;j<m.n_cols;j++) {
      if ( !Progress::check_abort())
      {
        p.increment(); // update progress
        arma::vec lj = m.col(j);
        d(i,j) = (dist) ? 1-arma::as_scalar(arma::cor(li,lj,0)) : arma::as_scalar(arma::cor(li,lj,0));
        d(j,i) = d(i,j);
      }
    }
  }
  return(d);
}

// [[Rcpp::export]]
Rcpp::DataFrame armaColMeans(const arma::sp_mat& m, int ncores=1,bool verbose=true)
{
  arma::vec meanV(m.n_cols,arma::fill::zeros); arma::vec meanAll(m.n_cols,arma::fill::zeros); arma::vec nobsV(m.n_cols,arma::fill::zeros);
  typedef arma::sp_mat::const_col_iterator iter;
  Progress p(m.n_cols, verbose);
#pragma omp parallel for num_threads(ncores) shared(meanV,meanAll,nobsV)
  for(unsigned int i=0;i<m.n_cols;i++) {
    if ( !Progress::check_abort())
    {
      double tot=0; int n=0;
      for(iter i_iter = m.begin_col(i); i_iter != m.end_col(i); ++i_iter) {
        tot+=(*i_iter);
        n++;
      }
      if(tot>0) {
        meanV[i] = tot/n;
        meanAll[i] = tot/m.n_cols;
        nobsV[i] = n;
      }
      p.increment(); // update progress
    }
  }
  return Rcpp::DataFrame::create(Named("mu1")=meanV,Named("mu2")=meanAll,Named("nobs",nobsV));
}

// [[Rcpp::export]]
arma::sp_mat scaleUMI(const arma::sp_mat& m, int ncores=1,bool verbose=false)
{
  typedef arma::sp_mat::const_col_iterator iter;
  typedef arma::sp_mat::const_row_col_iterator iter2;
  
  arma::sp_mat d(m.n_rows,m.n_cols);
  d.zeros();
  arma::vec totV(m.n_cols,arma::fill::zeros);
  arma::vec medianV(m.n_cols,arma::fill::zeros);
  Progress p(m.n_cols, verbose);
  
#pragma omp parallel for num_threads(ncores) shared(totV)
  for(unsigned int i=0;i<m.n_cols;i++) {
    if ( !Progress::check_abort())
    {
      double tot=0; int n=0;
      for(iter i_iter = m.begin_col(i); i_iter != m.end_col(i); ++i_iter) {
        tot+=(*i_iter);
      }
      totV(i) = tot;
      p.increment(); // update progress
    }
  }
  
  
  //medianV = arma::median(arma::mat(m),1);
  
  Progress p2(m.n_nonzero, verbose);
#pragma omp parallel for num_threads(ncores) shared(d)
  for(unsigned int i=0;i<m.n_cols;i++) {
    if ( !Progress::check_abort())
    {
      for(iter i_iter = m.begin_col(i); i_iter != m.end_col(i); ++i_iter) {
        d(i_iter.row(),i) = log10( (*i_iter) * 1000000 / totV(i) + 1);
        p2.increment(); // update progress
      }
    }
  }
  
  return(d);
}

