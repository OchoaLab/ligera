#include <Rcpp.h>
#include <cerrno>
#include <vector>
#include "bed.h"

// this function calculates two terms from one pass of the genotype matrix, needed for unbiased implicit kinship estimates

// params:
// file, m_loci, n_ind: BED file path and matrix dimensions
// mean_kinship: a value to unbias estimates
// indexes_ind_rm: vector of indeces to remove (can be NULL)

// [[Rcpp::export]]
Rcpp::NumericMatrix popkin_prod_bed_cpp(
					const char* file,
					std::vector<int>::size_type m_loci,
					std::vector<int>::size_type n_ind,
					Rcpp::NumericMatrix P_R,
					double b,
					Rcpp::Nullable<Rcpp::LogicalVector> indexes_ind_R
					) {
  // file must be full path (no missing extensions)

  // will need an index for individuals right away
  std::vector<int>::size_type j;
  // and individuals kept (go up to n_ind_kept instead of n_ind; differ when individuals are removed)
  std::vector<int>::size_type j_kept;
  // and for covariates
  std::vector<int>::size_type l;
  // recurrent product of j_kept * k_covars, used tp access P and write to KP
  std::vector<int>::size_type jK;

  // sort out indexes_ind_rm, which in C++ is a mess (but we want reasonable R-like behavior on the outside)
  bool do_ind_filt = false;
  // a pure C++ version of the input (negated), for ease in the loops
  bool* indexes_ind_rm = new bool[ n_ind ];
  if ( indexes_ind_R.isNotNull() ) {
    // change this boolean, will be easier to test in loops
    do_ind_filt = true;
    // also extract non-null values of indexes, but this is still an R variable
    Rcpp::LogicalVector indexes_ind_R_good( indexes_ind_R );
    // negate and copy values over to C++ version
    for (j = 0; j < n_ind; j++) {
      // have to do it in this awkward way, since indexes_ind_R_good is not type bool
      indexes_ind_rm[ j ] = ( indexes_ind_R_good[ j ] == FALSE );
    }
  }

  // we won't validate things that came from R, let's just assume they're all consistent
  // so here P_R gives us dimensions
  // n_ind_kept equals the number of individuals kept according to indexes_ind_R
  std::vector<int>::size_type n_ind_kept = P_R.nrow();
  std::vector<int>::size_type k_covars = P_R.ncol();
  // the product of both
  std::vector<int>::size_type size_P = n_ind_kept * k_covars;

  // might as well compute column sums of P, needed at the end
  double* csP = new double[ k_covars ];
  // first initialize all to zeroes
  for ( l = 0; l < k_covars; l++ )
    csP[ l ] = 0.0;
  
  // let's make a pure C++ copy of P, since we access it so much
  double* P = new double[ size_P ];
  // copy individual values
  for ( j_kept = 0; j_kept < n_ind_kept; j_kept++ ) {
    jK = j_kept * k_covars;
    for ( l = 0; l < k_covars; l++ ) {
      // copy P over
      P[ jK + l ] = P_R( j_kept, l );
      // compute column sums of P
      csP[ l ] += P_R( j_kept, l );
    }
  }
  
  // the output matrix, KP, has the same dimensions as P
  // there's the pure C++ version, which we'll update often
  double* KP = new double[ size_P ];
  // initialize with zeroes
  for ( j = 0; j < size_P; j++ )
    KP[ j ] = 0.0;

  // and an intermediate vector, this one is small
  double* xP = new double[ k_covars ];
  // initialize with zeroes
  for ( l = 0; l < k_covars; l++ )
    xP[ l ] = 0.0;
  
  
  // open input file in "binary" mode
  FILE *file_stream = fopen( file, "rb" );
  // die right away if needed, before initializing buffers etc
  if ( file_stream == NULL ) {
    // send error message to R
    Rcpp::stop( "Could not open BED file `%s` for reading: %s", file, strerror( errno ) );
  }

  /////////////////////////
  // check magic numbers //
  /////////////////////////

  // for header only
  unsigned char buffer_header[3];
  // for extra sanity checks, keep track of bytes actually read (to recognize truncated files)
  // reuse this one for genotypes below
  std::vector<int>::size_type n_buf_read;
  
  // read header bytes (magic numbers)
  n_buf_read = fread( buffer_header, sizeof(unsigned char), 3, file_stream );
  // this might just indicate an empty file
  if ( n_buf_read != 3 ) {
    // wrap up everything properly
    fclose( file_stream ); // close file
    // now send error message to R
    Rcpp::stop("Input BED file did not have a complete header (3-byte magic numbers)!");
  }
  
  // require that they match our only supported specification of locus-major order and latest format
  // was using strcmp but there are funky issues (wants signed, but we don't really want order anyway, just test for equality)
  // use explicit loop instead
  int pos;
  for (pos = 0; pos < 3; pos++) {
    if ( plink_bed_byte_header[pos] != buffer_header[pos] ) {
      // wrap up everything properly
      fclose( file_stream ); // close file
      // now send error message to R
      Rcpp::stop("Input BED file is not in supported format.  Either magic numbers do not match, or requested sample-major format is not supported.  Only latest locus-major format is supported!");
    }
  }

  ////////////////////
  // read genotypes //
  ////////////////////
  
  // number of columns (bytes) in input (for buffer), after byte compression
  // size set for full row, but overloaded used first for this header comparison
  // chose std::vector<int>::size_type to have it match n_buf_read value returned by fread
  std::vector<int>::size_type n_buf = ( n_ind + 3 ) / 4;
  // initialize row buffer
  unsigned char* buffer = new unsigned char[ n_buf ];
  
  // navigate data and process
  std::vector<int>::size_type i;
  std::vector<int>::size_type k;
  unsigned char buf_k; // working of buffer at k'th position
  unsigned char xij; // copy of extracted genotype
  for (i = 0; i < m_loci; i++) {
    
    // reset this one after every locus
    for ( l = 0; l < k_covars; l++ )
      xP[ l ] = 0;
    
    // read whole row into buffer
    n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
    
    // always check that file was not done too early
    if ( n_buf_read != n_buf ) {
      // wrap up everything properly
      delete [] buffer; // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      Rcpp::stop( "Truncated file: row %ld terminated at %ld bytes, expected %ld.", i+1, n_buf_read, n_buf ); // convert to 1-based coordinates
    }

    // process buffer now!

    // always reset these at start of row
    j = 0; // individuals
    j_kept = 0; // kept inds (different because of skips)
    
    // navigate buffer positions k (not individuals j)
    for (k = 0; k < n_buf; k++) {
      
      // copy down this value, which will be getting edited
      buf_k = buffer[k];

      // navigate the four positions
      // pos is just a dummy counter not really used except to know when to stop
      // update j too, accordingly
      // lastly, always shift packed data, throwing away genotype we just processed
      for (pos = 0; pos < 4; pos++, j++, buf_k = buf_k >> 2) {

	if (j < n_ind) {

	  // skip individual from calculations if there's a filter and its removed in this filter
	  if ( do_ind_filt && indexes_ind_rm[ j ] )
	    continue;

	  // extract current genotype using this mask
	  // (3 == 00000011 in binary)
	  xij = buf_k & 3;
	
	  // handle cases (in original encoding; canonical encoding in comments, followed by x-1, setting NAs to zeroes too)
	  // Homozygotes are more common, so test for those first
	  if (xij == 0) {
	    // xij = 2; // 0 -> 2 -> 1
	    // so in this case all of the values of P, for this individual, get *added* onto xP
	    // update this index
	    jK = j_kept * k_covars;
	    for ( l = 0; l < k_covars; l++ )
	      xP[ l ] += P[ jK + l ];
	  } else if (xij == 3) {
	    // xij = 0; // 3 -> 0 -> -1
	    // and in this case all of the values of P, for this individual, get *subtracted* from xP
	    // update this index
	    jK = j_kept * k_covars;
	    for ( l = 0; l < k_covars; l++ )
	      xP[ l ] -= P[ jK + l ];
	  }
	  // else nothing, as the remaining cases get multiplied by zeroes
	  // else if (xij == 2) {
	  //   // xij = 1; // 2 -> 1 -> 0
	  // } else { // only case left is NA
	  //   // xij = 0; // 1 -> NA -> 0
	  // }
	  
	  // increment individuals that weren't skipped, and only after we were done processing the individual
	  j_kept++;
	} else {
	  // when j is out of range, we're in the padding data now
	  // as an extra sanity check, the remaining data should be all zero (that's how the encoding is supposed to work)
	  // non-zero values would strongly suggest that n_ind was not set correctly
	  if (buf_k != 0) {
	    // wrap up everything properly
	    delete [] buffer; // free buffer memory
	    fclose( file_stream ); // close file
	    // now send error message to R
	    Rcpp::stop( "Row %ld padding was non-zero.  Either the specified number of individuals is incorrect or the input file is corrupt!", i+1 ); // convert to 1-based coordinates
	  }
	}
      }
      // finished byte
      
    }
    // finished row

    // process buffer again!
    // this is for the second part of the matrix product

    // always reset these at start of row
    j = 0; // individuals
    j_kept = 0; // kept inds (different because of skips)

    // navigate buffer positions k (not individuals j)
    for (k = 0; k < n_buf; k++) {
      
      // copy down this value, which will be getting edited
      buf_k = buffer[k];

      // navigate the four positions
      // pos is just a dummy counter not really used except to know when to stop
      // update j too, accordingly
      // lastly, always shift packed data, throwing away genotype we just processed
      for (pos = 0; pos < 4; pos++, j++, buf_k = buf_k >> 2) {

	if (j < n_ind) {

	  // skip individual from calculations if there's a filter and its removed in this filter
	  if ( do_ind_filt && indexes_ind_rm[ j ] )
	    continue;
	  
	  // extract current genotype using this mask
	  // (3 == 00000011 in binary)
	  xij = buf_k & 3;
	
	  // handle cases (in original encoding; canonical encoding in comments, followed by x-1, setting NAs to zeroes too)
	  // Homozygotes are more common, so test for those first
	  if (xij == 0) {
	    // xij = 2; // 0 -> 2 -> 1
	    // so in this case all of the values of xP, for this individual, get *added* onto KP
	    // update this index
	    jK = j_kept * k_covars;
	    for ( l = 0; l < k_covars; l++ )
	      KP[ jK + l ] += xP[ l ];
	  } else if (xij == 3) {
	    // xij = 0; // 3 -> 0 -> -1
	    // and in this case all of the values of xP, for this individual, get *subtracted* from KP
	    // update this index
	    jK = j_kept * k_covars;
	    for ( l = 0; l < k_covars; l++ )
	      KP[ jK + l ] -= xP[ l ];
	  }
	  // else nothing, as the remaining cases get multiplied by zeroes
	  // else if (xij == 2) {
	  //   // xij = 1; // 2 -> 1 -> 0
	  // } else { // only case left is NA
	  //   // xij = 0; // 1 -> NA -> 0
	  // }
	  
	  // increment individuals that weren't skipped, and only after we were done processing the individual
	  j_kept++;
	} else {
	  // when j is out of range, we're in the padding data now
	  // as an extra sanity check, the remaining data should be all zero (that's how the encoding is supposed to work)
	  // non-zero values would strongly suggest that n_ind was not set correctly
	  if (buf_k != 0) {
	    // wrap up everything properly
	    delete [] buffer; // free buffer memory
	    fclose( file_stream ); // close file
	    // now send error message to R
	    Rcpp::stop( "Row %ld padding was non-zero.  Either the specified number of individuals is incorrect or the input file is corrupt!", i+1 ); // convert to 1-based coordinates
	  }
	}
      }
      // finished byte
      
    }
    // finished row

  }
  // finished matrix!

  // let's check that file was indeed done
  n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
  // wrap up regardless
  delete [] buffer; // free buffer memory
  // and more troubleshooting messages (for windows)
  if ( fclose( file_stream ) != 0 )
    Rcpp::stop("Input BED file stream close failed!");
  if ( n_buf_read != 0 ) {
    // now send error message to R
    Rcpp::stop("Input BED file continued after all requested rows were read!  Either the specified the number of loci was too low or the input file is corrupt!");
  }

  // the R version of KP to return
  Rcpp::NumericMatrix KP_R( n_ind_kept, k_covars );
  // copy individual values back
  // also complete final step in KP calculation
  // - normalize by m_loci
  // - subtract b * colSums( P ) along the rows
  // - normalize by ( 1 - b )
  for ( j_kept = 0; j_kept < n_ind_kept; j_kept++ ) {
    jK = j_kept * k_covars;
    for ( l = 0; l < k_covars; l++ )
      KP_R( j_kept, l ) = ( KP[ jK + l ] / m_loci - b * csP[ l ] ) / ( 1 - b );
  }
  
  // free intermediate objects now
  delete [] indexes_ind_rm;
  delete [] P;
  delete [] csP;
  delete [] xP;
  delete [] KP;
  
  // return R version of KP matrix
  return KP_R;
}
