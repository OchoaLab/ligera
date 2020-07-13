#include <Rcpp.h>
#include <cerrno>
#include <vector>
#include "bed.h"

// this function calculates two terms from one pass of the genotype matrix, needed for unbiased implicit kinship estimates

// params:
// file, m_loci, n_ind: BED file path and matrix dimensions
// mean_kinship: a value to unbias estimates
// indexes_ind_R: vector of indeces to keep (can be NULL)

// [[Rcpp::export]]
Rcpp::List get_b_inbr_bed_cpp(
			      const char* file,
			      std::vector<int>::size_type m_loci,
			      std::vector<int>::size_type n_ind,
			      double mean_kinship,
			      Rcpp::Nullable<Rcpp::LogicalVector> indexes_ind_R
			      ) {
  // file must be full path (no missing extensions)

  // will need an index for individuals right away
  std::vector<int>::size_type j;

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
    
  // to have output length right, count the number of individuals we're keeping
  std::vector<int>::size_type n_ind_kept = n_ind; // decrement
  if ( do_ind_filt )
    for (j = 0; j < n_ind; j++)
      if ( indexes_ind_rm[ j ] )
	n_ind_kept--;
  
  // open input file in "binary" mode
  FILE *file_stream = fopen( file, "rb" );
  // die right away if needed, before initializing buffers etc
  if ( file_stream == NULL ) {
    // send error message to R
    char msg[100];
    sprintf(msg, "Could not open BED file `%s` for reading: %s", file, strerror( errno ));
    Rcpp::stop(msg);
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

  // preallocate the things we want:
  
  // (1) a summary of the MAF variance across the genome
  double b = 0.0;
  // a denominator in case there are loci that are completely NA (after removal of individuals)
  std::vector<int>::size_type m_loci_obs = m_loci; // decrement as we see NAs
  // intermediates for b, to calculate allele frequencies
  int x_sum = 0;
  std::vector<int>::size_type x_num = n_ind; // decrement NAs, just as for inbr_num above (assuming NAs are rare, this is faster)
  double x_mean;

  // (2) intermeditates for inbr
  // an uncorrected inbreeding coefficient estimate
  std::vector<int>::size_type* inbr_sum = new std::vector<int>::size_type[ n_ind ];
  // the number of non-NA cases, per individual
  std::vector<int>::size_type* inbr_num = new std::vector<int>::size_type[ n_ind ];
  // initialize
  // NOTE: if there is an individual filter, we'll just skip calculating those values in the big genotype loop (not here)
  for (j = 0; j < n_ind; j++) {
    inbr_sum[ j ] = 0;
    // all start as m_loci, will decrement later as we encounter NAs
    inbr_num[ j ] = m_loci;
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
  std::vector<int>::size_type k; // to match n_buf type
  unsigned char buf_k; // working of buffer at k'th position
  unsigned char xij; // copy of extracted genotype
  for (i = 0; i < m_loci; i++) {
    
    // read whole row into buffer
    n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
    
    // always check that file was not done too early
    if ( n_buf_read != n_buf ) {
      // wrap up everything properly
      delete [] buffer; // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      char msg[100];
      sprintf(msg, "Truncated file: row %ld terminated at %ld bytes, expected %ld.", i+1, n_buf_read, n_buf); // convert to 1-based coordinates
      Rcpp::stop(msg);
    }

    // process buffer now!

    // always reset these at start of row
    j = 0; // individuals

    // reset things we need to get allele frequency
    // NOTE: when there are filtered individuals, should start from their total (only NAs decrement this count)
    x_sum = 0;
    x_num = n_ind_kept; // decrement NAs, just as for inbr_num above (assuming NAs are rare, this is faster)
      
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
	
	  // handle cases (in original encoding; canonical encoding in comments)
	  // Homozygotes are more common, so test for those first
	  // inbr_sum is incremented upon each homozygote only!
	  // x_sum increments by ordinary xij, no need to overwrite xij at all
	  if (xij == 0) {
	    // xij = 2; // 0 -> 2
	    x_sum += 2;
	    inbr_sum[ j ]++;
	  } else if (xij == 3) {
	    // xij = 0; // 3 -> 0
	    inbr_sum[ j ]++;
	  } else if (xij == 2) {
	    // xij = 1; // 2 -> 1
	    x_sum++;
	  } else { // only case left is NA
	    // xij = 0; // 1 -> NA
	    // decrement non-NA count for this individual and locus, respectively
	    inbr_num[ j ]--;
	    x_num--;
	  }
	  
	} else {
	  // when j is out of range, we're in the padding data now
	  // as an extra sanity check, the remaining data should be all zero (that's how the encoding is supposed to work)
	  // non-zero values would strongly suggest that n_ind was not set correctly
	  if (buf_k != 0) {
	    // wrap up everything properly
	    delete [] buffer; // free buffer memory
	    fclose( file_stream ); // close file
	    // now send error message to R
	    char msg[200];
	    sprintf(msg, "Row %ld padding was non-zero.  Either the specified number of individuals is incorrect or the input file is corrupt!", i+1); // convert to 1-based coordinates
	    Rcpp::stop(msg);
	  }
	}
      }
      // finished byte
      
    }
    // finished row

    // finish processing terms for b
    // avoid division by zero (don't do anything if there were no non-NA observations; extremely rare)
    if ( x_num != 0 ) {
      // calculate mean genotype now that we're done with this one locus
      x_mean = static_cast<double>(x_sum) / x_num;
      // update b running sum
      b += x_mean * (2.0 - x_mean);
    } else {
      // NOTE: no observations treats this locus' contribution to b as zero, only happens if every individual that wasn't removed was NA
      // in this case we decrement the denominator of b, for when we normalize it into a mean
      m_loci_obs--;
    }
    
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

  // do some final processing of b
  // this includes all final steps, including the mean_kinship correction
  b = ( 1.0 - ( b / m_loci_obs ) - mean_kinship ) / ( 1.0 - mean_kinship );
  // create R version
  Rcpp::NumericVector bR(1);
  bR[ 0 ] = b;

  // final inbreeding vector
  // NOTE: contains only individuals that were not removed!
  Rcpp::NumericVector inbr(n_ind_kept);
  // the inbreeding data also now gets normalized fully, which requires the b above
  // we reuse i to indicate the output index
  i = 0;
  for (j = 0; j < n_ind; j++) {
    // skip individual if needed
    if ( do_ind_filt && indexes_ind_rm[ j ] )
      continue;
    // otherwise add to vector
    inbr[ i ] = ( ( 2.0 * inbr_sum[ j ] ) / inbr_num[ j ] - 1.0 - b ) / ( 1.0 - b );
    // increment output counter (only happens if individual wasn't skipped)
    i++;
  }

  // free intermediate objects now
  delete [] inbr_sum;
  delete [] inbr_num;
  delete [] indexes_ind_rm;

  // should return b and inbr together in an R List
  Rcpp::List ret_list = Rcpp::List::create(
					   Rcpp::Named("b") = bR,
					   Rcpp::Named("inbr") = inbr
					   );
  
  // return genotype matrix
  return ret_list;
}
